"""Python/numba port of the simulation loop of sqpid_sim.m, used only for the global gain search.
All results reported in the paper come from the Octave scripts.

Plant      : y = G0(s) [ u(t-L) + d(t) ],  G0 = kg*6.63/(s^2+14(1+delta)s+6.63)
Feedback   : ym = y + n   (measurement noise is fed back)
Actuator   : |u| <= Usat for every controller
Controllers (law / predictor flags):
  law 0 : linear PID, filtered derivative, clamp anti-windup
  law 1 : basic QPID (pi-normalised error, constant 0.5 weight, inner tanh, Kr/Ks)
  law 2 : enhanced QPID (state-dependent cos^2 / sin^2 weights, adaptive e_max, optional weights on D)
  sp 0  : no predictor;  sp 1 : Smith predictor with mismatch filter Fr (Tf=0 -> conventional Smith)
"""
import math
import numpy as np
from numba import njit
from scipy.signal import cont2discrete, tf2ss

Ts = 0.01
PSLOW = 0.5


_DSS = {}


def dss(kg=1.0, delta=0.0):
    key = (round(kg, 9), round(delta, 9))
    if key in _DSS:
        return _DSS[key]
    _DSS[key] = _dss(kg, delta)
    return _DSS[key]


def _dss(kg=1.0, delta=0.0):
    A, B, C, D = tf2ss([kg * 6.63], [1, 14 * (1 + delta), 6.63])
    Ad, Bd, Cd, Dd, _ = cont2discrete((A, B, C, D), Ts, method='zoh')
    return np.array([Ad[0, 0], Ad[0, 1], Ad[1, 0], Ad[1, 1], Bd[0, 0], Bd[1, 0], Cd[0, 0], Cd[0, 1]])


@njit(cache=True)
def core(pl, md, nL, nLm, r, d, nz, Usat, law, sp, g, y, u, aux):
    # g: [Kp, Ki, Kd, Ks, Tf, beta, Nf, alpha, Umax, emax, Ilim, wq, adapt_emax, mmax, dweight, smooth]
    N = r.shape[0]
    Kp = g[0]; Ki = g[1]; Kd = g[2]; Ks = g[3]; Tf = g[4]; beta = g[5]; Nf = g[6]
    al = g[7]; Umax = g[8]; emax0 = g[9]; Ilim = g[10]; wq = g[11]; adapt = g[12]; mmax = g[13]
    dwt = g[14]; smooth = g[15]; npow = g[16]; tm = g[17]; ew = g[18]
    a11, a12, a21, a22, b1, b2, c1, c2 = pl
    m11, m12, m21, m22, n1, n2, k1, k2 = md
    ubuf = np.zeros(N + nL + 2); ymbuf = np.zeros(N + nLm + 2)
    x1 = 0.0; x2 = 0.0; z1 = 0.0; z2 = 0.0
    I = 0.0; D = 0.0; xp = 0.0
    a1 = (2 - Nf * Ts) / (2 + Nf * Ts); b0 = Nf * 2 / (2 + Nf * Ts)
    # mismatch filters (Tustin)
    if Tf > 0:
        den0 = Ts + 2 * Tf; fa = -(Ts - 2 * Tf) / den0
        fl0 = (Ts + 2 * beta) / den0; fl1 = (Ts - 2 * beta) / den0   # lead-lag (1+beta s)/(1+Tf s)
        lp0 = Ts / den0; lp1 = Ts / den0                             # low-pass 1/(1+Tf s)
    else:
        fa = 0.0; fl0 = 1.0; fl1 = 0.0; lp0 = 1.0; lp1 = 0.0
    mF = 0.0; mL = 0.0; mprev = 0.0
    emax = emax0; rprev = 0.0; pd = 1.0; rm = 0.0; dref = 1.0
    pi = math.pi
    for k in range(N):
        yt = c1 * x1 + c2 * x2
        y[k] = yt
        if abs(yt) > 1e6:
            for j in range(k, N):
                y[j] = np.nan
            return
        ym = yt + nz[k]
        rk = r[k]
        if abs(rk - rprev) > 0.02:
            dref = abs(rk - rm)
            if dref < 0.05: dref = 0.05
            if adapt > 0: emax = max(emax0, abs(rk - ym))
        rprev = rk
        if tm > 0:
            rat = abs(rk - rm) / dref
            if rat > 1: rat = 1.0
            p0r = math.cos(0.5 * pi * rat) ** 2
            rm += (Ts / tm) * (rk - rm)
        else:
            p0r = 1.0
        if sp == 0:
            yfb = ym
        else:
            yh = k1 * z1 + k2 * z2
            m = ym - ymbuf[k]
            mF = fa * mF + fl0 * m + fl1 * mprev
            mL = fa * mL + lp0 * m + lp1 * mprev
            mprev = m
            if mmax > 0:
                ph = abs(mL) / mmax
                if ph > 1: ph = 1.0
                pd = math.sin(0.5 * pi * ph) ** 2
                mc = mL + pd * (mF - mL)
            else:
                mc = mF
            yfb = yh + mc
            ymbuf[k + nLm] = yh
        e = rk - yfb
        if law == 0:
            D = a1 * D + b0 * (e - xp); xp = e
            v = Kp * e + Ki * I + Kd * D
            uk = v
            if uk > Usat: uk = Usat
            if uk < -Usat: uk = -Usat
            if uk == v or e * v < 0:
                I += e * Ts
        else:
            th = pi * e / emax
            if th > pi: th = pi
            if th < -pi: th = -pi
            if tm > 0 or ew > 0:
                p0 = p0r
                if ew > 0: p0 = p0 * math.cos(th / 2) ** 2
            elif npow > 0:
                p0 = (math.cos(th / 2) ** 2) ** npow
            elif wq > 0:
                thI = pi * e / (wq * emax)
                if thI > pi: thI = pi
                if thI < -pi: thI = -pi
                p0 = math.cos(thI / 2) ** 2
            else:
                p0 = 0.5
            I += Ki * Ts * p0 * th
            if I > Ilim: I = Ilim
            if I < -Ilim: I = -Ilim
            D = a1 * D + b0 * (th - xp); xp = th
            if dwt > 0:
                dterm = Kd * D * (1 - p0) ** dwt if dwt < 50 else Kd * D
            else:
                dterm = Kd * D
            thu = Kp * th + I + dterm
            if law == 1:
                v = rk + Ks * e + Umax * math.tanh(al * thu)
            else:
                v = rk + Ks * e + Umax * al * thu
            if smooth > 0:
                uk = Usat * math.tanh(v / Usat)
            else:
                uk = v
                if uk > Usat: uk = Usat
                if uk < -Usat: uk = -Usat
        if sp == 1:
            z1n = m11 * z1 + m12 * z2 + n1 * uk
            z2 = m21 * z1 + m22 * z2 + n2 * uk
            z1 = z1n
        u[k] = uk
        ubuf[k + nL] = uk
        ud = ubuf[k] + d[k]
        x1n = a11 * x1 + a12 * x2 + b1 * ud
        x2 = a21 * x1 + a22 * x2 + b2 * ud
        x1 = x1n
        aux[k] = p0r


def beta_of(Tf, L):
    """Lead time constant that cancels the slow plant pole in the disturbance response."""
    return (1.0 - (1.0 - PSLOW * Tf) * math.exp(-PSLOW * L)) / PSLOW


CTRL = {  # name: (law, sp, filtered)
    'PID': (0, 0, 0), 'QPID': (1, 0, 0), 'S-PID': (0, 1, 0), 'S-QPID': (1, 1, 0),
    'FS-PID': (0, 1, 1), 'FS-QPID': (1, 1, 1), 'EQ': (2, 1, 1),
}
MD = dss()


import os
OCT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '')   # seed-1 random sequences exported from Octave
NZ_OCT = np.loadtxt(OCT + 'noise_unit.txt').ravel()
WN_OCT = np.loadtxt(OCT + 'turb_unit.txt').ravel()


_G = {}


def gust(N, amp=0.5, t_on=8.0, t_off=14.0, sig=0.2, seed=0, tc=0.3):
    key = (N, amp, t_on, t_off, sig, seed, tc)
    if key not in _G:
        if len(_G) > 5000: _G.clear()
        _G[key] = _gust(N, amp, t_on, t_off, sig, seed, tc)
    return _G[key]


def _gust(N, amp=0.5, t_on=8.0, t_off=14.0, sig=0.2, seed=0, tc=0.3):
    """Torque-equivalent wind gust: mean step amp on [t_on,t_off) with coloured turbulence (std = sig*amp)."""
    t = np.arange(N) * Ts
    d = np.zeros(N)
    w = (t >= t_on) & (t < t_off)
    d[w] = amp
    if sig > 0:
        xi = np.zeros(N); a = math.exp(-Ts / tc); s = math.sqrt(1 - a * a)
        wn = WN_OCT[:N] if seed == 0 else np.random.RandomState(seed).randn(N)
        for k in range(1, N):
            xi[k] = a * xi[k - 1] + s * wn[k]
        d[w] += amp * sig * xi[w]
    return d


_NZ = {}


def _noise(N, noise, nseed):
    key = (N, noise, nseed)
    if key not in _NZ:
        if len(_NZ) > 5000: _NZ.clear()
        _NZ[key] = noise * (NZ_OCT[:N] if nseed == 0 else np.random.RandomState(nseed).randn(N))
    return _NZ[key]


def run(name, L, p, T=24.0, R=1.0, Lerr=1.0, kg=1.0, delta=0.0, gamp=0.5, gsig=0.2, gseed=0, noise=0.005, nseed=0,
        Usat=20.0, emax=1.0, opts=None, r=None, d=None, full=False, t_on=8.0, t_off=14.0):
    """p: PID [Kp,Ki,Kd]; QPID/S-QPID [Kp,Ki,Kd,Ks]; FS-* add Tf; EQ: [Kp,Ki,Kd,Ks,Tf] + opts."""
    law, sp, filt = CTRL[name]
    N = int(round(T / Ts)); t = np.arange(N) * Ts
    if r is None: r = R * np.ones(N)
    if d is None: d = gust(N, gamp, t_on, t_off, gsig, gseed)
    nz = _noise(N, noise, nseed)
    g = np.zeros(19)
    g[0:3] = p[0:3]
    g[6] = 20.0 if law == 0 else 10.0
    g[7] = 0.10; g[8] = 20.0; g[9] = emax; g[10] = 6.0
    if law > 0: g[3] = p[3]
    if filt:
        g[4] = p[4] if law > 0 else p[3]
        g[5] = beta_of(g[4], L)
    if opts:
        for kk, vv in opts.items():
            g[{'wq': 11, 'adapt': 12, 'mmax': 13, 'dwt': 14, 'smooth': 15, 'npow': 16, 'tm': 17, 'ew': 18, 'Nf': 6, 'alpha': 7, 'Ilim': 10, 'beta': 5, 'Tf': 4}[kk]] = vv
    y = np.zeros(N); u = np.zeros(N); aux = np.zeros(N)
    core(dss(kg, delta), MD, int(round(L * Lerr / Ts)), int(round(L / Ts)), r, d, nz, Usat, law, sp, g, y, u, aux)
    return metrics(t, r, y, u, R, t_on, full, aux, t_off)


def metrics(t, r, y, u, R, t_on, full=False, aux=None, t_off=14.0):
    o = {}
    if not np.all(np.isfinite(y)):
        o.update(ISE=np.inf, IAE=np.inf, ITAE=np.inf, OS=np.inf, Ts=np.inf, rise=np.nan, su=np.inf, sy=np.inf,
                 gpk=np.inf, gIAE=np.inf, trIAE=np.inf, umax=np.inf, ok=False, div=True, trec=np.inf, ess=np.inf, sITAE=np.inf)
        return o
    e = r - y
    o['ISE'] = float(np.sum(e ** 2) * Ts); o['IAE'] = float(np.sum(np.abs(e)) * Ts); o['ITAE'] = float(np.sum(t * np.abs(e)) * Ts)
    m = t <= (t_on - 0.5); o['T_end'] = float(t[-1]); ye = y[m] / R; te = t[m]
    i10 = np.argmax(ye >= 0.1) if np.any(ye >= 0.1) else -1
    i90 = np.argmax(ye >= 0.9) if np.any(ye >= 0.9) else -1
    o['rise'] = float(te[i90] - te[i10]) if i90 >= 0 else float('nan')
    idx = np.where(np.abs(ye - 1) > 0.02)[0]
    o['Ts'] = float(te[idx.max()]) if len(idx) else 0.0
    o['OS'] = float(max((ye.max() - 1) * 100, 0))
    o['umax'] = float(np.abs(u).max())
    g = t >= t_on
    o['gpk'] = float(np.abs(e[g]).max()); o['gIAE'] = float(np.sum(np.abs(e[g])) * Ts)
    o['trIAE'] = float(np.sum(np.abs(e[~g])) * Ts)
    o['ess'] = float(abs(np.mean(e[(t >= t_off - 1.0) & (t < t_off)])))
    tau = np.where(t < t_on, t, np.where(t < t_off, t - t_on, t - t_off))
    o['sITAE'] = float(np.sum(tau * np.abs(e)) * Ts)
    q = (t >= t_on - 2.5) & (t < t_on)      # quiet window before the gust
    o['su'] = float(np.sqrt(np.mean(np.diff(u[q]) ** 2))); o['sy'] = float(np.std(e[q]))
    tail = np.abs(e[t > t[-1] - 2]).max()
    o['ok'] = bool(tail < 0.05 * R + 0.02); o['div'] = bool(tail > 0.5 * R)
    o['tv'] = float(np.sum(np.abs(np.diff(u))))
    if full:
        o['y'] = y; o['u'] = u; o['t'] = t; o['aux'] = aux
    return o
