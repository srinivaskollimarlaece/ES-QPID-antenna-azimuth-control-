"""Manuscript figures drawn from the Octave results file results_option<opt>.mat."""
import sys
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.linewidth': 0.8,
                     'axes.grid': True, 'grid.color': '#d9d9d9', 'grid.linewidth': 0.6,
                     'legend.frameon': True, 'legend.framealpha': 0.95, 'legend.edgecolor': '#bbbbbb',
                     'legend.fontsize': 7.8, 'mathtext.fontset': 'dejavusans', 'axes.titlesize': 9,
                     'axes.axisbelow': True})
STY = {'PID': ('#d55e00', (0, (5, 2, 1, 2)), 1.3), 'S-PID': ('#0072b2', (0, (4, 2)), 1.3),
       'FS-PID': ('#009e73', (0, (1, 1)), 1.5), 'S-QPID': ('#a50f2d', '-', 1.9),
       'ES-QPID': ('#a50f2d', '-', 1.9), 'QPID': ('#777777', (0, (1, 1)), 1.2),
       'ES-noweight': ('#8e6bb1', '-', 1.2), 'ES-nofilt': ('#5aa2d1', '-', 1.2), 'ES-noDw': ('#9a9a3a', '-', 1.2)}
GUST = '#fff1cc'
LS = [0.3, 0.8, 1.5]


def load(opt, d):
    R = sio.loadmat(f'{d}/results_option{opt}.mat', squeeze_me=True, struct_as_record=False)['R']
    return R


def names(x):
    return [str(v) for v in np.atleast_1d(x)]


def save(fig, out, name):
    fig.savefig(f'{out}/{name}.png', dpi=300, bbox_inches='tight', pad_inches=0.04)
    plt.close(fig)


def main(opt, d, out):
    R = load(opt, d)
    mainc = names(R.main); allc = names(R.allc); prop = str(R.prop)
    idx = {n: i for i, n in enumerate(allc)}
    traj = R.traj
    ton, toff = 8.0, 14.0

    def tr(n, j):
        return traj[idx[n], j]

    # ---- Fig 3: step responses
    fig, axs = plt.subplots(3, 1, figsize=(6.9, 6.6), sharex=True)
    for j, ax in enumerate(axs):
        ax.axvspan(ton, toff, color=GUST, lw=0)
        ax.plot([0, 24], [1, 1], color='k', lw=0.8, ls=(0, (1, 2)), label='Reference')
        for n in mainc:
            c, ls, lw = STY[n]; t = tr(n, j)
            ax.plot(t.t, t.y, color=c, ls=ls, lw=lw, label=n)
        ax.set_xlim(0, 24); ax.set_ylim(-0.05, 1.42); ax.set_ylabel('Azimuth (deg)')
        ax.set_title(f'({"abc"[j]}) L = {LS[j]} s', loc='left')
        ax.text(ton + 0.15, 1.32, 'wind gust', fontsize=7.5, color='#7a5c00')
        if j == 0: ax.legend(loc='lower right', ncol=len(mainc) + 1)
    axs[-1].set_xlabel('Time (s)')
    save(fig, out, 'fig_step_response')

    # ---- Fig 4: control effort (left: start-up transient at full scale; right: regulation, zoomed)
    fig, axs = plt.subplots(3, 2, figsize=(6.9, 6.0), gridspec_kw=dict(width_ratios=[1, 2.3]))
    for j in range(3):
        ax = axs[j, 0]
        ax.plot([0, 4], [20, 20], color='k', lw=0.8, ls=(0, (1, 2)), label='Actuator limit')
        for n in mainc:
            c, ls, lw = STY[n]; t = tr(n, j)
            ax.plot(t.t, t.u, color=c, ls=ls, lw=lw, label=n)
        ax.set_xlim(0, 4); ax.set_ylim(-2, 22); ax.set_ylabel('u (V)')
        ax.set_title(f'({"abc"[j]}) L = {LS[j]} s, start-up', loc='left')
        if j == 0: ax.legend(loc='upper right', fontsize=7)
        ax = axs[j, 1]
        ax.axvspan(ton, toff, color=GUST, lw=0)
        for n in mainc:
            c, ls, lw = STY[n]; t = tr(n, j)
            ax.plot(t.t, t.u, color=c, ls='-', lw=(1.3 if n == prop else 0.7), label=n, zorder=(5 if n == prop else 3))
        ax.set_xlim(5, 24); ax.set_ylim(-0.6, 2.6)
        ax.set_title('regulation and wind gust (zoomed)', loc='left')
        ax.text(ton + 0.15, 2.3, 'wind gust', fontsize=7.5, color='#7a5c00')
    axs[2, 0].set_xlabel('Time (s)'); axs[2, 1].set_xlabel('Time (s)')
    fig.tight_layout(w_pad=0.8, h_pad=0.6)
    save(fig, out, 'fig_control_effort')

    # ---- Fig 5: pointing error around the gust
    fig, axs = plt.subplots(3, 1, figsize=(6.9, 6.2), sharex=True)
    for j, ax in enumerate(axs):
        ax.axvspan(ton, toff, color=GUST, lw=0)
        ax.axhspan(-0.02, 0.02, color='#e6e6e6', lw=0)
        for n in mainc:
            c, ls, lw = STY[n]; t = tr(n, j)
            ax.plot(t.t, 1 - t.y, color=c, ls=ls, lw=lw, label=n)
        ax.set_xlim(6, 24); ax.set_ylim(-0.33, 0.37); ax.set_ylabel('Pointing error (deg)')
        ax.set_title(f'({"abc"[j]}) L = {LS[j]} s', loc='left')
        if j == 0: ax.legend(loc='upper right', ncol=len(mainc))
    axs[-1].set_xlabel('Time (s)')
    save(fig, out, 'fig_gust_error')

    # ---- Fig 6: indices
    fig, axs = plt.subplots(1, 3, figsize=(6.9, 2.5))
    lab = ['ISE (deg$^2$ s)', 'IAE (deg s)', 'ITAE (deg s$^2$)']
    w = 0.8 / len(mainc)
    for q, ax in enumerate(axs):
        for i, n in enumerate(mainc):
            v = R.nom[idx[n], :, q]
            ax.bar(np.arange(3) + (i - (len(mainc) - 1) / 2) * w, v, w * 0.92, color=STY[n][0], label=n,
                   alpha=0.55 if n != prop else 1.0, edgecolor=STY[n][0], lw=0.8)
        ax.set_xticks(range(3)); ax.set_xticklabels(['0.3', '0.8', '1.5']); ax.set_xlabel('Delay L (s)')
        ax.set_ylabel(lab[q]); ax.grid(axis='x', visible=False)
        if q == 0: ax.legend(loc='upper left')
    fig.tight_layout(w_pad=1.0)
    save(fig, out, 'fig_indices')

    # ---- Fig 7: robustness
    fig, axs = plt.subplots(3, 3, figsize=(6.9, 6.0), sharex=True)
    par = ['Plant gain error', 'Damping error', 'Time-delay error']
    lev = 100 * np.asarray(R.lev)
    for j in range(3):
        for q in range(3):
            ax = axs[j, q]
            V = np.array(R.rob_IAE[:len(mainc), j, q, :], dtype=float); OK = np.array(R.rob_ok[:len(mainc), j, q, :]) > 0.5
            top = 1.22 * np.nanmax(np.where(OK, V, np.nan)); bot = 0.97 * np.nanmin(V)
            for i, n in enumerate(mainc):
                c, ls, lw = STY[n]
                v = np.minimum(V[i], top); ok = OK[i]
                ax.plot(lev, v, color=c, ls=ls, lw=lw, label=n)
                ax.plot(lev[ok], v[ok], 'o', color=c, ms=3)
                ax.plot(lev[~ok], v[~ok], 'x', color=c, ms=6, mew=1.5, clip_on=False, zorder=6)
            ax.set_ylim(bot, top)
            if j == 0: ax.set_title(par[q])
            if q == 0: ax.set_ylabel(f'IAE (deg s), L = {LS[j]} s')
            if j == 2: ax.set_xlabel('Error (%)')
            ax.set_xticks([-30, -15, 0, 15, 30])
    axs[0, 0].legend(loc='upper center', fontsize=7)
    fig.subplots_adjust(left=0.09, right=0.99, top=0.95, bottom=0.08, wspace=0.3, hspace=0.14)
    save(fig, out, 'fig_robustness')

    # ---- Fig 8: Monte Carlo box plots (IAE, ITAE, gust peak)
    fig, axs = plt.subplots(1, 3, figsize=(6.9, 2.7))
    col = [1, 2, 4]; lab = ['IAE (deg s)', 'ITAE (deg s$^2$)', 'Gust peak error (deg)']
    w = 0.8 / len(mainc)
    for q, ax in enumerate(axs):
        for i, n in enumerate(mainc):
            data = []
            for j in range(3):
                x = np.array(R.mc[i, j, :, col[q]], dtype=float); data.append(x[np.isfinite(x)])
            pos = np.arange(3) + (i - (len(mainc) - 1) / 2) * w
            bp = ax.boxplot(data, positions=pos, widths=w * 0.8, patch_artist=True, showfliers=False,
                            medianprops=dict(color='k', lw=1.0), whiskerprops=dict(color=STY[n][0], lw=0.9),
                            capprops=dict(color=STY[n][0], lw=0.9))
            for b in bp['boxes']:
                b.set(facecolor=STY[n][0], alpha=0.55 if n != prop else 1.0, edgecolor=STY[n][0], lw=0.8)
            ax.plot([], [], 's', color=STY[n][0], alpha=0.55 if n != prop else 1.0, label=n)
        ax.set_xticks(range(3)); ax.set_xticklabels(['0.3', '0.8', '1.5']); ax.set_xlabel('Delay L (s)')
        ax.set_ylabel(lab[q]); ax.grid(axis='x', visible=False); ax.set_xlim(-0.6, 2.6)
        if q == 0: ax.legend(loc='upper left')
    fig.tight_layout(w_pad=1.0)
    save(fig, out, 'fig_monte_carlo')

    # ---- Fig 9: step-size study of the proposed controller
    fig, axs = plt.subplots(1, 3, figsize=(6.9, 2.4), sharey=True)
    cs = ['#111111', '#a50f2d', '#0072b2', '#009e73', '#d55e00']; lss = ['-', '-', (0, (4, 2)), (0, (5, 2, 1, 2)), (0, (1, 1))]
    steps = np.atleast_1d(R.steps)
    for j, ax in enumerate(axs):
        for a, s in enumerate(steps):
            t = R.stptraj[j, a]
            ax.plot(t.t, t.y, color=cs[a], ls=lss[a], lw=1.3, label=f'{s:g}°')
        ax.set_xlim(0, 9); ax.set_ylim(0, 1.12); ax.set_xlabel('Time (s)'); ax.set_title(f'L = {LS[j]} s')
    axs[0].set_ylabel('y / step size'); axs[2].legend(loc='lower right', title='step size', title_fontsize=7.5)
    fig.tight_layout(w_pad=0.6)
    save(fig, out, 'fig_step_size')

    # ---- Fig 10: Popov plots
    fig, axs = plt.subplots(1, 3, figsize=(6.9, 2.5))
    for j, ax in enumerate(axs):
        s = R.pop[j]
        H = np.asarray(s.H); w_ = np.asarray(s.w); k = float(s.k); q = float(np.atleast_1d(s.q)[0])
        X = H.real; Y = w_ * H.imag
        ax.plot(X, Y, color='#a50f2d', lw=1.6, label='Popov locus')
        sel = X > -1 / k - 0.6
        xl = (-1 / k - 0.45, max(0.6, 1.08 * X.max())); yl = (max(1.12 * Y[sel].min() - 0.05, -3.0), 0.5)
        xx = np.array([xl[0] - 1, xl[1] + 1])
        if q > 0:
            yy = (xx + 1 / k) / q
            ax.fill_between(xx, yy, yl[1] + 5, color='#e9e9e9', lw=0, zorder=0)
            ax.plot(xx, yy, 'k--', lw=1.0, label=f'Popov line (q = {q:.3g})')
        else:
            ax.axvspan(xl[0] - 1, -1 / k, color='#e9e9e9', lw=0, zorder=0)
            ax.plot([-1 / k, -1 / k], yl, 'k--', lw=1.0, label='Popov line (q = 0)')
        ax.plot([-1 / k], [0], 'ko', ms=4, label='$-1/k$')
        ax.set_xlim(xl); ax.set_ylim(yl); ax.set_xlabel('Re $H(j\\omega)$'); ax.set_title(f'L = {LS[j]} s')
        ax.legend(loc='lower right', fontsize=6.8)
    axs[0].set_ylabel('$\\omega$ Im $H(j\\omega)$')
    fig.tight_layout(w_pad=0.8)
    save(fig, out, 'fig_popov')

    if opt == 2:
        # ---- Fig 11: ablation of the robust design cost
        nm = ['PID', 'S-PID', 'FS-PID', 'S-QPID', 'ES-noweight', 'ES-nofilt', 'ES-noDw', 'ES-QPID']
        fig, ax = plt.subplots(figsize=(6.9, 2.8))
        w = 0.84 / len(nm)
        for i, n in enumerate(nm):
            v = R.J[idx[n], :]
            c = '#e7a0ad' if n == 'S-QPID' else STY[n][0]
            ax.bar(np.arange(3) + (i - (len(nm) - 1) / 2) * w, v, w * 0.92, color=c, label=n,
                   alpha=1.0 if n in ('ES-QPID', 'S-QPID') else 0.6, edgecolor=c, lw=0.8)
        ax.set_xticks(range(3)); ax.set_xticklabels(['L = 0.3 s', 'L = 0.8 s', 'L = 1.5 s'])
        ax.set_ylabel('Robust design cost J'); ax.grid(axis='x', visible=False); ax.legend(loc='upper left', ncol=2)
        save(fig, out, 'fig_ablation')

        # ---- Fig 12: probability weight and its effect, L = 0.8 s
        j = 1; a = tr('ES-QPID', j); b = tr('ES-noweight', j)
        fig, axs = plt.subplots(2, 1, figsize=(6.9, 4.2), sharex=True, gridspec_kw=dict(height_ratios=[1, 1.5]))
        axs[0].axvspan(ton, toff, color=GUST, lw=0)
        axs[0].plot(a.t, a.p, color='#111111', lw=1.5)
        axs[0].set_ylim(-0.05, 1.08); axs[0].set_ylabel('$P_{reg}$')
        axs[0].set_title('(a) probability weight of the regulation state, L = 0.8 s', loc='left')
        axs[1].axvspan(ton, toff, color=GUST, lw=0)
        axs[1].plot([0, 24], [1, 1], color='k', lw=0.8, ls=(0, (1, 2)))
        axs[1].plot(b.t, b.y, color=STY['ES-noweight'][0], ls=(0, (4, 2)), lw=1.3, label='constant weight 1/2 (ES-noweight)')
        axs[1].plot(a.t, a.y, color=STY['ES-QPID'][0], lw=1.9, label='two-qubit weights (ES-QPID)')
        axs[1].set_xlim(0, 24); axs[1].set_ylim(-0.05, 1.3); axs[1].set_xlabel('Time (s)'); axs[1].set_ylabel('Azimuth (deg)')
        axs[1].set_title('(b) response with and without the probability weights', loc='left')
        axs[1].legend(loc='lower right')
        fig.tight_layout(h_pad=0.6)
        save(fig, out, 'fig_weights')


if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2], sys.argv[3])
