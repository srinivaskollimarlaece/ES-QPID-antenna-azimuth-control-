function o = sqpid_sim(c)
% SQPID_SIM  Closed-loop simulation of the antenna azimuth loop (corrected test bed).
%
%   Plant      : y = G0(s) [ u(t-L) + d(t) ],  G0 = kg*6.63/(s^2 + 14(1+delta)s + 6.63)
%                (true sample delay, no Pade approximation; d = torque-equivalent wind gust)
%   Feedback   : ym = y + n      (the measurement noise IS fed back to the controller)
%   Actuator   : |u| <= Usat     (same limit for every controller)
%
%   c.law : 0 linear PID (filtered derivative, clamp anti-windup)
%           1 quantum-inspired PID (QPID)
%   c.sp  : 0 no predictor
%           1 Smith predictor with mismatch filter Fr(s) = (1 + beta s)/(1 + Tf s)
%             (Tf = 0 gives the conventional Smith predictor, Fr = 1)
%   c.L, c.Lm     : true plant delay and model delay (s)
%   c.kg, c.delta : plant gain error factor and damping error
%   c.r, c.d, c.nz: reference (deg), disturbance (V), measurement noise (deg) - row vectors
%   c.g : gain structure
%         Kp Ki Kd       PID / QPID gains
%         Ks             output-feedback weight of the QPID law (Kr = Ks + 1)
%         Tf beta        mismatch filter constants
%         Nf             derivative filter corner (rad/s)
%         alpha Uq       tanh slope and bound of the QPID correction term (V)
%         emax Ilim      fine-pointing range (deg) and integral limit
%         pw             0 constant weight 0.5 (basic law)
%                        1 two-qubit probability weights (enhanced law)
%         tm             time constant of the set-point transition model (s), used when pw = 1
%         dw             1 weight the derivative term by (1 - Preg), used when pw = 1
%
%   Returns o.t, o.y, o.u, o.preg (probability weight of the regulation state).

  Ts = c.Ts;  r = c.r;  d = c.d;  nz = c.nz;  N = numel(r);
  [A, B, C]    = sqpid_zoh(c.kg*6.63, [1 14*(1+c.delta) 6.63], Ts);   % true plant
  [Am, Bm, Cm] = sqpid_zoh(6.63, [1 14 6.63], Ts);                    % model in the predictor
  a11=A(1,1); a12=A(1,2); a21=A(2,1); a22=A(2,2); b1=B(1); b2=B(2); c1=C(1); c2=C(2);
  m11=Am(1,1); m12=Am(1,2); m21=Am(2,1); m22=Am(2,2); n1=Bm(1); n2=Bm(2); k1=Cm(1); k2=Cm(2);
  nL = round(c.L/Ts);  nLm = round(c.Lm/Ts);
  ubuf = zeros(1, N+nL+2);  ymbuf = zeros(1, N+nLm+2);

  g = c.g;  law = c.law;  sp = c.sp;  Usat = c.Usat;
  Kp=g.Kp; Ki=g.Ki; Kd=g.Kd; Ks=g.Ks; Tf=g.Tf; beta=g.beta; Nf=g.Nf; al=g.alpha; Uq=g.Uq;
  emax=g.emax; Ilim=g.Ilim; pw=g.pw; tm=g.tm; dw=g.dw;

  a1 = (2 - Nf*Ts)/(2 + Nf*Ts);  b0 = Nf*2/(2 + Nf*Ts);       % Tustin form of Nf*s/(s+Nf)
  if Tf > 0                                                   % Tustin form of Fr(s)
    den0 = Ts + 2*Tf;  fa = -(Ts - 2*Tf)/den0;
    fl0 = (Ts + 2*beta)/den0;  fl1 = (Ts - 2*beta)/den0;
  else
    fa = 0; fl0 = 1; fl1 = 0;
  end

  x1=0; x2=0; z1=0; z2=0; I=0; D=0; xp=0; mF=0; mprev=0;
  rprev = 0;  rm = 0;  dref = 1;  p0r = 1;  p0 = 0.5;
  y = zeros(1,N);  u = zeros(1,N);  pv = ones(1,N);

  for k = 1:N
    yt = c1*x1 + c2*x2;
    y(k) = yt;
    if abs(yt) > 1e6, y(k:end) = NaN; break; end
    ym = yt + nz(k);
    rk = r(k);

    % ---- set-point transition qubit (driven by the reference only) ----
    if abs(rk - rprev) > 0.02
      dref = max(abs(rk - rm), 0.05);        % size of the commanded step
    end
    rprev = rk;
    if pw > 0
      rat = min(abs(rk - rm)/dref, 1);
      p0r = cos(0.5*pi*rat)^2;               % probability of the "set-point reached" state
      rm  = rm + (Ts/tm)*(rk - rm);
    end

    % ---- predictor ----
    if sp == 0
      yfb = ym;
    else
      yh = k1*z1 + k2*z2;                    % delay-free model output
      m  = ym - ymbuf(k);                    % plant minus delayed model: mismatch + disturbance
      mF = fa*mF + fl0*m + fl1*mprev;        % filtered by Fr
      mprev = m;
      yfb = yh + mF;
      ymbuf(k + nLm) = yh;
    end
    e = rk - yfb;

    if law == 0
      % ---------- linear PID ----------
      D = a1*D + b0*(e - xp);  xp = e;
      v = Kp*e + Ki*I + Kd*D;
      uk = max(min(v, Usat), -Usat);
      if uk == v || e*v < 0, I = I + e*Ts; end
    else
      % ---------- quantum-inspired PID ----------
      th = pi*e/emax;  th = max(min(th, pi), -pi);          % pointing-error qubit angle
      if pw > 0
        p0 = p0r*cos(th/2)^2;                               % P(regulation state)
      else
        p0 = 0.5;                                           % constant weight (basic law)
      end
      I = I + Ki*Ts*p0*th;  I = max(min(I, Ilim), -Ilim);
      D = a1*D + b0*(th - xp);  xp = th;
      if pw > 0 && dw > 0
        dterm = Kd*D*(1 - p0);
      else
        dterm = Kd*D;
      end
      thu = Kp*th + I + dterm;
      v   = rk + Ks*e + Uq*tanh(al*thu);                    % = Kr*r - Ks*yfb + Uq*tanh(.)
      uk  = max(min(v, Usat), -Usat);
    end

    if sp == 1
      z1n = m11*z1 + m12*z2 + n1*uk;
      z2  = m21*z1 + m22*z2 + n2*uk;
      z1  = z1n;
    end
    u(k) = uk;  pv(k) = p0;
    ubuf(k + nL) = uk;
    ud  = ubuf(k) + d(k);
    x1n = a11*x1 + a12*x2 + b1*ud;
    x2  = a21*x1 + a22*x2 + b2*ud;
    x1  = x1n;
  end
  o.t = (0:N-1)*Ts;  o.y = y;  o.u = u;  o.preg = pv;
end
