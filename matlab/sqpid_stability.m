function s = sqpid_stability(name, L)
% SQPID_STABILITY  Frequency-domain stability checks for a QPID-type controller (nominal model).
%
%   theta_u = (pi/e_max) C e,   C = Kp + cI*Ki*(integrator) + cD*Kd*(filtered derivative)
%       basic law : cI = 0.5, cD = 1 (constant)
%       enhanced law  : cI = P,   cD = 1 - P, with the probability weight P frozen at
%                       P = 1 (regulation), 0.75, 0.5, 0.25, 0 (set-point transient)
%   u = r + Ks e + w,   w = Uq*tanh(alpha*theta_u)   ->  sector (0, k],  k = alpha*Uq
%   With an exact model the Smith predictor gives e = r - G0 u - Fr G0 d, hence
%       theta_u = -H w + (exogenous terms),   H = (pi/e_max) C G0/(1 + Ks G0)
%
%   Fields (one entry per frozen weight in s.P; the first entry is the regulation mode):
%     s.popov, s.q : Popov criterion on the continuous-time equivalent,
%                    max_q min_w Re[(1+jwq)H(jw)] + 1/k   (positive = satisfied) and the multiplier
%     s.circle     : min_w Re H(jw) + 1/k  (circle criterion; applies when cI*Ki = 0)
%     s.GM_dB, s.PM_deg, s.wc : margins of the SAMPLED linearised loop Ct(z) G0(z),
%                    Ct = Ks + k (pi/e_max) C, G0(z) = zero-order-hold equivalent, Ts = 0.01 s
%     s.dL         : delay error guaranteed stable for the sampled linearised loop by the
%                    small-gain test |T(z) Fr(z)| |exp(-jw dL) - 1| < 1,  z = exp(jw Ts)
%   s.w, s.H, s.k : frequency grid, H(jw) of the regulation mode and sector bound (for plotting)

  c = sqpid_ctrl(name, L);  g = c.g;  Ts = 0.01;
  k  = g.alpha*g.Uq;
  if g.pw > 0, P = [1 0.75 0.5 0.25 0]; else, P = 1; end
  s.P = P;  s.k = k;
  qs = [0 logspace(-3, 2, 300)];

  % ---- continuous-time equivalent (Popov / circle criteria)
  w  = logspace(-4, 3, 40000);  jw = 1i*w;
  G0 = 6.63./(jw.^2 + 14*jw + 6.63);
  s.w = w;

  % ---- sampled loop (margins and delay-error bound)
  wd = logspace(-3, log10(pi/Ts), 30000);  z = exp(1i*wd*Ts);
  [A, B, C] = sqpid_zoh(6.63, [1 14 6.63], Ts);
  Gd = zeros(size(z));
  for n = 1:numel(z), Gd(n) = C*((z(n)*eye(2) - A)\B); end
  a1 = (2 - g.Nf*Ts)/(2 + g.Nf*Ts);  b0 = g.Nf*2/(2 + g.Nf*Ts);
  Dz = b0*(1 - 1./z)./(1 - a1./z);                 % Tustin filtered derivative
  Iz = Ts*z./(z - 1);                              % integrator as implemented
  if g.Tf > 0
    den0 = Ts + 2*g.Tf;
    Frd = ((Ts + 2*g.beta)/den0 + (Ts - 2*g.beta)/den0./z)./(1 + (Ts - 2*g.Tf)/den0./z);
  else
    Frd = ones(size(z));
  end

  for n = 1:numel(P)
    if g.pw > 0
      cI = P(n);
      if g.dw > 0, cD = 1 - P(n); else, cD = 1; end
    else
      cI = 0.5;  cD = 1;
    end
    % Popov / circle on the continuous-time equivalent
    Cc = g.Kp + cI*g.Ki./jw + cD*g.Kd*g.Nf*jw./(jw + g.Nf);
    H  = (pi/g.emax)*Cc.*G0./(1 + g.Ks*G0);
    if n == 1, s.H = H; end
    s.circle(n) = min(real(H)) + 1/k;
    best = -Inf;  bq = 0;
    for q = qs
      v = min(real((1 + jw*q).*H)) + 1/k;
      if v > best, best = v; bq = q; end
    end
    s.popov(n) = best;  s.q(n) = bq;

    % margins and delay-error bound on the sampled loop
    Cd = g.Kp + cI*g.Ki*Iz + cD*g.Kd*Dz;
    Lo = (g.Ks + k*(pi/g.emax)*Cd).*Gd;
    mag = abs(Lo);  ph = unwrap(angle(Lo))*180/pi;
    ic = find(mag(1:end-1) >= 1 & mag(2:end) < 1, 1, 'last');
    if isempty(ic), s.wc(n) = NaN; s.PM_deg(n) = Inf; else, s.wc(n) = wd(ic); s.PM_deg(n) = 180 + ph(ic); end
    ig = find(ph(1:end-1) > -180 & ph(2:end) <= -180, 1);
    if isempty(ig), s.GM_dB(n) = Inf; else, s.GM_dB(n) = -20*log10(mag(ig)); end
    TF = abs(Lo./(1 + Lo).*Frd);
    lo = 0;  hi = 10;
    for it = 1:60
      dL = 0.5*(lo + hi);
      if max(TF.*2.*abs(sin(wd*dL/2))) < 1, lo = dL; else, hi = dL; end
    end
    s.dL(n) = lo;
  end
end
