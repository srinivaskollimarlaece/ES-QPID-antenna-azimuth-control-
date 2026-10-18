% ES_QPID_OPTION2_COMPLETE  Complete MATLAB code for the article
%   "Smith Predictor-Based Quantum-Inspired PID Controller with Hyperbolic
%   Tangent Nonlinearity for Robust Azimuth Position Control of Earth Station Antennas under
%   Variable Time Delay" - enhanced Smith-QPID (ES-QPID): mismatch filter + two-qubit weights.
%
%   HOW TO RUN (MATLAB R2016b or later, no toolbox needed)
%     1. Keep this file and sqpid_random.mat in the same folder (the matlab folder of the repository).
%     2. Make that folder the current folder and type   ES_QPID_Option2_complete
%   Output: all tables in the Command Window and in results_option2.txt, the data in
%   results_option2.mat, and the figures as PNG files in the folder figures_option2.
%   sqpid_random.mat holds the exact random sequences of the paper, so the tables match the
%   manuscript to the last digit. Without it MATLAB's own generator is used (small differences).
%
%   One run, after the full script has been run once (functions are local to this file, so use
%   the separate function files in this folder for interactive single runs):
%     [m, o] = sqpid_run('ES-QPID', 0.8, 'Lerr', 1.2);  plot(o.t, o.y)

clc; clear; close all;
tic;
R = sqpid_experiments(2, 200);      % 2 = Option 2, 200 Monte Carlo runs per delay
sqpid_figures(2);
fprintf('Finished in %.0f s\n', toc);


%% =====================================================================================
function R = sqpid_experiments(opt, nmc)
% SQPID_EXPERIMENTS  Run every experiment of the article.
%
%   R = sqpid_experiments(opt, nmc)
%     opt = 1 : basic Smith-QPID (S-QPID) as main controller
%     opt = 2 : enhanced Smith-QPID (ES-QPID), as in the article
%     nmc     : number of Monte Carlo runs (default 200)
%
%   Experiments
%     E1 nominal comparison
%     E2 parametric robustness sweeps
%     E3 delay margins
%     E4 step-size / e_max study
%     E5 actuator-limit study
%     E6 Monte Carlo trials
%     E7 absolute-stability checks
%
%   Results are printed, returned in R and saved to results_option<opt>.mat

  if nargin < 2, nmc = 200; end
  Ls = [0.3 0.8 1.5];
  if opt == 1
    prop = 'S-QPID';
    main = {'PID', 'S-PID', 'S-QPID'};
    abl  = {'QPID'};
  else
    prop = 'ES-QPID';
    main = {'PID', 'S-PID', 'FS-PID', 'ES-QPID'};
    abl  = {'S-QPID', 'ES-noweight', 'ES-nofilt', 'ES-noDw'};
  end
  allc = [main abl];
  R.opt = opt;  R.prop = prop;  R.main = main;  R.abl = abl;  R.allc = allc;  R.Ls = Ls;
  fid = fopen(sprintf('results_option%d.txt', opt), 'w');
  out = @(varargin) dual_print(fid, varargin{:});

  out('==========================================================================================\n');
  out(' OPTION %d  -  proposed controller: %s\n', opt, prop);
  out(' Test bed: true delay, real plant feedback, torque gust 0.5 V (20%% turbulence) on [8,14) s,\n');
  out('           measurement noise 0.005 deg fed back, actuator limit +/-20 V for ALL controllers\n');
  out('==========================================================================================\n');

  %% ------------------------------------------------------------------ E1 nominal
  out('\nE1. NOMINAL PERFORMANCE (1 deg step, nominal plant)\n');
  out('%-12s %4s | %6s %6s %7s | %5s %5s %5s | %6s %6s %6s | %6s %5s\n', 'Controller', 'L', ...
      'ISE', 'IAE', 'ITAE', 'tr', 'ts', 'OS%', 'gPeak', 'gIAE', 'ess', 'chat', 'umax');
  fld = {'ISE','IAE','ITAE','rise','Ts','OS','gpk','gIAE','ess','du','umax'};
  R.nom = zeros(numel(allc), numel(Ls), numel(fld));
  R.traj = cell(numel(allc), numel(Ls));
  for j = 1:numel(Ls)
    for i = 1:numel(allc)
      [m, o] = sqpid_run(allc{i}, Ls(j));
      for f = 1:numel(fld), R.nom(i,j,f) = m.(fld{f}); end
      R.traj{i,j} = struct('t', o.t, 'y', o.y, 'u', o.u, 'p', o.preg, 'd', o.d);
      out('%-12s %4.1f | %6.3f %6.3f %7.3f | %5.2f %5.2f %5.2f | %6.3f %6.3f %6.3f | %6.3f %5.1f\n', ...
          allc{i}, Ls(j), m.ISE, m.IAE, m.ITAE, m.rise, m.Ts, m.OS, m.gpk, m.gIAE, m.ess, m.du, m.umax);
    end
    out('\n');
  end
  R.fld = fld;
  ip = find(strcmp(allc, prop));
  out('Improvement of %s (positive = better), nominal:\n', prop);
  out('%-10s %4s | %6s %6s %6s %6s %6s\n', 'versus', 'L', 'ISE%', 'IAE%', 'ITAE%', 'ts%', 'gIAE%');
  cmpc = main(~strcmp(main, prop));
  R.impr = zeros(numel(cmpc), numel(Ls), 5);
  for i = 1:numel(cmpc)
    ib = find(strcmp(allc, cmpc{i}));
    for j = 1:numel(Ls)
      k = [1 2 3 5 8];
      v = 100*(squeeze(R.nom(ib,j,k)) - squeeze(R.nom(ip,j,k)))./squeeze(R.nom(ib,j,k));
      R.impr(i,j,:) = v;
      out('%-10s %4.1f | %6.1f %6.1f %6.1f %6.1f %6.1f\n', cmpc{i}, Ls(j), v);
    end
  end
  R.cmpc = cmpc;

  %% ------------------------------------------------------------------ robust design cost
  out('\nROBUST DESIGN COST J (7 scenarios, same for every controller)\n');
  R.J = zeros(numel(allc), numel(Ls));
  for i = 1:numel(allc)
    for j = 1:numel(Ls)
      R.J(i,j) = sqpid_cost(allc{i}, Ls(j), []);
    end
    out('%-12s  J = %6.3f  %6.3f  %6.3f\n', allc{i}, R.J(i,:));
  end

  %% ------------------------------------------------------------------ E2 parametric sweeps
  out('\nE2. PARAMETRIC ROBUSTNESS (IAE; "x" = not settled at the end of the run)\n');
  lev = [-0.30 -0.20 -0.10 -0.05 0 0.05 0.10 0.20 0.30];
  par = {'kg', 'delta', 'Lerr'};  pname = {'plant gain', 'damping', 'time delay'};
  R.lev = lev;  R.par = pname;
  R.rob_IAE = zeros(numel(main), numel(Ls), 3, numel(lev));
  R.rob_OS  = R.rob_IAE;  R.rob_ok = R.rob_IAE;
  for q = 1:3
    out('\n  Varied parameter: %s\n', pname{q});
    out('  %-10s %4s |', 'Controller', 'L');  out(' %+5.0f%% ', 100*lev);  out('\n');
    for j = 1:numel(Ls)
      for i = 1:numel(main)
        out('  %-10s %4.1f |', main{i}, Ls(j));
        for a = 1:numel(lev)
          if q == 2, val = lev(a); else, val = 1 + lev(a); end
          m = sqpid_run(main{i}, Ls(j), par{q}, val);
          R.rob_IAE(i,j,q,a) = m.IAE;  R.rob_OS(i,j,q,a) = m.OS;  R.rob_ok(i,j,q,a) = m.ok;
          if m.ok, out(' %6.3f ', m.IAE); else, out(' %6.3fx', min(m.IAE, 99.999)); end
        end
        out('\n');
      end
    end
  end

  %% ------------------------------------------------------------------ E3 delay margins
  out('\nE3. DELAY MARGIN (simulation, no gust/noise): range of true delay for which the loop settles\n');
  out('  %-10s %4s | %9s %9s | %9s %9s\n', 'Controller', 'L', 'Lmin (s)', 'Lmax (s)', 'low %', 'high %');
  R.dm = zeros(numel(main), numel(Ls), 2);
  for j = 1:numel(Ls)
    for i = 1:numel(main)
      lo = 1;  hi = 1;
      for f = 1:-0.05:0.05
        m = sqpid_run(main{i}, Ls(j), 'Lerr', f, 'gamp', 0, 'noise', 0, 'T', 40);
        if m.ok, lo = f; else, break; end
      end
      for f = 1:0.05:3
        m = sqpid_run(main{i}, Ls(j), 'Lerr', f, 'gamp', 0, 'noise', 0, 'T', 40);
        if m.ok, hi = f; else, break; end
      end
      R.dm(i,j,:) = [lo hi];
      out('  %-10s %4.1f | %9.3f %9.3f | %+8.0f%% %+8.0f%%\n', main{i}, Ls(j), lo*Ls(j), hi*Ls(j), ...
          100*(lo-1), 100*(hi-1));
    end
  end

  %% ------------------------------------------------------------------ E4 step size / e_max
  out('\nE4. STEP-SIZE STUDY (no gust): overshoot, settling time, peak |u|\n');
  steps = [0.5 1 2 5 10];
  R.steps = steps;
  R.stp = zeros(numel(main), numel(Ls), numel(steps), 4);
  R.stptraj = cell(numel(Ls), numel(steps));
  out('  %-10s %4s |', 'Controller', 'L');
  for s = steps, out('  R=%4.1f: OS%%  ts   umax |', s); end
  out('\n');
  for j = 1:numel(Ls)
    for i = 1:numel(main)
      out('  %-10s %4.1f |', main{i}, Ls(j));
      for a = 1:numel(steps)
        [m, o] = sqpid_run(main{i}, Ls(j), 'R', steps(a), 'gamp', 0, 't_on', 20, 't_off', 21, 'T', 24);
        R.stp(i,j,a,:) = [m.OS m.Ts m.umax m.IAE/steps(a)];
        if strcmp(main{i}, prop), R.stptraj{j,a} = struct('t', o.t, 'y', o.y/steps(a), 'u', o.u); end
        out('       %5.2f %5.2f %5.1f |', m.OS, m.Ts, m.umax);
      end
      out('\n');
    end
  end
  out('\n  %s with the fine-pointing range e_max set equal to the step size R:\n', prop);
  R.stp_emax = zeros(numel(Ls), numel(steps), 4);
  for j = 1:numel(Ls)
    out('  %-10s %4.1f |', prop, Ls(j));
    for a = 1:numel(steps)
      m = sqpid_run(prop, Ls(j), 'R', steps(a), 'gamp', 0, 't_on', 20, 't_off', 21, 'T', 24, ...
                    'emax', max(steps(a), 1e-3));
      R.stp_emax(j,a,:) = [m.OS m.Ts m.umax m.IAE/steps(a)];
      out('       %5.2f %5.2f %5.1f |', m.OS, m.Ts, m.umax);
    end
    out('\n');
  end

  %% ------------------------------------------------------------------ E5 actuator limit
  out('\nE5. ACTUATOR-LIMIT STUDY (same limit for both controllers): IAE / settling time / OS%%\n');
  lims = [6 12 20 30];
  R.lims = lims;
  cl = {'PID', prop};
  R.lim = zeros(2, numel(Ls), numel(lims), 4);
  out('  %-10s %4s |', 'Controller', 'L');
  for s = lims, out('   Usat=%2d V: IAE    ts   OS%%  |', s); end
  out('\n');
  for j = 1:numel(Ls)
    for i = 1:2
      out('  %-10s %4.1f |', cl{i}, Ls(j));
      for a = 1:numel(lims)
        m = sqpid_run(cl{i}, Ls(j), 'Usat', lims(a));
        R.lim(i,j,a,:) = [m.IAE m.Ts m.OS m.ITAE];
        out('            %6.3f %5.2f %5.2f |', m.IAE, m.Ts, m.OS);
      end
      out('\n');
    end
  end

  %% ------------------------------------------------------------------ E6 Monte Carlo
  out('\nE6. MONTE CARLO (%d runs per delay): gain +/-20%%, damping +/-20%%, delay +/-15%%,\n', nmc);
  out('    gust mean 0.25-1.0 V, turbulence 10-40%%, gust onset 5-10 s, duration 3-8 s, new noise each run\n');
  P = sqpid_noise(-1, nmc);                         % uniform draws of the Monte Carlo conditions
  R.mc = zeros(numel(main), numel(Ls), nmc, 6);     % ISE IAE ITAE OS gpk ok
  for j = 1:numel(Ls)
    for n = 1:nmc
      kg = 0.8 + 0.4*P(n,1);  dl = -0.2 + 0.4*P(n,2);  le = 0.85 + 0.30*P(n,3);
      ga = 0.25 + 0.75*P(n,4);  gs = 0.1 + 0.3*P(n,5);  t1 = 5 + 5*P(n,6);  t2 = t1 + 3 + 5*P(n,7);
      for i = 1:numel(main)
        m = sqpid_run(main{i}, Ls(j), 'kg', kg, 'delta', dl, 'Lerr', le, 'gamp', ga, 'gsig', gs, ...
                      't_on', t1, 't_off', t2, 'seed', n, 'T', 28);
        R.mc(i,j,n,:) = [m.ISE m.IAE m.ITAE m.OS m.gpk m.ok];
      end
    end
    out('\n  L = %.1f s\n', Ls(j));
    out('  %-10s | %-15s %-15s %-17s %-13s %-15s %8s\n', 'Controller', 'ISE mean+-sd', ...
        'IAE mean+-sd', 'ITAE mean+-sd', 'OS%% mean+-sd', 'gPeak mean+-sd', 'settled');
    for i = 1:numel(main)
      X = squeeze(R.mc(i,j,:,:));
      fin = isfinite(X(:,2));
      out('  %-10s | %6.3f+-%5.3f   %6.3f+-%5.3f   %7.3f+-%6.3f   %5.2f+-%5.2f   %6.3f+-%5.3f   %5.1f%%\n', ...
          main{i}, mean(X(fin,1)), std(X(fin,1)), mean(X(fin,2)), std(X(fin,2)), mean(X(fin,3)), ...
          std(X(fin,3)), mean(X(fin,4)), std(X(fin,4)), mean(X(fin,5)), std(X(fin,5)), 100*mean(X(:,6)));
    end
    ipm = find(strcmp(main, prop));
    for i = 1:numel(main)
      if i == ipm, continue; end
      w = mean(squeeze(R.mc(ipm,j,:,2)) < squeeze(R.mc(i,j,:,2)));
      out('  %s has the lower IAE than %-7s in %5.1f%% of the runs\n', prop, main{i}, 100*w);
    end
  end
  R.mcP = P;

  %% ------------------------------------------------------------------ E7 stability
  [R, txt] = sqpid_e7(R);
  out('%s', txt);

  fclose(fid);
  save(sprintf('results_option%d.mat', opt), 'R', '-v7');
end

function dual_print(fid, varargin)
  fprintf(varargin{:});  fprintf(fid, varargin{:});
  drawnow;
end

%% =====================================================================================
function [m, o, c] = sqpid_run(name, L, varargin)
% SQPID_RUN  Run one closed-loop experiment and return its performance indices.
%
%   [m, o] = sqpid_run(name, L, 'option', value, ...)
%
%   name : controller (see sqpid_ctrl.m),  L : nominal (model) delay in s
%   Options (default):
%     'R'     (1)    step amplitude, deg (step applied at t = 0)
%     'T'     (24)   run time, s
%     'Lerr'  (1)    true delay = Lerr * L        (delay mismatch)
%     'kg'    (1)    plant gain error factor
%     'delta' (0)    damping error
%     'gamp'  (0.5)  mean gust torque, in equivalent volts (0.5 V -> 0.5 deg static offset)
%     'gsig'  (0.2)  turbulence intensity (std of the coloured component / mean)
%     't_on'  (8), 't_off' (14)   gust window, s
%     'noise' (0.005) measurement noise std, deg
%     'seed'  (0)    0 = nominal fixed sequences; >0 = another reproducible realisation
%     'Usat'  (20)   actuator limit, V
%     'emax'  ([])   override of the fine-pointing range e_max, deg
%     'gains' ([])   struct of gain fields to override (used by tune_gains.m)
%
%   m : indices from sqpid_metrics,  o : time histories,  c : configuration used

  p = struct('R',1,'T',24,'Lerr',1,'kg',1,'delta',0,'gamp',0.5,'gsig',0.2,'t_on',8,'t_off',14, ...
             'noise',0.005,'seed',0,'Usat',20,'emax',[],'gains',[]);
  for k = 1:2:numel(varargin), p.(varargin{k}) = varargin{k+1}; end

  Ts = 0.01;  N = round(p.T/Ts);  t = (0:N-1)*Ts;
  [nzu, wn] = sqpid_noise(p.seed, N);

  % torque-equivalent wind gust: mean step + first-order coloured turbulence (tc = 0.3 s)
  a = exp(-Ts/0.3);  xi = filter(sqrt(1 - a^2), [1 -a], [0 wn(2:N)]);
  w = (t >= p.t_on & t < p.t_off);
  d = zeros(1,N);  d(w) = p.gamp*(1 + p.gsig*xi(w));

  c = sqpid_ctrl(name, L);
  if ~isempty(p.emax), c.g.emax = p.emax; end
  if ~isempty(p.gains)
    f = fieldnames(p.gains);
    for k = 1:numel(f), c.g.(f{k}) = p.gains.(f{k}); end
    if c.g.Tf > 0, c.g.beta = (1 - (1 - 0.5*c.g.Tf)*exp(-0.5*L))/0.5; end
  end
  c.Ts = Ts;  c.L = L*p.Lerr;  c.Lm = L;  c.kg = p.kg;  c.delta = p.delta;  c.Usat = p.Usat;
  c.r = p.R*ones(1,N);  c.d = d;  c.nz = p.noise*nzu;

  o = sqpid_sim(c);
  o.r = c.r;  o.d = d;
  m = sqpid_metrics(t, c.r, o.y, o.u, p.R, p.t_on, p.t_off);
end

%% =====================================================================================
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

%% =====================================================================================
function [Ad, Bd, Cd] = sqpid_zoh(K, den, Ts)
% SQPID_ZOH  Zero-order-hold discretisation of G(s) = K/(s^2 + den(2) s + den(3)).
%   Base MATLAB only (expm); no Control System Toolbox needed.
%   Equivalent to [Ad,Bd,Cd] = ssdata(c2d(ss(tf(K,den)),Ts)) up to a change of state coordinates.
  A = [-den(2) -den(3); 1 0];  B = [1; 0];  C = [0 K];
  M = expm([A B; zeros(1,3)]*Ts);
  Ad = M(1:2,1:2);  Bd = M(1:2,3);  Cd = C;
end

%% =====================================================================================
function [nz, wn] = sqpid_noise(seed, N)
% SQPID_NOISE  Reproducible unit-variance sequences: measurement noise nz, turbulence driver wn.
%   seed = 0  : nominal sequences used for every table and figure of the nominal study
%   seed > 0  : independent realisation number 'seed' (Monte Carlo run)
%   seed = -1 : returns the N-by-7 matrix of uniform draws that defines the Monte Carlo conditions
%
%   The sequences are read from sqpid_random.mat. They are exactly the sequences used for the
%   manuscript, so MATLAB reproduces the published tables to the last digit. If the file is
%   missing, MATLAB's own generator is used instead (results then differ in the last digits).
  persistent D
  if isempty(D)
    if exist('sqpid_random.mat', 'file')
      D = load('sqpid_random.mat');
    else
      warning('sqpid_noise:nofile', ...
        'sqpid_random.mat not found: using rng() sequences; numbers will differ slightly from the paper.');
      s = rng;  rng(1);  D.NZ0 = randn(1, 4000);  D.WN0 = randn(1, 4000);
      D.MCNZ = zeros(200, 2800);  D.MCWN = zeros(200, 2800);
      for n = 1:200, rng(1000 + n);  D.MCNZ(n,:) = randn(1, 2800);  D.MCWN(n,:) = randn(1, 2800); end
      rng(2026);  D.P = rand(200, 7);  rng(s);
    end
  end
  if seed == -1
    if N > size(D.P, 1), error('sqpid_noise: at most %d Monte Carlo runs are stored', size(D.P, 1)); end
    nz = D.P(1:N, :);  wn = [];
  elseif seed == 0
    nz = D.NZ0(1:N);  wn = D.WN0(1:N);
  else
    if seed > size(D.MCNZ, 1) || N > size(D.MCNZ, 2)
      error('sqpid_noise: stored Monte Carlo sequences cover 200 runs of at most 28 s');
    end
    nz = D.MCNZ(seed, 1:N);  wn = D.MCWN(seed, 1:N);
  end
end

%% =====================================================================================
function m = sqpid_metrics(t, r, y, u, R, t_on, t_off)
% SQPID_METRICS  Performance indices of one run.
%   Tracking indices (rise, settling, overshoot) use the window before the gust.
  Ts = t(2) - t(1);
  if any(~isfinite(y))
    m = struct('ISE',Inf,'IAE',Inf,'ITAE',Inf,'sITAE',Inf,'rise',NaN,'Ts',Inf,'OS',Inf,'umax',Inf, ...
               'gpk',Inf,'gIAE',Inf,'ess',Inf,'du',Inf,'ok',false,'div',true);
    return;
  end
  e = r - y;
  m.ISE  = sum(e.^2)*Ts;
  m.IAE  = sum(abs(e))*Ts;
  m.ITAE = sum(t.*abs(e))*Ts;
  tau = t;  tau(t >= t_on) = t(t >= t_on) - t_on;  tau(t >= t_off) = t(t >= t_off) - t_off;
  m.sITAE = sum(tau.*abs(e))*Ts;                          % ITAE with the clock restarted at each event
  w   = (t <= t_on - 0.5);  ye = y(w)/R;  te = t(w);
  i10 = find(ye >= 0.1, 1);  i90 = find(ye >= 0.9, 1);
  if isempty(i10) || isempty(i90), m.rise = NaN; else, m.rise = te(i90) - te(i10); end
  idx = find(abs(ye - 1) > 0.02);
  if isempty(idx), m.Ts = 0; else, m.Ts = te(max(idx)); end
  m.OS   = max((max(ye) - 1)*100, 0);
  m.umax = max(abs(u));
  gw = (t >= t_on);
  m.gpk  = max(abs(e(gw)));                               % peak pointing error caused by the gust
  m.gIAE = sum(abs(e(gw)))*Ts;                            % IAE from gust onset to end of run
  m.ess  = abs(mean(e(t >= t_off - 1 & t < t_off)));      % residual offset at the end of the gust
  q = (t >= t_on - 2.5 & t < t_on);
  m.du   = sqrt(mean(diff(u(q)).^2));                     % control chatter (V per sample, rms)
  tail   = max(abs(e(t > t(end) - 2)));
  m.ok   = tail < 0.05*R + 0.02;                          % settled at the end of the run
  m.div  = tail > 0.5*R;
end

%% =====================================================================================
function c = sqpid_ctrl(name, L)
% SQPID_CTRL  Controller definition and tuned gains for delay L (0.3, 0.8 or 1.5 s).
%
%   All gains were obtained with the SAME procedure for every controller: minimisation of the
%   robust design cost of sqpid_cost.m (7 scenarios: nominal, delay -10%/+10%, gain -10%/+10%,
%   double gust, early gust) by differential evolution followed by a Nelder-Mead polish.
%   tune_gains.m repeats the local polish in MATLAB.
%
%   Rows of each table: L = 0.3, 0.8, 1.5 s.
%
%   Controllers:
%     PID          classical PID, no predictor
%     S-PID        conventional Smith predictor + PID
%     FS-PID       filtered Smith predictor + PID
%     QPID         basic QPID (constant weight 0.5), no predictor
%     S-QPID       basic Smith-QPID, constant weight 0.5
%     ES-QPID      enhanced Smith-QPID: mismatch filter + two-qubit probability weights (Option 2)
%     ES-noweight  ES-QPID without probability weights (ablation)
%     ES-nofilt    ES-QPID without mismatch filter (ablation)
%     ES-noDw      ES-QPID without derivative weighting (ablation)

  Ltab = [0.3 0.8 1.5];
  i = find(abs(Ltab - L) < 1e-9);
  if isempty(i), error('sqpid_ctrl: gains are tabulated for L = 0.3, 0.8 and 1.5 s only'); end

  g = struct('Kp',0,'Ki',0,'Kd',0,'Ks',0,'Tf',0,'beta',0,'Nf',10,'alpha',0.10,'Uq',20, ...
             'emax',1,'Ilim',6,'pw',0,'tm',1,'dw',0);
  switch name
    case 'PID'
      %               Kp         Ki         Kd
      G = [    3.5565     1.8266    0.53353 ; ...
               1.7825    0.78956    0.63293 ; ...
               1.1748     0.4673    0.61174 ];
      law = 0;  sp = 0;
      g.Kp = G(i,1);  g.Ki = G(i,2);  g.Kd = G(i,3);
      g.Nf = 20;
    case 'S-PID'
      %               Kp         Ki         Kd
      G = [    11.766     6.8994    0.98034 ; ...
               6.1312     3.2349     0.7571 ; ...
               5.2242       2.45     0.4546 ];
      law = 0;  sp = 1;
      g.Kp = G(i,1);  g.Ki = G(i,2);  g.Kd = G(i,3);
      g.Nf = 20;
    case 'FS-PID'
      %               Kp         Ki         Kd         Tf
      G = [    10.366     5.8395    0.84786    0.51313 ; ...
               4.5992     2.4711    0.36417    0.39955 ; ...
               2.8109     1.4755    0.28691     0.4876 ];
      law = 0;  sp = 1;
      g.Kp = G(i,1);  g.Ki = G(i,2);  g.Kd = G(i,3);  g.Tf = G(i,4);
      g.Nf = 20;
    case 'QPID'
      %               Kp         Ki         Kd         Ks
      G = [   0.44759      0.001   0.080679   0.010006 ; ...
              0.19351      0.001   0.080047   0.010008 ; ...
             0.096029      0.001   0.067459   0.010007 ];
      law = 1;  sp = 0;
      g.Kp = G(i,1);  g.Ki = G(i,2);  g.Kd = G(i,3);  g.Ks = G(i,4);
    case 'S-QPID'
      %               Kp         Ki         Kd         Ks
      G = [    2.1447      0.001    0.38762   0.010002 ; ...
              0.96716      0.001    0.16014   0.010015 ; ...
              0.50454      0.001    0.15695   0.010005 ];
      law = 1;  sp = 1;
      g.Kp = G(i,1);  g.Ki = G(i,2);  g.Kd = G(i,3);  g.Ks = G(i,4);
    case 'ES-QPID'
      %               Kp         Ki         Kd         Ks         Tf         tm
      G = [  0.011154      6.691    0.14883     8.8356    0.66495     2.9157 ; ...
              0.51373    0.68281    0.04893     0.0365    0.29327          5 ; ...
              0.19834    0.54655    0.04906    0.88553    0.50058          5 ];
      law = 1;  sp = 1;
      g.Kp = G(i,1);  g.Ki = G(i,2);  g.Kd = G(i,3);  g.Ks = G(i,4);  g.Tf = G(i,5);  g.tm = G(i,6);
      g.pw = 1;  g.dw = 1;
    case 'ES-noweight'
      %               Kp         Ki         Kd         Ks         Tf
      G = [    1.8361  0.0043835    0.31273       0.01    0.76246 ; ...
              0.83052  0.0010001   0.088897   0.010001     0.6545 ; ...
              0.39052      0.001   0.053666   0.024862    0.69485 ];
      law = 1;  sp = 1;
      g.Kp = G(i,1);  g.Ki = G(i,2);  g.Kd = G(i,3);  g.Ks = G(i,4);  g.Tf = G(i,5);
    case 'ES-nofilt'
      %               Kp         Ki         Kd         Ks         tm
      G = [  0.049383     8.6221    0.12762     8.9501      3.018 ; ...
              0.03786     2.7294    0.14593     4.5496          5 ; ...
             0.059744     1.3108    0.11096     2.7639          5 ];
      law = 1;  sp = 1;
      g.Kp = G(i,1);  g.Ki = G(i,2);  g.Kd = G(i,3);  g.Ks = G(i,4);  g.tm = G(i,5);
      g.pw = 1;  g.dw = 1;
    case 'ES-noDw'
      %               Kp         Ki         Kd         Ks         Tf         tm
      G = [  0.047283     8.4268   0.075276     7.6751    0.66785          5 ; ...
              0.45886    0.72776   0.039618    0.34475    0.28402          5 ; ...
              0.26846    0.60332   0.044002    0.46135    0.50324          5 ];
      law = 1;  sp = 1;
      g.Kp = G(i,1);  g.Ki = G(i,2);  g.Kd = G(i,3);  g.Ks = G(i,4);  g.Tf = G(i,5);  g.tm = G(i,6);
      g.pw = 1;  g.dw = 0;
    otherwise
      error('sqpid_ctrl: unknown controller %s', name);
  end
  if g.Tf > 0
    p = 0.5;                                   % slow plant pole (rad/s)
    g.beta = (1 - (1 - p*g.Tf)*exp(-p*L))/p;   % lead constant that cancels it in the gust response
  end
  c.name = name;  c.law = law;  c.sp = sp;  c.g = g;
end

%% =====================================================================================
function [J, parts] = sqpid_cost(name, L, gains)
% SQPID_COST  Robust design cost used to tune EVERY controller.
%
%   J = mean over 7 scenarios of
%         ITAEe + 0.5*max(OS - 2, 0) + rms(du) + 5*max(ess - 0.02, 0)   (+20 if the run diverges)
%   ITAEe = ITAE with the clock restarted at each event (set-point step, gust onset, gust end),
%   OS in %, rms(du) = control chatter in V/sample, ess = residual offset at the end of the gust.
%   Scenarios: nominal | true delay -10% | +10% | plant gain -10% | +10% | double gust | early gust
%   gains : [] to use the tabulated gains, or a struct of gain fields to override.

  S = { {}, {'Lerr', 0.9}, {'Lerr', 1.1}, {'kg', 0.9}, {'kg', 1.1}, {'gamp', 1.0}, ...
        {'t_on', 3, 't_off', 9} };
  parts = zeros(1, numel(S));
  for k = 1:numel(S)
    m = sqpid_run(name, L, S{k}{:}, 'gains', gains);
    if ~isfinite(m.IAE), parts(k) = 1e3; continue; end
    parts(k) = m.sITAE + 0.5*max(m.OS - 2, 0) + m.du + 5*max(m.ess - 0.02, 0) + 20*m.div;
  end
  J = mean(parts);
end

%% =====================================================================================
function [R, txt] = sqpid_e7(R)
% SQPID_E7  Stability analysis table (experiment E7) for the proposed controller.
  prop = R.prop;  Ls = R.Ls;  txt = '';
  txt = [txt sprintf('\nE7. STABILITY ANALYSIS OF %s (nominal model)\n', prop)];
  txt = [txt sprintf('    Loop seen by the tanh nonlinearity: H(s) = (pi/e_max) C(s) G0(s)/(1 + Ks G0(s)),\n')];
  txt = [txt sprintf('    sector of Uq*tanh(alpha*x): (0, k], k = alpha*Uq = 2.  P = frozen probability weight\n')];
  txt = [txt sprintf('    (P = 1 regulation mode ... P = 0 set-point transient; a single row for a constant weight)\n')];
  txt = [txt sprintf('  %4s %5s | %8s %8s %8s | %9s %9s %8s | %12s\n', 'L', 'P', 'GM (dB)', 'PM (deg)', ...
         'wc', 'circle', 'Popov', 'q', 'dL bound (s)')];
  R.pop = cell(1, numel(Ls));
  R.stab = zeros(numel(Ls), 4);        % worst case over P: PM, Popov, dL ; regulation-mode Popov
  for j = 1:numel(Ls)
    s = sqpid_stability(prop, Ls(j));
    R.pop{j} = s;
    for n = 1:numel(s.P)
      txt = [txt sprintf('  %4.1f %5.2f | %8.2f %8.2f %8.3f | %9.4f %9.4f %8.3f | %10.3f\n', Ls(j), s.P(n), ...
             s.GM_dB(n), s.PM_deg(n), s.wc(n), s.circle(n), s.popov(n), s.q(n), s.dL(n))];
    end
    R.stab(j,:) = [min(s.PM_deg) min(s.popov) min(s.dL) s.popov(1)];
    txt = [txt sprintf('  %4.1f worst | PM = %.2f deg, Popov margin = %.4f, guaranteed delay error = %.3f s (%.0f%% of L)\n', ...
           Ls(j), min(s.PM_deg), min(s.popov), min(s.dL), 100*min(s.dL)/Ls(j))];
  end
  txt = [txt sprintf('  circle / Popov: minimum over frequency of the criterion margin (positive = criterion satisfied)\n')];
  txt = [txt sprintf('  GM, PM and dL bound are computed for the sampled linearised loop (zero-order hold, Ts = 0.01 s);\n  dL bound: delay error guaranteed stable by the small-gain test |T Fr| |exp(-j w dL) - 1| < 1\n')];
end

%% =====================================================================================
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

%% =====================================================================================
function sqpid_figures(opt, dpi)
% SQPID_FIGURES  Draw and save every figure of the article from results_option<opt>.mat
%   sqpid_figures(opt)       saves PNG files at 300 dpi in the folder figures_option<opt>
%   sqpid_figures(opt, dpi)  chooses the resolution

  if nargin < 2, dpi = 300; end
  S = load(sprintf('results_option%d.mat', opt));  R = S.R;
  odir = sprintf('figures_option%d', opt);
  if ~exist(odir, 'dir'), mkdir(odir); end
  Ls = R.Ls;  main = R.main;  allc = R.allc;  prop = R.prop;
  nm = numel(main);

  % colour and line style per controller
  sty = struct();
  sty.PID      = {[0.90 0.45 0.00], '-.', 1.4};
  sty.S_PID    = {[0.00 0.45 0.74], '--', 1.4};
  sty.FS_PID   = {[0.20 0.63 0.17], '-',  1.2};
  sty.S_QPID   = {[0.70 0.05 0.15], '-',  2.0};
  sty.ES_QPID  = {[0.70 0.05 0.15], '-',  2.0};
  sty.QPID     = {[0.50 0.50 0.50], ':',  1.4};
  sty.ES_noweight = {[0.55 0.35 0.65], '--', 1.2};
  sty.ES_nofilt   = {[0.35 0.60 0.80], '-.', 1.2};
  sty.ES_noDw     = {[0.55 0.55 0.20], ':', 1.4};
  key = @(n) strrep(n, '-', '_');
  cgust = [1.00 0.95 0.80];
  ton = 8;  toff = 14;

  %% ---------------- step responses ----------------
  fh = newfig(7.2, 7.4);
  for j = 1:numel(Ls)
    subplot(3,1,j);  hold on;
    patch([ton toff toff ton], [-0.1 -0.1 1.45 1.45], cgust, 'EdgeColor', 'none');
    h = gobjects(1, nm+1);
    h(1) = plot([0 24], [1 1], 'k:', 'LineWidth', 1.0);
    for i = 1:nm
      ii = find(strcmp(allc, main{i}));  tr = R.traj{ii,j};  s = sty.(key(main{i}));
      h(i+1) = plot(tr.t, tr.y, s{2}, 'Color', s{1}, 'LineWidth', s{3});
    end
    box on;  grid on;  xlim([0 24]);  ylim([-0.05 1.45]);
    ylabel('Azimuth (deg)');
    title(sprintf('(%s) L = %.1f s', char('a'+j-1), Ls(j)));
    text(ton+0.15, 1.36, 'wind gust', 'FontSize', 8);
    if j == 1, legend(h, [{'Reference'} main], 'Location', 'southeast'); end
    if j == 3, xlabel('Time (s)'); end
  end
  save_png(fh, odir, 'fig_step_response', dpi);

  %% ---------------- control effort (left: start-up, full scale; right: regulation, zoomed) ----------------
  fh = newfig(7.2, 7.4);
  for j = 1:numel(Ls)
    subplot(3,2,2*j-1);  hold on;
    h = gobjects(1, nm+1);
    h(1) = plot([0 4], [20 20], 'k:', 'LineWidth', 1.0);
    for i = 1:nm
      ii = find(strcmp(allc, main{i}));  tr = R.traj{ii,j};  s = sty.(key(main{i}));
      h(i+1) = plot(tr.t, tr.u, s{2}, 'Color', s{1}, 'LineWidth', s{3});
    end
    box on;  grid on;  xlim([0 4]);  ylim([-2 22]);
    ylabel('u (V)');
    title(sprintf('(%s) L = %.1f s, start-up', char('a'+j-1), Ls(j)));
    if j == 1, legend(h, [{'Actuator limit'} main], 'Location', 'northeast'); end
    if j == 3, xlabel('Time (s)'); end
    subplot(3,2,2*j);  hold on;
    patch([ton toff toff ton], [-0.6 -0.6 2.6 2.6], cgust, 'EdgeColor', 'none');
    for i = 1:nm
      ii = find(strcmp(allc, main{i}));  tr = R.traj{ii,j};  s = sty.(key(main{i}));
      plot(tr.t, tr.u, '-', 'Color', s{1}, 'LineWidth', 0.8);
    end
    box on;  grid on;  xlim([5 24]);  ylim([-0.6 2.6]);
    title('regulation and wind gust (zoomed)');
    if j == 3, xlabel('Time (s)'); end
  end
  save_png(fh, odir, 'fig_control_effort', dpi);

  %% ---------------- pointing error around the gust ----------------
  fh = newfig(7.2, 7.4);
  for j = 1:numel(Ls)
    subplot(3,1,j);  hold on;
    yl = [-0.33 0.37];
    patch([ton toff toff ton], [yl(1) yl(1) yl(2) yl(2)], cgust, 'EdgeColor', 'none');
    plot([6 24], [0.02 0.02], 'k:', 'LineWidth', 0.8);  plot([6 24], [-0.02 -0.02], 'k:', 'LineWidth', 0.8);
    h = gobjects(1, nm);
    for i = 1:nm
      ii = find(strcmp(allc, main{i}));  tr = R.traj{ii,j};  s = sty.(key(main{i}));
      h(i) = plot(tr.t, 1 - tr.y, s{2}, 'Color', s{1}, 'LineWidth', s{3});
    end
    box on;  grid on;  xlim([6 24]);  ylim(yl);
    ylabel('Pointing error (deg)');
    title(sprintf('(%s) L = %.1f s', char('a'+j-1), Ls(j)));
    if j == 1, legend(h, main, 'Location', 'northeast'); end
    if j == 3, xlabel('Time (s)'); end
  end
  save_png(fh, odir, 'fig_gust_error', dpi);

  %% ---------------- performance indices ----------------
  fh = newfig(7.2, 3.4);
  idx = [1 2 3];  nmI = {'ISE (deg^2 s)', 'IAE (deg s)', 'ITAE (deg s^2)'};
  for q = 1:3
    subplot(1,3,q);
    V = zeros(numel(Ls), nm);
    for i = 1:nm, V(:,i) = squeeze(R.nom(find(strcmp(allc, main{i})), :, idx(q))); end
    hb = bar(V, 0.85);
    for i = 1:nm, s = sty.(key(main{i}));  set(hb(i), 'FaceColor', s{1}); end
    set(gca, 'XTick', 1:3, 'XTickLabel', {'0.3', '0.8', '1.5'});
    xlabel('Delay L (s)');  ylabel(nmI{q});  grid on;  box on;
    if q == 1, legend(main, 'Location', 'northwest'); end
  end
  save_png(fh, odir, 'fig_indices', dpi);

  %% ---------------- parametric robustness ----------------
  fh = newfig(7.2, 6.6);
  for j = 1:numel(Ls)
    for q = 1:3
      subplot(3,3,(j-1)*3+q);  hold on;
      h = gobjects(1, nm);
      for i = 1:nm
        s = sty.(key(main{i}));
        v = squeeze(R.rob_IAE(i,j,q,:));  v(squeeze(R.rob_ok(i,j,q,:)) < 0.5) = NaN;   % unsettled runs not drawn
        h(i) = plot(100*R.lev, v, [s{2} 'o'], 'Color', s{1}, 'LineWidth', s{3}, 'MarkerSize', 3);
      end
      box on;  grid on;  xlim([-32 32]);
      if j == 1, title(sprintf('%s error', R.par{q})); end
      if q == 1, ylabel(sprintf('IAE,  L = %.1f s', Ls(j))); end
      if j == 3, xlabel('Error (%)'); end
      if j == 1 && q == 1, legend(h, main, 'Location', 'north'); end
    end
  end
  save_png(fh, odir, 'fig_robustness', dpi);

  %% ---------------- Monte Carlo ----------------
  fh = newfig(7.2, 3.4);
  lab = {'IAE (deg s)', 'ITAE (deg s^2)', 'Gust peak error (deg)'};  col = [2 3 5];
  for q = 1:3
    subplot(1,3,q);  hold on;
    M = zeros(numel(Ls), nm);  Sd = M;
    for i = 1:nm
      for j = 1:numel(Ls)
        x = squeeze(R.mc(i,j,:,col(q)));  x = x(isfinite(x));
        M(j,i) = mean(x);  Sd(j,i) = std(x);
      end
    end
    hb = bar(M, 0.85);
    for i = 1:nm, s = sty.(key(main{i}));  set(hb(i), 'FaceColor', s{1}); end
    wgrp = 0.85;
    for i = 1:nm
      if isprop(hb(i), 'XEndPoints'), xo = hb(i).XEndPoints;          % MATLAB R2019b and later
      else, xo = (1:numel(Ls)) + (i - (nm+1)/2)*wgrp/nm; end
      for j = 1:numel(Ls)
        plot([xo(j) xo(j)], [M(j,i)-Sd(j,i) M(j,i)+Sd(j,i)], 'k-', 'LineWidth', 1.0);
      end
    end
    set(gca, 'XTick', 1:3, 'XTickLabel', {'0.3', '0.8', '1.5'});
    xlabel('Delay L (s)');  ylabel(lab{q});  grid on;  box on;
    if q == 1, legend(hb, main, 'Location', 'northwest'); end
  end
  save_png(fh, odir, 'fig_monte_carlo', dpi);

  %% ---------------- step-size study ----------------
  fh = newfig(7.2, 3.4);
  cs = [0.10 0.10 0.10; 0.70 0.05 0.15; 0.00 0.45 0.74; 0.20 0.63 0.17; 0.90 0.45 0.00];
  ls = {'-', '-', '--', '-.', ':'};
  for j = 1:numel(Ls)
    subplot(1,3,j);  hold on;
    h = gobjects(1, numel(R.steps));  lg = cell(1, numel(R.steps));
    for a = 1:numel(R.steps)
      tr = R.stptraj{j,a};
      h(a) = plot(tr.t, tr.y, ls{a}, 'Color', cs(a,:), 'LineWidth', 1.3);
      lg{a} = sprintf('%g deg', R.steps(a));
    end
    box on;  grid on;  xlim([0 10]);  ylim([0 1.15]);
    xlabel('Time (s)');  title(sprintf('L = %.1f s', Ls(j)));
    if j == 1, ylabel('y / step size'); end
    if j == 3, legend(h, lg, 'Location', 'southeast'); end
  end
  save_png(fh, odir, 'fig_step_size', dpi);

  %% ---------------- Popov plot ----------------
  fh = newfig(7.2, 3.4);
  for j = 1:numel(Ls)
    s = R.pop{j};
    subplot(1,3,j);  hold on;
    X = real(s.H);  Y = s.w.*imag(s.H);
    plot(X, Y, '-', 'Color', [0.70 0.05 0.15], 'LineWidth', 1.6);
    sel = X > -1/s.k - 0.6;
    xl = [-1/s.k - 0.45, max(0.6, 1.08*max(X))];
    yl = [max(1.12*min(Y(sel)) - 0.05, -3), 0.5];
    q = s.q(1);                                  % regulation mode (first entry)
    if q > 0
      xx = [xl(1), xl(2)];  plot(xx, (xx + 1/s.k)/q, 'k--', 'LineWidth', 1.0);
    else
      plot([-1/s.k -1/s.k], yl, 'k--', 'LineWidth', 1.0);
    end
    plot(-1/s.k, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    box on;  grid on;  xlim(xl);  ylim(yl);
    xlabel('Re H(j\omega)');  title(sprintf('L = %.1f s', Ls(j)));
    if j == 1, ylabel('\omega Im H(j\omega)'); end
    if j == 3, legend({'Popov locus', 'Popov line', '-1/k'}, 'Location', 'southeast'); end
  end
  save_png(fh, odir, 'fig_popov', dpi);

  %% ---------------- Option 2 only ----------------
  if opt == 2
    % ablation of the design cost
    fh = newfig(7.2, 3.0);
    names = {'PID', 'S-PID', 'FS-PID', 'S-QPID', 'ES-noweight', 'ES-nofilt', 'ES-noDw', 'ES-QPID'};
    V = zeros(numel(Ls), numel(names));
    for i = 1:numel(names), V(:,i) = R.J(find(strcmp(allc, names{i})), :)'; end
    hb = bar(V, 0.9);
    for i = 1:numel(names), s = sty.(key(names{i}));  set(hb(i), 'FaceColor', s{1}); end
    set(hb(4), 'FaceColor', [0.95 0.60 0.65]);
    set(gca, 'XTick', 1:3, 'XTickLabel', {'L = 0.3 s', 'L = 0.8 s', 'L = 1.5 s'});
    ylabel('Robust design cost J');  grid on;  box on;
    legend(names, 'Location', 'northwest');
    save_png(fh, odir, 'fig_ablation', dpi);

    % probability weight and its effect (L = 0.8 s)
    fh = newfig(7.2, 4.6);
    j = 2;
    a = R.traj{find(strcmp(allc, 'ES-QPID')), j};  b = R.traj{find(strcmp(allc, 'ES-noweight')), j};
    subplot(2,1,1);  hold on;
    patch([ton toff toff ton], [-0.05 -0.05 1.3 1.3], cgust, 'EdgeColor', 'none');
    h1 = plot(a.t, a.p, '-', 'Color', [0.10 0.10 0.10], 'LineWidth', 1.6);
    box on;  grid on;  xlim([0 24]);  ylim([-0.05 1.1]);
    ylabel('P_{reg}');  title('(a) probability weight of the regulation state, L = 0.8 s');
    subplot(2,1,2);  hold on;
    patch([ton toff toff ton], [-0.1 -0.1 1.3 1.3], cgust, 'EdgeColor', 'none');
    plot([0 24], [1 1], 'k:', 'LineWidth', 1.0);
    h2 = plot(b.t, b.y, '--', 'Color', sty.ES_noweight{1}, 'LineWidth', 1.4);
    h3 = plot(a.t, a.y, '-',  'Color', sty.ES_QPID{1}, 'LineWidth', 2.0);
    box on;  grid on;  xlim([0 24]);  ylim([-0.05 1.3]);
    xlabel('Time (s)');  ylabel('Azimuth (deg)');  title('(b) response with and without the probability weights');
    legend([h2 h3], {'constant weight 0.5 (ES-noweight)', 'two-qubit weights (ES-QPID)'}, 'Location', 'southeast');
    save_png(fh, odir, 'fig_weights', dpi);
  end
  fprintf('Figures saved in %s\n', odir);
end

function fh = newfig(w, h)
  fh = figure('Visible', 'off', 'PaperUnits', 'inches', 'PaperPosition', [0 0 w h], ...
              'PaperSize', [w h], 'Position', [50 50 round(w*100) round(h*100)]);
  set(fh, 'DefaultAxesFontSize', 9, 'DefaultTextFontSize', 9);
end

function save_png(fh, odir, name, dpi)
  print(fh, fullfile(odir, [name '.png']), '-dpng', sprintf('-r%d', dpi));
  close(fh);
end

%% =====================================================================================
function g = tune_gains(name, L, maxit)
% TUNE_GAINS  Local Nelder-Mead polish of the gains of one controller on the robust design cost.
%
%   g = tune_gains('ES-QPID', 0.8)        % starts from the tabulated gains (sqpid_ctrl.m)
%   g = tune_gains('PID', 1.5, 200)       % at most 200 iterations (default 150)
%
%   The same cost (sqpid_cost.m) and the same procedure were used for every controller.
%   The tabulated gains were found by a global search (differential evolution, population
%   16 x number of gains, 200 generations, logarithmic bounds) followed by this polish.
%   One cost evaluation runs 7 simulations, so a polish takes several minutes.

  if nargin < 3, maxit = 150; end
  c = sqpid_ctrl(name, L);
  f = {'Kp', 'Ki', 'Kd'};
  if c.law == 1, f{end+1} = 'Ks'; end
  if c.g.Tf > 0, f{end+1} = 'Tf'; end
  if c.g.pw > 0, f{end+1} = 'tm'; end
  x0 = zeros(1, numel(f));
  for k = 1:numel(f), x0(k) = log10(c.g.(f{k})); end
  J0 = sqpid_cost(name, L, []);
  cost = @(x) sqpid_cost(name, L, cell2struct(num2cell(10.^x(:)), f(:), 1));
  [x, J] = fminsearch(cost, x0, optimset('MaxIter', maxit, 'TolX', 1e-4, 'TolFun', 1e-6, 'Display', 'iter'));
  g = cell2struct(num2cell(10.^x(:)), f(:), 1);
  fprintf('\n%s, L = %.1f s:  J = %.4f (tabulated gains)  ->  %.4f (after polish)\n', name, L, J0, J);
  disp(g);
end
