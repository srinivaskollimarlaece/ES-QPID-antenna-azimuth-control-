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
