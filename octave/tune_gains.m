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
  printf('\n%s, L = %.1f s:  J = %.4f (tabulated gains)  ->  %.4f (after polish)\n', name, L, J0, J);
  disp(g);
end
