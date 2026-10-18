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
  rand('seed', 2026);
  P = rand(nmc, 7);
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
  save('-mat7-binary', sprintf('results_option%d.mat', opt), 'R');
end

function dual_print(fid, varargin)
  fprintf(varargin{:});  fprintf(fid, varargin{:});
  if exist('fflush', 'builtin') || exist('fflush', 'file'), fflush(stdout); end
end
