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
