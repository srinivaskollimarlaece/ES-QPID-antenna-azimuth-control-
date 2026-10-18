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
