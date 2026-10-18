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
