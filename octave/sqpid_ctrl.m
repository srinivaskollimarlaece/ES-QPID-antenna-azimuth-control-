function c = sqpid_ctrl(name, L)
% SQPID_CTRL  Controller definition and tuned gains for delay L (0.3, 0.8 or 1.5 s).
%
%   All gains were obtained with the SAME procedure for every controller: minimisation of the
%   robust design cost of sqpid_cost.m (7 scenarios: nominal, delay -10%/+10%, gain -10%/+10%,
%   double gust, early gust) by differential evolution followed by a Nelder-Mead polish.
%   tune_gains.m repeats the local polish in Octave.
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
