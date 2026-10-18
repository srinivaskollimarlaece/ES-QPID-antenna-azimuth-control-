function [Ad, Bd, Cd] = sqpid_zoh(K, den, Ts)
% SQPID_ZOH  Zero-order-hold discretisation of G(s) = K/(s^2 + den(2) s + den(3)).
%   Base MATLAB only (expm); no Control System Toolbox needed.
%   Equivalent to [Ad,Bd,Cd] = ssdata(c2d(ss(tf(K,den)),Ts)) up to a change of state coordinates.
  A = [-den(2) -den(3); 1 0];  B = [1; 0];  C = [0 K];
  M = expm([A B; zeros(1,3)]*Ts);
  Ad = M(1:2,1:2);  Bd = M(1:2,3);  Cd = C;
end
