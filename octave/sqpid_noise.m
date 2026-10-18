function [nz, wn] = sqpid_noise(seed, N)
% SQPID_NOISE  Reproducible unit-variance sequences: measurement noise nz, turbulence driver wn.
%   seed = 0 : nominal sequences used for every table and figure of the nominal study
%   seed > 0 : independent realisation (Monte Carlo)
%   Uses Octave's legacy generator, as in the original scripts: randn('seed', .)
  persistent NZ0 WN0
  if seed == 0
    if isempty(NZ0)
      randn('seed', 1);  NZ0 = randn(1, 4000);  WN0 = randn(1, 4000);
    end
    nz = NZ0(1:N);  wn = WN0(1:N);
  else
    randn('seed', 1000 + seed);  nz = randn(1, N);  wn = randn(1, N);
  end
end
