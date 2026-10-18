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
