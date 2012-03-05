function mu = ildPrior(itds, freqs, b, model)

% mu = ildPrior(itds, freqs[, b[, model]])
%
% Compute the ILD prior for a given set of ITDs and frequencies.
% The prior is based on a regression that was run on ITD vs
% frequency.  ITDs are measured in milliseconds, and frequencies
% are measured in kHz.  Returned value is the mean of the ILD prior,
% measured in dB.

if ~exist('b', 'var')
  b = [-0.06945 -16.01 2.636 -0.1702 0.0036 0.431 ...
       -0.07007 0.002534 -0.4577 0.09647 -0.003942];
end
if ~exist('model', 'var')
  model = [0 0; 1 1; 1 2; 1 3; 1 4; 2 2; 2 3; 2 4; 4 2; 4 3; 4 4];
end

itd_mat = repmat(itds(:), 1, length(freqs));
frq_mat = repmat(freqs(:)', length(itds), 1);

X = x2fx([itd_mat(:) frq_mat(:)], model);

mu = reshape(b * X', length(itds), [])';
