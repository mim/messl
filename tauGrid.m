function tau = tauGrid(d_m, fs, nTau, c_mps)

% Compute a grid of ITDs
%
%   tau = tauGrid(d_m, nTau, fs, c_mps)
%
% d_m    distance between mics in meters
% fs     sampling rate in Hz
% nTau   number of tau values in grid (default: 31)
% c_mps  speed of sound in meters per second (default: 340)

if ~exist('c_mps', 'var') || isempty(c_mps), c_mps = 340; end
if ~exist('nTau', 'var') || isempty(nTau), nTau = 31; end

maxItd = 2 * fs * d_m / c_mps;  % maximum ITD [samples]
tau = linspace(-maxItd, maxItd, nTau);
