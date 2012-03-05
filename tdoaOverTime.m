function T = tdoaOverTime(lr, t, frame, method)

% T = tdoaOverTime(lr, t, method)
%
% Find the probability of each time-delay of arrival listed in the
% vector t from the left and right channels in the signal lr.  T is
% a matrix with one column per column of the spectrogram for LR's
% channels and one row per entry in t specifying the generalized
% cross-correlation at each of those delays for each FFT frame.

if(nargin < 3) frame = 512; end
if(nargin < 4) method = 'phat'; end

E = probCC(lr, t, frame);

switch lower(method)
 case {'e','entropy','ent'}
  s = size(E);
  ee = reshape(E, [s(1) s(2)*s(3)]);
  h = npEnt(angle(ee), 1);
  T = reshape(h, s(2), s(3))';
 case {'gcc', 'phat', 'p', 'g'}
  T = squeeze(real(sum(E ./ abs(E), 1)))';
 case 'a'
  T = squeeze(abs(sum(E ./ abs(E), 1)))'; 
 case 'i'
  T = squeeze(imag(sum(E ./ abs(E), 1)))'; 
 case {'like', 'likelihood', 'l'}
  T = squeeze(sum(log(probDp(angle(E))), 1))';
 otherwise
  error(['unrecognized method: ' method])
end
