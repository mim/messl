function m = prob2mask(p, mode, b, a)
%
% m = prob2mask(p, mode, b, a)
%
% Converts a spectrogram's worth of probabilities into a soft mask
% using a point-wise non-linearity mapping [0,1] to [0,1].  When mode
% = 0, uses the sigmoid function 1 ./ (1 + exp(-(x-a)*b*5)), where a
% sets the offset and b sets the sharpness or slope.  When mode = 1, a
% piecewise linear function is used which is 0.5 at a, with slope b.
% Reasonable values of a are in [0,1], and reasonable values of b are
% in [1,200] for both modes.

if ~exist('a', 'var'), a = 0.5; end
if ~exist('mode', 'var'), mode = 0; end
if ~exist('b', 'var')
  if mode == 0
    b = 2.89;
  else 
    b = 1.7;
  end
end

if mode == 0
  m = 1 ./ (1 + exp(-b*(p - a)*5));
else
  % y = bx + c
  % .5 = ba + c
  % c = .5 - ba
  m = max(0, min(1, b*(p-a) + .5));
end
