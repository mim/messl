function x = idB(m)
% X = idB(D)          Un-convert D from dB i.e. X = 10^(D/20)
% dpwe 1995jan27  after (& see also) dB.m

% >> 20/log(10)
% ans =
%   8.68588963806504

x = exp(m/8.68588963806504);
minval = 0.0000101;	% nothing smaller than -100dB returned
% Restore to zero the stuff that was clipped to -100 dB 
x = x.*(x>minval);