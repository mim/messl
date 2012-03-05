function d = dB(x)
% D = dB(X)          Convert X to dB i.e. D = 20 * log10(X).  Clips at -100 dB.
% dpwe 1994apr14  see also idB.m
ax = real(x).^2 + imag(x).^2;
minval = 1.0201e-10; % nothing smaller than -100dB returned
mm = ax < minval;
% >> 10/log(10)
% ans =
%   4.3429448190325175005455
d = 4.3429448190325175005455*log(ax.*(1-mm)+mm*minval);
