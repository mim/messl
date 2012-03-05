function [E,w] = probCC(lr, t, frame, doProduct)

% [E,w] = probCC(lr, t, frame, doProduct)
%
% Do some sort of probabilistic cross correlation.  LR is a 2xN matrix
% with the left channel on top and the right channel on the bottom.
% It could also be an FxTx2 complex spectrogram matrix from the two
% ears, with the left channel in lr(:,:,1) and the right in
% lr(:,:,2). T is a vector of lags to evaluate the cross correlation
% over.  Returns the full distribution of the ratio of the two
% channels for each lag.  Frame is the size of the FFT window to use.

if(nargin < 3) frame = 1024; end
if(nargin < 4) doProduct = 0; end

[L,R] = binSpec(lr, frame);

frame = 2*(size(L,1)+1);
  
w = [1:frame/2-1]' * 2*pi / frame;

if doProduct
  % A real xcorr, but doesn't give ILD
  LonR = L .* conj(R);
else
  % Not the real xcorr, but does give ILD
  LonR = L ./ R;
end

E = zeros([size(L) length(t)]);
for i=1:length(t)
  predicted = sparse(diag(1./exp(j * w * t(i))));
  E(:,:,i) = predicted * LonR;
end
