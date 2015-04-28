function [A angE cc W T Nt L R x] = messlObsDerive(x, tau, nfft)

% Get cross correlation values semi-digested for probability
% measures
E = single(probCC(x, tau, nfft));
A = dB(E(:,:,1));
angE = angle(E);
cc = real(squeeze(sum(sum(E ./ abs(E), 1), 2)));

% Find appropriate dimension sizes
W = size(E,1); 
T = size(E,2); 
Nt = length(tau);

% SP observations
if ndims(x) == 3
  L = x(:,:,1); R = x(:,:,2);
else
  [L,R] = binSpec(x, nfft);
end
x = cat(3, L, R);
L = dB(L);  L = L - max(max(L));
R = dB(R);  R = R - max(max(R));
%clear lr
