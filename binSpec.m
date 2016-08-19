function [L,R,hop,nfft] = binSpec(lr, nfft, hop)

% [L,R,hop,nfft] = binSpec(lr, nfft, hop)
%
% Make a spectrogram of the left and right channels of the waveform
% LR.  LR(1,:) is the left channel, LR(2,:) is the right channel.
% NFFT is the size of the window to use and the number of samples
% in each window, HOP is the number of samples between window
% starts.  Note that this also chops off the DC and nyquist
% components.

if(nargin < 2) nfft = 1024; end
if(nargin < 3) hop = nfft/4; end


if(ndims(lr) == 3)
  % Already specgrammed
  L = lr(:,:,1);
  R = lr(:,:,2);
  nfft = 2*(size(lr,1)+1);
else
  % NB: stft gives negative phase of specgram, so need to reverse the
  % left and right signals
  L = stft(lr(2,:), nfft, nfft, hop);
  R = stft(lr(1,:), nfft, nfft, hop);

  % Probably the right thing to do, but not consistent with experiments
  %L = stft(lr(2,:), 2*nfft, nfft, hop);
  %R = stft(lr(1,:), 2*nfft, nfft, hop);

  L = L(2:end-1, :);
  R = R(2:end-1, :);
end
