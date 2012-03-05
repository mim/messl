function t = reconstruct(mask, varargin)

% t = reconstruct(mask, L, R, target)
% t = reconstruct(mask, LR, target)
% t = reconstruct(mask, LR)
%
% Reconstruct a stereo waveform from the left and right channel
% mixtures and a mask.  Assumes the hop size of the stft was 1/4 the
% window size.  Mask, L, and R are all WxTxI, with the third dimension
% being the source number.  Can also use LR, which can either be a
% stereo time-domain waveform or a WxTx2 tensor, the third dimension
% being the ear.  Target selects which signal to reconstruct and
% defaults to 1.

switch length(varargin)
 case 3
  [L R target] = deal(varargin{:});
  LR = cat(3, sum(L,3), sum(R,3));
 case 2
  if numel(varargin{2}) == 1
    [LR target] = deal(varargin{:});
  else
    [L R] = deal(varargin{:});
    LR = cat(3, sum(L,3), sum(R,3));
    target = 1;
  end
 case 1
  LR = varargin{1};
  target = 1;
 otherwise
  error('Too many arguments')
end

if ndims(mask) == 4
  mask = squeeze(mean(mask,1));
end

[F,T,I] = size(mask);
z       = zeros(1, T);
nfft    = (F+1)*2;

[L R] = binSpec(LR, nfft);

masked_L = [z; L .* mask(:,:,target); z];
masked_R = [z; R .* mask(:,:,target); z];

wfL = istft(masked_L, nfft, nfft, nfft/4);
wfR = istft(masked_R, nfft, nfft, nfft/4);

t = [wfR; wfL]';
