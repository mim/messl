function [p_lr_iwt params hardMasks] = messl(lr, tau, I, varargin)

% [p_lr_iwt params hardMasks] = messl(lr, tau, I, [name1, value1] ...)
%
% Perform MESSL algorithm (including ILD, frequency-dependent ITD,
% source prior GMMs) on a stereo mixture of I spatially separated
% sources.
%
% LR is the stereo mixture in the time domain, a 2xN matrix.  TAU is
% the grid of times over which to evaluate the probability of the
% various samples.  The three "mode" options behave in similar
% ways, but control different variables.  A mode of 0 indicates
% that that feature should not be used.  A positive mode specifies
% the number of bands to break the frequency axis into, a separate
% parameter being used for each one.  A negative number specifies
% how many frequencies to use per band.
%
% The returned arguments are the probability of each point in the
% spectrogram coming from each source and the parameters for each
% of the source models.
% 
% Valid named arguments are:
%
% ::Initialization::
% tauPosInit     ([]) init source pos's (in samples), otherwise from xcorr
% pTauIInit      ([]) initial p(tau|i) distributions, gaussians
%                     around tauPosInit if not specified here
% sigmaInit      ([]) initial stddev of IPD residual
% xiInit         ([]) initial IPD mean
% ildInit         (0) WxI, initial mean ILD, in dB, can be matrix for freq dep
% ildStdInit     (10) WxI, initial std ILD, in dB, can be matrix for freq dep
% maskInit       ([]) to specify a permanent prior on the WxTxI mask
% maskHold        (0) number of iterations to hold the prior mask
%
% ::Mode selection::
% modes          ([]) vector of [ipd ild sp xi sigma dct] modes all together
% ipdMode         (1) 0 => no IPD, 1 => use IPD as defined by
%                     xiMode and sigmaMode
% ildMode        (-1) 0 => no ILD, 1 => freq indep, -1 => freq dep
%                    -2
% spMode          (0) 0 => no SP, 1 => model L and R responses separately,
%                    -1 => model L+R
% xiMode         (-1) 0 => ipd mean=0, 1 => freq indep, -1 => freq dep
% sigmaMode      (-1) 0 => varies by source, 1 => varies by src,tau,
%                    -1 => varies by src,tau,freq
% dctMode         (0) 0 => use canonical basis for ILD and SP parameters,
%                     n => use n DCT bases for ILD and SP params
%
% ::Extended features::
% garbageSrc      (0) add a garbage source with special init and updates
% sourcePriors   ([]) list of I GMM structures.  Not used if empty
% reliability    ([]) relative weights for each spectrogram point
%                     for use in parameter estimation
% ildPriorPrec    (0) precision (inverse variance) of ILD prior
% sr          (16000) sampling rate, only used by ILD prior
%
% ::Nuts and bolts::
% vis             (0) plot informational displays
% nfft         (1024) window size of FFT
% Nrep           (16) number of EM iterations

[tauPosInit, pTauIInit, ildInit, ildStdInit, maskInit, garbageSrc, ...
 ipdMode, ildMode, xiMode, sigmaMode, dctMode, spMode, nfft, ...
 vis, Nrep, modes, sigmaInit, xiInit, sourcePriors, maskHold, ...
 reliability, ildPriorPrec, sr, mrfHardCompatExp, mrfCompatFile ...
 mrfCompatExpSched fixIPriors mrfLbpIter] = ...
    process_options(varargin, 'tauPosInit', [], 'pTauIInit', [], ...
    'ildInit', 0, 'ildStdInit', 10, 'maskInit', [], ...
    'garbageSrc', 0, 'ipdMode', 1, 'ildMode', -1, 'xiMode', -1, ...
    'sigmaMode', -1, 'dctMode', 0, 'spMode', 0, 'nfft', 1024, 'vis', 0, ...
    'Nrep', 16, 'modes', [], 'sigmaInit', [], 'xiInit', [], ...
    'sourcePriors', [], 'maskHold', 0, 'reliability', [], ...
    'ildPriorPrec', 0, 'sr', 16000, 'mrfHardCompatExp', 0, ...
    'mrfCompatFile', '', 'mrfCompatExpSched', [0 0 0 0 .02 .02 .02 .02 .05], ...
    'fixIPriors', 0, 'mrfLbpIter', 8);

if ~isempty(modes)
  ipdMode   = modes(1);
  ildMode   = modes(2);
  spMode    = modes(3);
  xiMode    = modes(4);
  sigmaMode = modes(5);
  dctMode   = modes(6);
end

if spMode && isempty(sourcePriors)
  spMode = 0;
  warning('spMode set to 1, but ''sourcePriors'' option was not set.');
end

[A angE cc W T Nt L R lr] = messlObsDerive(lr, tau, nfft);

% Giant matrix holding the responsibility of each part of the model
% for each point in the spectrogram.  Summing over tau and i
% gives a matrix of 1s.  Nu(w, t, i, tau)
% nuBin = 1/(I*Nt) * ones(W,T, I,Nt);

% Initialize the probability distributions and other parameters
[ipdParams itds] = messlIpdInit(I, W, Nrep, tau, sr, lr, tauPosInit, pTauIInit, ...
                           sigmaInit, xiInit, ipdMode, xiMode, sigmaMode, ...
                           garbageSrc, vis, fixIPriors);
clear lr
%FIXME
ildParams = messlIldInit(I, W, sr, Nrep, ildInit, ildStdInit, ...
                    ildPriorPrec*T/100, ildMode, itds, ...
                    false&dctMode,  garbageSrc);
                    %dctMode,  garbageSrc);
[spParams C] = messlSpInit(I, W, L, R, sourcePriors, ildStdInit, dctMode, ...
                      spMode, garbageSrc);
                  
mrfCompatPot = messlMrfLoadCompat(mrfCompatFile, I, garbageSrc);
                  
% The rest of the code should act like the garbage source is like
% any other source, except for the M step, which has to know about
% it so it doesn't update it
I = I + garbageSrc;

% Inialize the mask prior
logMaskPrior = [];
if ~isempty(maskInit)
  logMaskPrior = single(log(maskInit));
end

% Keep track of the total log likelihood
ll = [];

% Start EM
for rep=1:Nrep
  clear nu*

  % Turn on SP mode if we've reached the appropriate iteration.
  if rep == spParams.spStartRep
    spParams.spMode = spParams.origSpMode;

    % The GMMs could have been passed to this function in an
    % arbitrary order.  Fix this by permuting the gmms so that they
    % match the binaural sources.
    if spParams.spMode && ~isfield(spParams, 'ev_params')
      maskBin = maskIpd .* maskIld;
      maskBin = maskBin ./ repmat(sum(maskBin,3), [1 1 I]);
      spParams.sourcePriors(1:I-garbageSrc) = messlSpPermuteGmms( ...
          spParams.sourcePriors,  L, R, maskBin, I-garbageSrc);
    end
  end

 
  %%%% E step: calculate nu matrix
  lpIpd = 0;  lpIld = 0;  lpSp = 0;
  if ipdParams.ipdMode
    lpIpd = messlIpdLogLikelihood(W,T,I,Nt,C,rep,Nrep, ipdParams, angE);
    if any(~isfinite(lpIpd(:))), warning('IPD liklihood is not finite'); end
  end
  if ildParams.ildMode
    lpIld = messlIldLogLikelihood(W,T,I,Nt,C,rep,Nrep, ildParams, A);
    if any(~isfinite(lpIld(:))), warning('ILD liklihood is not finite'); end
  end
  if spParams.spMode
    lpSp = messlSpLogLikelihood(W,T,I,Nt,C,rep,Nrep, spParams, ...
        ildParams, L, R);
    if any(~isfinite(lpSp(:))), warning('SP liklihood is not finite'); end
  end

  % Combine binaural and GMM likelihoods and normalize:
  [ll(rep) p_lr_iwt nuIpd maskIpd nuIld maskIld nuSp maskSp] = ...
      messlPosterior(W, T, I, Nt, C, logMaskPrior, ...
      ipdParams.ipdMode, lpIpd, ildParams.ildMode, lpIld, spParams.spMode, lpSp, ...
      vis || rep == Nrep, reliability, ...
      mrfCompatPot, mrfCompatExpSched(min(end,rep)), mrfLbpIter);
  clear lp*
  if rep >= maskHold
    logMaskPrior = [];
  end

  % ll should be non-decreasing
  fprintf('ll(%02d) = %e\n', rep, ll(rep));

  
  %%%% M step: use nu matrix to calcuate parameters
  nuIpd = messlIpdEnforcePriors(nuIpd, ipdParams);
      
  if ipdParams.ipdMode
    ipdParams = messlIpdUpdateParams(W, T, I, Nt, C, rep, ipdParams, ...
                                nuIpd, angE);
  end
  if ildParams.ildMode
    ildParams = messlIldUpdateParams(W, T, I, Nt, C, rep, ildParams, nuIld, ...
        A, Nrep);
  end
  if spParams.spMode
    spParams = messlSpUpdateParams(W, T, I, Nt, C, rep, spParams, nuSp, ...
        L, R, Nrep);
  end
  
  if vis
    messlUtilVisualizeParams(W, T, I, tau, sr, ipdParams, ildParams, spParams, ...
        p_lr_iwt, maskIpd, maskIld, maskSp, L, R, reliability);
  end
end

params = struct('p_tauI', ipdParams.p_tauI, ...
    'xi_wit', ipdParams.xi_wit, 's_wit', sqrt(ipdParams.s2_wit), ...
    'mu_wi', ildParams.mu_wi, 'h_wi', sqrt(ildParams.h2_wi), ...
    'ipdParams', ipdParams, 'maskIpd', maskIpd, ...
    'ildParams', ildParams, 'maskIld', maskIld, ...
    'spParams', spParams, 'maskSp', maskSp, 'C', C, ...
    'garbageSrc', garbageSrc);

% Compute hard masks, potentially using the MRF model
[~,~,~,hardSrcs] = messlMrfApply(nuIld, nuIpd, p_lr_iwt, mrfCompatPot, mrfHardCompatExp, mrfLbpIter, 'max');
hardMasks = zeros(size(p_lr_iwt));
for i = 1:max(hardSrcs(:))
    hardMasks(:,:,:,i) = repmat(permute(hardSrcs == i, [3 1 2]), [2 1 1]);
end

