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
 mrfCompatExpSched fixIPriors] = ...
    process_options(varargin, 'tauPosInit', [], 'pTauIInit', [], ...
    'ildInit', 0, 'ildStdInit', 10, 'maskInit', [], ...
    'garbageSrc', 0, 'ipdMode', 1, 'ildMode', -1, 'xiMode', -1, ...
    'sigmaMode', -1, 'dctMode', 0, 'spMode', 0, 'nfft', 1024, 'vis', 0, ...
    'Nrep', 16, 'modes', [], 'sigmaInit', [], 'xiInit', [], ...
    'sourcePriors', [], 'maskHold', 0, 'reliability', [], ...
    'ildPriorPrec', 0, 'sr', 16000, 'mrfHardCompatExp', 0, ...
    'mrfCompatFile', '', 'mrfCompatExpSched', [0 0 0 0 .02 .02 .02 .02 .05], ...
    'fixIPriors', 0);

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

% Get cross correlation values semi-digested for probability
% measures
E = single(probCC(lr, tau, nfft));
A = dB(E(:,:,1));
angE = angle(E);
cc = real(squeeze(sum(sum(E ./ abs(E), 1), 2)));

% Find appropriate dimension sizes
W = size(E,1); T = size(E,2); Nt = length(tau);
clear E

% SP observations
if ndims(lr) == 3
  L = lr(:,:,1); R = lr(:,:,2);
else
  [L,R] = binSpec(lr, nfft);
end
lr = cat(3, L, R);
L = dB(L);  L = L - max(max(L));
R = dB(R);  R = R - max(max(R));
%clear lr

% DCT basis vectors
if dctMode
  B = zeros(W, dctMode);
  B(:,1) = 1/sqrt(W);
  B(:,2:end) = sqrt(2/W) * cos(pi/W * [0.5:W-0.5]'*[1:dctMode-1]);
else
  B = eye(W);
end
  
% Giant matrix holding the responsibility of each part of the model
% for each point in the spectrogram.  Summing over tau and i
% gives a matrix of 1s.  Nu(w, t, i, tau)
% nuBin = 1/(I*Nt) * ones(W,T, I,Nt);

% Initialize the probability distributions and other parameters
[ipdParams itds] = initIpd(I, W, Nrep, tau, sr, lr, tauPosInit, pTauIInit, ...
                           sigmaInit, xiInit, ipdMode, xiMode, sigmaMode, ...
                           garbageSrc, vis, fixIPriors);
clear lr
%FIXME
ildParams = initIld(I, W, sr, Nrep, ildInit, ildStdInit, ...
                    ildPriorPrec*T/100, ildMode, itds, ...
                    false&dctMode,  garbageSrc, B);
                    %dctMode,  garbageSrc, B);
[spParams C] = initSp(I, W, L, R, sourcePriors, ildStdInit, dctMode, ...
                      spMode, garbageSrc, B);
                  
mrfCompatPot = loadMrfCompat(mrfCompatFile, I, garbageSrc);
                  
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
      spParams.sourcePriors(1:I-garbageSrc) = permuteGmms( ...
          spParams.sourcePriors,  L, R, maskBin, I-garbageSrc);
    end
  end

 
  %%%% E step: calculate nu matrix
  lpIpd = 0;  lpIld = 0;  lpSp = 0;
  if ipdParams.ipdMode
    lpIpd = computeIpdLogLikelihood(W,T,I,Nt,C,rep,Nrep, ipdParams, angE);
    if any(~isfinite(lpIpd(:))), warning('IPD liklihood is not finite'); end
  end
  if ildParams.ildMode
    lpIld = computeIldLogLikelihood(W,T,I,Nt,C,rep,Nrep, ildParams, A);
    if any(~isfinite(lpIld(:))), warning('ILD liklihood is not finite'); end
  end
  if spParams.spMode
    lpSp = computeSpLogLikelihood(W,T,I,Nt,C,rep,Nrep, spParams, ...
        ildParams, L, R);
    if any(~isfinite(lpSp(:))), warning('SP liklihood is not finite'); end
  end

  % Combine binaural and GMM likelihoods and normalize:
  [ll(rep) p_lr_iwt nuIpd maskIpd nuIld maskIld nuSp maskSp] = ...
      computePosterior(W, T, I, Nt, C, logMaskPrior, ...
      ipdParams, lpIpd, ildParams, lpIld, spParams, lpSp, ...
      vis || rep == Nrep, reliability, ...
      mrfCompatPot, mrfCompatExpSched(min(end,rep)));
  clear lp*
  if rep >= maskHold
    logMaskPrior = [];
  end

  % ll should be non-decreasing
  fprintf('ll(%02d) = %e\n', rep, ll(rep));

  
  %%%% M step: use nu matrix to calcuate parameters
  nuIpd = enforceIPriors(nuIpd, ipdParams);
      
  if ipdParams.ipdMode
    ipdParams = updateIpdParams(W, T, I, Nt, C, rep, ipdParams, ...
                                nuIpd, angE);
  end
  if ildParams.ildMode
    ildParams = updateIldParams(W, T, I, Nt, C, rep, ildParams, nuIld, ...
        A, Nrep);
  end
  if spParams.spMode
    spParams = updateSpParams(W, T, I, Nt, C, rep, spParams, nuSp, ...
        L, R, Nrep);
  end
  
  if vis
    visualizeParams(W, T, I, tau, sr, ipdParams, ildParams, spParams, ...
        p_lr_iwt, maskIpd, maskIld, maskSp, L, R, reliability);
  end
end

params = struct('p_tauI', ipdParams.p_tauI, ...
    'xi_wit', ipdParams.xi_wit, 's_wit', sqrt(ipdParams.s2_wit), ...
    'mu_wi', ildParams.mu_wi, 'h_wi', sqrt(ildParams.h2_wi), ...
    'ipdParams', ipdParams, 'maskIpd', maskIpd, ...
    'ildParams', ildParams, 'maskIld', maskIld, ...
    'spParams', spParams, 'maskSp', maskSp);

% Compute hard masks, potentially using the MRF model
mrfLbpIter = 8;
[~,~,~,hardSrcs] = applyMrf(nuIld, nuIpd, p_lr_iwt, mrfCompatPot, mrfHardCompatExp, mrfLbpIter, 'max');
hardMasks = zeros(size(p_lr_iwt));
for i = 1:max(hardSrcs(:))
    hardMasks(:,:,:,i) = repmat(permute(hardSrcs == i, [3 1 2]), [2 1 1]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IPD functions
function [ipdParams, itds] = ...
    initIpd(I, W, Nrep, tau, sr, lr, tauPosInit, pTauIInit, sigmaInit, ...
            xiInit, ipdMode, xiMode, sigmaMode, garbageSrc, vis, fixIPriors)
% Initialize p(tau | i) with smooth peaks centered on the peaks of
% the cross-correlation, or at locations specified by the user.
% Could also just supply the full probability distributions.

% IPD parameters
Nt = length(tau);

if isempty(xiInit)
  xi_wit = zeros(W,I+garbageSrc,Nt);
  xiBands = makeBands(xiMode, W, 1, Nrep);
else
  xi_wit = xiInit;
  xiBands = makeBands(xiMode, W, W, Nrep);
end

if isempty(sigmaInit)
  s2_wit = 1^2*ones(W,I,Nt);
  if garbageSrc, s2_wit = cat(2, s2_wit, 3^2*ones(W,1,Nt)); end
  sigmaBands = makeBands(sigmaMode, W, 1, Nrep);
else
  s2_wit = sigmaInit.^2;
  sigmaBands = makeBands(sigmaMode, W, W, Nrep);
end

% Make p(tau|i)
if(~isempty(pTauIInit))
  [I_t Tau_t] = size(pTauIInit);
  if (I_t ~= I+garbageSrc)
    error(['pTauIInit wrong number of sources.  ' ...
           'Expected %d, found %d'], I+garbageSrc, I_t);
  end
  if (Tau_t ~= length(tau))
    error(['pTauIInit wrong number of Taus.  ' ...
           'Expected %d, found %d'], length(tau), Tau_t);
  end
  pTauI = pTauIInit;
  itds = zeros(I+garbageSrc,1);
%  return
%end
else
  if(isempty(tauPosInit))
    % assert: ndims(lr) == 3
    tauPosInit = phatLoc(lr, tau, I, 0, vis);
%     warning('You probably want to use phatLoc to initialize tauPosInit')
%     tauPosInit = pickPeaks(tau, cc, I);
%     if(vis), plot(tau, cc), drawnow, end
    disp(['tauPosInit = ' num2str(tauPosInit)]);
  end

  % ITD measured in ms
  itds = tauPosInit * 1000 / sr;
  if garbageSrc, itds = [itds 0]; tauPosInit = [tauPosInit 0]; end

  tau = repmat(tau, I+garbageSrc, 1);
  pTauI = zeros(size(tau));

  sig_max = 1;      % Biggest bandwidth
  dist_frac = .20;  % set bandwidth to this fraction of distance to
                    % closest source
  
  % Figure out appropriate bandwidths of Gaussians
  dt = abs(diff(tauPosInit));
  sigma = sig_max*ones(size(tauPosInit));
  for i=1:I-1
    if dt(i) < sig_max/dist_frac
      sigma(i) = dt(i) * dist_frac;
      sigma(i+1) = dt(i) * dist_frac;
    end
  end

  if garbageSrc, sigma(I+garbageSrc) = max(abs(tau(i,:)))/.02; end

  % Put Gaussians at the specified locations
  mu = tauPosInit;

  for i=1:I+garbageSrc
    pTauI(i,:) = exp(-1/(2*sigma(i)^2) * (tau(i,:) - mu(i)).^2) + 1e-6;
  end

  % Normalize all of the sources to be equally weighted to begin with
  pTauI = diag(1./sum(pTauI,2)) * pTauI;
  if garbageSrc, pTauI(end,:) = pTauI(end,:) * 1; end
  pTauI = pTauI ./ sum(sum(pTauI));
end

ipdParams = struct('p_tauI', pTauI, 'xi_wit', xi_wit, 's2_wit', ...
    s2_wit, 'xiBands', xiBands, 'sigmaBands', sigmaBands, ...
    'ipdMode', ipdMode, 'xiMode', xiMode, 'sigmaMode', sigmaMode, ...
    'garbageSrc', garbageSrc, 'fixIPriors', fixIPriors);


function lpBin = computeIpdLogLikelihood(W, T, I, Nt, C, rep, ...
    Nrep, ipdParams, angE)

% Erep is a 5th order tensor with size [W T I Nt] and dimensions
% [freq time source delay gaussian]
Erep = repmat(permute(angE, [1 2 4 3]), [1 1 I 1]);

lpBin = single(logProbGmm(Erep, ipdParams.xi_wit, ...
    ipdParams.s2_wit, log(ipdParams.s2_wit)));
clear Erep
lpt = repmat(permute(single(log(ipdParams.p_tauI)), [3 4 1 2]), [W T 1 1]);
lpBin = lpBin + lpt;


function ipdParams = updateIpdParams(W, T, I, Nt, C, rep, ipdParams, ...
    nuIpd, angE);

% Erep is a 5th order tensor with size [W T I Nt] and dimensions
% [freq time source delay gaussian]
Erep = repmat(permute(angE, [1 2 4 3]), [1 1 I 1]);

% start with the easy ones
ipdParams.p_tauI = permute(mean(mean(nuIpd,1),2), [3 4 1 2]);

% Pre-digest nu for denominator of IPD estimation, dims: [W Nt I 1]
nuSum2 = permute(sum(nuIpd, 2), [1 4 3 2]);

nonzero = 1e-12;
% Compute the IPD means as regular gaussian means, not circular
if ipdParams.xiMode
  xiAvg = makeMuAvg(W, ipdParams.xiBands(rep));
  ErepNuSum2 = permute(sum(nuIpd .* Erep,2), [1 4 3 2]);
  for i=1:I
    ipdParams.xi_wit(:,i,:) = (xiAvg * ErepNuSum2(:,:,i) + nonzero) ./ ...
        (xiAvg * nuSum2(:,:,i) + nonzero);
  end
  clear ErepNuSum2
  
  if ipdParams.garbageSrc && 0, ipdParams.xi_wit(:,end,:) = 0; end

  % Sufficient statistics for s2_wit
  xi_rep = repmat(permute(single(ipdParams.xi_wit), [1 4 2 3]), [1 T 1 1]);

  weightedDiffSum2 = sum(nuIpd .* (Erep - xi_rep).^2, 2);
  clear xi_rep  % save memory
else
  % Sufficient statistics for s2_wit    
  weightedDiffSum2 = sum(nuIpd .* Erep.^2, 2);
end
  
% Compute the variances s2_wit
if ipdParams.sigmaMode
  oldS2 = ipdParams.s2_wit;

  sAvg = makeMuAvg(W, ipdParams.sigmaBands(rep));
  weightedDiffSum2 = permute(weightedDiffSum2, [1 4 3 2]);

  for i=1:I
    ipdParams.s2_wit(:,i,:) = ...
        (sAvg * weightedDiffSum2(:,:,i) + nonzero) ./ ...
        (sAvg * nuSum2(:,:,i) + nonzero);
  end

  if ipdParams.garbageSrc && 0
    ipdParams.s2_wit(:,end,:) = oldS2(:,end,:);
  end
else
  ipdParams.s2_wit = permute(sum(sum(weightedDiffSum2, 1),4) ...
      ./ sum(sum(sum(nuIpd,1),2),4), [3 1 2 4]);
  ipdParams.s2_wit = repmat(ipdParams.s2_wit, [1 W Nt]);
  ipdParams.s2_wit = double(permute(ipdParams.s2_wit, [2 1 3]));
end

% Keep parameters of garbage source from changing
if ipdParams.garbageSrc && 0
  ipdParams.xi_wit(:,end,:) = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ILD functions
function ildParams = initIld(I, W, sr, Nrep, ildInit, ildStdInit, ...
                             priorPrec, ildMode, itds, dctMode, ...
                             garbageSrc, B);
ildParams = struct('ildMode', ildMode, 'dctMode', dctMode, ...
    'garbageSrc', garbageSrc, 'B', B);

[rm,cm] = size(ildInit);
[rs,cs] = size(ildStdInit);

% % ITDs measured in ms should go from -0.75 to 0.75.  They're
% % negatively correlated with ILD for some reason.  Maximum ILD is
% % 15 dB, so use that when the ild is at 0.75ms.  The ILD looks like
% % it goes kind of linearly down to 15 dB at 4kHz and then levels
% % out.
% freqs = [1:W]'/(W+1) * sr/2 / 1000;
% pieceWise = min(freqs, 4000);
% ildParams.priorMean = pieceWise/4000 * itds(:)'*(-15/.75);
% %ildParams.priorMean = (1-cos([0:W-1]'*pi/W))/2 * itds(:)'*(-15);

% Generate ILD prior based on regression on anechoic ILDs.
% Frequency regressor is in kHz, ITD regressor is in ms
freqs = [1:W]'/(W+1) * sr/2 / 1000;
ildParams.priorMean = ildPrior(itds, freqs);
ildParams.priorPrec = priorPrec*ones(W,I);

if (numel(ildInit) == 1) && (ildInit == 0) && (priorPrec > 0)
  % Initialize ILD from prior, not for garbage source
  ildParams.mu_wi = ildParams.priorMean;
  ildParams.h2_wi = repmat(ildStdInit.^2, [W I] ./ [rs cs]);
  %ildParams.h2_wi = 1./ildParams.priorPrec(:,1:end-garbageSrc);

%   plot(1:W, ildParams.mu_wi)
%   hold on
%   plot(1:W, ildParams.mu_wi + sqrt(ildParams.h2_wi), ':');
%   plot(1:W, ildParams.mu_wi - sqrt(ildParams.h2_wi), ':');
%   hold off
%   pause

  % Set up fake bands for ILD parameter averaging, should enforce no
  % averaging.
  ildParams.ildBands = makeBands(ildMode, W, W, Nrep);  
  %ildParams.ildBands = makeBands(ildMode, W, max(rm,rs), Nrep);
else
  % Make both mu and eta WxI, error if any dimension is not 1 or the
  % full size
  ildParams.mu_wi = repmat(ildInit, [W I] ./ [rm cm]);
  ildParams.h2_wi = repmat(ildStdInit.^2, [W I] ./ [rs cs]);

  % Set up the bands for ILD parameter averaging.  ildMode > 0 specifies
  % the final number of bands, ildMode < 0 specifies the final number of
  % frequencies per band.
  ildParams.ildBands = makeBands(ildMode, W, max(rm,rs), Nrep);
  %ildParams.ildBands = makeBands(ildMode, W, W, Nrep);
end

if garbageSrc
  ildParams.mu_wi = [ildParams.mu_wi zeros(W,1)];
  %ildParams.h2_wi = [ildParams.h2_wi max(ildStdInit).^2*ones(W,1)];
  ildParams.h2_wi = [ildParams.h2_wi 10.^2*ones(W,1)];
  ildParams.priorPrec = [ildParams.priorPrec zeros(W,1)];
end

Binv = pinv(B);
if dctMode
  ildParams.zeta_di = Binv*ildParams.mu_wi;
end


function lpIld = computeIldLogLikelihood(W,T,I,Nt,C,rep,Nrep, ildParams, A)
lpIld = single(zeros(W,T,I));
for i=1:I
  h2 = repmat(ildParams.h2_wi(:,i), 1, T);
  lpIld(:,:,i) = -.5*log(h2) ...
      - (A-repmat(ildParams.mu_wi(:,i),1,T)).^2 ./ (2*h2);
end


function ildParams = updateIldParams(W, T, I, Nt, C, rep, ildParams, ...
    nuIld, A, Nrep)
% Compute the ILD terms
if ~ildParams.dctMode
%   if rep <= 3 && ildParams.garbageSrc
%     % Wait a couple iterations if using garbage source so that the
%     % IPD can catch on
%     return
%   end

  muAvg = makeMuAvg(W, ildParams.ildBands(rep));
  h2Avg = muAvg;

  % Update means, including prior
  ildParams.mu_wi = ...
      (muAvg * (squeeze(sum(nuIld .* repmat(A, [1 1 I]),2)) + ...
                ildParams.priorMean .* ildParams.priorPrec)) ./ ...
      (muAvg * (squeeze(sum(nuIld,2)) + ildParams.priorPrec));
  if ildParams.garbageSrc && 0, ildParams.mu_wi(:,end) = 0; end

  lastH2 = ildParams.h2_wi(:,end);

  for i=1:I
    ildParams.h2_wi(:,i) = ...
        sum(nuIld(:,:,i) .* (A-repmat(ildParams.mu_wi(:,i),1,T)).^2, 2);
  end
  
  % Add in prior
  ildParams.h2_wi = ildParams.h2_wi + ...
      + ildParams.priorPrec .* (ildParams.priorMean - ildParams.mu_wi).^2;

  ildParams.h2_wi = (h2Avg * ildParams.h2_wi) ./ ...
      (h2Avg * (squeeze(sum(nuIld,2)) + ildParams.priorPrec));

  if ildParams.garbageSrc && 0
    ildParams.h2_wi(:,end) = lastH2;
  end
else
  ndct = ildParams.dctMode;
  for i=1:I
    tmpA = zeros(ndct, ndct);
    tmpb = zeros(ndct, 1);
    for t=1:T
      tmp = diag(nuIld(:,t,i) ./ ildParams.h2_wi(:,i));
      tmpA = tmpA + ildParams.B'*tmp*ildParams.B;
      tmpb = tmpb + ildParams.B'*tmp*A(:,t); 
    end
    ildParams.zeta_di(:,i) = inv(tmpA) * tmpb;
  end
  
  if rep < 0.75*Nrep
    start = 2+fix(0.3*rep);
    ildParams.zeta_di(start:end,:) = 0;
  end
  
  ildParams.mu_wi = real(ildParams.B * ildParams.zeta_di);

  % default to contant variance across frequency
  h2Avg = makeMuAvg(W, 1);
  for i=1:I
    ildParams.h2_wi(:,i) = ...
        sum(nuIld(:,:,i) .* (A-repmat(ildParams.mu_wi(:,i),1,T)).^2, 2);
  end
  ildParams.h2_wi = ...
      (h2Avg * ildParams.h2_wi) ./ (h2Avg * squeeze(sum(nuIld,2)));
  
  if ildParams.ildMode == -2
    for i=1:I
      tmpA = zeros(ndct, ndct);
      tmpb = zeros(ndct, 1);
      for t=1:T
        tmp = diag(nuIld(:,t,i));
        dzm = A(:,t) - ildParams.mu_wi(:,i);
        tmpA = tmpA + ildParams.B'*tmp*ildParams.B;
        tmpb = tmpb + ildParams.B'*tmp*dzm.^2;
      end
      ildParams.varzeta_di(:,i) = inv(tmpA) * tmpb;
    end
    
    if rep < 0.75*Nrep
      start = 2+fix(0.3*rep);
      ildParams.varzeta_di(start:end,:) = 0;
    end
    
    ildParams.h2_wi = real(ildParams.B * ildParams.varzeta_di);
  end
end

% Keep parameters of garbage source from changing
if ildParams.garbageSrc && 0
  ildParams.mu_wi(:,end)    = 0;
end


function ildBands = makeBandsOld(ildMode, W, startBands, Nrep)
% Make a vector that determines the number of bands to use in ild
% each round of EM.  Nonlin is between 0 and 1 and controls the
% amount of non-linearity in the progression of frequencies per
% band.  The closer it is to 1, the more repetitions ildBands stays
% close to startBands and the fewer it stays close to endBands.
if ildMode
  nonlin = 0.99;

  if ildMode >= 0
    endBands = ildMode;
  else
    endBands = -W/ildMode;
  end

  % Never go from more bands to fewer bands
  startBands = min(endBands, startBands);

  ildBands = round(ls10(startBands-nonlin, endBands-nonlin, Nrep)+nonlin);
else
  ildBands = [];
end


function muAvg = makeMuAvg(W, bands)
% Make a matrix that will average the bands of mu_wi when it
% multiplies it.  W is the number of frequencies.  If bands is a
% scalar it specifies the number of bands, if it's a vector it
% specifies the highest freq in each band.

muAvg = eye(W);
if isscalar(bands)
  B = round(linspace(0,W,bands+1));
else
  B = bands;
end

for i=2:length(B)
  muAvg(B(i-1)+1:B(i),B(i-1)+1:B(i)) = 1 / (B(i)-B(i-1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SP functions
function [spParams C] = initSp(I, W, L, R, sourcePriors, stdInit, ...
    dctMode, spMode, garbageSrc, B)
spParams = struct('spMode', spMode, 'dctMode', dctMode, ...
    'garbageSrc', garbageSrc, 'B', B);

% Set spMode based on spStartRep.  spMode will remain turned off
% until rep == spStartRep.
if spMode
  % FIXME
  %spParams.spStartRep = 10;
  %spParams.origSpMode = spMode;
  %spParams.spMode = 0;

  %spParams.spStartRep = 2;
  spParams.spStartRep = 4;
  spParams.origSpMode = spMode;
  spParams.spMode = 0;
else
  spParams.spStartRep = Inf;
end

C = 0;
% Set up source prior parameters
if spMode
  if isfield(sourcePriors, 'ev_params')
    spParams.ev_params = sourcePriors.ev_params;
    spParams.gmm = sourcePriors.gmm;
    sourcePriors = spParams.gmm;
  end

  if length(sourcePriors) == 1
    tmp = sourcePriors;
    for i = 1:I
      sourcePriors(i) = tmp;
    end
    clear tmp
  elseif length(sourcePriors) ~= I
    error('Length of sourcePriors needs to be either 1 or I');
  end
  
  C = max([sourcePriors.nmix]);
  % Make sure all GMMs have the same number of components
  for i=1:I
    if sourcePriors(i).nmix < C
      sourcePriors(i).priors = ...
          [sourcePriors(i).priors repmat(-Inf, [1 C-sourcePriors(i).nmix])];
      sourcePriors(i).means = ...
          [sourcePriors(i).means zeros(W, C-sourcePriors(i).nmix)];
      sourcePriors(i).covars = ...
          [sourcePriors(i).covars ones(W, C-sourcePriors(i).nmix)];      
      sourcePriors(i).nmix = C;
    end
  end
  
  if garbageSrc
    garbageSourcePrior.nmix = C;
    garbageSourcePrior.priors = [0 repmat(-Inf, [1 C-1])];
    garbageSourcePrior.means = zeros(W, C);
    garbageSourcePrior.covars = ones(W, C);
    for i=1:I
      gmm = sourcePriors(i);
      for c=1:C
        p = exp(gmm.priors(c))/I;
        garbageSourcePrior.means(:,1) = garbageSourcePrior.means(:,1) ...
            + p*gmm.means(:,c);
        garbageSourcePrior.covars(:,1) = garbageSourcePrior.covars(:,1) ...
            + p*gmm.covars(:,c);
      end
    end
    sourcePriors(I+1) = garbageSourcePrior;
  end
  spParams.sourcePriors = sourcePriors;

  I = I + garbageSrc;
  switch spMode
    case 1
      spParams.channel_response = zeros(2,I,W);
    case {-1,-2}
      spParams.channel_response = zeros(I,W);
      spParams.channel_response_var = repmat(eps, [I W]);
 end
end


function pgmm = computeSpLogLikelihood(W,T,I,Nt,C,rep,Nrep, spParams, ...
    ildParams, L, R)
% Compute the likelihood of each time frequency cell of the left and
% right channels under each component of the given GMMs.
pgmm = single(zeros(2, W, T, I, C));

spMode = spParams.spMode;
if false && spMode < 0 && ildParams.ildMode && rep == Nrep
  % Differentiate L and R masks
  spMode = abs(spMode);
  tmp_chan_resp = spParams.channel_response;
  spParams.channel_response(1,:,:) = 0.5*(tmp_chan_resp + ildParams.mu_wi');
  spParams.channel_response(2,:,:) = 0.5*(tmp_chan_resp - ildParams.mu_wi');
end

if spMode > 0
  for i=1:I
    gmm = spParams.sourcePriors(i);
    for c = 1:C
      mmu_L = repmat(gmm.means(:,c) ...
          + squeeze(spParams.channel_response(1,i,:)), [1 T]);
      mcv_L = repmat(sqrt(gmm.covars(:,c)), [1 T]);
      mmu_R = repmat(gmm.means(:,c) ...
          + squeeze(spParams.channel_response(2,i,:)), [1 T]);
      mcv_R = repmat(sqrt(gmm.covars(:,c)), [1 T]);
      
      % we assume the priors are the same between ears, because
      % both ears must be in the same component of the GMM

      pl = normpdf(L, mmu_L, mcv_L);
      pr = normpdf(R, mmu_R, mcv_R);
      prior = exp(gmm.priors(c));
      % FIXME - do we really want separate masks?
      %pgmm(1,:,:,i,c) = prior * pl;
      %pgmm(2,:,:,i,c) = prior * pr;
      p = prior * (pl .* pr);
      pgmm(1,:,:,i,c) = p;
      pgmm(2,:,:,i,c) = p;
    end
  end
elseif spMode < 0
  obs = L + R;
  for i=1:I
    gmm = spParams.sourcePriors(i);
    gmm.means = 2*gmm.means + ...
        repmat(squeeze(spParams.channel_response(i,:))', [1 C]);
    gmm.covars = 4*gmm.covars + ...
        repmat(squeeze(spParams.channel_response_var(i,:))', [1 C]);
    priors = exp(gmm.priors);
    for c=1:C
      mmu = repmat(gmm.means(:,c), [1 T]);
      mcv = repmat(sqrt(gmm.covars(:,c)), [1 T]);
      pgmm(1,:,:,i,c) = normpdf(obs, mmu, mcv) * priors(c);
      pgmm(2,:,:,i,c) = pgmm(1,:,:,i,c);
    end
  end
end

pgmm = log(pgmm + 1e-45);

% Don't use SP for lowest frequencies.
%FIXME
pgmm(:,1:32,:,:,:) = 0;



function gmms = permuteGmms(gmms, L, R, bin_mask, I)
% Fix arbitrary permutations in the ordering of gmms so it agrees
% with bin_mask.

obs = logsum(cat(3, L, R), 3);
T = size(obs, 2);

gmlik = zeros(I);
for i = 1:I
  mu = gmms(1,i).means;
  cv = sqrt(gmms(1,i).covars);
  priors = gmms(1,i).priors(:);
  for j = 1:I
    mask = squeeze(bin_mask(:,:,j));
    
    for c = 1:gmms(1,i).nmix
      muc = repmat(mu(:,c), [1 T]);
      cvc = repmat(cv(:,c), [1 T]);

      % Soft mask missing data likelihood:
      lpr(c,:) = sum(log(mask.*normpdf(obs, muc, cvc) + (1-mask)), 1) ...
          + priors(c);
    end
    
    gmlik(i,j) = sum(logsum(lpr, 1), 2);
  end
end

ngmm = 0;
while ngmm < I
  [tmp ind] = max(gmlik(:));
  [i j] = ind2sub(size(gmlik), ind); 
  idx(i) = j;
  ngmm = ngmm + 1;

  % Make sure we don't pick the same gmm or mask again
  gmlik(i,:) = -Inf;
  gmlik(:,j) = -Inf;
end
gmms = gmms(idx);


function channel_response = channelAdaptGmms(initial_gmm, obs, p_wtc, B)
% Update left and right channel responses for each source GMM.
% Returns the updated GMM structure and the magnitude response of the
% channel.

[W T] = size(obs);
ndct = size(B, 2);
% posterior weighted projection of observation onto DCT bases
op = zeros(ndct, 1);
% posterior weighted projection of DCT bases onto themselves
Bp = zeros(ndct);
mu = initial_gmm.means;
for t = 1:T
  % + 1e-200 to avoid dividing by zero
  p_wc = double(squeeze(p_wtc(:,t,:))) + 1e-200;
  pcv = initial_gmm.covars ./ p_wc;
  
  tobs = repmat(obs(:,t), [1 initial_gmm.nmix]);
  opt = B' * sum((tobs - mu)./pcv, 2);
  Bpt = B' * diag(sum(1./pcv, 2)) * B;
    
  op = op + opt;
  Bp = Bp + Bpt;
end
chan = inv(Bp)*op;
channel_response = B*chan;



function [w hl hr] = eigenvoiceAndChannelAdaptGmms(spParams, L, R, p_lr_wtc)
[W T] = size(L);
C = size(spParams.ev_params.mean, 2);
ndct = size(spParams.B, 2);
nev = length(spParams.ev_params.means);

icvl = zeros(W, C);
icvr = zeros(W, C);
icvlobs = zeros(W, C);
icvrobs = zeros(W, C);
for t = 1:T
  Sl = double(squeeze(p_lr_wtc(1,:,t,:))) ./ spParams.gmm.covars;
  Sr = double(squeeze(p_lr_wtc(2,:,t,:))) ./ spParams.gmm.covars;
  icvl = icvl + Sl;
  icvr = icvr + Sr;
  icvlobs = icvlobs + diag(sparse(L(:,t))) * Sl;
  icvrobs = icvrobs + diag(sparse(R(:,t))) * Sr;
end

Aww = zeros(nev, nev);
Awl = zeros(nev, ndct);
Awr = zeros(nev, ndct);
All = zeros(ndct, ndct);
Arr = zeros(ndct, ndct);
bw = zeros(nev, 1);
bl = zeros(ndct, 1);
br = zeros(ndct, 1);

u = spParams.ev_params.mean;
Uct = zeros(nev, W);
B = spParams.B;
for c = 1:C
  for j = 1:nev
    Uct(j,:) = spParams.ev_params.means{j}(:,c);
  end

  Sl = diag(sparse(icvl(:,c)));
  Sr = diag(sparse(icvr(:,c)));
  BSl = B' * Sl;
  BSr = B' * Sr;
  Slobs = icvlobs(:,c);
  Srobs = icvrobs(:,c);

  Aww = Aww + Uct * (Sl + Sr) * Uct';
  Awl = Awl + Uct * BSl';
  Awr = Awr + Uct * BSr';
  All = All + BSl * B;
  Arr = Arr + BSr * B;
  bw = bw + Uct * (Slobs - Sl * u(:,c) +  Srobs - Sr * u(:,c));
  bl = bl + B' * (Slobs - Sl * u(:,c));
  br = br + B' * (Srobs - Sr * u(:,c));
end
z = zeros(ndct, ndct);
tmp = inv([Aww Awl Awr; Awl' All z; Awr' z Arr]) * [bw; bl; br];

w = tmp(1:nev);
hl = tmp(nev+[1:ndct]);
hr = tmp(nev+ndct+[1:ndct]);



function spParams = updateSpParams(W, T, I, Nt, C, rep, spParams, ...
    nuSp, L, R, Nrep);
for i = 1:I
  if spParams.garbageSrc && i == I
    continue
  elseif spParams.spMode > 0
    if isfield(spParams, 'ev_params')
      [w hl hr] = eigenvoiceAndChannelAdaptGmms(spParams, L, R, ...
          squeeze(nuSp(:,:,:,i,:)));
      spParams.w(i,:) = w;
      spParams.hl(i,:) = hl;
      spParams.hr(i,:) = hr;
      spParams.channel_response(1,i,:) = spParams.B * hl;
      spParams.channel_response(2,i,:) = spParams.B * hr;
      spParams.sourcePriors(i) = ...
          construct_adapted_model(spParams.gmm, spParams.ev_params, w);
    else
      spParams.channel_response(1,i,:) = channelAdaptGmms( ...
	  spParams.sourcePriors(i), L, squeeze(nuSp(1,:,:,i,:)), spParams.B);
      spParams.channel_response(2,i,:) = channelAdaptGmms(  ...
	  spParams.sourcePriors(i), R, squeeze(nuSp(2,:,:,i,:)), spParams.B);
    end
  elseif spParams.spMode < 0
    obs = L + R;

    p_wtc = squeeze(mean(nuSp(:,:,:,i,:), 1));
    gmm = spParams.sourcePriors(i);
    ndct = size(spParams.B, 2);
    tmpA = zeros(ndct, ndct);
    tmpb = zeros(ndct, 1);
    for t=1:T
      tmp = squeeze(p_wtc(:,t,:)) ./ (4*gmm.covars + ...
          repmat(spParams.channel_response_var(i,:)', [1 C]));
      tobs = repmat(obs(:,t), [1 C]);
      tmpA = tmpA + spParams.B'*diag(sum(tmp, 2))*spParams.B;
      tmpb = tmpb + spParams.B'*sum(tmp.*(tobs - 2*gmm.means), 2); 
    end
    dct_channel_response_i = inv(tmpA) * tmpb;
    
    %if rep < 0.66*Nrep
    %  start = 3+fix(0.5*(rep-1));
    %  dct_channel_response_i(start:end,:) = 0;
    %end
    
    spParams.dct_channel_response(i,:) = dct_channel_response_i;
    spParams.channel_response(i,:) = spParams.B * dct_channel_response_i;

    % FIXME - the update here is incorrect.  It leads to negative
    % variances (b/c we just subtract off 4*gmm.covars)
    if spParams.spMode == -2
      warning(['messl: the SP-ILS variance update is known to be' ...
            ' incorrect.  Don''t expect this to work properly.'])
      % variance
      tmpA = zeros(ndct, ndct);
      tmpb = zeros(ndct, 1);
      for t=1:T
        tmp = squeeze(p_wtc(:,t,:));
        dzm = repmat(obs(:,t), [1 C]) - 2*gmm.means ...
            - repmat(squeeze(spParams.channel_response(i,:))', [1 C]);
        tmpA = tmpA + spParams.B'*diag(sum(tmp, 2))*spParams.B;
        tmpb = tmpb + spParams.B'*sum(tmp.*(dzm.^2 - 4*gmm.covars), 2); 
      end
      dct_channel_response_var_i = inv(tmpA) * tmpb;
    
      %if rep < 0.66*Nrep
      %  start = 3+fix(0.5*(rep-1));
      %  dct_channel_response_var_i(start:end,:) = 0;
      %end
    
      spParams.dct_channel_response_var(i,:) = dct_channel_response_var_i;
      spParams.channel_response_var(i,:) = spParams.B * dct_channel_response_var_i;
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MRF Functions
function mrfCompatPot = loadMrfCompat(mrfCompatFile, I, garbageSrc)
% Hard coded for 2 sources and a garbage source for now...

if isempty(mrfCompatFile) || ~exist(mrfCompatFile, 'file') || (I ~= 2)
    mrfCompatPot = [];
else
    load(mrfCompatFile);
    mrfCompatPot = counts;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ll p_lr_iwt nuIpd maskIpd nuIld maskIld nuSp maskSp] = ...
    computePosterior(W, T, I, Nt, C, logMaskPrior, ...
    ipdParams, lpIpd, ildParams, lpIld, spParams, lpSp, vis, reliability, ...
    mrfCompatPot, mrfCompatExp)
% Defaults
nuIpd = single(lpIpd); maskIpd = 0;
nuIld = single(lpIld); maskIld = 0;
nuSp  = single(lpSp);  maskSp  = 0;

mrfLbpIter = 8;

if vis || ~spParams.spMode
  % Normalize each term separated to demonstrate the contribution of
  % each component to the overall posterior mask.
  if ipdParams.ipdMode
    maskIpd = sum(exp(lpIpd), 4);
    maskIpd = maskIpd ./ repmat(sum(maskIpd,3), [1 1 I]);
  end
  if ildParams.ildMode
    maskIld = exp(lpIld);
    maskIld = maskIld ./ repmat(sum(maskIld,3), [1 1 I]);
  end
  if spParams.spMode
    maskSp = sum(squeeze(mean(exp(lpSp), 1)), 4);
    maskSp = maskSp ./ repmat(sum(maskSp,3), [1 1 I]);
  end
end

% This actually normalizes the joint distribution of i,tau,j,c not
% the conditionals
if ipdParams.ipdMode
  pBin = lpIpd;
  clear lpIpd
  if ~isempty(logMaskPrior)
    pBin = pBin + repmat(logMaskPrior, [1 1 1 Nt]);
  end
  if ildParams.ildMode
    pBin = pBin + repmat(lpIld, [1 1 1 Nt]);
  end
  clear lpIld

%   ls = linspace(-30,0,1000);
%   figure(5), bar(ls, histc(pBin(:), ls)); drawnow

  pBin = exp(pBin) + eps;

%   % Only use high relative probability points
%   threshold = quantile(pBin(:), 0);

  if any(pBin(:) <= 0), warning('pBin <= 0'), end
end

if ~spParams.spMode
  norm = sum(sum(pBin,3),4);
  if any(norm(:) <= 0), warning('norm <= 0 (1)'), end
  nuIpd = pBin ./ repmat(norm, [1 1 I Nt]);
  %nuIpd = (pBin .* (pBin > threshold) + eps) ./ repmat(norm, [1 1 I Nt]);
  nuIld = sum(nuIpd, 4);
  %p_lr_iwt = permute(repmat(nuIld, [1 1 1 2]), [4 1 2 3]);
  p_lr_iwt = repmat(permute(nuIld, [4 1 2 3]), [2 1 1 1]);
  if any(norm(:) <= 0), warning('norm <= 0 (2)'), end
  ll = sum(sum(log(norm)));
else
  if ~ipdParams.ipdMode && ~ildParams.ildMode
    if ~isempty(logMaskPrior)
      lpSp = lpSp + permute(repmat(logMaskPrior, [1 1 1 2 C]), [4 1 2 3 5]);
    end
    nuSp = exp(lpSp);
    normSp = sum(sum(nuSp,4),5);
    nuSp = nuSp ./ repmat(normSp, [1 1 1 I C]);
    norm = squeeze(mean(normSp, 1));
    p_lr_iwt = sum(nuSp, 5);
    ll = sum(sum(log(norm)));
  else

    if true
%    if false
    lpSp = exp(lpSp);
    
    normSp = single(zeros(2, W, T));
    increment = max(1, fix(W/C));
    %if C > 32, increment = 1; end
    curr_s = 1;
    ll = 0;
    while curr_s < W
      try 
        for s = curr_s:increment:W
          curr_s = s;
          w = s:min(s+increment-1, W);
          pw = repmat(pBin(w,:,:,:), [1 1 1 1 C 2]) ...
              .* repmat(permute(lpSp(:,w,:,:,:), [2 3 4 6 5 1]), [1 1 1 Nt 1 1]);

          normSp(:,w,:) = permute(sum(sum(sum(pw,3),4),5), [6 1 2 3 4 5]);
          lpSp(:,w,:,:,:) = permute(sum(pw,4), [6 1 2 3 5 4]);
          nuIpd(w,:,:,:) = mean(sum(pw,5),6);
          % this doesn't work for some reason and uses more memory?
          % tmp = permute(sum(sum(sum(pw,3),4),5), [6 1 2 3 4 5]);
          % nuSp(:,w,:,:,:) = permute(sum(pw,4), [6 1 2 3 5 4]) ...
          %     ./ repmat(tmp, [1 1 1 I C]);
          % %tmp = squeeze(geomean(tmp, 1));
          % tmp = squeeze(mean(tmp, 1));
          % nuIpd(w,:,:,:) = mean(sum(pw,5),6) ...
          %     ./ repmat(tmp, [1 1 I Nt]);
          % ll = ll + sum(log(tmp(:)));
        end
        curr_s = W;
      catch
        increment = max(1, fix(increment/2));
        warning('Subband %d: backing off to %d.', curr_s, increment);
      end
    end
    %norm = squeeze(geomean(normSp, 1));
    norm = squeeze(mean(normSp, 1));
    ll = sum(sum(log(norm)));
    nuIpd = nuIpd ./ repmat(norm, [1 1 I Nt]);
    nuSp = lpSp ./ repmat(normSp, [1 1 1 I C]);
    nuIld = sum(nuIpd, 4);
    p_lr_iwt = sum(nuSp, 5);
    else
      % save memory by not storing a temporary variable.
      lpSp = exp(lpSp);
      [nuIpd nuSp normSp] = normalizeBinAndSp(pBin, lpSp);
      %norm = squeeze(geomean(normSp, 1));
      norm = squeeze(mean(normSp, 1)); 
      ll = sum(sum(log(norm)));
      nuIpd = nuIpd ./ repmat(norm, [1 1 I Nt]);
      nuSp = lpSp ./ repmat(normSp, [1 1 1 I C]);
      nuIld = sum(nuIpd, 4);
      p_lr_iwt = sum(nuSp, 5);
    end
%    dbstop
  end
end

[p_lr_iwt nuIld nuIpd] = applyMrf(nuIld, nuIpd, p_lr_iwt, mrfCompatPot, mrfCompatExp, mrfLbpIter, 'sum');

if ~isempty(reliability)
  nuIpd = nuIpd .* repmat(reliability, [1 1 I Nt]);
  nuIld = nuIld .* repmat(reliability, [1 1 I]);
  nuSp  = nuSp  .* repmat(permute(reliability, [3 1 2]), [2 1 1 I C]);
end


function [p_lr_iwt nuIld nuIpd hardSrcs] = applyMrf(nuIld, nuIpd, p_lr_iwt, mrfCompatPot, mrfCompatExp, mrfLbpIter, bpType, doPlot)
if ~exist('doPlot', 'var') || isempty(doPlot), doPlot = 0; end

if ~isempty(mrfCompatPot) && (mrfCompatExp ~= 0)
    [newNuIld hardSrcs] = mrfGridLbp(nuIld, mrfCompatPot.^mrfCompatExp, mrfLbpIter, bpType, doPlot);
    nuIpd = bsxfun(@times, nuIpd, newNuIld ./ nuIld);
    nuIld = newNuIld;
    p_lr_iwt = repmat(permute(nuIld, [4 1 2 3]), [2 1 1 1]);
else
    [~,hardSrcs] = max(squeeze(p_lr_iwt(1,:,:,:)), [], 3);
end


function nuIpd = enforceIPriors(nuIpd, ipdParams)
if ipdParams.fixIPriors
    newSrcPriors = mean(mean(sum(nuIpd, 4), 1), 2);
    rescaling = permute(sum(ipdParams.p_tauI, 2), [2 3 1]) ./ newSrcPriors;
    nuIpd = bsxfun(@times, rescaling, nuIpd);
    assert(all(abs(squeeze(mean(mean(sum(nuIpd,4), 1), 2)) - squeeze(sum(ipdParams.p_tauI, 2))) < 1e-6));
end


function visualizeParams(W, T, I, tau, sr, ipdParams, ildParams, spParams, ...
    p_lr_iwt, maskIpd, maskIld, maskSp, L, R, reliability)
%%%% Visualize what's happening
% 4 figures:
%   1. Factored masks
%   2. Recovered signals (i.e. mask .* signal)
%   3. Learned model parameters
%   4. IPD vs frequency distribution

sp_cols = max(I,2);

modes = [ipdParams.ipdMode, ildParams.ildMode, spParams.spMode];
params_on = sum(modes ~= 0);
names = {'IPD', 'ILD', 'SP'};
masks = {maskIpd, maskIld, maskSp};
obs = idB(L) + idB(R);

% Factored masks
%
% IPD masks
% ILD masks
% SP masks
% total masks (if more than one other is present)
mask_plots = {};
mask_titles = {};
mask_layout = [params_on+(params_on>1) I];


% Recovered signals (i.e. mask .* signal)
%
% IPD reconstructions
% ILD reconstructions
% SP reconstructions
% total reconstructions
% original mixture
recovered_plots = {};
recovered_titles = {};
recovered_layout = [params_on+1+(params_on>1) I];

for m=1:length(modes)
  if(modes(m))
    for i=1:I
      mask_plots = {mask_plots{:}, masks{m}(:,:,i)};
      t = sprintf('%s probability of being in source %d', names{m}, i);
      mask_titles = {mask_titles{:}, t};
      
      recovered_plots = {recovered_plots{:}, dB(masks{m}(:,:,i) .* obs)};
      t = sprintf('%s reconstruction of source %d', names{m}, i);
      recovered_titles = {recovered_titles{:}, t};
    end
  end
end

if params_on > 1
  for i=1:I
    mask_plots = {mask_plots{:}, squeeze(mean(p_lr_iwt(:,:,:,i), 1))};
    t = sprintf('Overall probability of being in source %d',i);
    mask_titles = {mask_titles{:}, t};
    
    recovered_plots = {recovered_plots{:}, ...
                       dB(squeeze(mean(p_lr_iwt(:,:,:,i),1)) .* obs)};
    t = sprintf('Overall reconstruction of source %d', i);
    recovered_titles = {recovered_titles{:}, t};
  end
end

recovered_plots = {recovered_plots{:}, dB(obs)};
recovered_titles = {recovered_titles{:}, 'Mixed signal'};

if ~isempty(reliability)
  recovered_plots = {recovered_plots{:}, 50*(reliability-1)};
  recovered_titles = {recovered_titles{:}, 'Reliability'};
end

% plotall(mask_plots, 'title', mask_titles, 'figure', 1, 'subplot', ...
%         mask_layout, 'CLim', [0 1]);
% 
% plotall(recovered_plots, 'title', recovered_titles, 'figure', 2, ...
%         'subplot', recovered_layout, 'CLim', [-60 0]);

fig_no_focus(1)
subplots(mask_plots, mask_layout, mask_titles, @(r,c,i) caxis([0 1]));
drawnow

fig_no_focus(2)
subplots(recovered_plots, recovered_layout, recovered_titles, @(r,c,i) caxis([-60 0]));
drawnow

% imgsc(mask_plots, 'title', mask_titles, 'figure', 1, 'subplot', ...
%       mask_layout, 'caxis', [0 1]);

% imgsc(recovered_plots, 'title', recovered_titles, 'figure', 2, ...
%       'subplot', recovered_layout, 'caxis', [-60 0]);


% Learned model parameters
%
% p(tau | i), ILD means and stds
% SP parameters
fig_no_focus(3)

if modes(3) == 1
  sp_rows = 2;
  sp_cols = I;
else
  sp_rows = 1;
  sp_cols = sum(modes ~= 0);
end

subplot(sp_rows,sp_cols,1)
plot(tau, ipdParams.p_tauI')
title('p(tau, i)');
grid on

if ildParams.ildMode
  subplot(sp_rows,sp_cols,2)
  plot(1:W, ildParams.mu_wi)
  hold on
  plot(1:W, ildParams.mu_wi + sqrt(ildParams.h2_wi), ':');
  plot(1:W, ildParams.mu_wi - sqrt(ildParams.h2_wi), ':');
  hold off
  grid on
  title('ILD');
end

if spParams.spMode == 1
  for i = 1:I
    subplot(sp_rows,sp_cols,sp_cols+i)
    plot(squeeze(spParams.channel_response(:,i,:))')
    if i == 1
      ax = axis();
    else
      axis(ax);
    end
    title(sprintf('Room+HRTF channel for source %d', i))
    grid on;
  end
  if isfield(spParams, 'ev_params')
    subplot(sp_rows,sp_cols,3)
    plot(spParams.w')
    title('Eigenvoice parameters')
    grid on;
  end
elseif spParams.spMode
  subplot(sp_rows,sp_cols,3)
  mu = squeeze(spParams.channel_response');
  plot(1:W, mu)
  if spParams.spMode == -2
    hold on
    plot(1:W, mu + sqrt(spParams.channel_response_var'), ':');
    plot(1:W, mu - sqrt(spParams.channel_response_var'), ':');
    hold off
  end
  grid on
  title('L+R channel')
end
subplot 111
  

% IPD vs Frequency distribution
%
fig_no_focus(4);
visParams(ipdParams, tau, sr);


drawnow
