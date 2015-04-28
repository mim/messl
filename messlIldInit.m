function ildParams = messlIldInit(I, W, sr, Nrep, ildInit, ildStdInit, ...
                             priorPrec, ildMode, itds, dctMode, ...
                             garbageSrc, B)
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
