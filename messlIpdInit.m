function [ipdParams, itds] = ...
    messlIpdInit(I, W, Nrep, tau, sr, lr, tauPosInit, pTauIInit, sigmaInit, ...
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
