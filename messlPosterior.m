function [ll p_lr_iwt nuIpd maskIpd nuIld maskIld nuSp maskSp likeIwt] = ...
    messlPosterior(W, T, I, Nt, C, logMaskPrior, ...
    ipdMode, lpIpd, ildMode, lpIld, spMode, lpSp, vis, reliability, ...
    mrfCompatPot, mrfCompatExp, mrfLbpIter)
% Defaults
nuIpd = single(lpIpd); maskIpd = 0;
nuIld = single(lpIld); maskIld = 0;
nuSp  = single(lpSp);  maskSp  = 0;

if vis || ~spMode
  % Normalize each term separated to demonstrate the contribution of
  % each component to the overall posterior mask.
  if ipdMode
    maskIpd = sum(approxExp(lpIpd), 4);
    maskIpd = maskIpd ./ repmat(sum(maskIpd,3), [1 1 I]);
  end
  if ildMode
    maskIld = exp(lpIld);
    maskIld = maskIld ./ repmat(sum(maskIld,3), [1 1 I]);
  end
  if spMode
    maskSp = sum(squeeze(mean(exp(lpSp), 1)), 4);
    maskSp = maskSp ./ repmat(sum(maskSp,3), [1 1 I]);
  end
end

% This actually normalizes the joint distribution of i,tau,j,c not
% the conditionals
pBin = 0;
if ipdMode
  pBin = pBin + lpIpd;
  clear lpIpd
end

if ~isempty(logMaskPrior)
    pBin = bsxfun(@plus, pBin, logMaskPrior);
end

if ildMode
    pBin = pBin + repmat(lpIld, [1 1 1 Nt]);
end
clear lpIld

%   ls = linspace(-30,0,1000);
%   figure(5), bar(ls, histc(pBin(:), ls)); drawnow

pBin = approxExp(pBin) + eps;

%   % Only use high relative probability points
%   threshold = quantile(pBin(:), 0);

if any(pBin(:) <= 0), warning('pBin <= 0'), end

likeIwt = sum(pBin, 4);

if ~spMode
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
  if ~ipdMode && ~ildMode
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

[p_lr_iwt nuIld nuIpd] = messlMrfApply(nuIld, nuIpd, p_lr_iwt, mrfCompatPot, mrfCompatExp, mrfLbpIter, 'sum');

if ~isempty(reliability)
  nuIpd = nuIpd .* repmat(reliability, [1 1 I Nt]);
  nuIld = nuIld .* repmat(reliability, [1 1 I]);
  nuSp  = nuSp  .* repmat(permute(reliability, [3 1 2]), [2 1 1 I C]);
end


function Y = approxExp(X, thresh)
if nargin < 2, thresh = log(eps); end

Y = exp(X);
return

keep = X > thresh;
if mean(keep(:)) > 0.05  % seems to need about 5% non-zero sparsity for this to be worth it.
    Y = exp(X);
else
    % Only compute exp on large-enough entries
    Y = zeros(size(X), 'single');
    Y(keep) = exp(X(keep));
end
