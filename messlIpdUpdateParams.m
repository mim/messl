function ipdParams = messlIpdUpdateParams(W, T, I, Nt, C, rep, ipdParams, ...
    nuIpd, angE)

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
  xiAvg = messlUtilMakeMuAvg(W, ipdParams.xiBands(rep));
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

  sAvg = messlUtilMakeMuAvg(W, ipdParams.sigmaBands(rep));
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
