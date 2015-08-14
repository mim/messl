function ildParams = messlIldUpdateParams(W, T, I, Nt, C, rep, ildParams, ...
    nuIld, A, Nrep)
% Compute the ILD terms

minVar = 1e-8;

if ~ildParams.dctMode
%   if rep <= 3 && ildParams.garbageSrc
%     % Wait a couple iterations if using garbage source so that the
%     % IPD can catch on
%     return
%   end

  muAvg = messlUtilMakeMuAvg(W, ildParams.ildBands(rep));
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

  % Prevent h2 going to 0
  ildParams.h2_wi = max(ildParams.h2_wi, minVar);
  
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
  h2Avg = messlUtilMakeMuAvg(W, 1);
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
