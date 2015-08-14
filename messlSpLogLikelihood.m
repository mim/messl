function pgmm = messlSpLogLikelihood(W,T,I,Nt,C,rep,Nrep, spParams, ...
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
