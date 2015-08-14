function gmms = messlSpPermuteGmms(gmms, L, R, bin_mask, I)
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
