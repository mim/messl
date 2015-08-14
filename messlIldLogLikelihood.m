function lpIld = messlIldLogLikelihood(W,T,I,Nt,C,rep,Nrep, ildParams, A)
lpIld = single(zeros(W,T,I));
for i=1:I
  h2 = repmat(ildParams.h2_wi(:,i), 1, T);
  lpIld(:,:,i) = -.5*log(h2) ...
      - (A-repmat(ildParams.mu_wi(:,i),1,T)).^2 ./ (2*h2);
end
