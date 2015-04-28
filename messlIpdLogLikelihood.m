function lpBin = messlIpdLogLikelihood(W, T, I, Nt, C, rep, ...
    Nrep, ipdParams, angE)

% Erep is a 5th order tensor with size [W T I Nt] and dimensions
% [freq time source delay gaussian]
Erep = repmat(permute(angE, [1 2 4 3]), [1 1 I 1]);

lpBin = single(logProbGmm(Erep, ipdParams.xi_wit, ...
    ipdParams.s2_wit, log(ipdParams.s2_wit)));
clear Erep
lpt = repmat(permute(single(log(ipdParams.p_tauI)), [3 4 1 2]), [W T 1 1]);
lpBin = lpBin + lpt;
