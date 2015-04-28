function nuIpd = messlIpdEnforcePriors(nuIpd, ipdParams)
if ipdParams.fixIPriors
    newSrcPriors = mean(mean(sum(nuIpd, 4), 1), 2);
    rescaling = permute(sum(ipdParams.p_tauI, 2), [2 3 1]) ./ newSrcPriors;
    nuIpd = bsxfun(@times, rescaling, nuIpd);
    assert(all(abs(squeeze(mean(mean(sum(nuIpd,4), 1), 2)) - squeeze(sum(ipdParams.p_tauI, 2))) < 1e-6));
end
