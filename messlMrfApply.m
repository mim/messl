function [p_lr_iwt nuIld nuIpd hardSrcs] = messlMrfApply(nuIld, nuIpd, p_lr_iwt, mrfCompatPot, mrfCompatExp, mrfLbpIter, bpType, doPlot)
if ~exist('doPlot', 'var') || isempty(doPlot), doPlot = 0; end

if ~isempty(mrfCompatPot) && (mrfCompatExp ~= 0)
    [newNuIld hardSrcs] = mrfGridLbp(nuIld, mrfCompatPot.^mrfCompatExp, mrfLbpIter, bpType, doPlot);
    nuIpd = bsxfun(@times, nuIpd, newNuIld ./ nuIld);
    nuIld = newNuIld;
    p_lr_iwt = repmat(permute(nuIld, [4 1 2 3]), [2 1 1 1]);
else
    [~,hardSrcs] = max(squeeze(p_lr_iwt(1,:,:,:)), [], 3);
end
