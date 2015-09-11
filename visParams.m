function [prob P1 w] = visParams(params, tau, sr, cb, plots, T)
%
% visParams(params, tau, sr, cb, plots)
%
% Visualize the parameters of up to three sources in IPD space.
% Params is the structure that messl produces.  If cb is true or is
% not supplied, then include a colorbar on all plots.  Plots is a
% vector of indices indicating which parameters to plot.

if ~exist('sr', 'var') || isempty(sr), sr = pi; end
if ~exist('cb', 'var') || isempty(plots), cb = 1;  end
if ~exist('plots', 'var'), plots = 1:size(params.xi_wit,2); end
if ~exist('T', 'var') || isempty(T), T = 100; end

xi_wit = params.xi_wit;
p_tauI = params.p_tauI;
if isfield(params, 's_wit')
  s2_wit = params.s_wit.^2;
else
  s2_wit = params.s2_wit;
end

W  = size(xi_wit, 1);
I  = size(xi_wit, 2);
P1 = exp(1j*linspace(-pi,pi,T+1));
P = repmat(P1(1:end-1), W, 1);

[E,w] = probCC(cat(3,P,ones(size(P))), tau);
w     = w * sr/(2*pi*1000);
Erep  = repmat(permute(single(angle(E)), [1 2 4 3]), [1 1 I 1]);
clear E  % Save memory

lp   = logProbGmm(Erep, xi_wit, s2_wit, log(s2_wit));
lpt  = permute(log(p_tauI), [3 4 1 2]);
lp   = bsxfun(@plus, lp, lpt);
prob = sum(exp(lp), 4);
prob = permute(prob, [2 1 3]);

for i=1:length(plots)
  subplot(length(plots),1,i)
%   subplot(I,2,2*(i-1)+1)
  cmax = quantile(reshape(prob(:,:,plots(i)),[],1), .98);
  imagesc(w, angle(P1), prob(:,:,plots(i)), [0 cmax]), axis xy
  xlabel('Frequency (kHz)')
  ylabel('IPD (rad)')
  if cb, colorbar, end

%   subplot(I,2,2*i)
%   imagesc(w, angle(P1), log(prob(:,:,i))-log(sum(prob,3)-prob(:,:,i))), axis xy
%   xlabel('Frequency (kHz)')
%   ylabel('IPD (rad)')
%   caxis([-2 2])
%   if cb, colorbar, end

%   hold on
%   contour(repmat(w',T,1), -angle(P'), ...
%           2*prob(:,:,i)-sum(prob,3), [0 0], 'w')
%   hold off

end

subplot 111
