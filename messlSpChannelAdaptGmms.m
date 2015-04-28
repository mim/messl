function channel_response = messlSpChannelAdaptGmms(initial_gmm, obs, p_wtc, B)
% Update left and right channel responses for each source GMM.
% Returns the updated GMM structure and the magnitude response of the
% channel.

[W T] = size(obs);
ndct = size(B, 2);
% posterior weighted projection of observation onto DCT bases
op = zeros(ndct, 1);
% posterior weighted projection of DCT bases onto themselves
Bp = zeros(ndct);
mu = initial_gmm.means;
for t = 1:T
  % + 1e-200 to avoid dividing by zero
  p_wc = double(squeeze(p_wtc(:,t,:))) + 1e-200;
  pcv = initial_gmm.covars ./ p_wc;
  
  tobs = repmat(obs(:,t), [1 initial_gmm.nmix]);
  opt = B' * sum((tobs - mu)./pcv, 2);
  Bpt = B' * diag(sum(1./pcv, 2)) * B;
    
  op = op + opt;
  Bp = Bp + Bpt;
end
chan = inv(Bp)*op;
channel_response = B*chan;

