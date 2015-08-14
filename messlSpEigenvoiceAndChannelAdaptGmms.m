function [w hl hr] = messlSpEigenvoiceAndChannelAdaptGmms(spParams, L, R, p_lr_wtc)
[W T] = size(L);
C = size(spParams.ev_params.mean, 2);
ndct = size(spParams.B, 2);
nev = length(spParams.ev_params.means);

icvl = zeros(W, C);
icvr = zeros(W, C);
icvlobs = zeros(W, C);
icvrobs = zeros(W, C);
for t = 1:T
  Sl = double(squeeze(p_lr_wtc(1,:,t,:))) ./ spParams.gmm.covars;
  Sr = double(squeeze(p_lr_wtc(2,:,t,:))) ./ spParams.gmm.covars;
  icvl = icvl + Sl;
  icvr = icvr + Sr;
  icvlobs = icvlobs + diag(sparse(L(:,t))) * Sl;
  icvrobs = icvrobs + diag(sparse(R(:,t))) * Sr;
end

Aww = zeros(nev, nev);
Awl = zeros(nev, ndct);
Awr = zeros(nev, ndct);
All = zeros(ndct, ndct);
Arr = zeros(ndct, ndct);
bw = zeros(nev, 1);
bl = zeros(ndct, 1);
br = zeros(ndct, 1);

u = spParams.ev_params.mean;
Uct = zeros(nev, W);
B = spParams.B;
for c = 1:C
  for j = 1:nev
    Uct(j,:) = spParams.ev_params.means{j}(:,c);
  end

  Sl = diag(sparse(icvl(:,c)));
  Sr = diag(sparse(icvr(:,c)));
  BSl = B' * Sl;
  BSr = B' * Sr;
  Slobs = icvlobs(:,c);
  Srobs = icvrobs(:,c);

  Aww = Aww + Uct * (Sl + Sr) * Uct';
  Awl = Awl + Uct * BSl';
  Awr = Awr + Uct * BSr';
  All = All + BSl * B;
  Arr = Arr + BSr * B;
  bw = bw + Uct * (Slobs - Sl * u(:,c) +  Srobs - Sr * u(:,c));
  bl = bl + B' * (Slobs - Sl * u(:,c));
  br = br + B' * (Srobs - Sr * u(:,c));
end
z = zeros(ndct, ndct);
tmp = inv([Aww Awl Awr; Awl' All z; Awr' z Arr]) * [bw; bl; br];

w = tmp(1:nev);
hl = tmp(nev+[1:ndct]);
hr = tmp(nev+ndct+[1:ndct]);

