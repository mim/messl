function messlUtilVisualizeParams(W, T, I, tau, sr, ipdParams, ildParams, spParams, ...
    p_lr_iwt, maskIpd, maskIld, maskSp, L, R, reliability)
%%%% Visualize what's happening
% 4 figures:
%   1. Factored masks
%   2. Recovered signals (i.e. mask .* signal)
%   3. Learned model parameters
%   4. IPD vs frequency distribution

sp_cols = max(I,2);

modes = [ipdParams.ipdMode, ildParams.ildMode, spParams.spMode];
params_on = sum(modes ~= 0);
names = {'IPD', 'ILD', 'SP'};
masks = {maskIpd, maskIld, maskSp};
obs = idB(L) + idB(R);

% Factored masks
%
% IPD masks
% ILD masks
% SP masks
% total masks (if more than one other is present)
mask_plots = {};
mask_titles = {};
mask_layout = [params_on+(params_on>1) I];


% Recovered signals (i.e. mask .* signal)
%
% IPD reconstructions
% ILD reconstructions
% SP reconstructions
% total reconstructions
% original mixture
recovered_plots = {};
recovered_titles = {};
recovered_layout = [params_on+1+(params_on>1) I];

for m=1:length(modes)
  if(modes(m))
    for i=1:I
      mask_plots = {mask_plots{:}, masks{m}(:,:,i)};
      t = sprintf('%s probability of being in source %d', names{m}, i);
      mask_titles = {mask_titles{:}, t};
      
      recovered_plots = {recovered_plots{:}, dB(masks{m}(:,:,i) .* obs)};
      t = sprintf('%s reconstruction of source %d', names{m}, i);
      recovered_titles = {recovered_titles{:}, t};
    end
  end
end

if params_on > 1
  for i=1:I
    mask_plots = {mask_plots{:}, squeeze(mean(p_lr_iwt(:,:,:,i), 1))};
    t = sprintf('Overall probability of being in source %d',i);
    mask_titles = {mask_titles{:}, t};
    
    recovered_plots = {recovered_plots{:}, ...
                       dB(squeeze(mean(p_lr_iwt(:,:,:,i),1)) .* obs)};
    t = sprintf('Overall reconstruction of source %d', i);
    recovered_titles = {recovered_titles{:}, t};
  end
end

recovered_plots = {recovered_plots{:}, dB(obs)};
recovered_titles = {recovered_titles{:}, 'Mixed signal'};

if ~isempty(reliability)
  recovered_plots = {recovered_plots{:}, 50*(reliability-1)};
  recovered_titles = {recovered_titles{:}, 'Reliability'};
end

% plotall(mask_plots, 'title', mask_titles, 'figure', 1, 'subplot', ...
%         mask_layout, 'CLim', [0 1]);
% 
% plotall(recovered_plots, 'title', recovered_titles, 'figure', 2, ...
%         'subplot', recovered_layout, 'CLim', [-60 0]);

fig_no_focus(1)
subplots(mask_plots, mask_layout, mask_titles, @(r,c,i) caxis([0 1]));
drawnow

fig_no_focus(2)
subplots(recovered_plots, recovered_layout, recovered_titles, @(r,c,i) caxis([-60 0]));
drawnow

% imgsc(mask_plots, 'title', mask_titles, 'figure', 1, 'subplot', ...
%       mask_layout, 'caxis', [0 1]);

% imgsc(recovered_plots, 'title', recovered_titles, 'figure', 2, ...
%       'subplot', recovered_layout, 'caxis', [-60 0]);


% Learned model parameters
%
% p(tau | i), ILD means and stds
% SP parameters
fig_no_focus(3)

if modes(3) == 1
  sp_rows = 2;
  sp_cols = I;
else
  sp_rows = 1;
  sp_cols = sum(modes ~= 0);
end

curPlot = 1;
if ipdParams.ipdMode
    subplot(sp_rows,sp_cols,curPlot)
    curPlot = curPlot + 1;
    plot(tau, ipdParams.p_tauI')
    title('p(tau, i)');
    grid on
end

if ildParams.ildMode
  subplot(sp_rows,sp_cols,curPlot)
  curPlot = curPlot + 1;
  plot(1:W, ildParams.mu_wi)
  hold on
  plot(1:W, ildParams.mu_wi + sqrt(ildParams.h2_wi), ':');
  plot(1:W, ildParams.mu_wi - sqrt(ildParams.h2_wi), ':');
  hold off
  grid on
  title('ILD');
end

if spParams.spMode == 1
  for i = 1:I
    subplot(sp_rows,sp_cols,sp_cols+i)
    plot(squeeze(spParams.channel_response(:,i,:))')
    if i == 1
      ax = axis();
    else
      axis(ax);
    end
    title(sprintf('Room+HRTF channel for source %d', i))
    grid on;
  end
  if isfield(spParams, 'ev_params')
    subplot(sp_rows,sp_cols,3)
    plot(spParams.w')
    title('Eigenvoice parameters')
    grid on;
  end
elseif spParams.spMode
  subplot(sp_rows,sp_cols,3)
  mu = squeeze(spParams.channel_response');
  plot(1:W, mu)
  if spParams.spMode == -2
    hold on
    plot(1:W, mu + sqrt(spParams.channel_response_var'), ':');
    plot(1:W, mu - sqrt(spParams.channel_response_var'), ':');
    hold off
  end
  grid on
  title('L+R channel')
end
subplot 111
  

% IPD vs Frequency distribution
%
fig_no_focus(4);
visParams(ipdParams, tau, sr);


drawnow
