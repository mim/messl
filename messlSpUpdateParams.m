function spParams = messlSpUpdateParams(W, T, I, Nt, C, rep, spParams, ...
    nuSp, L, R, Nrep);
for i = 1:I
  if spParams.garbageSrc && i == I
    continue
  elseif spParams.spMode > 0
    if isfield(spParams, 'ev_params')
      [w hl hr] = messlSpEigenvoiceAndChannelAdaptGmms(spParams, L, R, ...
          squeeze(nuSp(:,:,:,i,:)));
      spParams.w(i,:) = w;
      spParams.hl(i,:) = hl;
      spParams.hr(i,:) = hr;
      spParams.channel_response(1,i,:) = spParams.B * hl;
      spParams.channel_response(2,i,:) = spParams.B * hr;
      spParams.sourcePriors(i) = ...
          construct_adapted_model(spParams.gmm, spParams.ev_params, w);
    else
      spParams.channel_response(1,i,:) = messlSpChannelAdaptGmms( ...
	  spParams.sourcePriors(i), L, squeeze(nuSp(1,:,:,i,:)), spParams.B);
      spParams.channel_response(2,i,:) = messlSpChannelAdaptGmms(  ...
	  spParams.sourcePriors(i), R, squeeze(nuSp(2,:,:,i,:)), spParams.B);
    end
  elseif spParams.spMode < 0
    obs = L + R;

    p_wtc = squeeze(mean(nuSp(:,:,:,i,:), 1));
    gmm = spParams.sourcePriors(i);
    ndct = size(spParams.B, 2);
    tmpA = zeros(ndct, ndct);
    tmpb = zeros(ndct, 1);
    for t=1:T
      tmp = squeeze(p_wtc(:,t,:)) ./ (4*gmm.covars + ...
          repmat(spParams.channel_response_var(i,:)', [1 C]));
      tobs = repmat(obs(:,t), [1 C]);
      tmpA = tmpA + spParams.B'*diag(sum(tmp, 2))*spParams.B;
      tmpb = tmpb + spParams.B'*sum(tmp.*(tobs - 2*gmm.means), 2); 
    end
    dct_channel_response_i = inv(tmpA) * tmpb;
    
    %if rep < 0.66*Nrep
    %  start = 3+fix(0.5*(rep-1));
    %  dct_channel_response_i(start:end,:) = 0;
    %end
    
    spParams.dct_channel_response(i,:) = dct_channel_response_i;
    spParams.channel_response(i,:) = spParams.B * dct_channel_response_i;

    % FIXME - the update here is incorrect.  It leads to negative
    % variances (b/c we just subtract off 4*gmm.covars)
    if spParams.spMode == -2
      warning(['messl: the SP-ILS variance update is known to be' ...
            ' incorrect.  Don''t expect this to work properly.'])
      % variance
      tmpA = zeros(ndct, ndct);
      tmpb = zeros(ndct, 1);
      for t=1:T
        tmp = squeeze(p_wtc(:,t,:));
        dzm = repmat(obs(:,t), [1 C]) - 2*gmm.means ...
            - repmat(squeeze(spParams.channel_response(i,:))', [1 C]);
        tmpA = tmpA + spParams.B'*diag(sum(tmp, 2))*spParams.B;
        tmpb = tmpb + spParams.B'*sum(tmp.*(dzm.^2 - 4*gmm.covars), 2); 
      end
      dct_channel_response_var_i = inv(tmpA) * tmpb;
    
      %if rep < 0.66*Nrep
      %  start = 3+fix(0.5*(rep-1));
      %  dct_channel_response_var_i(start:end,:) = 0;
      %end
    
      spParams.dct_channel_response_var(i,:) = dct_channel_response_var_i;
      spParams.channel_response_var(i,:) = spParams.B * dct_channel_response_var_i;
    end
  end
end
