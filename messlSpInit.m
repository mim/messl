function [spParams C] = messlSpInit(I, W, L, R, sourcePriors, stdInit, ...
    dctMode, spMode, garbageSrc)

B = messlUtilDctInit(dctMode, W);
spParams = struct('spMode', spMode, 'dctMode', dctMode, ...
    'garbageSrc', garbageSrc, 'B', B);

% Set spMode based on spStartRep.  spMode will remain turned off
% until rep == spStartRep.
if spMode
  % FIXME
  %spParams.spStartRep = 10;
  %spParams.origSpMode = spMode;
  %spParams.spMode = 0;

  %spParams.spStartRep = 2;
  spParams.spStartRep = 4;
  spParams.origSpMode = spMode;
  spParams.spMode = 0;
else
  spParams.spStartRep = Inf;
end

C = 0;
% Set up source prior parameters
if spMode
  if isfield(sourcePriors, 'ev_params')
    spParams.ev_params = sourcePriors.ev_params;
    spParams.gmm = sourcePriors.gmm;
    sourcePriors = spParams.gmm;
  end

  if length(sourcePriors) == 1
    tmp = sourcePriors;
    for i = 1:I
      sourcePriors(i) = tmp;
    end
    clear tmp
  elseif length(sourcePriors) ~= I
    error('Length of sourcePriors needs to be either 1 or I');
  end
  
  C = max([sourcePriors.nmix]);
  % Make sure all GMMs have the same number of components
  for i=1:I
    if sourcePriors(i).nmix < C
      sourcePriors(i).priors = ...
          [sourcePriors(i).priors repmat(-Inf, [1 C-sourcePriors(i).nmix])];
      sourcePriors(i).means = ...
          [sourcePriors(i).means zeros(W, C-sourcePriors(i).nmix)];
      sourcePriors(i).covars = ...
          [sourcePriors(i).covars ones(W, C-sourcePriors(i).nmix)];      
      sourcePriors(i).nmix = C;
    end
  end
  
  if garbageSrc
    garbageSourcePrior.nmix = C;
    garbageSourcePrior.priors = [0 repmat(-Inf, [1 C-1])];
    garbageSourcePrior.means = zeros(W, C);
    garbageSourcePrior.covars = ones(W, C);
    for i=1:I
      gmm = sourcePriors(i);
      for c=1:C
        p = exp(gmm.priors(c))/I;
        garbageSourcePrior.means(:,1) = garbageSourcePrior.means(:,1) ...
            + p*gmm.means(:,c);
        garbageSourcePrior.covars(:,1) = garbageSourcePrior.covars(:,1) ...
            + p*gmm.covars(:,c);
      end
    end
    sourcePriors(I+1) = garbageSourcePrior;
  end
  spParams.sourcePriors = sourcePriors;

  I = I + garbageSrc;
  switch spMode
    case 1
      spParams.channel_response = zeros(2,I,W);
    case {-1,-2}
      spParams.channel_response = zeros(I,W);
      spParams.channel_response_var = repmat(eps, [I W]);
 end
end
