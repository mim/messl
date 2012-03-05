function bands = makeBands(mode, W, startBands, Nrep)
% Make a vector that determines the number of bands to use in
% each round of EM.  Nonlin is between 0 and 1 and controls the
% amount of non-linearity in the progression of frequencies per
% band.  The closer it is to 1, the more repetitions bands stays
% close to startBands and the fewer it stays close to endBands.
if mode
  nonlin = 0.99;

  if mode >= 0
    endBands = mode;
  else
    endBands = -W/mode;
  end

  % Never go from more bands to fewer bands
  startBands = min(endBands, startBands);

  leadIn = 4;
  leadOut = 4;
  bands = [startBands*ones(1,leadIn-1) ...
           round(ls10(startBands, endBands, Nrep-leadIn-leadOut+2)) ...
           endBands*ones(1,leadOut-1)];

  %bands = round(ls10(startBands-nonlin, endBands-nonlin, Nrep)+nonlin);
else
  bands = [];
end
