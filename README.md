# Model-Based EM Source Separation and Localization

Copyright 2006-2009 Michael I Mandel and Ron Weiss, all rights reserved
<mim@ee.columbia.ed> and <ronw@ee.columbia.edu>
Last updated 2009-08-20


## Basic usage to separate two sources:

```matlab
addpath('./plottools/')

% % Load stereo wav files of the same length and mix them
% [y1 fs] = wavread('data/src1.wav');
% [y2 fs] = wavread('data/src2.wav');
% lr = y1' + y2';

% Load pre-mixed version of those two files (don't forget the transpose)
[lr fs] = wavread('data/mix.wav');
lr = lr';

% Derive grid of tau values (use your numbers here)
c    = 340;  % speed of sound in air [M/s]
d    = 0.3;  % slightly more than the acoustic distance between microphones [M]
nTau = 31;   % number of tau samples to use

maxItd = fs * d / c;  % maximum ITD [samples]
tau = linspace(-maxItd, maxItd, nTau);

% Run MESSL
[m,p] = messl(lr, tau, 2, 'vis', 1);

% Reconstruct wavforms from masks
yhat1 = reconstruct(m, lr, 1);
yhat2 = reconstruct(m, lr, 2);
```


## Fancier usage

Initialized from PHAT-histogram:
```matlab
% Localize and then run MESSL
tdoa = phatLoc(lr, tau, 2, 1024, 1);
[m,p] = messl(lr, tau, 2, 'vis', 1, 'tauPosInit', tdoa);
```


Even fancier usage, garbage source and ILD prior (better in reverb,
but only when using dummy-head recordings):

```matlab
[m,p] = messl(lr, tau, 2, 'vis', 1, 'ildPriorPrec', 3, ...
              'GarbageSrc', 1, 'sr', 16000);
```


Can also use prob2mask to make the mask more definitive, i.e. closer
to binary, but not binary.

```matlab
m2 = prob2mask(m);
yhat1 = reconstruct(m2, lr, 1);
yhat2 = reconstruct(m2, lr, 2);
```
