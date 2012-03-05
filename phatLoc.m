function [tdoa,masks,failure] = phatLoc(lr, tau, I, frame, vis)
 
% [tdoa,masks, failure] = phatLoc(lr, tau, I, frame, vis)
% 
% Localize I sources in the signal LR, where candidate locations are
% specified by TAU.  If I sources did not stand out in the signal, the
% remainder will be chosen arbitrarily and failure will be set to 1.
% Everything else will still proceed as normal.

if ~exist('vis', 'var'), vis = 0; end

% Run PHAT on chunks of this many taus
tau_chunk = 20;

% Do the first two taus to know the size to make T, allocate T
T = tdoaOverTime(lr, tau(1:2), frame, 'phat');
T(length(tau),end) = 0;

% Do the rest of the taus
for i=3:tau_chunk:length(tau)
  tau_end   = min(i+tau_chunk-1, length(tau));
  T(i:tau_end,:) = tdoaOverTime(lr, tau(i:tau_end), frame, 'phat');
end

% Look at histogram of maxima
aT = argmax(T,1);
ploc = hist(tau(aT), tau);

% Find I peaks in the histogram
if vis, plot(tau, ploc), drawnow, end
for i=1:I
  %plot(tau, ploc), pause
  [mxs(i), loc(i)] = max(ploc);
  %z = min(max(loc(i), 2), length(ploc)-1);
  to_cancel = abs(tau - tau(loc(i))) <= 1;
  mxs(i) = to_cancel * ploc';
  ploc = ploc .* (1-to_cancel);
  %ploc(z + [-1:1]) = 0;

  % Find frames associated with that point in the specgram
  masks(:,:,i) = repmat(abs(aT - loc(i)) <= 1, frame/2-1, 1);
end
tdoa = sort(tau(loc));
mx_left = max(ploc);

%% If the last peak is less than 10% the size of the first, it's
%% probably not real
failure = mxs(end) < mxs(1)/12;

% If the last peak is less than 5% of the first or it's less than
% twice the next peak down, it's probably not real
%failure = (mxs(end) < mxs(1)/15) || (mxs(end) < mx_left*1.5);
