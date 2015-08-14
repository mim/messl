function [perMic estPerPair err failed] = perMicTdoaLs(perPair, channelPairs, maxDelayErr, beRobust)

% Convert pair-wise ITDs measured in samples into approximate per-mic TDOAs using least squares

if ~exist('maxDelayErr', 'var') || isempty(maxDelayErr), maxDelayErr = 4; end
if ~exist('beRobust', 'var') || isempty(beRobust), beRobust = 1; end
minKeepPairFrac = 0.5;

Ns = size(perPair, 2);
Np = size(channelPairs, 1);
Ch = max(channelPairs(:));

% Matrix representing pairs of mics, plus mic 1 by itself at the bottom
A = full(sparse([Np+1; (1:Np)'; (1:Np)'], ...
    [1; channelPairs(:,1); channelPairs(:,2)], ...
    [1; ones(Np,1); -ones(Np,1)], ...
    Np+1, Ch));

% Augment with 0 TDOA for mic 1 at the bottom
perPair = [perPair; zeros(1, Ns)];

lastKeep = false(size(perPair));
keep = true(size(perPair));
unRobustObs = true(1, size(perPair,2));

% Loop until stable, unless too many pairs are excluded
while any(unRobustObs) % ~all(keep(:) == lastKeep(:)) || any(mean(keep,1) < minKeepPairFrac) || any(sum(keep,1) < Ch)
    for s = find(unRobustObs) %1:Ns
        % Robust least squares
        perMic(:,s) = A(keep(:,s),:) \ perPair(keep(:,s),s);
    end
    
    % Error per mic pair in the model
    estPerPair = A(1:end-1,:) * perMic;
    err = estPerPair - perPair(1:end-1,:);
 
    lastKeep = keep;
    keep = [abs(err) < maxDelayErr; true(1, Ns)];
    
    if ~beRobust
        failed = false(1, Ns);
        break
    end    

    failed = (mean(keep,1) < minKeepPairFrac) | (sum(keep,1) <= Ch+1);
    unRobustObs = any(keep ~= lastKeep, 1) & ~failed;
end
%plot(estPerPair, perPair(1:end-1,:), '.')
%drawnow
