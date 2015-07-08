function [perMic estPerPair err] = perMicTdoaLs(perPair, channelPairs)

% Convert pair-wise ITDs into approximate per-mic TDOAs using least squares

Np = size(channelPairs,1);
Ch = max(channelPairs(:));

% Matrix representing pairs of mics, plus mic 1 by itself at the bottom
A = full(sparse([Np+1; (1:Np)'; (1:Np)'], ...
    [1; channelPairs(:,1); channelPairs(:,2)], ...
    [1; ones(Np,1); -ones(Np,1)], ...
    Np+1, Ch));

% Augment with 0 TDOA for mic 1 at the bottom
perPair = [perPair; zeros(1, size(perPair,2))];

% Least squares
perMic = A \ perPair;

% Error per mic pair in the model
estPerPair = A(1:end-1,:) * perMic;
err = estPerPair - perPair(1:end-1,:);
