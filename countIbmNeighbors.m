function counts = countIbmNeighbors(outFile, cleanDir, thresh_db, nfft, fs)

% Count co-ocurrences of neighbors in ideal binary masks
%
% For use as compatPot in mrfGridLbp().  Neighbor order is up, right, down,
% left.

if ~exist('thresh_db', 'var') || isempty(thresh_db), thresh_db = 0; end
if ~exist('nfft', 'var') || isempty(nfft), nfft = 1024; end
if ~exist('fs', 'var') || isempty(fs), fs = 16000; end
nFiles = 20;

[cf nf] = getFileNames(cleanDir, nFiles);

F = nfft/2 + 1;  % number of frequencies
I = 2;           % number of classes (same vs different)
nNeigh = 4;      % number of neighbors in grid MRF
counts = ones([F I I nNeigh]);
for f = 1:length(cf)
    ibm = computeIbm(cf{f}, nf{f}, thresh_db, nfft, fs);
    ibm = trimIbm(ibm);
    counts = updateCounts(counts, ibm, nNeigh);
end
counts = counts(2:end-1,:,:,:);
clear ibm

save(outFile)


function [cf nf] = getFileNames(cleanDir, nFiles)
[~,cf] = findFiles(cleanDir, '-tgt.wav');
ord = randperm(length(cf));
cf = cf(ord(1:nFiles));

nf = cell(size(cf));
for i = 1:length(cf)
    nf{i} = mapCleanToNoisyFile(cf{i});
end

function nf = mapCleanToNoisyFile(cf)
nf = strrep(strrep(cf, '-tgt', ''), 'anech', 'reverb');


function ibm = computeIbm(cleanFile, noisyFile, thresh_db, nfft, fs)

[clr cfs] = wavReadBetter(cleanFile);
[nlr nfs] = wavReadBetter(noisyFile);

clr = resample(clr, fs, cfs);
nlr = resample(nlr, fs, nfs);

[cL,cR] = binSpec(clr', nfft);
[nL,nR] = binSpec(nlr', nfft);

% Use geometric average magnitude of the two channels
target = 0.5 * (db(abs(cL)) + db(abs(cR)));
noise  = 0.5 * (db(abs(nL - cL)) + db(abs(nR - cR)));

ibm = target - noise > thresh_db;


function ibm = trimIbm(ibm)
% Trim off beginning and ending frames with no target
nTarget = sum(ibm,1);
start = 1;
while nTarget(start) == 0
    start = start + 1;
end
stop = length(nTarget);
while nTarget(stop) == 0
    stop = stop - 1;
end
ibm = ibm(:,start:stop);


function counts = updateCounts(counts, ibm, nNeigh)
% Counts is FxIxIx4, where dimensions are (frequency, targetClass,
% neighborClass, neighborDirection).
%
% Neighbors:     1      [(-1,0), (0,1), (1,0), (0,-1)]
%             4  X  2   df = [-1  0  1  0] = mod(n,  2).*(n-2)
%                3      dt = [ 0  1  0 -1] = mod(n+1,2).*(3-n)

[F T] = size(ibm);
ibm = 1 + ibm;   % {0,1} -> {1,2}
for n = 1:nNeigh
    df = mod(n,   2).*(n - 2);  % [-1  0  1  0];
    dt = mod(n+1, 2).*(3 - n);  % [ 0  1  0 -1];
    
    fi = max(1,1-df):min(F,F-df);
    target = ibm(fi, max(1,1-dt):min(T,T-dt));
    neighbor = ibm(max(1,1+df):min(F,F+df), max(1,1+dt):min(T,T+dt));

    for i1 = 1:size(counts,2)
        for i2 = 1:size(counts,3)
            counts(fi,i1,i2,n) = counts(fi,i1,i2,n) ...
                + sum((target == i1) .* (neighbor == i2), 2);
        end
    end
end
