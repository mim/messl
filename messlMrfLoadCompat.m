function mrfCompatPot = messlMrfLoadCompat(mrfCompatFile, I, garbageSrc)
% Hard coded for 2 sources and a garbage source for now...

if isempty(mrfCompatFile) || ~exist(mrfCompatFile, 'file') % || (I ~= 2)
    mrfCompatPot = [];
else
    load(mrfCompatFile);
    mrfCompatPot = counts;
end
