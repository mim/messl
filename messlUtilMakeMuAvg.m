function muAvg = messlUtilMakeMuAvg(W, bands)
% Make a matrix that will average the bands of mu_wi when it
% multiplies it.  W is the number of frequencies.  If bands is a
% scalar it specifies the number of bands, if it's a vector it
% specifies the highest freq in each band.

muAvg = eye(W);
if isscalar(bands)
  B = round(linspace(0,W,bands+1));
else
  B = bands;
end

for i=2:length(B)
  muAvg(B(i-1)+1:B(i),B(i-1)+1:B(i)) = 1 / (B(i)-B(i-1));
end
