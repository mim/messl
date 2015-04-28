function B = messlUtilDctInit(dctMode, W)

% DCT basis vectors
if dctMode
  B = zeros(W, dctMode);
  B(:,1) = 1/sqrt(W);
  B(:,2:end) = sqrt(2/W) * cos(pi/W * [0.5:W-0.5]'*[1:dctMode-1]);
else
  B = eye(W);
end
