function peaks = pickPeaks(x, y, n)

% peaks = pickPeaks(x, y, n)
%
% Find the N highest peaks in a signal with range x and values y.

% Make sure y>=0
y = y - min(y);

% Use peaks of the cross-correlation
dt = floor(1 ./ median(diff(x)));
for i=1:n
  %plot(x,y), pause(.1)
  peaks(i) = min(max(1+dt, argmax(y)), length(y)-dt);

  % Find the parameters of a parabola around the peak
  xn = x(peaks(i) + [-dt:dt]);
  yn = y(peaks(i) + [-dt:dt]);
  A = [xn(:).^2 xn(:) ones(numel(xn),1)];
  F = A \ yn(:);
  
  % Fit the parabola to all of the x data
  A = [x(:).^2 x(:) ones(numel(x),1)];
  p = A * F;
  p = p .* (p > 0) * y(peaks(i))/p(peaks(i));
  
  % Subtract it from y
  y = y(:) - p;
end
%plot(x,y), pause(.1)
peaks = sort(x(peaks));
