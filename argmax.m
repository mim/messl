function am = argmax(X, dim)

% am = argmax(X, dim)
%
% Find the index of the maximum value of the matrix X.  If dime is
% supplied, find the maximum index along the dimension dim.

% Copyright (C) 2005 Michael Mandel, mim at ee columbia edu;
% distributable under GPL

if(nargin < 2)
  [dummy, am] = max(X);
else
  [dummy, am] = max(X, [], dim);
end
