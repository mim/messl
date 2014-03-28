function D = x2fx2(X, model)

% Convert matrix of predictors, X, into a design matrix for
% regression, D
%
% If X is NxM, then model must be PxM and D will be NxP. Each
% column of D will be a combination of the columns of X raised to
% the corresponding powers in model.  So if X has two columns and
% model is [0 1; 1 0; 1 2], then D will have four columns:
% [ones(N,1) X(:,1).^0.*X(:,2).^1 X(:,1).^1.*X(:,2).^0 X(:,1).^1.*X(:,2).^2] 
% The columns of ones is always added.

[N M] = size(X);
[P M2] = size(model);

D = zeros(N, P);
for p = 1:P
    D(:,p) = prod(bsxfun(@power, X, model(p,:)), 2);
end
