function [X, n, Sigma2] = Pro2TraceNorm(Z, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min: 1/2*||Z-X||^2 + ||X||_tr - alpha*||X||_F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% new
[m, n] = size(Z);
if 2*m < n
    AAT = Z*Z';
    [S, Sigma2, D] = svd(AAT);
    Sigma2 = diag(Sigma2);
    V = sqrt(Sigma2);
    tol = max(size(Z)) * eps(max(V));
    n = sum(V > max(tol, tau));
    mid = max(V(1:n)-tau, 0)  ;
    alpha = 0.6; % We fixed alpha here. Please see the paper
    coef = (norm(mid(1:n)) + alpha*tau)  / norm(mid(1:n)) ;
    mid = mid.*coef./ V(1:n);
    X = S(:, 1:n) * diag(mid) * S(:, 1:n)' * Z;
    return;
end
if m > 2*n
    [X, n, Sigma2] = Pro2TraceNorm(Z', tau);
    X = X';
    return;
end
[S,V,D] = svd(Z);
Sigma2 = diag(V).^2;
n = sum(diag(V) > tau);
X = S(:, 1:n) * max(V(1:n,1:n)-tau, 0) * D(:, 1:n)';




