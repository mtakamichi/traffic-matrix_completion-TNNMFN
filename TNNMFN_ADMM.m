%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADM algorithm: tensor completion
% paper: Tensor completion for estimating missing values in visual data
% date: 05-22-2011
% min_X: \sum_i \alpha_i \|X_{i(i)}\|_*
% s.t.:  X_\Omega = T_\Omega
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U, errList] = WLRTC_truetrueADMM(T, Omega, alpha, gamma, maxIter, epsilon, T_org)


errList = zeros(maxIter, 1);
dim = size(T);
W = cell(ndims(T), 1);


normT = norm(T(:));

U = T;
V = cell(4,1);
V{1}= zeros(size(Unfold(T,size(T),1)));
V{2}= zeros(size(Unfold(T,size(T),2)));
V{3}= zeros(size(Unfold(T,size(T),3)));
V{4}= zeros(size(T));
D = V;
 
X = T;
X(~Omega) = mean(T(Omega));

for i = 1:ndims(T)
    W{i} = alpha(i)./(svd(Unfold(X, dim, i))+1e-5);
end

 
for k = 1: maxIter
    if mod(k, 20) == 0
        fprintf('ttWLRTC: iterations = %d   difference=%f NMAE=%f\n', k, errList(k-1), calc_NMAE(T_org, Omega, U));
    end
    gamma = gamma*1.05;
    
    % update U
    lastU = U;
    U = (Fold(V{1}-D{1},dim,1)+Fold(V{2}-D{2},dim,2)+Fold(V{3}-D{3},dim,3)+V{4}-D{4})./4;
    
    % update V1~V3
    for i = 1:3
       tau = W{i}/gamma;
       % 重みあり
       %V{i} = Pro2WTraceNorm_minus_F(Unfold(U, dim, i)+D{i}, tau);

       % 重みなし
       V{i} = Pro2TraceNorm_minus_F(Unfold(U, dim, i)+D{i}, alpha(i)/gamma);
       %Pro2WSPTraceNorm
    end
    
    % update V4
    temp = U+D{4};
    temp(Omega) = T(Omega);
    %temp = (temp>0).*temp; % non-negativity constraint
    V{4} = temp;
    
    % update D1~D3
    for i=1:3
        D{i} = D{i} + Unfold(U,dim,i) - V{i};
    end
    
    % update D4
    D{4} = D{4} + U - V{4};
    
    % compute the error
    errList(k) = norm(U(:)-lastU(:)) / normT;
    if errList(k) < epsilon
       break;
    end
end

errList = errList(1:k);
fprintf('ttWLRTC ends: total iterations = %d   difference=%f\n\n', k, errList(k));

