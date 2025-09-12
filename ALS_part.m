function [X err errobs]=ALS_simple(A, B, r, Xtrue, lambda, maxiter); 

[n,m]=size(Xtrue);

Ir=eye(r,r);
L=ones(r,n);
R=ones(r,m);

for j=1:m
    rowind{j} = find(A(:,j));
end
for j=1:n
    colind{j} = find(A(j,:));
end

err=[];
errobs=[];
for i=1:maxiter
    
    %% vector derivative
    %R = (L*L'+lambda*Ir)\(L*B);
    %L = (R*R'+lambda*Ir)\(R*B');
    
    %% vector by vector (slow)
    % ‰ð‚«•û‚Í‚±‚±‚©‚çhttp://stanford.edu/~rezab/classes/cme323/S15/notes/lec14.pdf
    for j=1:m
        L2 = L(:,rowind{j});
        R(:,j)=(L2*L2'+lambda*Ir)\(L2*B(rowind{j},j));
    end
        
    for j=1:n
        R2 = R(:,colind{j});
        L(:,j)=(R2*R2'+lambda*Ir)\(R2*B(j,colind{j})');
    end
    err(i)=norm(Xtrue-L'*R,'fro');
    errobs(i)=norm(A.*Xtrue- A.*(L'*R), 'fro');
end

X=L'*R;

