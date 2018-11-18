function [Q, sumS, isFullRank]=unt(A)

% If A is normal Q=eigVec * sgn(eigVal) * eigVec' (normal).
% If A is hermittian then Q is hermitian.
% If A is hpd then Q is identity.

[U,S,V]=svd(A,'econ');
r=sum(diag(S)>1e-8);
if true == true
    Q=U(:,1:r)*V(:,1:r)';
else
    sz=min(size(U,1),size(V,1));
    Q=U(:,1:sz)*V(:,1:sz)';
end
sumS=sum(diag(S));
isFullRank= r==min(size(A));
