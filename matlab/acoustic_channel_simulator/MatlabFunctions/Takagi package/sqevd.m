function [d,ifail,Q] = sqevd(a,b)
% [d,ifail,Q] = sqevd(a,b)
%
% Symmetric pentadiagonal eigenvalue decomposition using
% divide-and-conquer method. The symmetric quindiagonal
% matrix is given by T*T', where
%       T = diag(b,-1) + diag(a) + diag(b,1)
% is complex-symmetric and tridiagonal.
% Assume no zeros in the subdiagonal b.
%
% Inputs:
%   a -- main diagonal of a complex-symmetric and tridiagonal matrix
%   b -- subdiagonal of a complex-symmetric and tridiagonal matrix
% Outputs:
%     d -- eigenvalue vector of T*T'
% ifail -- 0 for success and 1 for failure.
%     Q -- eigenvector matrix of T*T'
%
% On successful return
%          T*T' = Q*diag(d)*Q'
%
% Dependency:
%   CSSVD - complex-symmetric SVD
%   sevr  - symmetric eigenproblem rank-one modification 
%
% W. Xu and S. Qiao  McMaster Univ.  Nov. 2004
%  Revised  Feb., Apr. 2005

% References
%  W. Xu and S. Qiao.
%  A divide-and-conquer method for the Takagi factorization.
%  Technical Report No. CAS 05-01-SQ. Department of Computing and
%  Software, McMaster University, Hamilton, Ontario L8S 4K1, Canada.
%  February 2005.
%
%  J.W.Demmel.
%  Applied Numerical Linear Algebra.
%  SIAM, 1997, pp 216-226.

n = length(a);
d = zeros(n,1);
Q = eye(n,n);

% if the matrix is small, apply the implicit QR method
if n < 5
    [d,Q] = CSSVD(a,b);
    d = d.*d;
    ifail = 0;
    return
end
% if the matrix is large, apply the divide-and-conquer method
% for symmetric eigenvalue problem
%
% divide the matrix
mid = floor(n/2);
% first block T1 = diag(b1,-1) + diag(a1) + diag(b1,1)
a1 = a(1:mid);          % diagonal of the first block
b1 = b(1:mid-1);        % subdiagonal of first block
% second block T2 = diag(b2,-1) + diag(a2) + diag(b2,1)
a2 = a(mid+1:n);        % diagonal of second block
b2 = b(mid+1:n-1);      % subdiagonal of second block
% mid entry
bm = b(mid);
%
% four rank-one symmetric modifications 
% T*T' =  [T1*T1' - (T1*em)*(T1*em)'              0             ]
%         [           0                T2*T2' - (T2*e1)*(T2*e1)']
%       + z1*z1' + z2*z2'
% where the two rank-one modification vectors are
% z1 = [zeros(mid - 1,1);
%             bm;
%            a2(1);
%            b2(1);
%       zeros(n - mid - 2,1)]
% z2 = [zeros(mid - 2,1);
%         b1(mid - 1);
%           a1(mid);
%             bm;
%       zeros(n - mid - 1,1)]
%
% apply the divide-and-conquer to the first block T1
[d1,ifail,Q1] = sqevd(a1,b1);
if ifail == 1
    return
end
% T1 rank-one modification vector Q1'*(T1*em)
rhow1 = b1(mid -1)'*b1(mid-1) + a1(mid)'*a1(mid);
nrm = sqrt(rhow1);
w1 = (b1(mid -1)/nrm)*Q1(mid-1,:)' + (a1(mid)/nrm)*Q1(mid,:)';
% symmetric rank-one modification diag(d1) - rhow*w1*w1'
[d1,ifail,Q1] = sevr(d1,w1,-rhow1,Q1);
if ifail == 1
    return
end
%
% apply the divide-and-conquer to the second block T2
[d2,ifail,Q2] = sqevd(a2,b2);
if ifail == 1
    return
end
% T2 rank-one modification vector Q2'*(T2*e1)
rhow2 = a2(1)'*a2(1) + b2(1)'*b2(1);
nrm = sqrt(rhow2);
w2 = (a2(1)/nrm)*Q2(1,:)' + (b2(1)/nrm)*Q2(2,:)';
% symmetric rank-one modification diag(d2) - rhow2*w2*w2'
[d2,ifail,Q2] = sevr(d2,w2,-rhow2,Q2);
if ifail == 1
    return
end
%
% combine T1 and T2
d = [d1; d2];
% symmetric rank-one update z1*z1' 
rhoz1 = bm'*bm + a2(1)'*a2(1) + b2(1)'*b2(1);
nrm = sqrt(rhoz1);
z1 = [(bm/nrm)*Q1(mid,:)';
      (a2(1)/nrm)*Q2(1,:)' + (b2(1)/nrm)*Q2(2,:)'];
[d,ifail,QQ] = sevr(d,z1,rhoz1);
if ifail == 1
    return
end
% update the eigenvectors
Q(1:mid,:) = Q1*QQ(1:mid,:);
Q(mid+1:n,:) = Q2*QQ(mid+1:n,:);
%
% symmetric rank-one update z2*z2'
rhoz2 = b1(mid-1)'*b1(mid-1) + a1(mid)'*a1(mid) + bm'*bm;
nrm = sqrt(rhoz2);
z2 = (b1(mid-1)/nrm)*Q(mid-1,:)' + (a1(mid)/nrm)*Q(mid,:)' + ...
      (bm/nrm)*Q(mid+1,:)';
[d,ifail,Q] = sevr(d,z2,rhoz2,Q);
if ifail == 1 
    return
end
