function [d,Q] = cstsvdt(a,b)
% [d,ifail,Q] = cstsvdd(a,b)
%
% Symmetric SVD of a complex-symmetric and tridiagonal matrix
%            T = diag(b,-1) + diag(a) + diag(b,1)
% using twisted factorization method.
%
% Inputs: 
%     a -- main diagonal of a complex-symmetric and tridiagonal matrix
%     b -- subdiagonal of a complex-symmetric and tridiagonal matrix
%
% Outputs:
%     d -- singular (Takagi) value vector
%     Q -- left singular vector (Takagi vector) matrix
%          so that
%              diag(b,-1) + diag(a) + diag(b,1) = Q*diag(d)*Q.'
%
% Dependency:
%    cssing -- singular values of a complex symmetric tridiagonal T
%              given by its diagonals
%    cstlqd -- LQ decomposition of a complex-symmetric and 
%              tridiagonal matrix
%    cssvdt -- compute the corresponding eigenvector using the 
%              twisted factorization.
%   eig2tak -- convert eigenvalues and eigenvectors into Takagi
%              values and Takagi vectors
%
% Reference:
%  W. Xu and S. Qiao.
%  A divide-and-conquer method for the Takagi factorization.
%  Technical Report No. CAS 05-01-SQ. Department of Computing and
%  Software, McMaster University, Hamilton, Ontario L8S 4K1, Canada.
%  February 2005.
%
% W. Xu and S. Qiao  McMaster Univ. May 2007

% get dimension of the matrix Q
n = length(a);

% get the singular values of a complex symmetric tridiagonal T, squared to get
% eigenvalues of T*T'
e = (cssing(a,b)).^2;

% get the LQ decomposition of a complex-symmetric and tridiagonal matrix
[d,l1,l2] = cstlqd(a,b);

% get the corresponding eigenvector for the eigenvalue, LQ decomposition
for k = 1:n
    Q(:,k) = twist(d,l1,l2,e(k));
end

% convert eigenvalues and eigenvectors of T*T' into
% Takagi (singular) values and Takagi (left singular)
% vectors of T
[d, Q] = eig2tak(a,b,e,Q);