function [d,ifail,Q] = cstsvdd(a,b)
% [d,ifail,Q] = cstsvdd(a,b)
%
% Symmetric SVD of a complex-symmetric and tridiagonal matrix
%            T = diag(b,-1) + diag(a) + diag(b,1)
% using divide-and-conquer method.
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
%     ifail -- 0 for success and 1 for fail
%
% Dependency:
%     sqevd -- eigenvalue decomposition of a symmetric quindiagonal
%              matrix using divide-and-conquer method
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
% W. Xu and S. Qiao  McMaster Univ. Nov. 2004
% Revised
%     January, 2005
%     February, 2005
%     April, 2005

SQRTEPS = sqrt(eps);	% tolerance for numerical zero
n = length(a);
d = zeros(n,1);
Q = eye(n,n);
% 
% partition T = diag(b,-1) + diag(a) + diag(b,1) into
% block diagonal with tridiagonal blocks.
% 
i = 1;
while i <= n
    % find a block with no zeros on the subdiagonal
    j = i;
    while ((j<n) & (b(j)>SQRTEPS))
        j = j + 1;
    end
    % block size is j-i+1
    if j > i
        % eigendecomposition of the block
        [d(i:j),ifail,Q(i:j,i:j)] = sqevd(a(i:j),b(i:j-1));
        if ifail == 1
            display('Failed'),
            return
        end
    else % bs=1
        d(i) = a(i)'*a(i);
        ifail = 0;
    end
    i = j + 1;
end

% convert eigenvalues and eigenvectors of T*T' into
% Takagi (singular) values and Takagi (left singular)
% vectors of T
[d, Q] = eig2tak(a,b,d,Q);
