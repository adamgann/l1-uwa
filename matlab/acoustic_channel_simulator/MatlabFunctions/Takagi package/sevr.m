function [d, ifail, Q] = sevr(din, z, rho, Qin)
% [d, ifail, Q] = sevr(din, z, rho, Qin)
%
% Symmetric eigenproblem rank-one modification
%      Q*diag(d)*Q' = Qin*(diag(din) + rho*z*z')*Qin'
%
% inputs
%   din -- input eigenvalues
%     z -- updating vector, normalized in 2-norm
%   rho -- scalar factor in updating
%   Qin -- input unitary matrix, optional, if not provided
%          the outputs are the same as if Qin = eye
% outputs
%     d -- updated eigenvalues
% ifail -- 0 for success, 1 for fail
%     Q -- updated unitary (eigenvector) matrix
%
% dependency
%  deflate -- introduce zeros in z for equal diagonal elements
%             in din
%      slv -- secular equation solver for eigenvalues and eigenvectors
%
% Reference:
%  W. Xu and S. Qiao.
%  A divide-and-conquer method for the Takagi factorization.
%  Technical Report No. CAS 05-01-SQ. Department of Computing and
%  Software, McMaster University, Hamilton, Ontario L8S 4K1, Canada.
%  February 2005.
%
% W. Xu and S. Qiao   McMaster Univ.   April 2005

if nargin < 4		      % Qin is not provided
    n = length(din);      % set Qin to identity
    Qin = eye(n,n);
end

% deflation
[d, z, inzeros, indx, Q] = deflate(din, z, rho, Qin);
dd = d(inzeros);        % entries are distinct and sorted in decreasing order
zz = z(inzeros);        % nonzero entries
arho = abs(rho);

% compute eigenvalues and eigenvectors
k = length(inzeros);
lam = zeros(k,1);	% eigenvalues
G = zeros(k,k);		% eigenvectors
for i=1:k
    [lam(i), G(:,i), ifail] = slv(i, dd, zz, arho);
end

% update eigenvalues
d(inzeros) = lam;
d = sign(rho)*d;
% update eigenvectors
Q(:,indx(inzeros)) = Q(:,indx(inzeros))*G;
Q = Q(:,indx);
