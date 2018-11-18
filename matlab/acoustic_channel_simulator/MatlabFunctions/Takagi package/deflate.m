function [d, z, inzeros, indx, Q] = deflate(din, zin, rho, Qin)
% [d, z, inzeros, indx, Q] = deflate(din, zin, rho, Qin)
%
% deflation stage in the symmetric eigenvalue rank-one
% modification problem: diag(din) + rho*zin*zin'.
% sort sign(rho)*din and introduce zeros into zin for equal entries
% in din
%
% inputs
%    din -- eigenvalue vector to be updated
%    zin -- symmetric rank-one update vector, normalized
%    rho -- scalar factor in update
%    Qin -- unitary matrix
% outputs
%       d -- eigenvalue vector after deflation
%            d sorted in decreasing order
%       z -- symmetric rank-one update vector after deflation
% inzeros -- indices to the nonzero entries of z
%            d(inzeros) entries are distinct
%            z(inzeros) entries are nonzeros
%    indx -- permutation vector for Q
%    Q    -- updated unitary matrix
% on return
% Q(:indx)*sign(rho)*(sign(rho)*diag(d) + abs(rho)*z*z')*Q(:,indx)'
% = Qin*(diag(din) + rho*zin*zin')*Qin'
%
% Reference:
%  W. Xu and S. Qiao.
%  A divide-and-conquer method for the Takagi factorization.
%  Technical Report No. CAS 05-01-SQ. Department of Computing and
%  Software, McMaster University, Hamilton, Ontario L8S 4K1, Canada.
%  February 2005.
%
% W. Xu and S. Qiao  McMaster Univ.   April 2005

n = length(din);
d = din; z = zin;
Q = Qin;

inzeros = [];    	% indices of nonzeros in z

signrho = sign(rho);
% sort signrho*din in decreasing order
[d,indx] = sort(-signrho*d);
d = -d;
z = z(indx);		% permute z correspondingly

dmax = max(abs(d(1)), abs(d(n)));
tolz = (n*eps*dmax)/abs(rho);	% tolerance for numerical zeros in z

j = 1;
while j < n
    if abs(z(j)) > tolz		% z(j) not very small
        normz = norm([z(j);z(j+1)]);
        told = (dmax*eps)/(abs(z(j))/normz);
			% tolerance for d(j) = d(j+1) numerically
        if abs(d(j) - d(j+1)) < told  % d(j) ~ d(j+1)
            % find a rotation to eliminate z(j) using z(j+1)
            [c,s] = rotate(z(j+1),z(j));
            z(j) = 0; z(j+1) = normz;
            % update Q
            Q(:,indx(j:j+1)) = Q(:,indx(j:j+1))*[c' s; -s' c];
         else % d(j) does not equal d(j+1)
            inzeros = [inzeros; j];
         end
     end
     j = j+1;
end % while
if (abs(z(n)) > tolz)  		% check the last z(n)
    inzeros = [inzeros; n];
end
