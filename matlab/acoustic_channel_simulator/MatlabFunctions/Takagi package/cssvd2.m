function [s,Q] = cssvd2(a11,a22,a12)
% [s,Q] = cssvd2(a11,a22,a12)
%
% Compute Takagi factorization of a 2x2 complex symmetric matrix
%
%        [a11 a12] = Q*diag(s)*conj(Q')
%        [a12 a22]
%    where Q is unitary
%
% Inputs
%   a11,a22,a12    entries of the 2x2 symmetric matrix
% Outputs
%   s              singular values, nonnegative and ordered
%   Q              unitary, optional
%
% Dependency: csvd2.m

% References
% Sanzheng Qiao and Xiaohong Wang.
% Computing the Singular Vaues of 2-by-2 Complex Matrices.
% Software Quality Research Laboratory Report No. 5,
% Department of Computing and Software, McMaster University,
% Hamilton, Ontario, Canada. L8S 4L7. June 2002.
%
% S. Qiao         McMaster Univ., January 2001
% Revised June, August 2001
%         March 2002 (made 2x2 eigen decomp explicit)
%         May 2002, use complex 2x2 svd csvd2.m in place of
%                   eigen decomp

A = [a11 a12; a12 a22];

% svd of the 2x2 complex-symmetric matrix
[V,S,VV] = csvd2(a11, a12, a12, a22);
s = diag(S);		% singular values
 
if nargout>1
   % find y such that A*conj(y) = alpha*y
   x = A*conj(V(:,1));

   if 1/cond([x V(:,1)]) < eps*norm([x V(:,1)])
           % x and eigenvector V(:,1) are linearly dependent
      y = V(:,1);
   else    % x and eigenvector V(:,1) are linearly independent
      y = x + s(1)*V(:,1);
      y = y/norm(y);
   end
   % construct unitary Q such that Q(:,1) = y
   Q = [y(1) -y(2)'; y(2) y(1)'];
   s = diag(Q'*A*conj(Q));
end % if nargout>1

% now we have Takagi factorization A = Q*diag(s)*conj(Q')
% the following steps make s nonnegative and ordered so that
% the factorization is consistent with SVD

if nargout>1
   % make the singular values nonnegative
   for i=1:2
      if s(i)~=0
         Q(:,i) = sqrt(sign(s(i)))*Q(:,i);
      end
   end
end
s = abs(s);

% order the singular values
if s(2)>s(1)
   s = [s(2); s(1)];
   if nargout>1 Q = [Q(:,2) Q(:,1)]; end
end

