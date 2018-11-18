function [u,beta,alpha]= house(x)
% [u,beta,alpha] = house(x)
%
% Given a nonzero vector x, this function computes vector u
% and scalars beta and alpha such that the Householder matrix
%    H = eye(length(x)) - beta*u*u'
% satisfies
%          [-alpha]
%          [  0   ]
%    H*x = [ ...  ]
%          [  0   ]
% and
% u is scaled so that u(1)=1, so u(2:n) is the essential part.

% References:
% G.W. Stewart, Introduction to Matrix
% Computations, Academic Press, Inc., 1973.
% Theorem 3.6 (p.234) and Ex.19 (p.262).
%
% G.H. Golub and C.F. Van Loan,
% Matrix Computations, Third Edition,
% The John Hopkins University Press, 1996.
% pp. 209.
%
%
% Implemented by
% S. Qiao	McMaster Univ.  July 2000
% revised March 2001
%         April 2001 (added quick return)
%         May 2002 (removed quick return, caused inaccuracy)

[nr,nc] = size(x);
if nc~=1
   error('Input x to house is not a column vector.')
end
% scaling to avoid over(under)flow
m = max(abs(x));
x = x/m;
%
normx = norm(x);	% normx >= 1.0
%
% get sign(x1)
if x(1)==0		% in matlab sign(0)=0
    signx1 = 1.0;
else
    signx1 = sign(x(1));
end
%
alpha = signx1*normx;   % |alpha| >= 1.0
%
u = x;
u(1)= u(1) + alpha;     % |u(1)| >= 1.0
beta = 1/(alpha'*u(1));
%
% make u(1)=1, so u(2:n) is the essential part of u
beta = beta*(u(1)*u(1)');	% update beta
u = u/u(1);
% scale back
alpha = m*alpha;
