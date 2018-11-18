function s = cssing(a,b)
% function s = cssing(a,b)
%
% singular values of a complex symmetric tridiagonal T
% given by its diagonals
%
% Inputs
%    a  main diagonal
%    b  subdiagonal
% Outputs
%    s  singular values, nonnegative, ordered
%
% Dependency
%   ./csssingstep.m (one step of implicit QR, singular values only)

% S.Qiao	McMaster Univ	January 2001
% Revised June 2001
%         July 2002 (condition for setting b(i) to 0)
%         June 2003 (add the trivial case when length(a) is one)
%         Nov  2005 (modified for computing singular values only)

% Reference
% G.H. Golub and C.F. Van Loan,
% Matrix Computations, 2nd Ed,
% The Johns Hopkins University Press, 1989,
% pp.423--424
%
% ~/papers/takagi/pgms/cssvd.m

n = length(a);
%
if n == 1		% quick return
    s = abs(a(1));
    return
end
%
% initial s, and Q
s = a;
%
q = 0;			% size of the diagonalized block
%
while q < n-1,		% while not diagonal
    for i=1:n-q-1,
	% set small subdiagonal elements to zero
        tmp = abs(s(i)) + abs(s(i+1));
        if i>1 tmp = tmp + abs(b(i-1)); end
        if i<n-1 tmp = tmp + abs(b(i+1)); end
        if abs(b(i)) <= tmp*eps
	    b(i) = 0.0;
	end
    end
%
%   update the diagonal block
    while q < n-1
	if b(n-q-1) == 0.0
	    q = q+1;
	else
	    break
	end
    end
%
%   find p so that S(p+1:n-q,p+1:n-q) is the largest block
%   with no zero subdiagonal elements
    p = n - q - 1;     % starting value of p
    while p > 0
	if b(p) ~= 0.0
	    p = p-1;
	else
	    break
	end
    end
%
    if q < n-1
%	do one step of the implicit QR method with shifting
        [s(p+1:n-q), b(p+1:n-q-1)] = cssingstep(s(p+1:n-q), b(p+1:n-q-1));
    end
 end
 
% make s nonnegative and in descending order
s = sort(-abs(s));
s = -s;
