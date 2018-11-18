function [s,Q] = CSSVD(a,b)
% function [s,Q] = CSSVD(a,b)
%
% SVD of a complex symmetric tridiagonal T given by its diagonals
%
% Inputs
%    a  main diagonal
%    b  subdiagonal
% Outputs
%    s  singular values, nonnegative, ordered
%    Q  unitary, optional
% so that
%        Q*diag(s)*conj(Q') = diag(a) + diag(b,1) + diag(b,-1)	
%
% Dependency
%   ./cssvdstep.m (one step of implicit QR)

% S.Qiao	McMaster Univ	January 2001
% Revised June 2001
%         July 2002 (condition for setting b(i) to 0)
%         June 2003 (add the trivial case when length(a) is one)
%         August 2005 (updating Q is moved inside cssvdstep for
%                      for efficiency)

% References
% G.H. Golub and C.F. Van Loan,
% Matrix Computations, 2nd Ed,
% The Johns Hopkins University Press, 1989,
% pp.423--424
%
% F.T. Luk and S. Qiao.
% A fast singular value algorithm for Hankel matrices.
% Fast Algorithms for Structured Matrices:
% Theory and Applications, Contemporary Mathematics 323,
% Editor V. Olshevsky,
% American Mathematical Society. 2003. 169--177.

n = length(a);
% check the length of b
if n ~= length(b) + 1
    error('Error: Incorrect input b to CSSVD.')
end
%
if n == 1		% quick return
    s = abs(a(1));
    Q = 1.0;
    if s ~= 0
        Q = sqrt(a(1)/s);
    end
    return
end
%
% initial s, and Q
s = a;
if (nargout>1) Q = eye(n); end;
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
	if (nargout>1)
	    [tmps,tmpb,Q(:,p+1:n-q)] = ...
            cssvdstep(s(p+1:n-q), b(p+1:n-q-1),Q(:,p+1:n-q));
	else	% nargout=1
	    [tmps, tmpb] = cssvdstep(s(p+1:n-q), b(p+1:n-q-1));
	end;
%	update s, b
	s(p+1:n-q) = tmps; b(p+1:n-q-1) = tmpb;
    end
 end
 
% now we have Takagi factorization T = Q*diag(s)*conj(Q')
% the following steps make s nonnegative and in descending order
% so that the factorization is consistent with SVD
tmp = s;
s = abs(tmp);
if (nargout>1)
    for i=1:n
       if s(i)~=0.0
          Q(:,i) = sqrt(tmp(i)/s(i))*Q(:,i);
       end
    end
end
for i=1:n
    for k=i+1:n
        if s(k)>s(i)
            tmp = s(i);
            s(i) = s(k);
            s(k) = tmp;
            Q(:,[i,k]) = Q(:,[k,i]);
        end
    end
end

