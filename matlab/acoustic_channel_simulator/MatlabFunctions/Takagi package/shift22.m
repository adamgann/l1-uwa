function s = shift22(a, b)
% s = shift22(a, b)
%
% Compute the second order class two shift in the Takagi factorization
% of a tridiagonal and complex-symmetric matrix.
% Inputs
%   a    main diagonal of the input tridiagonal matrix
%   b    subdiagonal of the input tridiagonal matrix
% Output
%   s    second order class two shift
%
% Dependency
%   ./rotate.m (rotation transformation)

% Reference
%   Angelika Bunse-Gerstner and William B. Gragg,
%   Singular value decomposition of complex symmetric matrices,
%   Journal of Computational and Applied Mathematics 21(1988) 41-54,
%   Page 49 (second order class two shift)
%
%   G.H. Golub and C.F. Van Loan,
%   Matrix Computations, 3rd Ed,
%   The Johns Hopkins University Press, 1996.
%   Page 418 (Wilkinson shift)

% S. Qiao    McMaster Univ.  June 2002
% Nov. 2005, modified for efficiency

n = length(a);

% trailing 3x3 submatrix T3'*T3 of the input matrix, where
% if n==3
%    T3 = [a(n-2) b(n-2)   0;
%          b(n-2) a(n-1) b(n-1);
%            0    b(n-1)  a(n)  ]
% n>3
%    T3 = [b(n-3)    0      0;
%          a(n-2) b(n-2)    0;
%          b(n-2) a(n-1) b(n-1);
%             0   b(n-1)  a(n) ];
T11 = real(conj(a(n-2))*a(n-2)) + real(conj(b(n-2))*b(n-2));
if n>3
    T11 = T11 + real(conj(b(n-3))*b(n-3));
end
T12 = conj(a(n-2))*b(n-2) + conj(b(n-2))*a(n-1);
T13 = conj(b(n-2))*b(n-1);
T22 = real(conj(b(n-2))*b(n-2)) + real(conj(a(n-1))*a(n-1)) ...
          + real(conj(b(n-1))*b(n-1));
T23 = conj(a(n-1))*b(n-1) + conj(b(n-1))*a(n);
T33 = real(conj(b(n-1))*b(n-1)) + real(conj(a(n))*a(n));

% determine a Givens rotation
% Q = [c -s; s' c'] eliminating T13 using T23:
% Q*[T13] = [0]
%   [T23] = [*]
[c,s] = rotate(T23, T13);
% premultiply Q, unnecessary to update T11, T12, T13
T21 = conj(s)*T11 + conj(c*T12);
T22 = conj(s)*T12 + conj(c)*T22;
T23 = conj(s)*T13 + conj(c)*T23;
% postmultiply Q'
T22 = real(T21*s + T22*c);
% calculate Wilkinson shift using
% [T22  T23]
% [T32  T33]
t = (T22 - T33)/2;
if (t~=0.0)
    s = T33 + t - sign(t)*sqrt(t*t + real(conj(T23)*T23));
else % t=0.0
    if (T33>=0.0)
        s = T33 + abs(T23);
    else
        s = T33 - abs(T23);
    end
end

