function [d,l1,l2] = cstlqd(a,b)
% [d,l1,l2] = cstlqd(a,b)
%
% LQ decomposition of a complex-symmetric and tridiagonal matrix.
% That is
%     diag(a) + diag(b,+1) + diag(b,-1) = L*Q,
% where Q is unitary and L = diag(d) + diag(l1,-1) + diag(l2,-2)
% the lower triangular Cholesky factor.
% This function computes L only.
%
% inputs
%   a,b     diagonal and subdiagonal of the complex-symmetric
%           and tridiagonal matrix to be lower triangularized.
%           No zero entries in b.
% outputs
%   d       diagonal of L, positive
%   l1,l2   first and second subdiagonals of L
% dependency
%   rotate.m   Givens rotation

% S. Qiao   McMaster Univ.    Nov. 2005

n = length(a);
% variable declarations
d = zeros(n,1); l1 = zeros(n-1,1); l2 = zeros(n-2,1);

d(1) = a(1); l1(1) = b(1);
for i=1:n-2
    % eliminate b(i) using d(i)
    [c,s] = rotate(d(i),b(i));
    d(i) = sqrt(real(d(i)*d(i)' + b(i)*b(i)'));
    % update columns i and i+1
    d(i+1) = c*a(i+1) - s*l1(i);
    l1(i) = c'*l1(i) + s'*a(i+1);
    l2(i) = s'*b(i+1);
    l1(i+1) = c*b(i+1);
end
[c,s] = rotate(d(n-1),b(n-1));
d(n-1) = sqrt(real(d(n-1)*d(n-1)' + b(n-1)*b(n-1)'));
d(n) = abs(c*a(n) - s*l1(n-1));
l1(n-1) = c'*l1(n-1) + s'*a(n);
