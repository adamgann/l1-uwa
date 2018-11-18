function [A,U] = csgen(s)
% [A,U] = csgen(s)
% generate a random complex symmetric matrix
% with specified singular values
%
% Input
%    s		given singular values
% Output
%    A          random complex symmetric matrix with singular values s
%    U          unitary matrix such that
%                    A = U*diag(s)*conj(U')
%
% S. Qiao   McMaster University   March 2001

% dependency:
%    ./unitarand (randome unitary matrix generator)

% singular values must be real.
% Matlab relational operator "<" does compare complex numbers,
% but only compares their real parts
if (any(imag(s) ~= 0))
    error('Singular values must be real.')
end

if (any(s < 0))
    error('Singular values must be nonnegative.')
end

n = length(s);
% generate a random unitary matrix
U = unitarand(n,n);
% form the complex symmetric matrix using s as singular values
A = U*diag(s)*conj(U');
% force symmetry
A = triu(A) + tril(conj(A'),-1);
