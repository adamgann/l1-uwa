function y = fhmvmul(c,r,w)
% y = fhmvmul(c,r,w)
%
% Fast Hankel matrix-vector multiplication
% 
% Inputs
%    c  first column of the Hankel
%    r  last row of the Hankel
%    w  vector to be multiplied by the Hankel
% Output 
%    y  product of the Hankel and w
% Note. When c and r conflict, the last element of c wins
%       over the first element of r.
%
% This version embeds the Hankel in a Toeplitz circulant
% it also inserts a zero so the circulant is 2nx2n
% it runs faster when n is a power of 2

% S. Qiao       McMaster Univ.  Aug. 1997
% Revised March 2004, simplified code
%
% Reference
% F.T. Luk and S. Qiao.
% A fast singular value algorithm for Hankel matrices.
% Fast Algorithms for Structured Matrices:
% Theory and Applications, Contemporary Mathematics 323,
% Editor V. Olshevsky,
% American Mathematical Society. 2003. 169--177.
%
n = length(r);
%
if (length(c) ~= n)
    error('Error: Col dim does not equal row dim in FHMVMul.m')
end
if (length(w) ~= n)
    error('Error: Inconsistent dim in matrix-vector multiplication FHMVMul')
end
if (c(n) ~= r(1))
    disp('Col wins over row in FHMVMul.m')
end
%
% reverse w
wrev = w(n:-1:1);
% change Hankel to Toeplitz by reversing columns
tc = zeros(n,1);		% initialize two vectors for the
tr = zeros(n,1);		% first col and first row of Toeplitz
tc(1) = c(n); tr(1) = c(n);
for i = 2:n
    tc(i) = r(i);		% first col of Toeplitz
    tr(i) = c(n-i+1);		% first row of Toeplitz
end
%
% embed Toeplitz in a circulant
cc = [tc; 0; tr(n:-1:2)]; 	% first column of circulant
% expand wrev
ww = [wrev; zeros(n,1)];
%
% use fft to multiply circulant matrix-vector
yy = ifft(fft(cc).*fft(ww));

% extract the product
y = yy(1:n);

