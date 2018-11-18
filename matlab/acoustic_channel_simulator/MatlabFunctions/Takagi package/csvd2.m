function [U,S,V] = csvd2(a11, a12, a21, a22)
% [U,S,V] = csvd2(a11, a12, a21, a22)
%
% SVD of a 2-by-2 complex matrix:
%   [a11  a12] = U*S*V'
%   [a21  a22]
%
% Dependencies: ./rotate.m, ./slasv2.m
%
% Reference
% Sanzheng Qiao and Xiaohong Wang.
% Computing the Singular Vaues of 2-by-2 Complex Matrices.
% Software Quality Research Laboratory Report No. 5,
% Department of Computing and Software, McMaster University,
% Hamilton, Ontario, Canada. L8S 4L7. June 2002.
%
% S. Qiao   McMaster Univ.  May 2002

% reduce the matrix into real upper triangular
%   [a11  a12] = U*[f  g]*V'
%   [a21  a22]     [0  h]
[c,s] = rotate(a11, a21); 	% eliminate a21
U = [c -s'; s c'];
f = c'*a11 + s'*a21;
tmp = U'*[a12; a22];		% update the second column

V = eye(2);			% make a12 real
if (tmp(1)~=0.0)
    V(2,2) = tmp(1)'/abs(tmp(1));
    tmp(2) = tmp(2)*(tmp(1)'/abs(tmp(1)));
end
g = abs(tmp(1));

if (tmp(2)~=0.0)		% make a22 real
   U(:,2) = U(:,2)*(tmp(2)/abs(tmp(2)));
end
h = abs(tmp(2));

% SVD of a real upper triangular matrix
%   [f  g] = [cl -sl]*[smax  0  ]*[cr -sr]'
%   [0  h]   [sl  cl] [ 0   smin] [sr  cr]
[smin,smax,sr,cr,sl,cl] = slasv2(f,g,h);
S = [smax 0; 0 smin];
U = U*[cl -sl; sl cl];		% update U and V
V = V*[cr -sr; sr cr];
