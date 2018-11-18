function [c,s]= rotate(a, b)
% [c,s] = rotate(a,b)
%
% Generate cosine and sine of Givens rotation
% such that
%	( c' s')(a) = (*)
% 	(-s  c )(b)   (0)
% INPUT: 	a, b  complex numbers
% OUTPUT:	c, s  cosine and sine of the
%		rotation which eliminates b

% Reference:
% G.Golub and C.Van Loan, Matrix Computations,
% 2ed Edition, The John Hopkins University Press,
% 1989, p.202.
%
% Implemented by S. Qiao, McMaster University

if b==0,
    if (a==0),
	c = 1.0;
    else
	c = sign(a);
    end;
    s = 0;
    return;
end;
if a==0,
    c = 0;
    s = sign(b);
    return;
end;
if (abs(a)>abs(b)),
    t= b/a;
    tt= sqrt(1 + t*t');
    c= sign(a)/tt;
    s= t*c;
else
    t= a/b;
    tt= sqrt(1 + t'*t);
    s= sign(b)/tt;
    c= t*s;
end;
