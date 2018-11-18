function [smin,smax,sr,cr,sl,cl]=slasv2(f,g,h)
% [smin,smax,sr,cr,sl,cl] = slasv2(f,g,h)
%
% Compute singular value decomposition of a 2-by-2 
% real upper triangular matrix:
%	[f g] = [cl -sl]*[smax  0  ]*[ cr sr]
%	[0 h]   [sl  cl] [ 0   smin] [-sr cr]
% |smax| is larger singular value and |smin| is
% smaller singular value.

% Reference
% Z.Bai and J.Demmel,
% Computing the Generalized Singular Value Decomposition,
% SIAM J. Sci. Comput., Vol. 14, No. 6, pp. 1464-1486, November 1993
%
% Implemented by S. Qiao	McMaster Univ.	Nov. 1993

ft = f; fa = abs(f);
gt = g; ga = abs(g);
ht = h; ha = abs(h);
pmax = 1;		%pmax points to max abs entry
swap = 0;
if (fa<ha),
    pmax = 3;
    temp = ft;	%swap(ft,ht)
    ft = ht;
    ht = temp;
    temp = fa;	%swap(fa,ha)
    fa = ha;
    ha = temp;
    swap = 1;
end; %if fa<ha

if (ga==0),		%diagonal matrix
    smin = ha;
    smax = fa;
    clt = 1.0; slt = 0.0;
    crt = 1.0; srt = 0.0;
else	%not diagonal
    glarge = 0;	%g is not very large
    if (ga>fa),		%g is the largest entry
	pmax = 2;
	if ((fa/ga)<eps),	%g is very large
	    glarge = 1;
	    smax = ga;			%1 ulp
	    if (ha>1.0)
		smin = fa/(ga/ha);	%2 ulps
	    else
		smin = (fa/ga)*ha;	%2 ulps
	    end; %if ha>1
	    clt = 1.0; slt = ht/gt;
	    crt = 1.0; srt = ft/gt;
	end; %if g very large
    end; %if ga>fa
    if (glarge==0),	%normal case
	d = fa - ha;			%1 ulp
	if (d==fa),	%cope with infinite f or h
	    l = 1.0;			%0 ulp
	else
	    l = d/fa;			%2 ulps
	end; %if d	%note that 0<=l<=1
	m = gt/ft;	%note |m|<1/eps, 1 ulp
	t = 2.0 - l;	%note t>=1,	 3 ulps
	mm = m*m; tt = t*t;
	s = sqrt(tt + mm);		%5 ulps
			%note 1<=s<=1+1/eps
	if (l==0.0),
	    r = abs(m);			%0 ulps
	else
	    r = sqrt(l*l + mm);		%3.5 ulps
	end; %if l	%note 0<=r<=1+1/eps
	a = 0.5*(s + r);		%6 ulps
			%note 1<=a<=1+|m|
	smin = ha/a;			%7 ulps
	smax = fa*a;			%7 ulps

	if (mm==0.0),	%m*m underflow
	    if (l==0.0),
		t = sign(ft)*2*sign(gt);	%0 ulps
	    else
		t = gt/(sign(ft)*d) + m/t;	%6 ulps
	    end; %if l
	else
	    t = (m/(s+t) + m/(r+l))*(1.0 + a);	%17 ulps
	end; %if mm
	l = sqrt(t*t + 4.0);		%18.5 ulps
	crt = 2.0/l;			%19.5 ulps
	srt = t/l;			%36.5 ulps
	clt = (crt + srt*m)/a;		%46.5 ulps
	slt = (ht/ft)*srt/a;		%45.5 ulps
    end; %if glarge
end; %if ga

if (swap==1),
    cl = srt; sl = crt;
    cr = slt; sr = clt;
else
    cl = clt; sl = slt;
    cr = crt; sr = srt;
end; %if swap
%correct the signs of smax and smin
if (pmax==1), tsign= sign(cr)*sign(cl)*sign(f); end;
if (pmax==2), tsign= sign(sr)*sign(cl)*sign(g); end;
if (pmax==3), tsign= sign(sr)*sign(sl)*sign(h); end;
smax = sign(tsign)*smax;
smin = sign(tsign*sign(f)*sign(h))*smin;
