function [a,b] = cssingstep(a,b)
% [a,b] = cssingstep(a,b)
%
% One QR step in complex symmetric SVD, computes singular values only.
% Inputs
%   a    main diagonal of the input tridiagonal matrix Tin
%   b    subdiagonal of the input tridiagonal matrix Tin
% Outputs
%   a    updated main diagonal
%   b    updated subdiagonal
% Tout = diag(a) + diag(b,1) + diag(b,-1) is closer to diagonal than Tin

% Dependency
%   cssvd2 (2x2 complex symmetric svd)
%   shift22 (compute second degree class two shift)

% S. Qiao     McMaster Univ.  July 2000
% Revised     Juanuary 2001
%             March 2001
%             June 2002 (call shift22)
%             Nov 2005, modified to compute singular values only
%                       improved efficiency

n = length(a);
%
if n==2 		% trivial case: n=2
    s = cssvd2(a(1),a(2),b(1));
    a = s; b = 0.0;
    return
end
%  
c = zeros(n-2,1); d = zeros(n-3,1);  % two temporary subdiagonals
%
% find second degree class two shift
shift = shift22(a, b);
%
% determine the first Householder matrix in the QRD of Tin'*Tin - shift*I
% u = [u1; u2; u3], first column of Tin'*Tin - shift*I
u1 = conj(a(1))*a(1) + conj(b(1))*b(1) - shift;
u2 = conj(b(1))*a(1) + conj(a(2))*b(1);
u3 = conj(b(2))*b(1);
% [u,beta,alpha] = house(u), Householder matrix to eliminate u2, u3
alpha = sqrt(conj(u1)*u1 + conj(u2)*u2 + conj(u3)*u3);
if (u1~=0.0)
    alpha = sign(u1)*alpha;
end
u1 = u1 + alpha;
beta = 1/(conj(alpha)*u1);
% premultiply
%        [ 0  ]
%        [ 0  ]
%        [b(3)]
% with (I - beta*u*u').' when n>3
if n>3
    ts = -beta*u3*b(3);  % -beta*u.'*[0;0;b(3)]
    d(1) = ts*conj(u1);
    c(2) = ts*conj(u2);
    b(3) = b(3) + ts*conj(u3);
end
% premultiply and postmultiply
% B = [a(1) b(1)  0;
%      b(1) a(2) b(2);
%       0   b(2) a(3)]
% with (I - beta*u*u').' and I - beta*u*u' respectively

% tv = u.'*B
tv1 = u1*a(1) + u2*b(1);
tv2 = u1*b(1) + u2*a(2) + u3*b(2);
tv3 = u2*b(2) + u3*a(3);
% ts = beta*beta*tv.'*u
ts = beta*beta*(tv1*u1 + tv2*u2 + tv3*u3);
% B = B - beta*(conj(u)*tv.' + tv*u') + ts*conj(u)*u'
u1 = conj(u1); u2 = conj(u2); u3 = conj(u3);
a(1) = a(1) - 2*beta*u1*tv1 + ts*u1*u1;
b(1) = b(1) - beta*(u1*tv2 + tv1*u2) + ts*u1*u2;
c(1) = -beta*(u1*tv3 + tv1*u3) + ts*u1*u3;
a(2) = a(2) - 2*beta*tv2*u2 + ts*u2*u2;
b(2) = b(2) - beta*(u2*tv3 + tv2*u3) + ts*u2*u3;
a(3) = a(3) - 2*beta*tv3*u3 + ts*u3*u3;

% restore tridiagonal structure
for k=2:n-3
    % determine a Householder matrix to eliminate c(k-1) and d(k-1)
    % [u,beta,alpha] = house([b(k-1); c(k-1); d(k-1)]);
    u1 = b(k-1); u2 = c(k-1); u3 = d(k-1);
    alpha = sqrt(conj(u1)*u1 + conj(u2)*u2 + conj(u3)*u3);
    if (u1~=0.0)
        alpha = sign(u1)*alpha;
    end
    u1 = u1 + alpha;
    beta = 1/(conj(alpha)*u1);
    b(k-1) = -alpha; c(k-1) = 0; d(k-1) = 0;
    % premultiply
    %        [ d(k);
    %        [c(k+1);
    %        [b(k+2)]
    % with I - beta*u*u'
    
    % ts = beta*u'*[d(k);c(k+1);b(k+2)]
    ts = beta*(conj(u1)*d(k) + conj(u2)*c(k+1) + conj(u3)*b(k+2));
    d(k) = d(k) - ts*u1;
    c(k+1) = c(k+1) - ts*u2;
    b(k+2) = b(k+2) - ts*u3;
    % premultiply and postmultiply
    %    B = [a(k)  b(k)   c(k);
    %        [b(k) a(k+1) b(k+1);
    %        [c(k) b(k+1) a(k+2)]
    % with I - beta*u*u' and (I - beta*u*u').' respectively

    % tv = u'*B
    tv1 = conj(u1)*a(k) + conj(u2)*b(k) + conj(u3)*c(k);
    tv2 = conj(u1)*b(k) + conj(u2)*a(k+1) + conj(u3)*b(k+1);
    tv3 = conj(u1)*c(k) + conj(u2)*b(k+1) + conj(u3)*a(k+2);
    % ts = beta*beta*tv.'*conj(u);
    ts = beta*beta*(tv1*conj(u1) + tv2*conj(u2) + tv3*conj(u3));
    % B = B - beta*(u*tv.' + tv*u.') + ts*u*u.'
    a(k) = a(k) - 2*beta*(u1*tv1) + ts*u1*u1;
    b(k) = b(k) - beta*(u1*tv2 + tv1*u2) + ts*u1*u2;
    c(k) = c(k) - beta*(u1*tv3 + tv1*u3) + ts*u1*u3;
    a(k+1) = a(k+1) - 2*beta*(u2*tv2) + ts*u2*u2;
    b(k+1) = b(k+1) - beta*(u2*tv3 + tv2*u3) + ts*u2*u3;
    a(k+2) = a(k+2) - 2*beta*(u3*tv3) + ts*u3*u3;
end % for k=2:n-3
%
% step n-2
if n > 3
    % determine a Householder matrix to eliminate c(n-3) and d(n-3)
    u1 = b(n-3); u2 = c(n-3); u3 = d(n-3);
    alpha = sqrt(conj(u1)*u1 + conj(u2)*u2 + conj(u3)*u3);
    if (u1~=0.0)
        alpha = sign(u1)*alpha;
    end
    u1 = u1 + alpha;
    beta = 1/(conj(alpha)*u1);
    b(n-3) = -alpha; c(n-3) = 0; d(n-3) = 0;
    % premultiply and postmultiply
    %        [a(n-2) b(n-2) c(n-2)]
    %        [b(n-2) a(n-1) b(n-1)]
    %        [c(n-2) b(n-1)  a(n) ]
    % with I - beta*u*u' and (I - beta*u*u').' respectively

    % tv = u'*B
    tv1 = conj(u1)*a(n-2) + conj(u2)*b(n-2) + conj(u3)*c(n-2);
    tv2 = conj(u1)*b(n-2) + conj(u2)*a(n-1) + conj(u3)*b(n-1);
    tv3 = conj(u1)*c(n-2) + conj(u2)*b(n-1) + conj(u3)*a(n);
    % ts = beta*beta*tv.'*conj(u)
    ts = beta*beta*(tv1*conj(u1) + tv2*conj(u2) + tv3*conj(u3));
    % B = B - beta*(u*tv.' + tv*u.') + ts*u*u.'
    a(n-2) = a(n-2) - 2*beta*(u1*tv1) + ts*u1*u1;
    b(n-2) = b(n-2) - beta*(u1*tv2 + tv1*u2) + ts*u1*u2;
    c(n-2) = c(n-2) - beta*(u1*tv3 + tv1*u3) + ts*u1*u3;
    a(n-1) = a(n-1) - 2*beta*(u2*tv2) + ts*u2*u2;
    b(n-1) = b(n-1) - beta*(u2*tv3 + tv2*u3) + ts*u2*u3;
    a(n) = a(n) - 2*beta*(u3*tv3) + ts*u3*u3;
end % if n>3
% last step in tridiagonalization
% determine a Householder matrix to eliminate c(n-2)
u1 = b(n-2); u2 = c(n-2);
% [u,beta,alpha] = house(u);
alpha = sqrt(conj(u1)*u1 + conj(u2)*u2);
if (u1~=0.0)
    alpha = sign(u1)*alpha;
end
u1 = u1 + alpha;
beta = 1/(conj(alpha)*u1);
b(n-2) = -alpha; c(n-2) = 0;;
% premultiply and postmultiply
%    B = [a(n-1) b(n-1)]
%        [b(n-1)  a(n) ]
% with I - beta*u*u' and (I - beta*u*u').' respectively

% tv = u'*B
tv1 = conj(u1)*a(n-1) + conj(u2)*b(n-1);
tv2 = conj(u1)*b(n-1) + conj(u2)*a(n);
% ts = beta*beta*tv.'*conj(u)
ts = beta*beta*(tv1*conj(u1) + tv2*conj(u2));
% B = B - beta*(u*tv.' + tv*u.') + ts*u*u.'
a(n-1) = a(n-1) - 2*beta*(u1*tv1) + ts*u1*u1;
b(n-1) = b(n-1) - beta*(u1*tv2 + tv1*u2) + ts*u1*u2;
a(n) = a(n) - 2*beta*(u2*tv2) + ts*u2*u2;

