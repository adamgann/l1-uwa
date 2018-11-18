function [a,b,Q] = cssvdstep(aIn,bIn,Q)
% [a,b,Q] = cssvdstep(aIn,bIn,Q)
%
% One QR step in complex symmetric SVD using second degree
% class two shift.
%
% Inputs
%   aIn    main diagonal of the input tridiagonal matrix
%   bIn    subdiagonal of the input tridiagonal matrix
%   Q      unitary matrix to be updated, optional
% Outputs
%   a      main diagonal of the output matrix
%   b      subdiagonal of the output matrix
%   Q      updated unitary matrix
% Inputs and outputs satisfy
%   Qin*Tin*conj(Qin') = Qout*T*conj(Qout')
%   where TIn = diag(aIn) + diag(bIn,1) + diag(bIn,-1)
%   and     T = diag(a) + diag(b,1) + diag(b,-1)
%           is closer to diagonal

% Dependency
%   cssvd2 (2x2 complex symmetric svd)
%   house (Householder transformation)
%   shift22 (compute second degree class two shift)

% S. Qiao     McMaster Univ.  July 2000
% Revised     Juanuary 2001
%             March 2001
%             June 2002 (call shift22)
%
% References
% F.T. Luk and S. Qiao.
% A fast singular value algorithm for Hankel matrices.
% Fast Algorithms for Structured Matrices:
% Theory and Applications, Contemporary Mathematics 323,
% Editor V. Olshevsky,
% American Mathematical Society. 2003. 169--177.
%
% A. Bunse-Gerstner and W.B. Gragg.
% Singular value decompositions of complex symmetric matrices.
% Journal of Computational and Applied Mathematics, 21(1988) 41--54.

a = aIn; b = bIn;
n = length(a);
if (n ~= length(b) + 1)
   error('Error: Incorrect length(b) input to cssvdstep');
end
%
if n==2 		% trivial case: n=2
    if (nargin>2)
        [s,QQ] = cssvd2(a(1),a(2),b(1));
        Q = Q*QQ;
    end
    if (nargin<3)
        s = cssvd2(a(1),a(2),b(1));
    end
    a = s; b = 0.0;
    return
end
%  
c = zeros(n-2,1); d = zeros(n-3,1);  % two temporary subdiagonals
%
shift = shift22(a, b);
%
for k=1:n-2
   if k==1            % first column of (Tin'*Tin - shift*eye)
      temp = [a(1)*a(1)' + b(1)*b(1)' - shift;
              a(1)*b(1)' + a(2)'*b(1);
              b(2)'*b(1)];
   else               % for restoring tridiagonal
      temp = [b(k-1); c(k-1); d(k-1)];
   end
   [u,beta,alpha] = house(temp);	% Householder matrix
   				        % to eliminate temp(2:3)
   if k==1
      u = conj(u); beta = conj(beta);
   end
   H3 = eye(3) - (beta*u)*u';
%
   if k>1				% restore tridiagonal
      b(k-1) = -alpha; c(k-1) = 0; d(k-1) = 0;
   end
   % update [a(k)  b(k)  c(k)    d(k) ]
   %        [b(k) a(k+1) b(k+1) c(k+1)]
   %        [c(k) b(k+1) a(k+2) b(k+2)]
   %        [d(k) c(k+1) b(k+2)  ...  ]
   if k<n-2
      temp = H3*[d(k); c(k+1); b(k+2)];
      d(k) = temp(1); c(k+1) = temp(2); b(k+2) = temp(3);
   end
   					% update a block
   block = [a(k) b(k)   c(k);
            b(k) a(k+1) b(k+1);
            c(k) b(k+1) a(k+2)];
   block = H3*block*conj(H3');
   a(k) = block(1,1); b(k) = block(1,2); c(k) = block(1,3);
   a(k+1) = block(2,2); b(k+1) = block(2,3);
   a(k+2) = block(3,3);
   if (nargin>2)		% update Q
      Q(:,k:k+2) = Q(:,k:k+2) - (Q(:,k:k+2)*u)*(beta*u');
   end
end
%
% last step in tridiagonalization
% eliminate c(n-2) in
%       [a(n-2) b(n-2) c(n-2)]
%       [b(n-2) a(n-1) b(n-1)]
%       [c(n-2) b(n-1)  a(n) ]
[u,beta,alpha] = house([b(n-2); c(n-2)]);
H2 = eye(2) - (beta*u)*u';
b(n-2) = -alpha; c(n-2) = 0;
%
% update [a(n-1) b(n-1)]
%        [b(n-1)  a(n) ]
block = [a(n-1) b(n-1); b(n-1) a(n)];
block = H2*block*conj(H2); 
a(n-1) = block(1,1); b(n-1) = block(1,2); a(n) = block(2,2);
%
if (nargin>2)
   Q(:,n-1:n) = Q(:,n-1:n) - (Q(:,n-1:n)*u)*(beta*u');
end

