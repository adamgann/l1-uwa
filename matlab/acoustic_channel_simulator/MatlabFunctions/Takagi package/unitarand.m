function U = unitarand(m,n)
% U = unitarand(m,n)
% Generate an m-by-n random unitary matrix
%
% S. Qiao	McMaster University
% revised March 2001

% Reference:
% G.W.Stewart, "The efficient generation of random orthogonal
%		matrices with an application to condition estimators",
%		SIAM J. Numer. Anal. vol. 17, No. 3, June 1980,
%		pp.403-409.
%
% Dependency: house (Householder transformation)

% initialize a random diagonal unitary matrix
x= (ones(m,1) - 2*rand(m,1)) + sqrt(-1)*(ones(m,1) - 2*rand(m,1));
U= diag(sign(x));
%
for i=m:-1:max(2,m-n+1),
    % generate an i-dim. random vector with
    % uniform distribution on [-1, 1]
    x= (ones(i,1) - 2*rand(i,1)) + sqrt(-1)*(ones(i,1) - 2*rand(i,1));
    % generate Householder matrix
    [u,beta,alpha]= house(x);
    % accumulate Householder transformations
    U(:,m-i+1:m)= U(:,m-i+1:m) - (U(:,m-i+1:m)*u)*(beta*u');
end;
U= U(:,1:n);
