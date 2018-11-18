function z = twist(d,l1,l2,mu)
% z = twist(d,l1,l2,mu)
%
% Given the Cholesky factorization of a Hermitian positive definite
% pentadiagonal matrix and a computed eigenvalue,
% compute the corresponding eigenvector using the twisted factorization. 
% inputs
%     d -- diagonal of Cholesky factor
% l1,l2 -- first and second subdiagonals of Cholesky factor 
%    mu -- computed eigenvalue
% output
%     z -- eigenvector corresponding to mu
%  dependency
%   cqds.m
%
% Nov. 2005

n = length(d);

% constant
tiny = 1.0e-12;    % tolerance for small gamma
% variable declaration
z = zeros(n,1);	   % eigenvector
%
% find quotient-differences with shift mu
[dl,ll1,ll2,du,u1,u2] = cqds(d,l1,l2,mu);
% 
% for j=1,..,n, compute the twisted factorization
%        L*L' - mu*I = Nj*Dj*Nj'
% where L = diag(d) + diag(l1,-1) + diag(l2,-2),
% Dj = diag([dl(1:j-2);x;gamma;du(j+1:n)]), and
% Nj = eye(n) + diag([ll1(1:j-2);yj;zeros(n-j,1)],-1) ...
%             + diag([ll2(1:j-2);zeros(n-j,1)],-2) ...
%             + diag([zeros(j-1,1);u1(j:n-1)],+1) ...
%             + diag([zeros(j-2,1);u2(j-1:n-2)],+2),
% and find the index of the min(abs(gamma)) and the corresponding yj. 
% j=1
gmin = abs(du(1));		   % abs(gamma(1))
indx = 1;                          % index of minimum abs(gamma)
% j=2
x = dl(1) - real(conj(u2(1))*u2(1))*du(3);
y = (ll1(1)*dl(1) - u1(2)*conj(u2(1))*du(3))/x;
gamma = du(2) - x*real(conj(y)*y);  % gamma(2)
agamma = abs(gamma);
if gmin > agamma 
    gmin = agamma;
    indx = 2;
end
for j=3:n-1
    x = dl(j-1) - real(u2(j-1)*conj(u2(j-1)))*du(j+1);
    yj = (ll1(j-1)*dl(j-1) - u1(j)*conj(u2(j-1))*du(j+1))/x;
    gamma = du(j) - real(conj(ll2(j-2))*ll2(j-2))*dl(j-2) ...
            - real(conj(yj)*yj)*x;   % gamma(j)
    agamma = abs(gamma);
    if gmin > agamma 
        gmin = agamma;
        indx = j;
        y = yj;
        if gmin < tiny              % found a small gamma
            break;
        end
    end
end
% j=n
gamma = dl(n);                      % gamma(n)
agamma = abs(gamma);
if gmin > agamma
    gmin = agamma;
    indx = n;
    y = ll1(n-1);
end
%
% the eigenvector z corresponding to mu is the solution for
%                     N_{indx}'*z = e_{indx}.
z(indx) = 1;
if indx > 1
    z(indx-1) = -conj(y);
end
if indx < n
    z(indx+1) = -conj(u1(indx));
    if indx > 1
        z(indx+1) = conj(u2(indx-1)*y) + z(indx+1);
    end
end
for j = indx-2:-1:1
    z(j) = -conj(ll1(j))*z(j+1) - conj(ll2(j))*z(j+2);
end
for j = indx+2:n
    z(j) = -conj(u2(j-2))*z(j-2) - conj(u1(j-1))*z(j-1);
end
z = z/norm(z);       % normalize
