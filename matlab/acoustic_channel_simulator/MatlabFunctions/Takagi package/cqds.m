function [dl,ll1,ll2,du,u1,u2] = cqds(d,l1,l2,mu)
% [dl,ll1,ll2,du,u1,u2] = cqds(d,l1,l2,mu)
%
% Given Cholesky factorization, compute complex quotient-differences
% with shift, that is,
%    L*L' - mu*I = LL*DL*LL' = U*DU*U'
% where
%     L = diag(d) + diag(l1,-1) + diag(l2,-2)
% is the Cholesky factor and
%    LL = eye(n) + diag(ll1,-1) + diag(ll2,-2)
%    DL = diag(dl)
%     U = eye(n) + diag(u1,+1) + diag(u2,+2)
%    DU = diag(du).
%
%  Inputs
%      d -- diagonal of the Cholesky factor L, positive
%  l1,l2 -- first and second subdiagonals of L
%     mu -- shift, a real scalar
%  Outputs
%      dl -- diagonal of DL, real
% ll1,ll2 -- first and second subdiagonals of LL
%      du -- diagonal of DU, real
%   u1,u2 -- first and second superdiagonals of U

% Nov. 2005

n = length(d);
%%%%%%%%%%% compute the lower triangular part LL %%%%%%%%%%%%
% variable declarations
dl = zeros(n,1);
ll1 = zeros(n-1,1);
ll2 = zeros(n-2,1);
% temporary variables
lt = zeros(n-1,1);     % lt(i) is the (i+1,i)-entry of L*L - mu*I
dt = zeros(n,1);       % dt(i) is the (i,i)-entry of L*L - mu*I
if mu == 0  % trivial case when mu=0
    dl = d.*d;
    ll1 = l1./d(1:n-1);
    ll2 = l2./d(1:n-2);
else % mu ~= 0
    dt(1) = d(1)*d(1) - mu;      % (1,1)-entry of L*L' - mu*I
    dl(1) = dt(1);               % equating (1,1)-entries of
                                 % L*L - mu*I and LL*DL*LL'
    lt(1) = l1(1)*d(1);          % (2,1)-entry of L*L' - mu*I
    ll1(1) = lt(1)/dl(1);        % equating (2,1)-entries
    dt(2) = real(l1(1)*conj(l1(1))) + d(2)*d(2) - mu;
                                 % (2,2)-entry of L*L' - mu*I
    dl(2) = dt(2) - real(ll1(1)*conj(ll1(1)))*dl(1);
                                 % equating (2,2)-entries
    for i = 1:n-2
        ll2(i) = l2(i)*d(i)/dl(i);
                                 % equating (i+2,i)-entries
        lt(i+1) = l2(i)*conj(l1(i)) + l1(i+1)*d(i+1);
                                 % (i+2,i+1)-entry of L*L' - mu*I
        ll1(i+1) = (lt(i+1) - ll2(i)*dl(i)*conj(ll1(i)))/dl(i+1);
                                 % equating (i+2,i+1)-entries
        dt(i+2) = real(l2(i)*conj(l2(i))) + real(l1(i+1)*conj(l1(i+1))) ...
                  + d(i+2)*d(i+2) - mu;
                                 % (i+2,i+2)-entry of L*L' - mu*I
        dl(i+2) = dt(i+2) - real(ll2(i)*conj(ll2(i)))*dl(i) ...
                  - real(ll1(i+1)*conj(ll1(i+1)))*dl(i+1);
                                 % equating (i+2,i+2)-entries
    end
end

%%%%%%%%%% compute the upper triangular part U %%%%%%%%%%%
% variable declarations
du = zeros(n+2,1); % pad zeros to the ends to simplify the program
u1 = zeros(n,1);
u2 = zeros(n,1);
for i = n-2:-1:1
    du(i+2) = dt(i+2) - real(u1(i+2)*conj(u1(i+2)))*du(i+3) ...
              - real(u2(i+2)*conj(u2(i+2)))*du(i+4);
                                  % equating (i+2,i+2)-entries of
                                  % L*L' - mu*I and U*DU*U'
    u1(i+1) = (conj(lt(i+1)) - u2(i+1)*du(i+3)*conj(u1(i+2)))/du(i+2);
                                  % equating (i+1,i+2)-entries
    u2(i) = d(i)*conj(l2(i))/du(i+2);
                                  % equating (i,i+2)-entries
end
du(2) = dt(2) - real(u1(2)*conj(u1(2)))*du(3) ...
        - real(u2(2)*conj(u2(2)))*du(4);
                                  % equating (2,2)-entries
u1(1) = (conj(lt(1)) - u2(1)*du(3)*conj(u1(2)))/du(2);
                                  % equating (1,2)-entries
du(1) = dt(1) - real(u1(1)*conj(u1(1)))*du(2) ...
        - real(u2(1)*conj(u2(1)))*du(3);
                                  % equating (1,1)-entries

% discard the padded zero entries
du = du(1:n);
u1 = u1(1:n-1);
u2 = u2(1:n-2);
