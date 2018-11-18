% script file: TakagiDemo
% demonstrate Takagi factorization
%
% Dependency
%    ./csgen.m     generate complex-symmetric matrix
%    ./LanMPOR.m   Lanczos using modified partial orthogonalization
%                   with restart  
%    ./CSSVD.m     SVD of tridiagonal complex-symmetric matrix (QR)
%    ./cstsvdd.m   SVD of tridiagonal complex-symmetric matrix (D & C)
%    ./cstsvdt.m   SVD of tridiagonal complex-symmetric matrix (Twisted)

% S. Qiao       McMaster Univ.  March 2004

% get size of matrix
n = input('Enter the size of the matrix: ');

% generate singular values
sv = rand(n,1);
sv = sort(-sv); sv = -sv;

% generate complex symmetric matrix with singular values sv
A = csgen(sv); A= Wp_cov; n=length(Wp_cov);
nSteps = input('Enter the number of Lanczos iterations: ');


% Lanczos tridiagonalization using partial orthogonalization
%[a,b,Q1,nSteps] = LanPO(A,rand(n,1),nSteps);
%
% Lanczos tridiagonalization using modified partial orthogonalization
%[a,b,Q1,nSteps] = LanMPO(A,rand(n,1),nSteps);
%
% Lanczos tridiagonalization using modified partial orthogonalization
% and restart
[a,b,Q1,nSteps,nVec] = LanMPOR(A,rand(n,1),nSteps);

% get number of iterations to run
fprintf('\nNumber of iterations actually run: %d', nSteps);

% calculate and report errors
tmp = norm(Q1'*Q1 - eye(nSteps), 'fro')/(nSteps*nSteps);
fprintf('\nError in orthogonality in tridiagonalization: %E', tmp);
tmp = norm(Q1'*A*conj(Q1) - (diag(a)+diag(b,1)+diag(b,-1)), 'fro');
tmp = tmp/(nSteps*nSteps);
fprintf('\nError in tridiagonalization: %E', tmp);

% prompt the user as to which SVD algorithm to use
svd_selection = 1;
fprintf('\n\n(1) implicit QR');
fprintf('\n(2) divide-and-conquer');
fprintf('\n(3) twisted factorization');
svd_selection = input('\nSelect SVD method: ');

% svd of tridiagonal 
if (svd_selection == 1)
    [s,Q2] = CSSVD(a, b);		% pure QR 
elseif (svd_selection == 2)
    [s,ifail,Q2] = cstsvdd(a, b);	% divide-and-conquer
else
    [s,Q2] = cstsvdt(a, b);             % twisted factorization
end

Q = Q1*Q2;

% check results
tmp = norm(Q'*Q - eye(nSteps), 'fro')/(nSteps*nSteps);
fprintf('\nError in orthogonality: %E', tmp);

% calculate and report errors
if nSteps==n
    fprintf('\nError in singular values: %E', norm(s - sv)/n);

    tmp = norm(A - Q*diag(s)*conj(Q'), 'fro')/(n*n);
    fprintf('\nError in Takagi factorization: %E\n', tmp);
end
