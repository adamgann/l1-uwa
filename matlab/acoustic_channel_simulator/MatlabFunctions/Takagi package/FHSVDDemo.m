% script file: FHSVDDemo.m
% Demonstrate fast Hankel SVD
%
% Dependency
%    ./FHLanMPOR.m fast Lanczos for Hankel using modified partial
%                  orthogonalization and restart
%    ./CSSVD.m     SVD of tridiagonal complex-symmetric matrix (QR)
%    ./cstsvdd.m   SVD of tridiagonal complex-symmetric matrix (D & C)
%    ./cstsvdt.m   SVD of tridiagonal complex-symmetric matrix (Twisted)

% S. Qiao       McMaster Univ.  March 2004
% Revised May 2007

n = input('Enter the size: ');
% generate the first col and last row of Hankel
col = (2*rand(n,1) - ones(n,1)) + j*(2*rand(n,1) - ones(n,1));
row = (2*rand(n,1) - ones(n,1)) + j*(2*rand(n,1) - ones(n,1));
row(1) = col(n);	% last element of col wins over
			% first element of row

nSteps = input('Enter the number of Lanczos iterations: ');
%
% fast Lanczos tridiagonalization of Hankel using modified
% partial orthogonalization
%[a,b,Q1,nSteps] = FHLanMPO(col, row, ones(n,1), nSteps);
%
% fast Lanczos tridiagonalization of Hankel using modified
% partial orthogonalization with restart
[a,b,Q1,nSteps] = FHLanMPOR(col, row, ones(n,1), nSteps);

fprintf('\nNumber of iterations actually run: %d', nSteps);

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

if nSteps==n
    H = hankel(col, row);
    fprintf('\nError in singular values: %E', norm((s - svd(H))/n));

    tmp = norm(H - Q*diag(s)*conj(Q'), 'fro')/(n*n);
    fprintf('\nError in fast Hankel SVD: %E\n', tmp);
end

