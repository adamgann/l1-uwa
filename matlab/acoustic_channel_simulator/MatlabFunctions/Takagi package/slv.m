function [lam, g, ifail] = slv(i, d, z, rho)
% [lam, g, ifail] = slv(i, d, z, rho)
%
% This function finds the root lam of the rational function
%          1 + rho*sum(|z(1:n)|^2/(d(1:n) - lam))
% in the interval
%      (d(i),d(i-1)) for i>1  or  (d(1),d(1) + rho) for i=1.
% It is used in the rank-one modification of symmetric eigenproblem:
%         diag(d) + rho*z*z'
% So, this function computes the ith eigenvalue lam of the rank-one
% modified matrix and its corresponding eigenvector.
%
% inputs:
%   i -- interval (d(i),d(i-1)) where the eigenvalue to be computed
%   d -- the original eigenvalues, the entries are sorted
%        in decreasing order
%   z -- the updating vector, normalized (2-norm)
% rho -- a positive scalar factor in the updating
%
% outputs:
%   lam -- ith eigenvalue after updating
%     g -- the corresponding eigenvector
% ifail -- 0 when successful or 1 when failed

% References
% J.R. Bunch, C.P. Nielsen and D.C. Sorensen,
% Rank-one modification of the symmetric eigenproblem,
% Numer. Math. 31, 31-48 (1978).
%
% W. Xu and S. Qiao, McMaster Univ.   Nov. 2004
% revised  April 2005

n = length(z);
z2 = conj(z).*z; % |z|^2
im1 = i - 1;

maxiter = 50; 	% maximum number of iterations
ifail = 1;
niter = 0; 	% number of iterations
lambda = 0;
del = d(i);
delta = (d - del*ones(n,1))/rho;

% caculate an initial guess
if i > 1
    del = d(im1);
else
    del =  d(1) + rho;
end
a = rho*sum(z2(i+1:n)./(d(i+1:n) - del*ones(n-i,1)));
b = rho*sum(z2(1:i-2)./(d(1:i-2) - del*ones(i-2,1)));
a = a + b + 1;
if i == 1
    t = z2(i)/a;
else
    t = a*delta(im1);
    b = t + z2(i) + z2(im1);
    if b > 0
        t = 2*z2(i)*delta(im1)/(b + sqrt(abs(b*b - 4*t*z2(i))));
    else
        t = (b - sqrt(b*b - 4*t*z2(i)))/(2*a);
    end
end

while (niter < maxiter)
    niter = niter + 1;
    if (i > 1)&(t > 0.9*delta(im1)) % t is too close to endpoint
        t = 0.9*delta(im1);         % back off a little
    end
    delta = delta - t*ones(n,1);
    lambda = lambda + t;
    lam = d(i) + rho*lambda;        % adjust the root
    
    % evaluate psi and its derivative dpsi
    psi = 0;
    dpsi = 0;
    for j = i : n        
        t = z2(j)/delta(j);
        psi = psi + t;
        dpsi = dpsi + t/delta(j);
    end
    % evaluate phi and its derivative dphi
    phi = 0;
    dphi = 0;
    if i ~= 1
        for j = 1 : im1
            t = z2(j)/delta(j);
            phi = phi + t;
            dphi = dphi + t/delta(j);
        end
    end
    
    % test for convergence
    w = 1 + phi + psi;
    eps1 = eps*n*(abs(phi) + abs(psi) + 1);
    if (abs(w) < eps1)    % converged
        g = z./delta;     % eigenvector
        g = g/norm(g);
       % norm(delta - (d - lam*ones(n,1))/rho),
        ifail = 0;  % success
        return
    end
    
    % calculate the new estimate
    if (i == 1)
        t = (w*psi)/dpsi;
    else
        del = delta(im1);
        temp = psi/dpsi;        
        a = 1 + phi - del*dphi;
        if a == 0
            return  % fail
        end
        b = (del*(1 + phi) + psi*temp)/a + temp;
        c = (2*temp*del*w)/a;
        t = c/(b + sqrt(abs(b*b - 2*c)));
    end
end % while
