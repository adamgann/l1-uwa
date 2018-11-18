function [a,b,Q,steps] = LanPO(A,r,steps)
% [a,b,Q] = LanPO(A,r,steps)
%
% Lanczos tridiagonalization with partial orthogonalization.
% Given a (complex) symmetric matrix A, this computes
% a, b, and unitary Q such that
%       A*conj(Q) = Q*(diag(a) + diag(b,1) + diag(b,-1))
%
% Input
%   A     complex-symmetric matrix
%   r     starting vector
%   steps number of iterations
% Outputs
%   a     main diagonal of the tridiagonal
%   b     subdiagonal of the tridiagonal
%   Q     unitary
%   steps number of iterations actaully run


% References
% Chengshu Guo and Sanzheng Qiao.
% A Stable Lanczos Tridiagonalization of Complex Symmetric Matrices.
% Technical Report No. CAS 03-08-SQ,
% Department of Computing and Software, McMaster University,
% Hamilton, Ontario, Canada. L8S 4K1. July 2003.
%
% Horst D.~Simon.
% The Lanczos algorithm with partial reorthogonalization.
% Mathematics of Computation 42 (1984) 115--142.
%
% Demmel, James W.
% Applied Numerical Linear Algebra.
% SIAM, Philadelphia, PA. 1997.
% 382--383.
%

% S. Qiao       McMaster Univ.  Sept 03
%
n = length(r);
%
% constants
IM = sqrt(-1);
SQRTEPS = sqrt(eps);
%
n = min(steps,n);
% initialize two column vectors for diagonals
a = zeros(n,1); b = zeros(n-1,1);
%
wOld = zeros(n,1);		% orthogonality estimates
wCur = zeros(n,1);
wOld(1) = 1.0;
%
second = 0;			% if this is the second partial orthog
%
nVec = 0;			% number of vectors selected for
				% orthogonalization
%
qCur = r/norm(r);		% a unit 2-norm starting vector
%
Q(:,1) = qCur;
%
for j=1:n
    tmp = A*conj(qCur);
    a(j) = qCur'*tmp;
    if j == 1
        r = tmp - a(j)*qCur;
    else
        r = tmp - a(j)*qCur - b(j-1)*qOld;
    end
%
    if (j < n)
        b(j) = norm(r);
%
        if (j > 2) 		% compute orthogonality estimates
            wOld(1) = (b(1)*conj(wCur(2)) + a(1)*conj(wCur(1)) ...
                       - a(j)*wCur(1) - b(j-1)*wOld(1))/b(j) ...
                      + eps*(b(1)+b(j))*0.3*(randn + IM*randn);
            wOld(2:j-1) = (b(2:j-1).*conj(wCur(3:j)) ...
                           + a(2:j-1).*conj(wCur(2:j-1)) ...
                           - a(j)*wCur(2:j-1) ...
                           + b(1:j-2).*conj(wCur(1:j-2)) ...
                           - b(j-1)*wOld(2:j-1))/b(j) ...
  + eps*0.3*(b(2:j-1)+b(j)*ones(j-2,1)).*(randn(j-2,1) + IM*randn(j-2,1));
%
            % swap wOld and wCur
            tmp = wOld(1:j-1);
            wOld(1:j-1) = wCur(1:j-1);
            wCur(1:j-1) = tmp;
            wOld(j) = 1.0;
        end % if j>2
        wCur(j) = eps*n*(b(1)/b(j))*0.6*(randn + IM*randn);
        wCur(j+1) = 1.0;
%
        if ((max(abs(wCur(1:j))) > SQRTEPS) | (second==1))
            			% reorthogonalize
% fprintf('\northogonalization in iteration %d', j)
            if (max(abs(wCur(1:j))) > SQRTEPS)
                                % if orthogonalization is required
                second = 1;     % always do it again in next iteration
            else                % if this is the second time
                second = 0;     % reset
            end
%
	    for k = 1:j		% orthogonalization
                r = r - (Q(:,k)'*r)*Q(:,k);
                wCur(k) = eps*1.5*(randn + IM*randn);
				% set wCur(1:j) after orthogonalization
            end
%
            nVec = nVec + j;	% count the number of vectors selected
%
            b(j) = norm(r);	% recalculate b(j)
        end % if ((max ...

        if (abs(j) < eps)
            a = a(1:j); b = b(1:j-1);
            steps = j;
            return
        else
            qOld = qCur;
            qCur = r/b(j);
            Q(:,j+1) = qCur;
	end
    end % if j<n
end % for j
steps = n;
%
% fprintf('\nNumber of vectors selected for orthogonalization: %d', nVec)
