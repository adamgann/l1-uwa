function [a,b,Q,nsteps] = FHLanMPOR(col,row,r,nsteps)
% [a,b,Q,nsteps] = FHLanMPOR(A,r,nsteps)
%
% Fast Lanczos tridiagonalization of a Hankel matrix
% with modified partial orthogonalization and restart.
% Given the first column and last row of a Hankel matrix,
% this computes a, b, and unitary Q such that
%       hankel(col,row)*conj(Q) = Q*(diag(a) + diag(b,1) + diag(b,-1))
%
% Inputs
%   col    first column of Hankel
%   row    last row of Hankel
%   r      starting vector
%   nsteps number of iterations
% Outputs
%   a      main diagonal of the tridiagonal
%   b      subdiagonal of the tridiagonal
%   Q      unitary
%   nsteps number of iterations actually run
%
% Dependency
%     ./fhmvmul.m   fast Hankel matrix-vector multiplication

% References
% Chengshu Guo and Sanzheng Qiao.
% A Stable Lanczos Tridiagonalization of Complex Symmetric Matrices.
% Technical Report No. CAS 03-08-SQ,
% Department of Computing and Software, McMaster University,
% Hamilton, Ontario, Canada. L8S 4K1. July 2003.
%
% Horst D. Simon.
% The Lanczos algorithm with partial reorthogonalization.
% Mathematics of Computation 42 (1984) 115--142.
%

% S. Qiao       McMaster Univ.  May 2007
%
% is A complex?
A = hankel(col, row);
if (imag(A) ~= 0)
    isComplex = 1;
else
    isComplex = 0;
end
%
n = length(r);
%
% constants
IM = sqrt(-1);
SQRTEPS = sqrt(eps);
MIDEPS = sqrt(SQRTEPS)^3;	% between sqrt(eps) and eps
RSTOL = SQRTEPS*(norm(A,'fro')/(n*n));
                                % tolerance for restart
MAXRS = 5;                      % maximal number of restarts
%
if nsteps > n nsteps = n; end
% initialize two column vectors for diagonals
a = zeros(nsteps,1); b = zeros(nsteps-1,1);
%
wOld = zeros(nsteps,1);		% orthogonality estimates
wCur = zeros(nsteps,1);
wOld(1) = 1.0;
%
up = ones(nsteps,1);		% upper and lower bounds for
low = ones(nsteps,1);           % orthogonalization intervals
interNum = 0;			% orthogonalization interval number
%
second = 0;			% if this is the second partial orthog
%

%
if (norm(r) == 0)
    r = ones(n,1) - 2*rand(n,1);
    if (isComplex == 1)
        r = r + IM*(ones(n,1) - 2*rand(n,1));
    end
end
qCur = r/norm(r);		% a unit 2-norm starting vector
%
Q(:,1) = qCur;
%
% Lanczos tridiagonalization
for j=1:nsteps
    tmp = fhmvmul(col, row, conj(qCur));
    				% fast Hankel matrix-vector multiplication
    a(j) = qCur'*tmp;
    if j == 1
        r = tmp - a(j)*qCur;
    else
        r = tmp - a(j)*qCur - b(j-1)*qOld;
    end
%
    if (j < nsteps)
        normr = norm(r);
        b(j) = normr;

        % determine orthogonalization intervals
        nRS = 0;
        nEntry = 0;         % number of entries into the following loop
        while ((nEntry < 1) | ((normr < RSTOL) & (nRS < MAXRS)))
            nEntry = nEntry + 1;
            if (normr < RSTOL)  % small subdiagonal detected
fprintf('\nSmall subdiagonal encountered. Restart ...');
                % orthogonalize r against all previous Q(:,1:j)
                interNum = 1;   % one ortho interval [1:j]
                low(interNum) = 1;
                up(interNum) = j;
                r = ones(n,1) - 2*rand(n,1);  % restart
                if (isComplex == 1)
                    r = r + IM*(ones(n,1) - 2*rand(n,1));
                end
                nRS = nRS + 1;
                % update wOld and wCur
                if (j >=2)
                    wOld(1:j-1) = wCur(1:j-1);
                    wOld(j) = 1.0;
                end
                wCur(j+1) = 1.0;
            else % subdiagonal is not small, no restart
                % compute orthogonality estimates
                if (j >= 2)
                    wOld(1) = (b(1)*conj(wCur(2)) + a(1)*conj(wCur(1)) ...
                            - a(j)*wCur(1) - b(j-1)*wOld(1))/b(j) ...
                            + eps*(b(1)+b(j))*0.3*(randn + IM*randn);
                    if (j > 2)
                        wOld(2:j-1) = (b(2:j-1).*conj(wCur(3:j)) ...
                                    + a(2:j-1).*conj(wCur(2:j-1)) ...
                                    - a(j)*wCur(2:j-1) ...
                                    + b(1:j-2).*conj(wCur(1:j-2)) ...
                                    - b(j-1)*wOld(2:j-1))/b(j) ...
  + eps*0.3*(b(2:j-1)+b(j)*ones(j-2,1)).*(randn(j-2,1) + IM*randn(j-2,1));
                    end
%
                    % update wOld and wCur
                    tmp = wOld(1:j-1);
                    wOld(1:j-1) = wCur(1:j-1);
                    wCur(1:j-1) = tmp;
                    wOld(j) = 1.0;
                end % if j>=2
                wCur(j) = eps*n*(b(1)/b(j))*0.6*(randn + IM*randn);
                wCur(j+1) = 1.0;
%
                if (second == 0) % not the second time, determine intervals
	            interNum = 0;   % reset
                    % search for orthogonalization intervals
	            k = 1;  % starting point
	            while k <= j
		        if (abs(wCur(k)) >= SQRTEPS)	% lost orthogonality
		            interNum = interNum + 1;
		            % find the upper bound
		            p = k + 1;
                            while ((p < (j + 1)) & (abs(wCur(p)) >= MIDEPS))
			        p = p + 1;	% nearly lost orthogonality
                            end % while
		            up(interNum) = p - 1;
                            % find the lower bound
		            p = k - 1;
                            while ((p > 0) & (abs(wCur(p)) >= MIDEPS))
			        p = p - 1;	% nearly lost orthogonality
		            end % while
		            low(interNum) = p + 1;
%
		            k = up(interNum) + 1;	% continue search
		        else
		            k = k + 1;
                        end % if lost orthogonality
	            end % while k <= j
	        end % if not second time
            end % if-else normr < RSTOL
            % Now that we have found interNum of intervals for
            % orthogonalization and their upper and lower bounds
            % are up(i) and low(i), carry out orthogonalization.
            if ((interNum > 0) | (second == 1))
% fprintf('\northogonalization in iteration %d', j);
	        for (k = 1:interNum)		% for each interval
		    for (i = low(k):up(k))
		        r = r - (Q(:,i)'*r)*Q(:,i); % do orthogonalization
		        wCur(i) = eps*1.5*(randn + IM*randn);
						% reset ortho estimates
		    end % for i
%
		    if (second == 1)	% this is the second time
% fprintf ('\nSecond orthogonalization. Vector #%d against Vectors #%d--#%d', j+1, low(k), up(k));                            
%
                        second = 0;		% reset
		        low(k) = 0; up(k) = 0; 
                    else % this is the first time
% fprintf ('\nOrthogonalization.        Vector #%d against Vectors #%d--#%d', j+1, low(k), up(k));                            
%
                        second = 1;	% do second time in the next iteration
		        % adjust orthogonalization intervals
		        low(k) = max(1, low(k) - 1);
                        up(k) = min(j + 1, up(k) + 1);
		    end % if-else
	        end % for k
	        normr = norm(r);
	    end % if ((interNum > 0) ... )
        end % while (nEntry | ... )
        
        if (normr < RSTOL)
            fprintf('\nRestart %d times', nRS);
            a = a(1:j); b = b(1:j-1);
            nsteps = j;
            return
        end
%
        % if no restart, we may have performed orthogonalization,
        % recalculate b(j); if restart, we use the old small b(j)
        % instead of normr.
        if (nRS == 0)       % no restart
            b(j) = normr;   % recalculate b(j)
        end
        % update qOld and qCur
        qOld = qCur;
        qCur = r/normr;
        Q(:,j+1) = qCur;
    end % if j < nsteps
end % for j
