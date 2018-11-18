function [a,b,Q,steps] = FHLanMPO(col,row,r,steps)
% [a,b,Q] = FHLanMPO(col,row,r,steps)
%
% Fast Lanczos tridiagonalization of a Hankel matrix
% with modified partial orthogonalization.
% Given the first column and last row of a Hankel matrix,
% this computes a, b, and unitary Q such that
%       hankel(col,row)*conj(Q) = Q*(diag(a) + diag(b,1) + diag(b,-1))
%
% Inputs
%   col   first column of Hankel
%   row   last row of Hankel
%   r     starting vector
%   steps number of iterations
% Outputs
%   a     main diagonal of the tridiagonal
%   b     subdiagonal of the tridiagonal
%   Q     unitary
%   steps number of iterations actually run
%
% Dependency
%     ./fhmvmul.m   fast Hankel matrix-vector multiplication

% Reference
% F.T. Luk and S. Qiao.
% A fast singular value algorithm for Hankel matrices.
% Fast Algorithms for Structured Matrices:
% Theory and Applications, Contemporary Mathematics 323,
% Editor V. Olshevsky,
% American Mathematical Society. 2003. 169--177.

% S. Qiao       McMaster Univ.  March 2004
%
n = length(r);
%
% constants
IM = sqrt(-1);
SQRTEPS = sqrt(eps);
MIDEPS = sqrt(SQRTEPS)^3;	% between sqrt(eps) and eps
%
if steps > n steps = n; end
% initialize two column vectors for diagonals
a = zeros(steps,1); b = zeros(steps-1,1);
%
wOld = zeros(steps,1);		% orthogonality estimates
wCur = zeros(steps,1);
wOld(1) = 1.0;
%
up = ones(steps,1);		% upper and lower bounds for
low = ones(steps,1);            % orthogonalization intervals
interNum = 0;			% orthogonalization interval number
%
doOrtho = 0;			% if do orthogonalization
second = 0;			% if this is the second partial orthog
%
nVec = 0;			% number of vectors selected for
				% orthogonalization
%
qCur = r/norm(r);		% a unit 2-norm starting vector
%
Q(:,1) = qCur;
%
for j=1:steps
    tmp = fhmvmul(col, row, conj(qCur));
				% fast Hankel matrix-vector multiplication
    a(j) = qCur'*tmp;
    if j == 1
        r = tmp - a(j)*qCur;
    else
        r = tmp - a(j)*qCur - b(j-1)*qOld;
    end
%
    if (j < steps)
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
        if (second == 0)	% not the second time, determine intervals
	    doOrtho = 0;	% initialization
	    interNum = 0;
	    k = 1;
	    while k <= j
		if (abs(wCur(k)) >= SQRTEPS)	% lost orthogonality
		    doOrtho = 1;
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
	    end % while k
	end % if not second time
%
        if (doOrtho | (second == 1))	% now we have intervals,
					% carry out orthogonalization
% fprintf('\northogonalization in iteration %d', j);
	    for (k = 1:interNum)		% for each interval
		for (i = low(k):up(k))
		    r = r - (Q(:,i)'*r)*Q(:,i);	% do orthogonalization
		    wCur(i) = eps*1.5*(randn + IM*randn);
						% reset ortho estimates
		end % for i
%
		nVec = nVec + up(k) - low(k) + 1;
				% count the number of vectors selected
		if (second == 1)	% this is the second time
		    second = 0;		% reset
		    low(k) = 0; up(k) = 0; 
                else
		    second = 1;		% do second time
		    doOrtho = 0;	% reset
		    % adjust orthogonalization intervals for the second time
		    low(k) = max(1, low(k) - 1);
                    up(k) = min(j + 1, up(k) + 1);
		end % if
	    end % for k
	    b(j) = norm(r);     % recalculate b(j)
	end % if
%
        if (abs(b(j)) < eps)
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
%
% fprintf('\nNumber of vectors selected for orthogonalization: %d', nVec);
