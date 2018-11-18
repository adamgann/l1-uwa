function [d,Q] = eig2tak(a,b,din,Qin)
% [d,Q] = eig2tak(a,b,din,Qin)
%
% Convert the eigenvalues din and eigenvectors Qin of T*T',
%       T = diag(b,-1) + diag(a) + diag(b,1)
% is symmetric tridiagonal, into the Takagi (singular) values
% d and Takagi (left singular) vectors of T,
%            T = Q*diag(d)*Q.'
%
% inputs
%    a -- main diagonal of the symmetric tridiagonal T
%    b -- subdiagonal of the symmetric tridiagonal T
%  din -- eigenvalues of T*T'
%  Qin -- eigenvectors of T*T'
% outputs
%    d -- Takagi (singular) values
%    Q -- Takagi (left singular) vectors
%       diag(b,-1) + diag(a) + diag(b,1) = Q*diag(d)*Q.'
%
% W. Xu and S. Qiao    McMaster Univ.   April 2005

n = length(din);
NEPS = n*eps;
SQRTEPS = sqrt(eps);

% Takagi values
d = sqrt(din);
[d,indx] = sort(-d);
d = -d;

% Takagi vectors
Q = Qin(:,indx);
i = 1;
while i<=n
    % find the multiplicity of d(i)
    m = 1;	% multiplicity
    while (i+m <= n) & (abs(d(i) - d(i+m)) < 10*NEPS)
        m = m+1;
    end

    if m==1	% d(i) is a single Takagi value
        % compute sign(Q(:,i)'*T*conj(Q(:,i)))
        tq = stmv(a,b,conj(Q(:,i)));
        tmp = Q(:,i)'*tq;
        if tmp==0
            sgn = 1;
        else
            sgn = tmp/abs(tmp);
        end
        % Takagi vector
        Q(:,i) = sqrt(sgn)*Q(:,i);
    else % d(i) is a multiple Takagi value, multiplicity m
        if d(i)>10*SQRTEPS	% d(i) not numerically zero
            for j=i:i+m-1 % all Q(:,j) corresponding to the multiple value
                tq = stmv(a,b,conj(Q(:,j))); % T*conj(Q(:,j))
                q = tq + d(j)*Q(:,j);
                % Modified Gram-Schmidt orthogonalization
                for k=i:j-1  % agaist j-1 previous vectors
                    q = q - (Q(:,k)'*q)*Q(:,k);
                end
                Q(:,j) = q/norm(q);
             end % for
        end % if
    end % if-else
    i = i + m;
end % while


function y = stmv(a,b,x)
% y = stmv(a,b,x)
%
% Symmetric tridiagonal matrix-vector multiplication.
%     T = diag(b,-1) + diag(a) + diag(b,1)
% is symmetric tridiagonal,
%              y = T*x
%
% inputs
%   a -- main diagonal of symmetric tridiagonal matrix
%   b -- subdiagonal of symmetric tridiagonal matrix
%   x -- column vector
% output
%   y -- (diag(b,-1) + diag(a) + diag(b,1))*x

n = length(x);
y = [0; b.*x(1:n-1)] + a.*x + [b.*x(2:n); 0];

