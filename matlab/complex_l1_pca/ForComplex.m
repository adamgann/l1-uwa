function [Q, B, metric, cnts, init, kept, freq]=ForComplex(X,K,M,both,init,tol,use_mex)

%
% DOC of this function: Local maximizer of || X'*Q ||_1 over Q unitary. 
% The algorithms used are two 
%           I)  Stoica Iteration     (Guaranteed convergence)
%           II) Tsagarakis Iteration (If convergent, fast)
% Many initializations.
%c
%
%
% OUTPUTS:
% -------------------------------------------------
% Q: Basis/Bases for local optimal subspaces.
% B: Associated optimal unimodular matrices.
% metric: Objective value (L1 norm) acheived for the above bases.
% init: The initialization used (either from input or generated herein).
% kept: The indeces of the converged bases that are non equivalent.
% frep: The number of occurness of each non-equivalent bases.
%
%
% INPUTS:
% -------------------------------------------------
% X: The data matrix.
% K: The dimension of the subspace we are looking for.
% both: If false, try both iterations. If false, only Stoica is emplyed.
% M: Number of initializations.
% init: Option to specify the initial points.
% tol: Converging constant.
% countFreqs: If true, keep track of the frequencies; if false, don't.
%

%%
[D,N]=size(X);

%%
if ~exist('K','var') || isempty(K),K=1;end % if isempty(K),K=1;end
if K>1,both=false;end % there is no alternative alg for K>1
if ~exist('both','var') || isempty(both),both=false;end % if isempty(both),both=false;end

if exist('init','var') && ~isempty(init) % 'init' has priority over 'M' % if ~isempty(init) % 'init' has priority over 'M' 
    assert(size(init,1)==N & size(init,2)==K)
    M=size(init,3);
else
    if ~exist('M','var') || isempty(M),M=1e2;end %     if isempty(M),M=1e2;end
    
    %init = scmat(N,K,M);
    init = sign(randn(N,K,M)+1i*randn(N,K,M));
    if 1 || M>1 % sign of L2-PC
        [~,S,V]=svd(X,'econ');
        init(:,:,1)=sign(V(:,1:K));
        if (0 && (K==1 && (1 || M>2)))
            diag_S=diag(S);
            init(:,:,2)=sign(V(:,1:length(diag_S))*diag_S);
        end
    end
end
if ~exist('tol','var') || isempty(tol),tol=1e-7;end % if isempty(tol),tol=1e-7;end
if ~exist('use_mex','var') || isempty(use_mex),use_mex=false;end % if isempty(use_mex),use_mex=false;end

%%
if both % K==1 is implied
    
    if use_mex 
        [Q,B,metric,cnts]=ForComplex_K1_Both_mex(X,D,N,init,M,tol);
    else
        Q=czeros(D,1,M,2);
        B=cat(4,init,init);
        metric=zeros(M,2);
        cnts=zeros(M,2);
        % First and second dimentions of B(and Q) is of size [N,K] 
        % (and [D,K]). If K==1 (this case) then, B(and Q) could be 
        % squeezed. I choose not to squeeze it for consistency.

%         A=X'*x;Ad=A-diag(diag(A));
        Ad=X'*X;Ad(1:N+1:end)=0;
        for m=1:M
            % Stoica
            mtr_prev=0;
            while true
                cnts(m,1)=cnts(m,1)+1;
                
                [Q(:,1,m,1),metric(m,1)]=unt(X*B(:,1,m,1));
                B(:,1,m,1)=sign(X'*Q(:,1,m,1));
                
                if metric(m,1)-mtr_prev<tol,break,end
                mtr_prev=metric(m,1);
            end
            Q(:,1,m,1)=Q(:,1,m,1)/B(1,1,m,1);
            B(:,1,m,1)=B(:,1,m,1)*conj(B(1,1,m,1));

            % tsaga
            mtr_prev=0;
            while true
                cnts(m,2)=cnts(m,2)+1;
                
                tmp=Ad*B(:,1,m,2);
                metric(m,2)=norm(tmp,1);
                B(:,1,m,2)=sign(tmp);
                
                if abs(metric(m,2)-mtr_prev)<(tol)^.675,break,end % SHOULD THERE BE abs()????????
%                 if abs(metric(m,2)-mtr_prev)<tol,break,end
                mtr_prev=metric(m,2);
            end
            Q(:,1,m,2)=orth(X*B(:,1,m,2));
            B(:,1,m,2)=B(:,1,m,2)*conj(B(1,1,m,2));
        end
    end

%     return
    
    %% ********************************************************************
    % Get rid of equivalent B's (O(M^2))
    [B, kept, freq]=removeDuplicates(B);
    Q={Q(:,1,kept{1},1), Q(:,1,kept{2},2)};
    metric={metric(kept{1},1), metric(kept{2},2)};
    
%     return
    
    %% ********************************************************************
    % Pick a representative for each value of metric
    % ('kept' becomes meaningless unless it turns length(inds)=length(I))
    for i=1:2
        [metric{i},I]=sort(metric{i});
        inds=[0;find(diff(metric{i})>1e-6)]+1;
        metric{i}=metric{i}(inds);
        Q{i}=Q{i}(:,:,I(inds));
        B{i}=B{i}(:,:,I(inds));
    end
    
%     return
    
    %% Optional
    % return only the best (w.r.t. the metric) basis and B
    for i=1:2
        [metric{i},mx_pos]=max(metric{i});
        Q{i}=Q{i}(:,mx_pos);
        B{i}=B{i}(:,mx_pos);
    end
    
else % both = false
    
    if use_mex
        if K==1
            [Q,B,metric,cnts]=ForComplex_K1_only_Kwak_mex(X,D,N,init,M,tol);
        else
            % Mex is slower than matlab if K is big. The reason 
            % for this is that, gsl does not have svd for complex.
            [Q,B,metric,cnts]=ForComplex_mex(X,D,N,K,init,M,tol);
            for m=1:M,for k=1:K,B(:,k,m)=B(:,k,m)/B(1,k,m);Q(:,k,m)=Q(:,k,m)*conj(B(1,k,m));end,end
        end
    else % ,X'*X
       % Q=czeros(D,K,M); % K might be 1 (no problem)
        Q=zeros(D,K,M) + 1i*zeros(D,K,M); % K might be 1 (no problem)
        B=init;
        metric=zeros(M,1);
        cnts=zeros(M,1);
        
        for m=1:M
            % Stoica
            mtr_prev=0;
            while true
                cnts(m)=cnts(m)+1;
                
                [Q(:,:,m),metric(m)]=unt(X*B(:,:,m));
                B(:,:,m)=sign(X'*Q(:,:,m));
                
                if metric(m)-mtr_prev<tol,break,end
                mtr_prev=metric(m);
            end
            for k=1:K,Q(:,k,m)=Q(:,k,m)/B(1,k,m);B(:,k,m)=B(:,k,m)*conj(B(1,k,m));end
        end  
    end % if use_mex
    
%     return

%     diff(sort(metric))
%     figure,plot(sort(metric))
    
    %% ********************************************************************
    % Get rid of equivalents B's (O(M^2))
    [B, kept, freq]=removeDuplicates(B);
    Q=Q(:,:,kept);
    metric=metric(kept);
    
%     return
    
    %% ********************************************************************
    % Pick a representative for each value of metric 
    % ('kept' becomes meaningless unless it turns length(inds)=length(I))
    [metric,I]=sort(metric);
    inds=[0;find(diff(metric)>1e-6)]+1;
    metric=metric(inds);
    Q=Q(:,:,I(inds));
    B=B(:,:,I(inds));
    
%     return
    
    %% Optional
    % return only the best (w.r.t. the metric) basis and B
    [metric,mx_pos]=max(metric);
    Q=Q(:,:,mx_pos);
    B=B(:,:,mx_pos);
    
end
