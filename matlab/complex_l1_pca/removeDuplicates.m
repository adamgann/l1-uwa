function [B, kept, freq]=removeDuplicates(B)

%
% This is for B unimodular (binary is special case).
%

assert(norm(B(:)-sign(B(:)))<1e-8)
B = sign(B);

[N,K,M,Na]=size(B); % Na represents the number of algorithms tested

Bnew=cell(1,Na);
kept=cell(1,Na);
freq=cell(1,Na);

toKeep=cell(Na,1);
for i=1:Na
    Bnew{i}=zeros(N,K,0);
    toKeep{i}=true(M,1);
    freq{i}=[];
    for j=1:M
        if toKeep{i}(j)
            % jth does not exist in Bnew yet
            Bnew{i}=cat(3,Bnew{i},B(:,:,j,i));
            freq{i}=[freq{i}, 1];
            for k=j+1:M
                if toKeep{i}(k) && areEquiv(B(:,:,j,i),B(:,:,k,i),N)
                    toKeep{i}(k)=false;
                end
            end
        else
            % jth already exists somewhere
            for k=1:size(Bnew{i},3)
                if areEquiv(B(:,:,j,i),Bnew{i}(:,:,k),N)
                    freq{i}(k)=freq{i}(k)+1;
                    break
                end
            end
        end
    end
    kept{i}=find(toKeep{i});
end
if Na==1
    B=Bnew{1};
    kept=kept{1};
    freq=freq{1};
else
    B=Bnew;
end

function yn=areEquiv(B1,B2,N)
yn = all(max(abs(B1'*B2))>N-1e-5);
