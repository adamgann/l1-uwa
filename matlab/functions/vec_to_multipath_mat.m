function [ dataMat ] = vec_to_multipath_mat( dataVec, L, M, N)
% Create a multipath-processed matrix from a vector using 'buffer'

    % Set params
    L_M = L+M-1;

    % Pad the vector with zeros for the buffer operation. 
    targetLength = N*L + L_M; 
    bufferLength = targetLength - length(dataVec);
    dataVecBuf = [dataVec;zeros(bufferLength,1)];
    
    % Buffer and trim the resulting matrix
    dataMatRaw = buffer(dataVecBuf,L_M,M-1,'nodelay');
    dataMat = dataMatRaw(:,1:N); 

end

