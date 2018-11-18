function [berVec] = run_one_trial(hVec,pnSeq,gT,Params)
% Run one trial of the simulation. Output the BER of different algorithms
%
% hVec:     Mx1 channel vector
% pnSeq:    PN sequence for sync
% gT:       Taps of RRC pulse-shape filter
% Params:   Struct containing simulation parameters
%
% berVec:   BER results in key-value pairs
%
%
% Adam Gannon - SUNY Buffalo - 2018.

        %% Generate TX Signal
        
        % Random sequence of bits
        bitsUser = sign(randn(Params.nPacket,1));

        % Random code sequence 
        sigUser = sign(randn(Params.L,1))./sqrt(Params.L);
        sigMat = toeplitz([sigUser;zeros(Params.Mrec-1,1)],[sigUser(1),zeros(1,Params.Mrec-1)]);

        % Generate the continuous time signal 
        bitAmpl = 10^(Params.currentSnr/20);
        
        % Get the CT transmitted signal, including the PN sequence.
        nGuard = 4000;
        sigVec = get_signal_static_chan(sigUser, bitsUser.', bitAmpl, hVec, Params.psfAlpha, Params.sps, true, true, pnSeq, nGuard); 


        %% Generate Noise and Interference 
        
        
        % Ambient noise: either from noise PSD or AWGN 
        if (Params.useAmbientNoise)
            ambientNoise = get_colored_noise(Params.L,Params.M,length(sigVec),Params);
        else
            ambientNoise = noiseAmpl.*(randn(length(sigVec),1) + 1j* randn(length(sigVec),1))./sqrt(2);
        end
        
        
        % Add shrimp noise
        shrimpNoise = zeros(size(ambientNoise));
        if (Params.useShrimpNoise)
            
            shrimpNoise = get_shrimp_noise(length(sigVec),Params.shrimpAlpha);

            if (Params.enableDebug)
                gen_noise_plots(ambientNoise,shrimpNoise);
                return
            end

        end
        noiseVec = ambientNoise + shrimpNoise; 
        
       
        

        
        % Generate interference 
        autocorrInterf = zeros(Params.L_M);
        interfVec = zeros(size(sigVec));
        for kk=1:Params.K
            
            sigVecInterf = sign(randn(Params.L,1))./sqrt(Params.L);
            sigMatInterf = toeplitz([sigVecInterf;zeros(Params.M-1,1)],[sigVecInterf(1),zeros(1,Params.M-1)]);
            bitsInterf = sign(randn(N_packet,1));
            amplInterf = 10^(interfSnr/10);
            hvecInterf = (randn(M,1) + 1j*randn(M,1))./sqrt(2);
            
            interfVec = interfVec + get_signal_static_chan(sigVecInterf,bitsInterf.', amplInterf, hvecInterf, alpha, sps, true, true);
            
            % Update the ideal autocorrelation matrix
            autocorrInterf = autocorrInterf + ((amplInterf.^2)*(sigMatInterf*(hvecInterf*hvecInterf')*sigMatInterf'));

        end
        
        
        % Combine the signal with noise and interference
        distVec = noiseVec + interfVec;
        
        % Only add noise over the signal, not the PN sequence. 
        % the PN sequence is just for sync not part of the results
        pnUp = conv(upsample(pnSeq, Params.sps), gT);
        distVecStart = length(pnUp)+(nGuard);
        distVec(1:distVecStart) = 0;
        
        
        rxVec = sigVec + distVec;

        
        %% Symbol Sync

        % Filter and sample 
        rxFilt = conv(rxVec,gT);

        % Get the index corresponding to the start of the bits
        sampInd = sync_to_pn(rxFilt,pnUp,bitAmpl,gT, nGuard);
        


%%      Create RX Matrix
        
        rxSamp = rxFilt(sampInd:Params.sps:end);


        % Create the RX matrix
        rxMat = vec_to_multipath_mat(rxSamp,Params.L,Params.Mrec,Params.nPacket).';
 
        
        % Channel estimation
        if (Params.estChan)

            % Estimate channel with pilots
            hEst = (1/Params.nPilot)*(pinv(sigMat))*rxMat((1:Params.nPilot),:).'*bitsUser(1:Params.nPilot);
            hEst = hEst./bitAmpl; %So the scatterplots and EVM work
            
            %h_estimate_m_final_store(:,iPkt) = hEst;

        else

            % Use the known channel
            hEst = hVec;

        end


        %% Filter Signals

        % L2-based channel estimation 
        [u,~,~] = svd(rxMat.','econ');
        w_l2pca = u(:,1);
        
        % 1-bit for phase ambiguity resolution
        phi_l2 = angle(w_l2pca'*rxMat(1:1,:).'*bitsUser(1:1));

        % L1-based channel estimation 
        nInit=5; %5 reinitializations, probably a bit overkill
        w_l1pca = ForComplex(rxMat.',1,nInit);
        
        % 1-bit for phase ambiguity resolution
        phi_l1 = angle(w_l1pca'*rxMat(1:1,:).'*bitsUser(1:1));


        %% Equalization
        
        % MF filter
        w_r_mf_norm = (sigMat*hEst)/norm(sigMat*hEst).^2;
        equalization = w_r_mf_norm'*rxMat.';   

 
        % PCA filters
        equalization_l2 = exp(-1i*phi_l2)*w_l2pca'*rxMat.';
        equalization_l1 = exp(1i*phi_l1)*w_l1pca'*rxMat.';

        
        % Calculate SNR and BERs
        
        % Cut the last few bits to account for filter responses and other
        % effects in simulation that produce errors at the end of the
        % packet. Also cut off any known bits (pilots). For the L1/L2 case,
        % that's the first bit only. For the MF that's the first N_pilot.
        nCut = 5;
        bitsCut = bitsUser(2:end-nCut).';
        bitsPayload = bitsUser(Params.nPilot:end-nCut).';
        
        
        % Cut the end of the RX signal vectors
        rxMf = equalization(Params.nPilot:end-nCut)./bitAmpl;
        rxL2 = equalization_l2(2:end-nCut)./bitAmpl;
        rxL1 = equalization_l1(2:end-nCut)./bitAmpl;
         

        % BER Calculation
        bRecMf = sign(real(rxMf));
        errVecMf = (bRecMf ~= bitsPayload);
        berMf = nnz(errVecMf)/length(errVecMf);
       
        bRecL2 = sign(real(rxL2));
        errVecL2 = (bRecL2  ~= bitsCut);
        berL2 = nnz(errVecL2)/length(errVecL2);
        
        bRecL1 = sign(real(rxL1));
        errVecL1 = (bRecL1  ~= bitsCut);
        berL1 = nnz(errVecL1)/length(errVecL1);
        


        %% Store the BER in a vector 
        keySet = {'MF_PILOT','L1','L2'};
        berSet = [berMf, berL1, berL2];
        berVec = containers.Map(keySet,berSet);

end

