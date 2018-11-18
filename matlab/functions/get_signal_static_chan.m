function [ symVec] = get_signal_static_chan( sigVec, bitVec, bitAmpl, chanVec, alpha, sps, resetFilt, padFilt,pnSeq,nGuard)
% Generate a continuous time signal and pack it into a matrix
% This will contain a PN sequence of very high amplitude, 
% followed by a guard interval, then the packet at bitAmpl. 
% The PN sequence is thrown away after sync
% 
% sigVec:       binary signature vector
% bitVec:       information bits
% bitAmpl:      amplitude of information bits
% chanVec:      taps of a FIR filter to simulate channel
% Params:       simulation parameters struct
% resetFilt:    reset filter if using the dsp.FIRFilter object,
%               used when processing a new 'chunk' of a TX signal
% padFilt:      pad the input signal with zeros to flush the filter
% pnSeq:        Sequence of +/- for sync
% nGuard:       Number of zeros between PN seq and packet
%
% This function *should* work on MATLAB 2011b+
% DSP toolbox is preferred but should run without it. 
% 
% Adam Gannon, SUNY Buffalo, 2018


    %% Setup 
    
    persistent firObj prevData useDsp;
    
    % Use the DSP toolbox 
    if isempty(useDsp)
        useDsp = logical(exist('dsp.FIRFilter','class'));
    end
    

    % Instantiate a FIR for the channel filter (if DSP toolbox is used)
    if isempty(firObj)
        if (useDsp)
            firObj = dsp.FIRFilter;
        else
            firObj = false;
        end
    end
    
    % Initialize the previous data to zeros 
    dataLen = length(sigVec)*length(bitVec)*sps;
    if isempty(prevData)
        prevData = zeros(dataLen,1);
    end
    

    % By default, we will reset and pad the filter for each new packet 
    if (nargin == 5)
        resetFilt = true;
        padFilt = true;
    elseif (nargin == 6)
        padFilt = true;
    elseif (nargin < 5)
        error('Not enough input parameters specified');
    end
    
    
    %% Processing

    
    % Update filter taps vector.
    % If using dsp.FIR then update the Numerator parameter 
    firTaps = upsample(chanVec,sps).'; 
    
    if (useDsp)
        firObj.Numerator = firTaps;
    end
    
    
    % Reset the FIR filter, if new packet
    if (resetFilt)
        if (useDsp)
            reset(firObj)
        else
            prevData = zeros(dataLen,1);
        end
        
    end

        
    % Spread and upsample bits 
    bitsSpreaded = kron((bitAmpl.*bitVec)',sigVec);
    bitsUpsamp = upsample(bitsSpreaded,sps);
    
    
    % Pulse shape filter bits, using rcosdesign if available.
    % Maintain backwards compatability with older versions (rcosine)
    if (exist('rcosdesign','file'))
        g_T = rcosdesign(alpha,6,sps,'sqrt');
    else
        g_T = rcosine(1,sps,'sqrt',alpha);
    end
    g_T = g_T./norm(g_T);

    
    % Pad the input data with zeros to flush the filter (number of zeros
    % should equal length of filter taps). 
    % Then pulse shape filter the bits. 
    if (padFilt)
        bitsUpsamp=[bitsUpsamp;zeros(length(firTaps),1)];
    end    
    bitsShaped = conv(g_T,bitsUpsamp);

    % Filte the PN sequence which is 100x higher ampl than packet
    pnUpsamp = upsample((100*bitAmpl*pnSeq),sps);
    pnShaped = conv(g_T,pnUpsamp);
    
    % PN Seq | Gaurd Interval | Packet
    bitsShaped = [pnShaped; zeros(nGuard,1); bitsShaped];

    
    % Run the TX signal through a channel 
    if (useDsp)
        symVec = step(firObj,bitsShaped);
    else
        
        % Get the initial state of the filter by filtering the OLD data
        % with the NEW filter taps. 
        [~,filtState] = filter(firTaps,1,prevData);
        prevData = bitsShaped; 
        
        % Filter NEW data with the NEW channel taps
        [symVec] = filter(firTaps,1,bitsShaped,filtState);
        
    end
    
    


end

