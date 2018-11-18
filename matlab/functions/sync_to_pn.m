function [sampInd] = sync_to_pn(rxFilt,pnUp,bitAmpl,gT, nGuard)
% Correlate with the PN sequence and determine the start index for sampling
%
% rxFilt:   Pulse-shape filtered received signal vector
% pnUp:     Upsampled PN Sequence
% bitAmpl:  Amplitude of data bits
% gT:       Filter taps of RRC pulse-shape filter
% nGuard:   Number of zeros between PN seq and data
%
% sampInd:  Index of first data bit in packet
%
% Adam Gannon, SUNY Buffalo, 2018.
%%    
    % Correlate 
    pnConv = abs(xcorr(pnUp,rxFilt));

    
    % Peak detection - trying to catch second arrival 
    mpd = 5;
    [peaks,locs] = findpeaks(pnConv,'MinPeakDistance',mpd,'MinPeakHeight',bitAmpl);

    % Get the two highest value peaks
    [peakSort,inds] = sort(peaks,'descend');
    maxPeaks = peakSort(1:2);
    maxInds = locs(inds(1:2)); 



    % If peak detection returned a second peak that is signficantly
    % weaker than the max one, ignore it. 
    maxAmplitudeDiff = 5;
    if ((maxPeaks(1)/maxPeaks(2)) > maxAmplitudeDiff)
        ind = maxInds(1);
    else
        % Then find out if the most prominent peak was the first or second
        % channel path. 
        if (maxInds(1) > maxInds(2))
            ind = maxInds(1);
        else
            ind = maxInds(2);
        end
    end

    if (0)
        figure;plot(pnConv);
        hold on
        plot(maxInds,pnConv(maxInds),'ro')
        xlim([min(maxInds)-100 max(maxInds)+100])
        return
    end

    % Get the start index 
    start1 = abs(length(rxFilt)-ind+1);
    start2 = start1 + floor(length(gT)/2);
    start3 = start2 + length(pnUp);% - floor(length(g_T)/2);
    start4 = start3 + nGuard;% - length(g_T) +1;
    sampInd = start4;
%         sampInd = sampInd - 88;
%         sampInd = 5148

    if (0)
        maxVec = zeros(size(pnConv));
        maxVec(maxInds) = maxPeaks;
        figure;plot(real(pnConv));
        hold on
        plot(maxVec,'r')
        xlim([maxInds(1)-200 maxInds(1)+200])

        plotdat(rxFilt)
        hold on
        plot(start3,rxFilt(start3),'r*')
        plot(start4,rxFilt(start4),'g*')
        return
    end
end

