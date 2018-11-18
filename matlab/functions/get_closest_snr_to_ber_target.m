function [snrVal] = get_closest_snr_to_ber_target(snrVec,berVec,berTarget)
% Given a BER vs SNR curve, return the SNR that produces a certain bit
% error rate

    % Perform linear interpolation to get a smaller SNR step
    deltaSnr = 0.5;
    interpSnr = snrVec(1):deltaSnr:snrVec(end);
    interpBer = interp1(snrVec,berVec,interpSnr);

    % Calculate the closest SNR point to target 
    diffVec = abs(interpBer-berTarget);
    [~,ind] = min(diffVec);
    snrVal = interpSnr(ind);

end

