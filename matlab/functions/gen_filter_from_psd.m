function [ filt_obj, norm_db ] = gen_filter_from_psd( psd_db, freq_vec, nfilt, fstart, fstop )
% Create a filter object from a PSD in dbHz
    
    bw = fstop-fstart;

    psd_ln=10.^(psd_db./10);           %Filter expects linear scale
    psd_norm = psd_ln./max(psd_ln);   %Normalize to 0dB to build filter
    norm_db = 10*log10(psd_norm);

    freq_norm = (freq_vec-fstart)./bw;                                %Normalize to [0,1]
    d = fdesign.arbmag('N,F,A',nfilt,freq_norm,psd_norm);            % Design filter
    filt_obj = design(d);  
    

end

