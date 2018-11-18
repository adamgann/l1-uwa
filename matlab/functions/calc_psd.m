function [ psd_vec, f_vec ] = calc_psd( y_mat, Nfft, bw )
%Calculate a single-sided PSD

    N = size(y_mat,2);

    % Split and take FFT
    y_split = reshape(y_mat(:,1:Nfft*floor(N/Nfft)),Nfft,[]);
    Y_split = fft(y_split,Nfft);
    
    % Calculate the PSD by taking mag squared of FFT
    Sk = mean( abs(Y_split).^2, 2) / Nfft;
    
    % Convert to Single-Sided Spectrum 
    Sk_ss = Sk(1:Nfft/2+1);
    Sk_ss(2:end-1) = 2*Sk_ss(2:end-1);
    Sk_wattsHz = Sk_ss;
    Sy_dbHz = 10*log10(Sk_wattsHz)/2; %FIXME: Make sure div by 2 is right
    
    psd_vec = Sy_dbHz;
    
    deltaf = 2*bw/Nfft; %Why factor of 2?
    f_vec = (0:Nfft/2)*deltaf; 

end

