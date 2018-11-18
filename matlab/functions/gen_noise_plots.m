function gen_noise_plots( ambientNoise, shrimpNoise )
% Create some time domain plots to double check the noise generation

      figure;
      subplot(2,1,1);
      plot(real(ambientNoise),'r');
      xlabel('Sample Index')
      ylabel('Amplitude (uPa)')
      title('Ambient Noise')

      subplot(2,1,2);
      plot(20*log10(abs(ambientNoise)),'r')
      xlabel('Sample Index')
      ylabel('dB re uPa')

      figure;
      subplot(2,1,1);
      plot(real(shrimpNoise),'b');
      xlabel('Sample Index')
      ylabel('Amplitude (uPa)')
      title('Shrimp Noise')

      subplot(2,1,2);
      plot(20*log10(abs(shrimpNoise)),'b')
      xlabel('Sample Index')
      ylabel('dB re uPa')

      figure;
      plot(real(shrimpNoise),'b');
      hold on
      plot(real(ambientNoise),'r')
      
end

