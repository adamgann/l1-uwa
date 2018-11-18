% Thorp: absorption loss in dB/kyd; f is in kHz.
% Corrected for km.

function alpha=absorption(f)


alpha=0.11*f.^2./(1+f.^2)+44*f.^2./(4100+f.^2)+2.75*10^(-4)*f.^2+0.003;
indvlf=find(f<0.3);
alphas=2*10^(-3);
alpha(indvlf)=alphas+0.11*f(indvlf).^2./(1+f(indvlf).^2)+0.011*f(indvlf).^2;
