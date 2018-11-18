function [ noiseVec ] = get_shrimp_noise( nSamps, alpha )
% Use a univariate SaS distribution to generate snapping shrimp noise.
%
% Input: number of samples (nSamps), and char param (alpha).
% Output: noise vector in uPa
%
% See shrimp_noise.m for more details. 
%
% Adam Gannon, SUNY Buffalo, 2018.

%% Parameters 

% Debugging, set false for running sims
verifyApd = false;

%% Generate noise

 % Using scale param from Chirtre2006 
gamma = 1.5e5;

% Use stblrnd to generate a symmetric (beta=delta=0)
noiseVec = stblrnd(alpha,0,gamma,0,nSamps,1);



if (verifyApd)
    
    [aX,pX] = apd(abs(noiseVec));
    
    figure;
    semilogy(aX,pX)
    ylim([1e-2 1])
    xlim([0 6e5])
    grid on

    xlabel('Amplitude (uPa)')
    ylabel('Probability( |data| > amplitude )')
end



end

