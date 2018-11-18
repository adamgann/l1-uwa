% shrimp_noise.m
%
% Generate plots from other research to verify that our algorithm is
% working correctly. Looking specifically at:
%
% (Mahmood2015) A. Mahmood and M. Chitre, "Modeling colored impulsive 
% noise by Markov chains and alpha-stable processes," OCEANS 2015. 
%
% (Chitre2006) M. A. Chitre, J. R. Potter, S-H. Ong, "Optimal and 
% Near-Optimal Signal Detection in Snapping Shrimp Dominated Ambient 
% Noise," JOE 2006.
%
% Requires: Statistics and Machine Learning Toolbox
%           Mark Veillette's stbl library (included in repo)
%           R. J. Achatz's amplitude probability dist function (in repo)
%
% Adam Gannon, SUNY Buffalo, 2018


clear all;
close all;
clc

addpath(genpath('functions/'))


%% Generate Fig. 6 of Mahood2015 

% Given parameters from Mahood2015
alpha = 1.82;
R = [1, 0.7; 0.7, 1];

nSamps = 1e5;


% Generate Gaussian vector G
mu = zeros(1,2);
G = mvnrnd(mu,R,nSamps);

% Generate the A parameter (Described in A. Mahmood, M. Chitre "Generating 
% Radom Variables for Stable Sub-Gaussian Processes with Memory", 2017.)
scaleParam = 2*cos((pi*alpha)/4)^(2/alpha);
A = stblrnd(alpha/2, 1, scaleParam, 0, nSamps, 1);


% Multiply to generate RV X, the bivariate SAS
pdfX = repmat(sqrt(A),1,2).*G;


% Replicate Fig. 6
figH = figure;

subplot(1,2,1)
scatter(pdfX(:,1),pdfX(:,2))

xlim([-15 15])
ylim([-15 15])
grid on
xlabel('x_1')
ylabel('x_2')

subplot(1,2,2)
scatter(G(:,1),G(:,2))

xlim([-15 15])
ylim([-15 15])
grid on
xlabel('g_1')
ylabel('g_2')


%% Generate Fig. 2b of Optimal.Near-Optimal Signal Detection (Chitre2006)

% Using parameters from Chirtre2006
alpha = 1.82;

%This is c in Chirtre's notation. You can verify this by comparing the
%characteristic function (1) in Chitre2006 with the charateristic function
% for the general case (a != 1) in the README of Veillette's code. 
gamma = 1.5e5;

% Use stblrnd to generate a symmetric (beta=delta=0)
pdfX = stblrnd(alpha,0,gamma,0,nSamps,1);
[aX,pX] = apd(abs(pdfX));

% Get best Gaussian fit
% Occasionally normfit totally goes off in the weeds. May need to run again
% if results make no sense.
[mu,sigma] = normfit(pdfX);


% Generate Gaussian APD
G = sigma.*randn(nSamps,1) + mu;
[aG,pG] = apd(abs(G));


figure
semilogy(aX,pX)
hold on
semilogy(aG,pG,'--')


ylim([1e-2 1])
xlim([0 6e5])
grid on

xlabel('Amplitude (uPa)')
ylabel('Probability( |data| > amplitude )')
legend('SaS Noise','Gaussian Noise')

%% Generate Fig. 1 of Chitre2006

figure;
normplot(pdfX./(1e6))
xlim([-1 1])
xlabel('Pressure (Pa)')
title('')


%% Generate a combined plot that shows APD and a time series for paper

figure;

subplot(2,1,1)
plot(pdfX(1:0.5e4))

xH = xlabel('Sample Index')
yH = ylabel('Amplitude ($\mu$Pa)')

set(xH,'Interpreter','latex','Fontsize',14)
set(yH,'Interpreter','latex','Fontsize',14)

% APD
subplot(2,1,2)
semilogy(aX,pX)

ylim([1e-2 1])
xlim([0 6e5])
grid on

xH = xlabel('Amplitude ($\mu$Pa)')
yH = ylabel('Pr($|$data$|>$amplitude)')

set(xH,'Interpreter','latex','Fontsize',14)
set(yH,'Interpreter','latex','Fontsize',14)
