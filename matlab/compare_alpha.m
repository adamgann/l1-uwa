% compare_alpha.m
%
% Generate plots showing alpha versus time from saved runs of the 
% script equalizer_comparison_simulation.m
%
% Adam Gannon, SUNY Buffalo, 2018.

clear all;
close all;
clc

addpath(genpath('functions/'))

%% Load Data

a180Name = 'data/april18/workspace_alpha18_npacket500.mat';
a200Name = 'data/april18/workspace_alpha20_npacket_500.mat';
a150Name = 'data/april18/workspace_alpha150_npacket_500.mat';
a170Name = 'data/april26/workspace_alpha170.mat';
a160Name = 'data/april26/workspace_alpha160.mat';
a190Name = 'data/april26/workspace_alpha190.mat';

%% Create cells 

dataFileNames = {};
dataFileNames{1} = a150Name;
dataFileNames{2} = a160Name;
dataFileNames{3} = a170Name;
dataFileNames{4} = a180Name;
dataFileNames{5} = a190Name;
dataFileNames{6} = a200Name;

nPoints = size(dataFileNames,2);


mfPilot = cell(1,nPoints);
l1Blind = cell(1,nPoints);
l2Blind = cell(1,nPoints);
snrVecStore = cell(1,nPoints);
alphaStore = cell(1,nPoints);

%%



for iData = 1:nPoints

    load(dataFileNames{iData})
    mfPilot{iData} = berMfStore;
    l1Blind{iData} = berL1Store;
    l2Blind{iData} = berL2Store;
    snrVecStore{iData} = snrVec;
    alphaStore{iData} = shrimpAlpha;

end


%% Calculate the SNR loss 

berTarget = (5e-3);

mfVec = zeros(1,nPoints);
l1Vec = zeros(1,nPoints);
l2Vec = zeros(1,nPoints);

for iData = 1:nPoints

mfVec(iData) = get_closest_snr_to_ber_target(snrVecStore{iData},mfPilot{iData},berTarget);
l1Vec(iData) = get_closest_snr_to_ber_target(snrVecStore{iData},l1Blind{iData},berTarget);
l2Vec(iData) = get_closest_snr_to_ber_target(snrVecStore{iData},l2Blind{iData},berTarget);
end


l1Diff = l1Vec - mfVec;
l2Diff = l2Vec - mfVec;

alphaVec = cell2mat(alphaStore);

%% Plot the SNR loss

figure;
hL(1) = plot(alphaVec,l1Diff,'s-');
hold on
hL(2) = plot(alphaVec,l2Diff,'o-');


hY=ylabel(('Transmit Energy Difference from Pilot ($\mathrm{dB\:re\ \mu Pa}$)'))
hX=xlabel('$\alpha$');
hLeg = legend('L1-PCA (Proposed)','L2-PCA');

% Format fonts
set(hLeg,'Interpreter','Latex','Fontsize',12);
set(hX,'Interpreter','Latex','Fontsize',12);
set(hY,'Interpreter','Latex','Fontsize',12);

for hh=1:length(hL)
    set(hL(hh),'LineWidth',1.5);
    set(hL(hh),'MarkerSize',10);
end


