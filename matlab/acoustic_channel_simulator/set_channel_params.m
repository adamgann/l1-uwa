% This code reads the channel parameters and writes them in data 
% files [file_name].prm and [file_name].dop which will be used 
% by the main simulator code channel_simulator.
%
% Change the following channel parameters by updating this code. 
% Alternatively, you can directly update the .prm and .dop data 
% files via a text editing software.
% 
% file_name: data file name
% h0: surface height (depth) [m]
% ht0: TX height [m]
% hr0: RX height [m]
% d0: channel distance [m]
% k: spreading factor
% c: speed of sound in water [m/s]
% c2: speed of sound in bottom [m/s] (>1500 for hard, < 1500 for soft)
% cut: do not consider arrivals whose strength is below that of direct arrival divided by cut
% fmin: minimum frequency [Hz]
% B: bandwidth [Hz]
% df: frequency resolution [Hz]
% dt: time resolution [seconds] 
% T_SS: coherence time of the Small-Scale (S-S) variations [seconds]
% sig2s: variance of S-S surface variations 
% sig2b: variance of S-S bottom variations 
% B_delp: 3-dB width of the p.s.d. of intra-path delays (assumed constant for all paths)
% Sp: number of intra-paths (assumed constant for all paths)
% mu_p: mean of intra-path amplitudes (assumed constant for all paths)
% nu_p: variance of intra-path amplitudes (assumed constant for all paths)
% T_tot: total duration of the simulated signal [seconds]
% h_bnd: range of surface height variation (Large-Scale (L-S) realizations are limited to h+h_band)
% ht_bnd: range of transmitter height variation
% hr_bnd: range of receiver height variation
% d_bnd: range of channel distance variation
% sig_h: standard deviation of L-S variations of surface height 
% sig_ht: standard deviation of L-S variations of transmitter height 
% sig_hr: standard deviation of L-S variations of receiver height 
% sig_d: standard deviation of L-S variations of distance height 
% a_AR: AR parameter for generating L-S variations (constant for variables h, ht, hr, d) 
% method: Large-Scale (L-S) & Small-Scale (S-S) Simulation Methods:
% vtd: TX drifting speed
% theta_td: TX drifting angle
% vrd: RX drifting speed
% theta_rd: RX drifting angle
% vtv: TX vehicular speed
% theta_tv: TX vehicular angle
% vrv: RX vehicular speed
% theta_rv: RX vehicular angle 
% Aw: surface variation amplitude
% fw: surface variation frequency
%
% Author: Parastoo Qarabaqi, Northeastern University, 2013
%
% Modifications for this project: Adam Gannon, UB, 2017.

clc, clear all;
format long;

addpath('channels/')

%% Data file name: 
file_name= 'channels/chan_wide_0';

%% Deterministic channel geometry:

% Geometry of our channel 
h0=20; % surface height (depth) [m]
ht0=12; % TX height [m]
hr0=12; % RX height [m]
d0=500; % channel distance [m]

% Constants describing the channel environment
k=1.7; % spreading factor
c=1500; % speed of sound in water [m/s]
c2=1300; % speed of sound in bottom [m/s] (>1500 for hard, < 1500 for soft)
cut=10; % do not consider arrivals whose strength is below that of direct arrival divided by cut

% Bandwidth of 97.6kHz at a center frequency of 100kHz. 
B=97.6e3; % bandwidth [Hz]
fmin=(100e3-(B/2)); % minimum frequency [Hz]

% Calculate frequency resolution from desired number of taps
num_taps = 200; % Number of taps to include 
df= ceil(B/num_taps); % frequency resolution [Hz], f_vec=fmin:df:fmax; 

% Make the coherence time exact multiples of dt
coherence_time = 200e-3;
dt = coherence_time/10; % Each vector is indepedent for BER curves
T_SS = floor(coherence_time/dt)*dt  % coherence time of the small-scale variations [seconds] 


total_time = 2000; % In seconds
T_tot= floor(total_time/T_SS)*T_SS

%% Small-Scale (S-S) parameters: 
sig2s= 1.125; % variance of S-S surface variations 
sig2b=sig2s/2; % variance of S-S bottom variations 
 
B_delp= 5e-4; % 3-dB width of the p.s.d. of intra-path delays (assumed constant for all paths)
Sp= 20; % number of intra-paths (assumed constant for all paths)
mu_p= .5/Sp; % mean of intra-path amplitudes (assumed constant for all paths)
nu_p= 1e-6; % variance of intra-path amplitudes (assumed constant for all paths)

%% Large-Scale (L-S) parameters:
%T_tot= 100*T_SS; % total duration of the simulated signal [seconds]
t_tot_vec=(0:dt:T_tot); Lt_tot=length(t_tot_vec);

h_bnd=[-0.0001 0.0001]; % range of surface height variation (L-S realizations are limited to h+h_band)
ht_bnd=[-0.0001 0.0001]; % range of transmitter height variation
hr_bnd=[-0.0001 0.0001]; % range of receiver height variation
d_bnd=[-10 100]; % range of channel distance variation
sig_h=1; % standard deviation of L-S variations of surface height 
sig_ht=1; % standard deviation of L-S variations of transmitter height 
sig_hr=1; % standard deviation of L-S variations of receiver height 
sig_d=1; % standard deviation of L-S variations of distance height 
a_AR= .9; % AR parameter for generating L-S variations (constant for variables h, ht, hr, d) 

%% Large-Scale (L-S) & Small-Scale (S-S) Simulation Methods:
% 1)'sim_dir': SIMplified L-S model & DIRect S-S model
% 2)'sim_seq': SIMplified L-S model & Statistically EQuivalent S-S model
% 3)'bel_dir': BELlhop L-S model & DIRect S-S model
% 4)'bel_seq': BELlhop L-S model & Statistically EQuivalent S-S model
method='sim_dir';%'sim_dir'; 

%% Doppler: 
% All parameters can be entered as scalars or vectors of size(t_tot_vec).


% drifting: 
vtd= 0.1*cos((1/20)*(2*pi*t_tot_vec)-rand*2*pi); 
theta_td= linspace(rand*2*pi,rand*2*pi,length(t_tot_vec)); 

vrd= 0.02*cos((1/5)*(2*pi*t_tot_vec)); 
theta_rd= rand*2*pi; 

% vehicular: 
% Assume no vehicular motion for these experiments 
vtv= 0;   
theta_tv= 0; 
vrv= 0;  
theta_rv= 0; 

% Assume no drifting for these experiments 
vtd = 0;
vrd = 0;

% surface: 
Aw= .05; fw= .01; 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of reading data. Start writing the data files: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Write parameters on prm file: 
f1 = fopen([file_name, '.prm'], 'w');

fprintf(f1, '%.32f surface_height_(depth)_[m]\n', h0); 
fprintf(f1, '%.32f TX_height_[m]\n', ht0); 
fprintf(f1, '%.32f RX_height_[m]\n', hr0); 
fprintf(f1, '%.32f channel_distance_[m]\n', d0); 

fprintf(f1, '%.32f spreading_factor\n', k); 
fprintf(f1, '%.32f speed_of_sound_in_water_[m/s]\n', c); 
fprintf(f1, '%.32f speed_of_sound_in_bottom_[m/s]\n', c2); 
fprintf(f1, '%.32f minimum_relative_path_strength\n', cut); 

fprintf(f1, '%.32f minimum_frequency_[Hz]\n', fmin); 
fprintf(f1, '%.32f bandwidth_[Hz]\n', B); 
fprintf(f1, '%.32f frequency_resolution_[Hz]\n', df); 
fprintf(f1, '%.32f time_resolution_[seconds]\n', dt); 
fprintf(f1, '%.32f coherence_time_of_the_small-scale_variations_[seconds]\n', T_SS); 

fprintf(f1, '%.32f variance_of_S_S_surface_variations_[m^2]\n', sig2s); 
fprintf(f1, '%.32f variance_of_S_S_bottom_variations_[m^2]\n', sig2b); 
fprintf(f1, '%.32f 3-dB_width_of_the_p.s.d._of_intra-path_delays_[Hz]\n', B_delp); 
fprintf(f1, '%d number_of_intra-paths\n', Sp); 
fprintf(f1, '%.32f mean_of_intra-path_amplitudes\n', mu_p); 
fprintf(f1, '%.32f variance_of_intra-path_amplitudes\n', nu_p); 

fprintf(f1, '%.32f total_duration_of_the_simulated_signal_[seconds]\n', T_tot); 
fprintf(f1, '%.32f min_of_surface_height\n', h_bnd(1));
fprintf(f1, '%.32f max_of_surface_height\n', h_bnd(2));
fprintf(f1, '%.32f min_of_transmitter_height\n', ht_bnd(1));
fprintf(f1, '%.32f max_of_transmitter_height\n', ht_bnd(2));
fprintf(f1, '%.32f min_of_receiver_height\n', hr_bnd(1));
fprintf(f1, '%.32f max_of_receiver_height\n', hr_bnd(2));
fprintf(f1, '%.32f min_of_channel_distance\n', d_bnd(1));
fprintf(f1, '%.32f max_of_channel_distance\n', d_bnd(2));

fprintf(f1, '%.32f L_S_standard_deviation_of_surface_height_[m]\n', sig_h); 
fprintf(f1, '%.32f L_S_standard_deviation_of_transmitter_height_[m]\n', sig_ht); 
fprintf(f1, '%.32f L_S_standard_deviation_of_receiver_height_[m]\n', sig_hr); 
fprintf(f1, '%d L_S_standard_deviation_of_receiver_height_[m]\n', sig_d); 
fprintf(f1, '%.32f AR_parameter_for_generating_L_S_variations\n', a_AR); 

switch method
    case {'sim_dir', 1}, fprintf(f1, '1 method\n');
    case {'sim_seq', 2}, fprintf(f1, '2 method\n');
    case {'bel_dir', 3}, fprintf(f1, '3 method\n');
    case {'bel_seq', 4}, fprintf(f1, '4 method\n');        
end

fclose(f1);

%% form Doppler matrix: 
Dopp_mat= zeros(Lt_tot, 10);
if length(vtd)==1; Dopp_mat(:,1)= repmat(vtd, Lt_tot, 1); 
else Dopp_mat(:,1)= vtd.';
end
if length(theta_td)==1; Dopp_mat(:,2)= repmat(theta_td, Lt_tot, 1); 
else Dopp_mat(:,2)= theta_td.';
end
if length(vrd)==1; Dopp_mat(:,3)= repmat(vrd, Lt_tot, 1); 
else Dopp_mat(:,3)= vrd.';
end
if length(theta_rd)==1; Dopp_mat(:,4)= repmat(theta_rd, Lt_tot, 1); 
else Dopp_mat(:,4)= theta_rd.';
end
if length(vtv)==1; Dopp_mat(:,5)= repmat(vtv, Lt_tot, 1); 
else Dopp_mat(:,5)= vtv.';
end
if length(theta_tv)==1; Dopp_mat(:,6)= repmat(theta_tv, Lt_tot, 1); 
else Dopp_mat(:,6)= theta_tv.';
end
if length(vrv)==1; Dopp_mat(:,7)= repmat(vrv, Lt_tot, 1); 
else Dopp_mat(:,7)= vrv.';
end
if length(theta_rv)==1; Dopp_mat(:,8)= repmat(theta_rv, Lt_tot, 1); 
else Dopp_mat(:,8)= theta_rv.';
end
if length(Aw)==1; Dopp_mat(:,9)= repmat(Aw, Lt_tot, 1); 
else Dopp_mat(:,9)= Aw.';
end
if length(fw)==1; Dopp_mat(:,10)= repmat(fw, Lt_tot, 1); 
else Dopp_mat(:,10)= fw.';
end


%% Write Doppler rates on dop file: 
f2 = fopen([file_name, '.dop'], 'w');

fprintf(f2, '%f\n', reshape(Dopp_mat, 10*Lt_tot, 1)); 

fclose(f2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the channel parameters:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the Doppler rates: 

% drifting:
figure; subplot(411); plot(t_tot_vec, (1/c)*Dopp_mat(:,1).*cos(Dopp_mat(:,2))); 
ylabel('TX drifting'), set(gca,'XTick', (0:T_SS:T_tot)), 
title('Doppler rates relative to horizontal')
subplot(412); plot(t_tot_vec, (1/c)*Dopp_mat(:,3).*cos(Dopp_mat(:,4))); 
ylabel('RX drifting'), set(gca,'XTick', (0:T_SS:T_tot)),

% vehicular:
subplot(413); plot(t_tot_vec, (1/c)*Dopp_mat(:,5).*cos(Dopp_mat(:,6))); 
ylabel('TX vehicular'), set(gca,'XTick', (0:T_SS:T_tot)), 
subplot(414); plot(t_tot_vec, (1/c)*Dopp_mat(:,7).*cos(Dopp_mat(:,8))); 
ylabel('RX vehicular'), set(gca,'XTick', (0:T_SS:T_tot)),
xlabel('time [s]')

%% Plot the nominal channel geometry: 

figure; 
plot([0,0], [0,ht0], 'r'), hold on,
plot(0,ht0, 'rx', 'linewidth', 2, 'markersize', 10), 
plot([d0,d0], [0,hr0], 'r'), 
plot(d0,hr0, 'rx', 'linewidth', 2, 'markersize', 10),
plot([-.1*d0, 1.1*d0], [0, 0], 'k', 'linewidth', 2),
plot([-.1*d0, 1.1*d0], [h0,h0], 'b', 'linewidth', 2),

axis([-.1*d0, 1.1*d0, -.1*h0, 1.1*h0]),
xlabel('range [m]'), ylabel('height above bottom [m]'),
text(.02*d0, ht0, 'TX'), text(1.02*d0, hr0, 'RX'), 
title('Nominal channel geometry')


