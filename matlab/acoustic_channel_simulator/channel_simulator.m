% This code generates an ensemble of channel impulse responses.
% Data files [file_name].prm and [file_name].dop provide the 
% channel parameters and Doppler information. 
% 
% Before running this code: 
% -If you wish to use Bellhop for Large-Scale (L-S) simulations, 
% edit the function runBellhop.m to suit to your desired channel 
% properties, e.g. bathymetry and sound speed profile. 
% Alternatively, a simplified simulator mpgeometry.m enclosed in 
% the Matlab functions folder can be used which assumes flat 
% bathymetry, surface shape and sound speed profile.
% -Set the .prm and .dop data file name which you wish to use. 
%
% The following outputs are generated:
% hmat: matrix of T_tot/dt channel impulse response columns of delay resolution 1/df
%       which is saved in a Matlab data file [file_name].mat.
% Gtilde: vector of T_tot/dt instantaneous channel gains. 
%
% Author: Parastoo Qarabaqi, Northeastern University, 2013
%
% Modified for this project: Adam Gannon at UB, 2017. 

close all
clc, clear all, warning off, 
tstart= tic;
cd(fileparts(which('channel_simulator')));
path(genpath([pwd,'/MatlabFunctions']), path);


%% read data from prm file: 
file_name= 'channels/chan_wide_0';
[parameters, descriptions] = textread([file_name, '.prm'], '%f %s');


%% Deterministic channel parameters:

h0= parameters(1); % surface height (depth) [m]
ht0= parameters(2); % TX height [m]
hr0= parameters(3); % RX height [m]
d0= parameters(4); % channel distance [m]

k= parameters(5); % spreading factor
c= parameters(6); % speed of sound [m/s]
c2= parameters(7); % speed of sound in bottom [m/s] (>1500 for hard, < 1500 for soft)
cut= parameters(8); % do not consider arrivals whose relative strength is below this level

fmin= parameters(9); % minimum frequency [Hz]
B= parameters(10); % bandwidth [Hz]
fmax=fmin+B; 

df= parameters(11); % frequency resolution [Hz]
f_vec=(fmin:df:fmax).'; Lf=length(f_vec);
fc=(fmax+fmin)/2; 
f0=fmin;
f_vec2=(f0-B/2:df:f0+B+B/2); dif_f=f_vec2-fc;

dt= parameters(12); % time resolution [seconds] 
T_SS= parameters(13); % coherence time of the small-scale variations [seconds]
t_vec=(0:dt:T_SS); Lt=length(t_vec);

%% Small-scale channel parameters: 
sig2s= parameters(14); % variance of SS surface variations (surfampl^2/2?)
sig2b= parameters(15); % variance of SS bottom variations 

B_delp= parameters(16); % 3-dB width of the p.s.d. of intra-path delays (assumed constant for all paths)
Sp= parameters(17); % number of intra-paths (assumed constant for all paths)
mu_p= parameters(18); % mean of intra-path amplitudes (assumed constant for all paths)
nu_p= parameters(19); % variance of intra-path amplitudes (assumed constant for all paths)

%% Large-scale channel parameters:
T_tot= parameters(20); % total duration of the simulated signal [seconds]
N_LS=round(T_tot/T_SS); 
t_tot_vec=(0:dt:T_tot); Lt_tot=length(t_tot_vec);


h_bnd= parameters(21:22).'; % range of surface height (LS realizations \in h+h_band)
ht_bnd= parameters(23:24).'; % range of transmitter height
hr_bnd= parameters(25:26).'; % range of receiver height
d_bnd= parameters(27:28).'; % range of channel distance

sig_h= parameters(29); % standard deviation of LS variations of surface height 
sig_ht= parameters(30); % standard deviation of LS variations of transmitter height 
sig_hr= parameters(31); % standard deviation of LS variations of receiver height 
sig_d= parameters(32); % standard deviation of LS variations of distance height 
a_AR= parameters(33); % AR parameter for generating L-S variations (constant for variables h, ht, hr, d) 

%% Large-Scale (L-S) & Small-Scale (S-S) Simulation Methods:
% 1)'sim_dir': SIMplified L-S model & DIRect S-S model
% 2)'sim_seq': SIMplified L-S model & Statistically EQuivalent S-S model
% 3)'bel_dir': BELlhop L-S model & DIRect S-S model
% 4)'bel_seq': BELlhop L-S model & Statistically EQuivalent S-S model
switch  parameters(34)
    case 1, method='sim_dir';
    case 2, method='sim_seq';
    case 3, method='bel_dir';
    case 4, method='bel_seq';
end
method_LS= lower(method(1:3));
method_SS= lower(method(5:7));

%% Doppler parameters: 
[Dopp_params] = textread([file_name, '.dop'], '%f');
Dopp_params= reshape(Dopp_params, Lt_tot, 10);

% drifting: 
vtd_tot= Dopp_params(:,1).'; 
theta_td_tot= Dopp_params(:,2).'; 
vrd_tot= Dopp_params(:,3).';
theta_rd_tot= Dopp_params(:,4).';

% vehicular: 
vtv_tot= Dopp_params(:,5).';
theta_tv_tot= Dopp_params(:,6).';
vrv_tot= Dopp_params(:,7).';
theta_rv_tot= Dopp_params(:,8).';
 
% surface: 
Aw_tot= Dopp_params(:,9).';
fw_tot= Dopp_params(:,10).'; 

%% Large-scale loop: 
H_LS= zeros(Lf, Lt*N_LS);
del_h=0; del_ht=0; del_hr=0; del_d=0; 
h=h0; ht=ht0; hr=hr0; d= d0; % initial values

%% initialize plotting of large-scale channel parameters: 
% figure(1), 
% subplot(411), hold on, plot(0,h0, 'r*'), box on, ylabel('h [m]'), ylim(h0+h_bnd),
% subplot(412), hold on, plot(0,ht0, 'r*'), box on, ylabel('h_t [m]'), ylim(ht0+ht_bnd), 
% subplot(413), hold on, plot(0,hr0, 'r*'), box on, ylabel('h_r [m]'), ylim(hr0+hr_bnd), 
% subplot(414), hold on, plot(0,d0, 'r*'), box on, ylabel('d [m]'), xlabel('Large-scale realization number'), ylim(d0+d_bnd), 

adopp0=zeros(1,50);

for LScount= 1:N_LS
rndvec= randn(1, 4); 

% figure(1), subplot(411), hold on,
del_h= a_AR * del_h + sqrt(1-a_AR^2)*sig_h*rndvec(1); 
if (del_h > h_bnd(2))||(del_h < h_bnd(1)), 
    del_h= del_h- 2*sqrt(1-a_AR^2)*sig_h*rndvec(1); 
end
htemp=h;
h= h0+del_h;
% plot([LScount-1, LScount], [htemp, h], 'b')
% ylim(h0+h_bnd);

% subplot(412), hold on,
del_ht= a_AR * del_ht + sqrt(1-a_AR^2)*sig_ht*rndvec(2); 
if (del_ht > ht_bnd(2))||(del_ht < ht_bnd(1)), 
    del_ht= del_ht- 2*sqrt(1-a_AR^2)*sig_ht*rndvec(2); 
end
httemp=ht;
ht= ht0+del_ht;
% plot([LScount-1, LScount], [httemp, ht], 'b')
% ylim(ht0+ht_bnd);

% subplot(413), hold on,
del_hr= a_AR * del_hr + sqrt(1-a_AR^2)*sig_hr*rndvec(3); 
if (del_hr > hr_bnd(2))||(del_hr < hr_bnd(1)), 
    del_hr= del_hr- 2*sqrt(1-a_AR^2)*sig_hr*rndvec(3); 
end
hrtemp=hr;
hr= hr0+del_hr;
% plot([LScount-1, LScount], [hrtemp, hr], 'b')
% ylim(hr0+hr_bnd);

% subplot(414), hold on,
del_d= a_AR * del_d + sqrt(1-a_AR^2)*sig_d*rndvec(4); 
if (del_d > d_bnd(2))||(del_d < d_bnd(1)), 
    del_d= del_d- 2*sqrt(1-a_AR^2)*sig_d*rndvec(4); 
end
dtemp=d;
d= d0+del_d;
% plot([LScount-1, LScount], [dtemp, d], 'b')
% ylim(d0+d_bnd);


% Find Large-scale model parameters: 

switch method_LS
case 'sim'
    [lmean,taumean,Gamma,theta,ns,nb,hp]= mpgeometry(h,h-ht,h-hr,d,fc,k,cut,c,c2);

case 'bel'
    [lmean,taumean,theta,ns,nb,hp,p0_ind]=runBellhop(h,h-ht,h-hr,d,fc,c,c2);
    lmean=[lmean(p0_ind), lmean(1:p0_ind-1), lmean(p0_ind+1:end)];
    taumean=[taumean(p0_ind), taumean(1:p0_ind-1), taumean(p0_ind+1:end)];
    theta=[theta(p0_ind), theta(1:p0_ind-1), theta(p0_ind+1:end)];
    ns=[ns(p0_ind), ns(1:p0_ind-1), ns(p0_ind+1:end)];
    nb=[nb(p0_ind), nb(1:p0_ind-1), nb(p0_ind+1:end)];
    hp=[hp(p0_ind), hp(1:p0_ind-1), hp(p0_ind+1:end)];

end

% ignore paths with delays longer than allowed by frequency resolution: 
lmean= lmean(taumean<1/df); 
theta= theta(taumean<1/df); 
ns= ns(taumean<1/df); 
nb= nb(taumean<1/df); 
hp= hp(taumean<1/df); 
taumean= taumean(taumean<1/df); 

P=length(lmean); % number of paths

% Reference path transfer function: 
H0= 1./sqrt(lmean(1)^k* (10.^(absorption(f_vec/1000)/10000)).^lmean(1) );
H= hp(1)*repmat( exp(-1j*2*pi*f_vec*taumean(1)) , 1, Lt);

%% Find small-scale model parameters:

sig_delp= sqrt(1/c^2*((2*sin(theta)).^2).*(ns*sig2s+nb*sig2b)); 
rho_p= exp(-((2*pi*f_vec).^2) * (sig_delp.^2/2));
rho_p_bb= exp(-(2*pi*(f_vec-fmin).^2) * (sig_delp.^2/2));


Bp= ((2*pi*f_vec*sig_delp).^2).*B_delp;
sig_p= sqrt(.5*(mu_p.^2.*Sp.*(1-rho_p.^2)+Sp.*nu_p.^2));

%% Find doppler rates:
 
% drifting: 
vtd= vtd_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1)); 
theta_td= theta_td_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
vrd= vrd_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
theta_rd= theta_rd_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));

% vehicular: 
vtv= vtv_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
theta_tv= theta_tv_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
vrv= vrv_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
theta_rv= theta_rv_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
 
% surface: 
Aw= Aw_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
fw= fw_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
vw= 2*pi*fw.*Aw;


% first path doppler:

vdrift= vtd.*cos(theta(1)-theta_td)-vrd.*cos(theta(1)+theta_rd); 
adrift= vdrift/c; 

vvhcl= 0; 
avhcl= vvhcl/c;

vsurf=  0; 
asurf= vsurf/c; 


adopp= adrift + avhcl + asurf*ns(1); 
eff_adopp = adopp0(1)+cumsum(adopp);
Dopp= exp(1j*2*pi*f_vec*(eff_adopp*dt)); 
adopp0(1)=eff_adopp(end);
H= H.*Dopp;


%% small-scale simulation:

for p=2:P; 

switch method_SS
case 'dir', % DIRectly generate gamma: 

    gamma=zeros(Lf, Lt);
    for counti=1:Sp;
        gamma_pi= mu_p+nu_p*randn(1,Lt); 
        gamma_pi=repmat(gamma_pi, Lf, 1)*Sp; 
        deltau_pi=zeros(Lf, Lt);
        w_delpi= sig_delp(p)*sqrt(1-exp(-1*pi*B_delp*dt)^2)*randn(1, 2*Lt);
        temp_deltau_pi=filter(1, [1, -exp(-pi*B_delp*dt)], w_delpi);
        for countf= 1:Lf
            deltau_pi(countf, :)= temp_deltau_pi(Lt+1:end);
        end
        gamma=gamma+ gamma_pi.*exp(-1j*2*pi*repmat(f_vec, 1,Lt).*deltau_pi);
    end
    
case 'seq', % use the Statistically EQuivalent model:
    if sig_delp(p)*f0 < .4,
    fprintf(['\n Warning! sig_delp*f0 < 0.4 for path ' num2str(p),' of L-S realization ', num2str(LScount),'. The statistically equivalent model is not accurate. \n'])
    end
    alpha_p= exp(-pi*Bp(:,p)*dt); 
    Alpha_p= diag(alpha_p);    
    
    Wp_cov= 2*(ones(Lf, Lf)-alpha_p*alpha_p.').*toeplitz(rho_p_bb(:,p)).*(sig_p(:,p)*sig_p(:,p).');

    error_takagi=-1;
    while sign(error_takagi)~=1;
        [Wp_Q, Wp_s, error_takagi]= find_takagi_factor(Wp_cov);  
    end
    
    Wp_factor= Wp_Q*sqrt(diag(Wp_s)); 
    wp= Wp_factor*(randn(Lf,Lt) + 1j*randn(Lf,Lt));
        
    gammabar=repmat(mu_p + mu_p*Sp*rho_p(:,p), 1, Lt);

    Delgamma=zeros(Lf,Lt+1);
    for count_t=1:Lt
        Delgamma(:,count_t+1)= Alpha_p*Delgamma(:,count_t)+wp(:,count_t);
    end
    Delgamma=Delgamma(:, 2:end);

   gamma= gammabar+Delgamma; 

end

    %% Doppler term:
    vdrift= vtd.*cos(theta(p)-theta_td)-vrd.*cos(theta(p)+theta_rd); 
    adrift= vdrift/c; 
    
    vvhcl= vtv.*cos(theta(p)-theta_tv)-vrv.*cos(theta(p)+theta_rv)-(vtv.*cos(theta(1)-theta_tv)-vrv.*cos(theta(1)+theta_rv));
    avhcl= vvhcl/c;
    
    phi_pj= 2*pi*rand(1, ns(p))-pi;
    sum_j= zeros(1, Lt); 
    for jcount= 1: ns(p);
        sum_j= sum_j + sin(phi_pj(jcount)+2*pi*fw.*t_vec);
    end
    vsurf= 2*vw.*sin(theta(p)).*sum_j;
    asurf= vsurf/c; 
    
    adopp= adrift + avhcl + asurf*ns(p);
    eff_adopp = adopp0(p)+cumsum(adopp);
    Dopp= exp(1j*2*pi*f_vec*eff_adopp*dt); 
    adopp0(p)=eff_adopp(end);
    % If adopp was constant, we would have: 
    % Dopp= exp(1j*2*pi*repmat(adopp, Lf, 1).*repmat(f0, Lf, Lt));     
    
    % Multiply gamma by hp:
    gamma_noDopp=gamma; 
    gamma=gamma.*Dopp;
    H= H+ hp(p)*repmat( exp(-1j*2*pi*f_vec*taumean(p)) , 1, Lt).*gamma;    

end

H= repmat(H0 ,1, Lt).*H;
H_LS(:, (LScount-1)*Lt+1:(LScount)*Lt)= H; 

disp(sprintf('Iteration %d/%d',LScount,N_LS))
end


%% find channel impulse response: 

Lt_tot= size(H_LS, 2);
hmat=zeros(Lf,Lt_tot);
% N=ceil(log10(size(H_LS,1))/log10(2)); 
for countt=1:Lt_tot
    hmat(:, countt)= ifft(H_LS(:, countt)); %hmat(:, countt)=fftshift(hmat(:, countt));
end

telapsed = toc(tstart);
fprintf(['Total runtime is ', num2str(telapsed), ' seconds.\n']),

shift=10; skip=10; 
figure; axes('fontsize', 16);
% hmat2plot= abs(hmat(1:end, 1:skip:end)).';

xVec = ((0:Lf-1)-shift)/B*1000;
yVec = (0:skip:Lt_tot-1)*dt;
cData = (circshift(abs(hmat(1:end, 1:skip:end)), shift)).';

image(xVec, yVec , cData ,  'CDataMapping','scaled');
xlabel('delay [ms]', 'fontsize', 16), ylabel('time [s]', 'fontsize', 16), axis('ij'); 
set(gca,'YTick', 0:T_SS:T_tot);
colorbar;

if sum((ns==0)&(nb==0))>0
text(min(1000*taumean((ns==0)&(nb==0))), -.05*T_tot, 0, 'p_0', 'fontsize', 14)
end
if sum((ns==1)&(nb==0))>0
text(min(1000*taumean((ns==1)&(nb==0))), -.05*T_tot, 'p_s', 'fontsize', 14)
end
if sum((ns==0)&(nb==1))>0
text(min(1000*taumean((ns==0)&(nb==1))), -.05*T_tot, 'p_b', 'fontsize', 14)
end

%% Plot for paper

yLimit = 60;
yInd = find(yVec==yLimit);

figure;
image(xVec, yVec(1:yInd) , cData(1:yInd,:) ,  'CDataMapping','scaled');
ylim([0 60]); xlim([-inf 1.85])

hY = ylabel('Time (s)');
hX = xlabel('Delay (ms)');
set(hX,'Interpreter','Latex','Fontsize',12);
set(hY,'Interpreter','Latex','Fontsize',12);


text(-.02,-1.5,'D','Interpreter','Latex','Fontsize',12)
text(0.15,-1.5,'S','Interpreter','Latex','Fontsize',12)
text(0.32,-1.5,'B','Interpreter','Latex','Fontsize',12)
text(0.9,-1.5,'B-S/S-B','Interpreter','Latex','Fontsize',12)
%% find channel gain: 

hmat2=hmat.'; 
Gtilde=zeros(length(hmat2),1);
for countt= 1:size(hmat2,1)
    Gtilde(countt)= ((hmat2(countt, :)*hmat2(countt, :)'));
end

figure; plot((0:length(Gtilde)-1)*dt, 10*log10(Gtilde)),
xlabel('time [s]'), ylabel('instantaneous channel gain [dB]')

%% Calculate prop time between transmitters 

SPEED_OF_SOUND = 1500; %m/s

dist_euc = sqrt( (d0).^2 + (ht0-hr0).^2 );
T_prop = dist_euc / SPEED_OF_SOUND; % Prop time in seconds

%% Cut down channel matrix to independent realizations for BER curves

nCut = (200e-3)/dt
hmat=hmat(:,1:nCut:end);
size(hmat)


save([file_name, '.mat'], 'hmat', 'dt', 'df', 'T_SS', 'dt', 'T_tot', 'T_prop')
