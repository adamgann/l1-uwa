% mpgeometry.m
% geometric multipath structure
function [l,tau,Gamma,theta,ns,nb,hp]=mpgeometry(h,ht,hr,d,f,k,cut,c,c2) 

% Input:
% h: depth of water [m]
% ht: depth of tx [m]
% hr: depth of rx [m]
% d: distance [m]
% c: speed of sound in water [m/s], used to calculate reflcoef 
% c2: speed of sound in bottom [m/s], used to calculate reflcoef (1300 for soft, 1800 for hard bottom)
% %reflloss [dB]: can use fixed per-bottom-bounce loss instead
% k: spreading factor
% f: frequency [kHz]: use lower band-edge freq. to get max. number of paths;
% cut: stop after a multipath arrival has strength below that of  direct arrival divided by cut
%
% Output:
% l: path lengths
% tau: path delays, relative to direct (add l(1)/c to get absolute values)
% Gamma: cumulative reflection coefficient 
% theta: angles of path arrivals
% ns/nb: number of surface/bottom  reflections
% hp: path gains
%
% Calls: absorption.m, reflcoeff.m
% Note: frequency is used only to determine the number of multipaths as that for which path gains remain above a threshold;
% one could have used a fixed number of paths instead (but not sure in general what this number should be)



a=10^(absorption(f/1000)/10);a=a^(1/1000);
%reflloss=10^(reflloss/10); 

nr=0; % direct path, no reflections
theta(1)=atan((ht-hr)/d);
l(1)=sqrt((ht-hr)^2+d^2);
del(1)=l(1)/c;
A(1)=(l(1)^k).*(a.^l(1)); 
ns(1)=0; nb(1)=0;
G(1)=1/sqrt(A(1));
Gamma(1)=1;
hp(1)=1;

path=[0]; % begin with surface reflection; 

while min(abs(G))>=G(1)/cut;
%p=0; while p<=10; % could do for a certain number of paths, but better by attenuation
    nr=nr+1;
   
    p=2*nr; 
    first=path(1); last=path(end); 
    nb(p)=sum(path); ns(p)=nr-nb(p);
    heff=(1-first)*ht + first*(h-ht) + (nr-1)*h + (1-last)*hr + last*(h-hr);
    l(p)=sqrt(heff^2+d^2);
    theta(p)=atan(heff/d); if first==1; theta(p)=-theta(p); end; % corrected by Parastoo to address the grazing angle (rather than the angle of arrival)
    del(p)=l(p)/c;
    tau(p)=del(p)-del(1);
    A(p)=(l(p)^k)*(a^l(p)); 
    %Gamma(p)=1/reflloss^nb(p)*(-1)^ns(p); % fixed reflection coefficient
    Gamma(p)=reflcoeff(abs(theta(p)),c, c2)^nb(p)*(-1)^ns(p); % refl. coeff. calulated as f. of grazing angle
    G(p)=Gamma(p)/sqrt(A(p)); 
    hp(p)=Gamma(p)/sqrt(( (l(p)/l(1))^k )*( a^(l(p)-l(1)) ));
    
    p=2*nr+1; path=not(path);
    first=path(1); last=path(end); 
    nb(p)=sum(path); ns(p)=nr-nb(p);
    heff=(1-first)*ht + first*(h-ht) + (nr-1)*h + (1-last)*hr + last*(h-hr);
    l(p)=sqrt(heff^2+d^2);
    theta(p)=atan(heff/d); if first==1; theta(p)=-theta(p); end; % corrected by Parastoo to address the grazing angle (rather than the angle of arrival)
    del(p)=l(p)/c;
    tau(p)=del(p)-del(1);
    A(p)=(l(p)^k)*(a^l(p));
    %Gamma(p)=1/reflloss^nb(p)*(-1)^ns(p); % fixed reflectoin coefficient
    Gamma(p)=reflcoeff(abs(theta(p)),c, c2)^nb(p)*(-1)^ns(p); % refl. coeff. calulated as f. of grazing angle
    G(p)=Gamma(p)/sqrt(A(p)); 
    hp(p)=Gamma(p)/sqrt(( (l(p)/l(1))^k )*( a^(l(p)-l(1)) ));
    
    path=[path not(path(end))];
    
end;


    
 

