function refl=reflcoeff(theta, c1, c2)

%c1=1500;
%c2=1300;
rho1=1000; % in kg/m3
rho2=1800; % in kg/m3

x1=rho2/c1*sin(theta);
x2=rho1/c2*sqrt(1-(c2/c1)^2*cos(theta)^2);
thetac=acos(c1/c2);  % for theta below critical, total reflection
thetac=real(thetac); % in case c1>c2 
if theta<thetac 
    if thetac==0 
        refl=-1; 
    else 
        refl=exp(1j*pi*(1-theta/thetac)); 
    end
end
if theta>=thetac
    refl=(x1-x2)/(x1+x2);
end
