function alpha=Absorp_FK(E,wavelength, PD_Str,zero_or_one)
% Use the maximum(regular alpha, Franz Keldsly alpha)
E(abs(E)<1) = 1; % to avoid 1/0 division
q=1.602e-19;
hbar = 1.05457e-34;
omga = 2*pi*2.99792458e+8/wavelength;
omga1= PD_Str.Eg*q/hbar;
C = PD_Str.alpha.*omga./sqrt(abs(omga-omga1));
% C(omga<=omga1) = 0;
mu = 1./(1./PD_Str.mes+1./PD_Str.mhs)*9.10938e-31;
thetaF=(q^2*E.^2./(2*mu*hbar)).^(1/3);
beta=(omga1-omga)./thetaF;
alpha1 = C./omga.*thetaF.^(1/2).*((abs(airy(1,beta))).^2-beta.*(abs(airy(beta))).^2);
alpha = zero_or_one*alpha1+PD_Str.alpha;
