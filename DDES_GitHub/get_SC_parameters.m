function [me, mh, Eg, CnAu, CpAu, Br, NC, NV, Ae, Be, ...
    Ah, Bh, betai, Ep, DeltaEA, DeltaED, Ne_ref, Np_ref, ...
    Vnsat, Vpsat, Ebslon, un, up, tor1, tor2, n_i, Xi, alpha, ...
    eta_n, eta_p, alpha_abs, resistivity, delta_Fermi, k_L] = get_SC_parameters(material,np_type, Nd, Q, Ni, wavelength,NP,donor_acceptor)
% some parameters are taken from
% https://www.iue.tuwien.ac.at/phd/palankovski/node37.html
% for GaInAs
% http://www.ioffe.ru/SVA/NSM/Semicond/GaInAs/bandstr.html

% icc = intrinsic carrier concentration
% int_car_cont_InP = 1e7;
% int_car_cont_GaAs = 2.1e6;
% int_car_cont_InGaAs = 6.3e11;
% int_car_cont_GaP = 2;
% int_car_cont_InAs = 1e15;
% int_car_cont_Si = 1e10;
% int_car_cont_Ge = 2e13;
% int_car_cont_InGaAsP = 10.^(8.6335+4.3733*(y0-0.27));

if nargin < 2
    kB=1.38e-23;                % Boltzmann constant,unit,J/K
    q=1.6e-19;                  % unit charge (coulomb)
    wavelength = 1.55e-6;
    np_type = 1;
    me = 9.1e-31;
    T= 300;     % Assume room temperature
    Nd = 0;     % No dopping
    Q = 1.2;    % for InGaAsP
    Ni = 1.4e16;
    epsilon_normalization = 1.21e-12;
    Diff_Coef = 8000;
    %
    VT=kB*T/q;                  % normalized voltage or potential
    NX=sqrt(epsilon_normalization/q*VT/Ni);     % normalized position x
    D0=Diff_Coef*VT;                 % normalized diffusion coeificient
    NT=NX^2/D0;                 % normalized time t
    k_L = 80*NT^3*T/(me*NX);                %k_L = 0.68/(me*L*NT^2*T);
    Nmu=D0/VT;                  % normalized mobility;
else
    T = NP.T;
    kB=NP.KB;
    q=NP.q;
    VT=NP.VT;                  % normalized voltage or potential
    NX = NP.NX;
    NT = NP.NT;
    Nmu = NP.Nmu;
    me = 9.1e-31;
    k_L = 80*NT^3*T/(me*NX);                 %k_L = 0.68/(me*L*NT^2*T);
end
E_inc_energy = 1.24e-6/wavelength;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electron and Hall Mobilities: GaAs, InAs, InP, InGaAs, InGaAsP
un_GaAs0=8500;
up_GaAs0=800;
un_InAs0=32500;
up_InAs0=510;
un_InP0=5300;
up_InP0=200;

un_GaAs=un_GaAs0*(300/T).^2.2;   %GaAs
up_GaAs=up_GaAs0*(300/T).^0.9;   %GaAs
un_InAs=un_InAs0*(300/T).^1.7;   % InAs
up_InAs=up_InAs0*(300/T).^2.3;   % InAs
un_InP= un_InP0*(300/T).^1.9;    %InP
up_InP= up_InP0*(300/T).^1.2;    %InP

un_InGaAs = 1/(0.53/un_GaAs+0.47/un_InAs);
up_InGaAs = 1/(0.53/up_GaAs+0.47/up_InAs);

% Electron and Hall Saturation Velocities: GaAs, InAs, InP
Vnsat_GaAs = 7.2e6/(0.44+0.56*T/300);
Vpsat_GaAs = 9e6/(0.43+0.57*T/300);
Vnsat_InAs = 9e6/(0.59+0.41*T/300);
Vpsat_InAs = 5e6/(0.7+0.3*T/300);
Vnsat_InP = 6.6e6/(0.7+0.3*T/300);
Vpsat_InP = 8e6/(0.7+0.3*T/300);
Vnsat_GaP = 1.25e7/(0.7+0.3*T/300);
Vpsat_GaP = 1e7/(0.7+0.3*T/300);

Vnsat_InGaAs = Vnsat_GaAs*0.47+0.53*Vnsat_InAs-0.47*0.53*0.196e7;
Vpsat_InGaAs = Vpsat_GaAs*0.47+0.53*Vpsat_InAs-0.47*0.53*0.196e7;
Vnsat_InGaP = Vnsat_InP*0.47+0.53*Vnsat_GaP-0.47*0.53*0.196e7;
Vpsat_InGaP = Vpsat_InP*0.47+0.53*Vpsat_GaP-0.47*0.53*0.196e7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% taun_GaAs = 1e-9;
% taup_GaAs = 2e-8;

switch material
    case 'Si'
        % icc = 1e10;
        [n, k, alpha] = get_Si_n_k(wavelength);
        CnAu=8.3e-32; % Auger for e, Ref: Atlas User Manual
        CpAu=1.8e-31; % Auger for h,Ref: Atlas User Manual
        Xi = 4.17;  % Vacuum-Ec difference (eV) Ref: Atlas User Manual
        NC = 5.4e15*T^1.5;  % Atlas uses 2.8e19 at T = 300
        NV = 2.0e15*T^1.5;  % Atlas uses 1.04e19 at T = 300
        Eg = 1.1255 - 4.73e-4*T^2/(T+636); % silvaco

        me = 0.98;  % effective electron mass (m0)
        mh = 0.49;   % effective hole mass (m0)
        Br=1.1e-14;
        Ae= 7.03e5;           % ref: http://impact-ionisation.group.shef.ac.uk/ionisation_coeff/Si/
        Be= 1.23e6;
        Ah=1.58e6;
        Bh=2.04e6;
% somehow these parameters fail        
%         betai = 8e-4;
%         Ep=1e32;
%         DeltaEA=0.045;
%         DeltaED=0.054;
%         switch donor_acceptor  % https://www.ioffe.ru/SVA/NSM/Semicond/Si/bandstr.html#Donors
%             case 'As'
%                 DeltaED=0.054;
%             case 'P'
%                 DeltaED=0.045;
%             case 'Sb'
%                 DeltaED=0.043;
%             case 'Al'
%                 DeltaEA=0.072;
%             case 'B'
%                 DeltaEA=0.045;
%             case 'Ga'
%                 DeltaEA=0.074;
%             case 'In'
%                 DeltaEA=0.157;
%         end

        betai = 8e-8;
        Ep=1e32;
        DeltaED=0.085; %
        DeltaEA=0.007; %

        Ne_ref = 8.5e16;    % https://en.wikipedia.org/wiki/Electron_mobility
        Np_ref = 6.3e16;
        eta_n = 0.72;
        eta_p = 0.76;
        Vnsat = 2.4e7/(1+0.8*exp(T/600));
        Vpsat = 1.67e7/(1+0.8*exp(T/600));    % REF: 10.1007/bf00885852
        Ebslon=11.8*8.854e-14;
        %         un=1400/(300/T).^2.4;
        %         up=500/(300/T).^2.2;    % Ref: Atlas User Manual

        un = get_Si_electron_mobility(Nd, T);
        up = get_Si_hole_mobility(Nd, T);

        %
        if np_type == 1
            Eg = Eg - (10.3e-9*Nd^(1/3)+4.43e-7*Nd^0.25+3.38e-12*Nd^0.5);
        else
            Eg = Eg -  22.5e-9*Nd^(1/3);
        end

%         KSRHTN = 2.5e-3;
%         KSRHCN = 3.0e-13;
%         KSRHGN =  1.77;
%         KSRHTP = 2.5e-3;
%         KSRHCP = 11.76e-13;
%         KSRHGP = 0.57;
% 
%         TL = 300;
%         Ndd = logspace(3,20,100);
% 
%         taun = 1./((1./KSRHTN+KSRHCN*Ndd)*(300/TL)^KSRHGN);
%         taup = 1./((1./KSRHTP+KSRHCP*Ndd)*(300/TL)^KSRHGP);
%         %
    case 'Ge'
        %  icc = 2e13;
        [n, k, alpha] = get_Ge_n_k(wavelength);
        me = 1.59;  % effective electron mass (m0)
        mh = 0.33;   % effective hole mass (m0)
        Xi = 4;  % Vacuum-Ec difference (eV)
        Eg = 0.742- 4.8e-4*T^2/(T+235) ;
        CnAu=1e-30; % Auger for e
        CpAu=1e-30; % Auger for h
        Br=6.41e-14;
        NC = 1.98e15*T^1.5;
        NV = 9.6e14*T^1.5;
        Ae= 4.9e5;
        Be= 7.9e5;
        Ah= 2.15e5;
        Bh= 7.1e5;
        
        betai = 8e-8;
        Ep=4e3;
        DeltaED=0.085; %
        DeltaEA=0.007; %

%         betai = 8e-4;
%         Ep=4e3;
%         DeltaEA=0.012;
%         DeltaED=0.013;
%         switch donor_acceptor  % https://www.ioffe.ru/SVA/NSM/Semicond/Ge/bandstr.html#Donors
%             case 'As'
%                 DeltaED=0.014;
%             case 'P'
%                 DeltaED=0.013;
%             case 'Sb'
%                 DeltaED=0.010;
%             case 'Bi'
%                 DeltaED=0.013;
%             case 'Li'
%                 DeltaED=0.093;
%             case 'Al'
%                 DeltaEA=0.011;
%             case 'B'
%                 DeltaEA=0.011;
%             case 'Ga'
%                 DeltaEA=0.011;
%             case 'In'
%                 DeltaEA=0.012;
%             case 'Tl'
%                 DeltaEA=0.013;
%         end
        Ne_ref = 2.8e16;   % https://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=1594&context=ecetr
        Np_ref = 1.37e17;
        eta_n = 0.42;
        eta_p = 0.467;
        Vnsat = 1.4e7/(1+0.8*exp(T/600));
        Vpsat = 1.67e7/(1+0.8*exp(T/600));    % https://sci-hub.do/10.1139/p79-172
        Ebslon= 16.2*8.854e-14;
        un=3900/(300/T).^1.7/(1+Nd/1e17)^0.5;   % ref: https://ecee.colorado.edu/~bart/book/book/chapter2/ch2_7.htm
        up=1900/(300/T).^2.3;
        %
        if np_type == 1
            Eg = Eg - 1e-3*(8.21*(Nd/1e18)^(1/3)+9.18*(Nd/1e18)^0.25+5.77*(Nd/1e18)^0.5);
        else
            Eg = Eg - 1e-3*(8.67*(Nd/1e18)^(1/3)+8.14*(Nd/1e18)^0.25+4.31*(Nd/1e18)^0.5);
        end
%         
%         taun = 1e-4;
%         taup = 1e-3;

    case 'InP'
        % icc = 1e7;
        alpha = get_SC_alpha('InP', wavelength*1e6, 0);
        me = 0.08;  % effective electron mass (m0)
        % mh = 0.64;   % effective hole mass (m0)
        mlh = 0.089;
        mhh = 0.6;
        mh = (mlh.^1.5+mhh.^1.5)^(2/3);
        k_L = 68; 
        Xi = 4.38;  % Vacuum-Ec difference (eV)
        Eg = 1.421 - 4.9e-4*T*T/(T+327);
        CnAu=1e-31; % Auger for e
        CpAu=1e-31; % Auger for h
        Br=1e-11;
        NC = 1.1e14*T^1.5;
        NV = 2.2e15*T^1.5;
        Ae=1.12e7;
        Be=6.22e6;
        Ah=4.79e6;
        Bh=2.55e6;
        betai = 8e-4;
        Ep=10e3;
        DeltaEA=0.035;
        DeltaED=0.0057;
        % https://www.ioffe.ru/SVA/NSM/Semicond/InP/bandstr.html#Donors
        switch donor_acceptor  
            case 'C'
                DeltaEA=0.04;
            case 'Hg'
                DeltaEA=0.098;
            case 'Zn'
                DeltaEA=0.035;
            case 'Cd'
                DeltaEA=0.057;
            case 'Si'
                DeltaEA=0.03;
            case 'Cu'
                DeltaEA=0.06;
            case 'Be'
                DeltaEA=0.03;
            case 'Mg'
                DeltaEA=0.03;
            case 'Ge'
                DeltaEA=0.021;
            case 'Mn'
                DeltaEA=0.27;
        end
        Ne_ref = 1e17;
        Np_ref = 6e17;
        eta_n = 0.34;
        eta_p = 0.64;
        Vnsat = Vnsat_InP;
        Vpsat = Vpsat_InP;
        Ebslon=12.5*8.854e-14;
        un=un_InP;
        up=up_InP;

        
        taup = 3e-6;
        taun = 2e-9;
        %
        if np_type == 1
            Eg = Eg - (10.3e-9*Nd^(1/3)+4.43e-7*Nd^0.25+3.38e-12*Nd^0.5);
        else
            Eg = Eg -  22.5e-9*Nd^(1/3);
        end
        %
        % mlp = 0.089;    % effective hole mass (m0)
        % EgammaX = 0.96 - 3.7e-4*T;
        % Ni0 = 1.3e7; % intrinsic carrier concentration (cm-3)
    case 'GaN'
        %<><><><><><><><><><><><><><><><><><><><><><><><><>
        % http://www.ioffe.ru/SVA/NSM/Semicond/GaN/basic.html
        % https://aip.scitation.org/doi/am-pdf/10.1063/1.4948794?class=chorus+notVisible
        [n, k, alpha] = get_GaN_n_k(wavelength);       
        % Wurtzite crystal structure
        me=0.2;   % effective mass electrons
        mh=0.8;    % effective mass for holes
%         Ebslon=8.9*8.854e-14;   % high frequency
       NC = 4.3e14*T^1.5;
       NV = 8.9e15*T^1.5;
%         % Zinc Blende crystal structure
%         me=0.13;   % effective mass electrons
%         mh=1.4;    % effective mass for holes
%         NC = 2.3e14*T^1.5;
%         NV = 8e15*T^1.5;
        Xi=4.1;     % eV affinity
        Ebslon=n^2*8.854e-14;   % high frequency
        Eg = 3.427-9.39e-4*T*T/(T+772);
        un = 400;
        up = 200;
        betai = 1e-10;
        Vnsat = 2.5e7;
        Vpsat = 2e7;        
        Ae=1.5e5;
        Be=1.41e7;
        Ah=6.4e5;
        Bh=1.46e7;
        Ne_ref = 1.0e17;
        Np_ref = 1.0e17;
        CnAu=3e-28;
        CpAu=3e-28;
        Br = 1.1e-10;
        DeltaEA=0.14; % could be anything in between 0.12-0.02 eV
        DeltaED=0.12; % could be anything in between 0.14-0.21 eV
        % I don't know what these parameters should be
        Ep=4e3;
        eta_n = 0.1;
        eta_p = 0.1;
        
        taun = 1e-9;    % in region-i
        taup = 2e-8;   % o/wise

        %<><><><><><><><><><><><><><><><><><><><><><><><><>
    case 'InGaAs'
        % icc = 6.3e11;
        % http://www.ioffe.ru/SVA/NSM/Semicond/GaInAs/bandstr.html
        x=0.47;
        alpha = get_SC_alpha('InGaAs', wavelength*1e6, 0);

        me=0.041;   % effective mass electrons
        % mh=0.59;    % effective mass for holes
        mlh = 0.052;
        mhh = 0.45;
        mh = (mlh.^1.5+mhh.^1.5)^(2/3);
        k_L = 4.5; 
        CnAu=1.8e-31;
        CpAu=1.2e-31;
        Br = 1e-10;
        NC = 4.82e15*(0.023+0.037*x+0.003*x*x)^1.5*T^1.5;
        NV = 4.82e15*(0.41-0.1*x)^1.5*T^1.5;
        Xi=4.9 - 0.83*x;     % eV affinity
        Ae=6.64e7;
        Be=4e6;
        Ah=9.34e7;
        Bh=2.26e6;
        betai = 6e-4;
        Ep=4e3;

        DeltaEA=0.025; %
        DeltaED=0.03; %
        switch donor_acceptor
            case 'Te'
                DeltaED=0.03; %
            case 'C'
                DeltaEA=0.02; %
            case 'Si'
                DeltaEA=0.03; %
            case 'Ge'
                DeltaEA=0.03; %
            case 'Zn'
                DeltaEA=0.025; %
            case 'Sn'
                DeltaEA=0.2; %
        end


        Ne_ref = 1e17;
        Np_ref = 1e18;
        Vnsat = Vnsat_InGaAs;
        Vpsat = Vpsat_InGaAs;
        Ebslon=14.2*8.854e-14;
        un = un_InGaAs;
        up = up_InGaAs;

        eta_n = 0.5;
        eta_p = 0.45;
        %
        Eg = 0.42+0.625*x-(5.8/(T+300)-4.19/(T+271))*1e-4*T*T*x-4.19e-4*T*T/(T+271)+0.475*x*x;
        if np_type == 1  %  p type
            Eg = Eg -1e-3*(9.2e-9*Nd^(1/3)+3.57e-7*Nd^0.25+3.65e-12*Nd^0.5);
        else             % n type
            Eg = Eg -1e-3*(15.5e-9*Nd^(1/3)+1.95e-7*Nd^0.25+159e-12*Nd^0.5);
        end
        taun = 1e-9;    % in region-i
        taup = 2e-8;   % o/wise

    case 'InGaAsP'
        % http://www.ioffe.ru/SVA/NSM/Semicond/GaInAsP/bandstr.html
        % http://www.ioffe.ru/SVA/NSM/Semicond/GaInAsP/bandstr.html
        [x, y, me, mh, un, up, Eg, permittivity, Auger] = get_InGaAsP_y(Q);
        % icc = 10.^(8.6335+4.3733*(y-0.27));
        alpha = get_SC_alpha('InGaAsP', wavelength*1e6, y);
        CnAu=Auger;
        CpAu=Auger;
        Br = 0.5e-10;
        NC = 7e13*T^1.5;
        NV = 5.04e15*(0.6-0.18*y)^1.5*T^1.5;
        % GaxIn(1-x)AsyP(1-y)
        % Example: x = 0.47, y =0.3 ==> Ga0.47In0.53As0.3P0.7
        % so it has to be 0.3GaAs 0.17 GaP 0.53InP
        Xi = 4.32;  % silvaco
        Ae=1e4;
        Be=1e6;
        Ah=1e4;
        Bh=4e6;
        betai = 1e-4;
        Ep = 6e3;
        k_L = 3;
        DeltaED=0.044; %
        DeltaEA=0.035; %
        switch donor_acceptor % https://www.ioffe.ru/SVA/NSM/Semicond/GaInAsP/bandstr.html
            case 'Mg'
                DeltaEA=0.035; %
            case 'Zn'
                DeltaEA=(37.5-(y-0.3)*(37.5-22)/0.6)/1000; %
            case 'Cd'
                DeltaEA=(60-(y-0.2)*(30)/0.7)/1000;
            case 'Be'
                DeltaEA=0.04; %
        end


        Ne_ref = (0.08-0.039*y)^1.5*2.5e19;
        Np_ref = (0.6-0.18*y)^1.5*2.5e19;
        Vnsat= (1-y-x)*Vnsat_GaP+Vnsat_GaAs*y+Vnsat_InP*x;
        % Vnsat= 1.5e7;
        Vpsat=(1-y-x)*Vpsat_GaP+Vpsat_GaAs*y+Vpsat_InP*x;
        Ebslon=permittivity*8.854e-14;
        eta_n = 0.42;
        eta_p = 0.54;
        % Temperature effect
        %Eg = Eg-4.3e-4*T*T/(T+224)-923.6e-12/permittivity.^1.5*sqrt(Nd)*sqrt(300/T);
        % Not Used
        % Ni0 = 2.5e19*(0.08-0.039*y)^1.5;
        % Eso = 0.11+0.24*y;
        taun = 1e-9;    % in region-i
        taup = 2e-8;   % o/wise
        % keyboard
end

%
if np_type == 1
    resistivity = 1/q/(NC*un+(NV+Nd)*up)*1e4; % 1e4 to convert cm to um
else
    resistivity = 1/q/((NC+Nd)*un+NV*up)*1e4; % 1e4 to convert cm to um
end
tor1 =  2e-8;   % low doping
tor2 =  1e-9;    % high doping
 %keyboard
% Normalizations
Ep=Ep*NX/VT;
NV = NV/Ni;
NC = NC/Ni;
n_i = sqrt(NC*NV)*exp(-Eg/(2*kB*T/q));
tor1 = tor1/NT;
tor2 = tor2/NT;
CnAu = CnAu*NT/(NX)^6;
CpAu = CpAu*NT/(NX)^6;
Br = Br*NT/(NX)^3;
k_L = k_L*NT^3*T/(9.1e-31*NX);
un = un/Nmu;
up = up/Nmu;
Vnsat = Vnsat*NX/VT/Nmu;
Vpsat = Vpsat*NX/VT/Nmu;
betai= betai/(NX/VT);
% kL=k_L*NX*NT^3;
alpha_abs = alpha; % it might not generate e-h pairs but still the field decays
if Eg>E_inc_energy
    alpha = 0;
end
delta_Fermi = -kB*T*log(Nd/n_i)/q;
if np_type ~=1
    delta_Fermi =-delta_Fermi;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha = get_SC_alpha(material, micro_m, y)
%  input: micro_m (wavelength in micro meter, e.g. micro_m = 1.55)
% y for InGaAs P
%             y 1-y
% Reference:
% Optical dispersion relations for GaP, GaAs, GaSb, InP, InAs, InSb,
% Al[x]Ga[1?x]As, and In[1?x]Ga[x]As[y]P[1?y]
% Sadao Adachi
% J. Appl. Phys. 66, 6030 (1989)
% doi: 10.1063/1.343580

if nargin<2
    y = 1;
end
%
EeV = 4.13566733e-1*2.99792458/micro_m;
% parameters for GaP, GaAs, GaSb, InP, InAs, InSb
params = [  2.74    1.42    0.72    1.35    0.36    0.18;
    2.84    1.77    1.46    1.45    0.76    0.99;
    3.70    2.90    2.05    3.10    2.50    1.80;
    0.00    3.13    2.50    3.25    2.78    2.30;
    5.00    4.7     4.0     4.7     4.45    3.9;
    2.26    1.73    0.76    2.05    1.07    0.93;
    13.76   3.45    0.71    6.57    0.61    0.19;
    6.35    6.37    6.68    4.93    6.59    6.37;
    9.49    13.08   14.29   10.43   13.76   12.26;
    0.06    0.10    0.09    0.10    0.2     0.16;
    2.08    2.39    5.69    1.49    1.78    5.37;
    0.132   0.146   0.290   0.094   0.108   0.318;
    4.6     24.2    7.4     60.4    20.8    19.5;
    0.0     1.6     1.0     1.6     2.8     3.1];
% parameters for AlGaAs, InGaAs[y]P[1-y]
params2 = [1.83 2.42 1.18 0.75
    2.15 2.73 1.34 1.04
    3.13 3.43 2.96 2.57
    0 0 3.14 2.83
    4.7 4.7 4.65 4.41
    1.92 2.03 1.83 1.20
    8.80 23.20 4.39 1.20
    6.05 5.41 4.30 3.84
    0 0 0.53 1.48
    11.05 9.55 8.76 7.57
    0 0 1.06 2.96
    0.11 0.12 0.12 0.14
    2.30 1.76 1.98 2.90
    0.135 0.103 0.145 0.225
    16.1 8.1 39.0 20.7
    0.6 -0.3 2.1 2.8];

switch material
    case 'InGaAsP'
        E0 =        params2(1,3)+(params2(1,4)-params2(1,3))/0.76*(y-0.24);
        delta0 =    params2(2,3)+(params2(2,4)-params2(2,3))/0.76*(y-0.24);
        E1 =        params2(3,3)+(params2(3,4)-params2(3,3))/0.76*(y-0.24);
        delta1 =    params2(4,3)+(params2(4,4)-params2(4,3))/0.76*(y-0.24);
        E2 =        params2(5,3)+(params2(5,4)-params2(5,3))/0.76*(y-0.24);
        Eg =        params2(6,3)+(params2(6,4)-params2(6,3))/0.76*(y-0.24);
        A =         params2(7,3)+(params2(7,4)-params2(7,3))/0.76*(y-0.24);
        B1 =        params2(8,3)+(params2(8,4)-params2(8,3))/0.76*(y-0.24);
        B2 =        params2(9,3)+(params2(9,4)-params2(9,3))/0.76*(y-0.24);
        B11 =       params2(10,3)+(params2(10,4)-params2(10,3))/0.76*(y-0.24);
        B21 =       params2(11,3)+(params2(11,4)-params2(11,3))/0.76*(y-0.24);
        Gamma =     params2(12,3)+(params2(12,4)-params2(12,3))/0.76*(y-0.24);
        C =         params2(13,3)+(params2(13,4)-params2(13,3))/0.76*(y-0.24);
        Y =         params2(14,3)+(params2(14,4)-params2(14,3))/0.76*(y-0.24);
        D =         params2(15,3)+(params2(15,4)-params2(15,3))/0.76*(y-0.24);
        eps_inf =   params2(16,3)+(params2(16,4)-params2(16,3))/0.76*(y-0.24);
    case 'InGaAs'
        eps_inf = 2.8;
        E0   = 0.75;    %EeV
        delta0   = 1.04-E0; %EeV
        E1   = 2.57;    %EeV
        delta1   = 2.83-E1; %EeV
        E2   = 4.41;    %EeV
        Eg   = 1.20;    %EeV
        A    = 1.20;    %EeV.^1.5
        B1   = 3.84;
        B2   = 1.48;
        B11  = 7.57;    %EeV.^-0.5
        B21  = 2.96;    %EeV.^-0.5
        Gamma = 0.14;    %EeV
        C    = 2.90;
        Y    = 0.225;
        D    = 20.7;
    case 'GaP'
        col = 1;
        E0 = params(1,col);
        delta0 = params(2,col);
        E1 = params(3,col);
        delta1 = params(4,col);
        E2 = params(5,col);
        Eg = params(6,col);
        A = params(7,col);
        B1 = params(8,col);
        B11 = params(9,col);
        Gamma = params(10,col);
        C = params(11,col);
        Y = params(12,col);
        D = params(13,col);
        eps_inf = params(14,col);
        B2 = 0;
        B21 = 0;
    case 'GaAs'
        col = 2;
        E0 = params(1,col);
        delta0 = params(2,col);
        E1 = params(3,col);
        delta1 = params(4,col);
        E2 = params(5,col);
        Eg = params(6,col);
        A = params(7,col);
        B1 = params(8,col);
        B11 = params(9,col);
        Gamma = params(10,col);
        C = params(11,col);
        Y = params(12,col);
        D = params(13,col);
        eps_inf = params(14,col);
        B2 = 0;
        B21 = 0;
    case 'GaSb'
        col = 3;
        E0 = params(1,col);
        delta0 = params(2,col);
        E1 = params(3,col);
        delta1 = params(4,col);
        E2 = params(5,col);
        Eg = params(6,col);
        A = params(7,col);
        B1 = params(8,col);
        B11 = params(9,col);
        Gamma = params(10,col);
        C = params(11,col);
        Y = params(12,col);
        D = params(13,col);
        eps_inf = params(14,col);
        B2 = 0;
        B21 = 0;
    case 'InP'
        col = 4;
        E0 = params(1,col);
        delta0 = params(2,col);
        E1 = params(3,col);
        delta1 = params(4,col);
        E2 = params(5,col);
        Eg = params(6,col);
        A = params(7,col);
        B1 = params(8,col);
        B11 = params(9,col);
        Gamma = params(10,col);
        C = params(11,col);
        Y = params(12,col);
        D = params(13,col);
        eps_inf = params(14,col);
        B2 = 0;
        B21 = 0;
    case 'InAs'
        col = 5;
        E0 = params(1,col);
        delta0 = params(2,col);
        E1 = params(3,col);
        delta1 = params(4,col);
        E2 = params(5,col);
        Eg = params(6,col);
        A = params(7,col);
        B1 = params(8,col);
        B11 = params(9,col);
        Gamma = params(10,col);
        C = params(11,col);
        Y = params(12,col);
        D = params(13,col);
        eps_inf = params(14,col);
        B2 = 0;
        B21 = 0;
    case 'InSb'
        col = 6;
        E0 = params(1,col);
        delta0 = params(2,col);
        E1 = params(3,col);
        delta1 = params(4,col);
        E2 = params(5,col);
        Eg = params(6,col);
        A = params(7,col);
        B1 = params(8,col);
        B11 = params(9,col);
        Gamma = params(10,col);
        C = params(11,col);
        Y = params(12,col);
        D = params(13,col);
        eps_inf = params(14,col);
        B2 = 0;
        B21 = 0;
end
eps_A  = Epsilon_A(EeV, eps_inf, E0, delta0, E1, delta1, E2, Eg, A, B1, B2, B11, B21, Gamma, C, Y, D);
eps_B  = Epsilon_B(EeV, eps_inf, E0, delta0, E1, delta1, E2, Eg, A, B1, B2, B11, B21, Gamma, C, Y, D);
eps_C  = Epsilon_C(EeV, eps_inf, E0, delta0, E1, delta1, E2, Eg, A, B1, B2, B11, B21, Gamma, C, Y, D);
eps_D  = Epsilon_D(EeV, eps_inf, E0, delta0, E1, delta1, E2, Eg, A, B1, B2, B11, B21, Gamma, C, Y, D);
eps_total = eps_A + eps_B + eps_C + eps_D + eps_inf;
n = real(sqrt(eps_total));
k = imag(sqrt(eps_total));
alpha = 4*pi*k/micro_m*1e4;  %(1/cm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Heaviside = H(x)
Heaviside = 0.5*(sign(x) + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function epsilon_complex = Epsilon_A(hbar_omega, eps_inf, E0, delta0, E1, delta1, E2, Eg, A, B1, B2, B11, B21, Gamma, C, Y, D)
x0 = hbar_omega/E0;
xso = hbar_omega / (E0+delta0);
H0 = H(1-x0);
Hso = H(1-xso);
fx0 = x0.^-2 * ( 2 -(1+x0).^0.5 - ((1-x0)*H0).^0.5 );
fxso = xso.^-2 * ( 2 - (1+xso).^0.5 - ((1-xso)*Hso).^0.5 );
H0 = H(x0-1);
Hso = H(xso-1);
eps_2 = A/(hbar_omega).^2 * ( ((hbar_omega-E0)*H0).^0.5 + 0.5*((hbar_omega-E0-delta0)*Hso).^0.5);
eps_1 = A*E0.^-1.5 * (fx0+0.5*(E0/(E0+delta0)).^1.5*fxso);
epsilon_complex = eps_1 + 1i*eps_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eps_complex = Epsilon_B(hbar_omega, eps_inf, E0, delta0, E1, delta1, E2, Eg, A, B1, B2, B11, B21, Gamma, C, Y, D) %E1
x1 = hbar_omega/E1;
x1s = hbar_omega/(E1+delta1);
H1 = H(1-x1);
H1s = H(1-x1s);

eps_2 = (pi*x1.^-2*(B1-B11*((E1-hbar_omega)*H1).^0.5)...
    + pi*x1s.^-2*(B2-B21*((E1+delta1-hbar_omega)*H1s).^0.5) );
eps_2 = eps_2*H(eps_2); % undocumented trick: ignore negative eps_2

x1 = (hbar_omega+1j*Gamma)/E1;
x1s = (hbar_omega+1j*Gamma)/(E1+delta1);
eps_1 = -B1*x1.^-2*log(1-x1.^2) - B2*x1s.^-2*log(1-x1s.^2);
eps_complex = real(eps_1) + 1i*real(eps_2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eps_complex = Epsilon_C(hbar_omega, eps_inf, E0, delta0, E1, delta1, E2, Eg, A, B1, B2, B11, B21, Gamma, C, Y, D)     %E2
x2 = hbar_omega/E2;
eps_2 = C*x2*Y / ((1-x2.^2).^2+(x2*Y).^2);
eps_1 = C*(1-x2.^2) / ((1-x2.^2).^2+(x2*Y).^2);
eps_complex = eps_1 + 1i*eps_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eps_complex = Epsilon_D(hbar_omega, eps_inf, E0, delta0, E1, delta1, E2, Eg, A, B1, B2, B11, B21, Gamma, C, Y, D)  %Eg
Ech = E1;
xg = Eg/hbar_omega;
xch = hbar_omega/Ech;
Hg = H(1-xg);
Hch = H(1-xch);
eps_2 = D/hbar_omega.^2 * (hbar_omega-Eg).^2*Hg*Hch;
eps_complex = 1i*eps_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% calculates Si hole mobility as a function of doping concentration (N)
% and temperature (T)
% ref: Atlas User Manual
% Ergun (UMBC, 11/22/22)
function mup0 = get_Si_hole_mobility(N, T)
MU1P_CAUG = 49.7;  % cm2/(V ⋅ s)
MU2P_CAUG = 479.37; % cm2/(V ⋅ s)
ALPHAP_CAUG = 0.0;  % arbitrary
BETAP_CAUG = -2.2; % arbitrary
GAMMAP_CAUG = -3.7; % arbitrary
DELTAP_CAUG = 0.70; % arbitrary
NCRITP_CAUG = 1.606e17; % cm-3

term1 = MU1P_CAUG.*(T/300).^ALPHAP_CAUG;
term2 = (MU2P_CAUG.*(T/300).^BETAP_CAUG - MU1P_CAUG.*(T/300).^ALPHAP_CAUG);
term3 = 1+(T/300).^GAMMAP_CAUG.*(N/NCRITP_CAUG).^DELTAP_CAUG;

mup0 = term1+term2./term3;


% calculates Si electron mobility as a function of doping concentration (N)
% and temperature (T)
% ref: Atlas User Manual
% Ergun (UMBC, 11/22/22)
function mun0 = get_Si_electron_mobility(N, T)
MU1N_CAUG = 55.24;  % cm2/(V ⋅ s)
MU2N_CAUG = 1429.23; % cm2/(V ⋅ s)
ALPHAN_CAUG = 0.0;  % arbitrary
BETAN_CAUG = -2.3; % arbitrary
GAMMAN_CAUG = -3.8; % arbitrary
DELTAN_CAUG = 0.73; % arbitrary
NCRITN_CAUG = 1.072e17; % cm-3

term1 = MU1N_CAUG.*(T/300).^ALPHAN_CAUG;
term2 = (MU2N_CAUG.*(T/300).^BETAN_CAUG - MU1N_CAUG.*(T/300).^ALPHAN_CAUG);
term3 = 1+(T/300).^GAMMAN_CAUG.*(N/NCRITN_CAUG).^DELTAN_CAUG;

mun0 = term1+term2./term3;


function [x, y0, me, mh, un, up, Eg, permittivity, Auger] = get_InGaAsP_y(Q)
% This function first finds the correct "y" value for a given Q value
% Then calculates mu_n, mu_p, Eg, me
% INPUT: Q (in microns), e.g. 1.15
% Reference: Appl. Phys. Lett. 33, 659 (1978)
Eg = 1.24/Q;
options = optimoptions('fsolve','Display','off','TolFun',1e-2,'MaxIter',200,'Tolx',1e-4...
    ,'Diagnostics','off','Jacobian','off','DerivativeCheck','off');
y0 = fsolve(@(y) 1.344-0.738*y+0.138*y*y-Eg, 0.5,options);
x = 0.137+(y0-0.3)/0.7*0.33;
% http://www.ioffe.ru/SVA/NSM/Semicond/GaInAsP/electric.html#Basic
% x = 0.47;
% Eg_check = 1.35 + x*(0.642 + 0.758*x)+ (0.101*y0-1.101 )*y0-(0.28*x-0.109*y0 + 0.159 )*x*y0;

% permittivity = 12.5 +1.44*y0;
% old formula: %me = 0.08-0.05*y0+0.017*y0*y0;
% mhh = 0.6-0.18*y0;  %  heavy holes
% mlh = 0.12 -0.099*y0 +0.03*y0^2; %  light holes
% un = 5400-7750*y0+14400*y0^2;
% up = 200-400*y0+500*y0^2;

% silvaco
permittivity =14.6*(1- x)*y0+ 12.5*(1-x)*(1 - y0)+ 13.18*x*y0+11.11*x*(1 - y0);
me= 0.08-0.116*y0+0.026*x-0.059*x*y0+(0.064-0.02*y0)*x^2+(0.06+0.032*x)*y0^2;
mhh = 0.46; % silvaco
mlh = 0.12-0.116*y0+0.03*x^2;   % silvaco
mh = (mlh.^1.5+mhh.^1.5)^(2/3);
mun1 = 33000 + (8500 - 33000) * x;
mup1 = 460+(400-460)*x;
mun2 = 4600 + (300 - 4600) * x;
mup2 = 150+(100-150)*x;
un = mun1 + (1-y0)*(mun2 - mun1);
up = mup1 + (1-y0)*(mup2 - mup1);
Auger = 4.7e-31*y0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n, k, alpha] = get_Si_n_k(wl)
wn = [2.5000e-01	1.6650e+00
    2.6000e-01	1.7570e+00
    2.7000e-01	2.0680e+00
    2.8000e-01	2.9590e+00
    2.9000e-01	4.3560e+00
    3.0000e-01	4.9760e+00
    3.1000e-01	5.1210e+00
    3.2000e-01	5.1120e+00
    3.3000e-01	5.1950e+00
    3.4000e-01	5.3010e+00
    3.5000e-01	5.4940e+00
    3.6000e-01	6.0260e+00
    3.7000e-01	6.8910e+00
    3.8000e-01	6.6160e+00
    3.9000e-01	6.0390e+00
    4.0000e-01	5.6130e+00
    4.1000e-01	5.3300e+00
    4.2000e-01	5.1190e+00
    4.3000e-01	4.9490e+00
    4.4000e-01	4.8120e+00
    4.5000e-01	4.6910e+00
    4.6000e-01	4.5870e+00
    4.7000e-01	4.4970e+00
    4.8000e-01	4.4190e+00
    4.9000e-01	4.3500e+00
    5.0000e-01	4.2940e+00
    5.1000e-01	4.2410e+00
    5.2000e-01	4.1930e+00
    5.3000e-01	4.1510e+00
    5.4000e-01	4.1120e+00
    5.5000e-01	4.0770e+00
    5.6000e-01	4.0450e+00
    5.7000e-01	4.0150e+00
    5.8000e-01	3.9880e+00
    5.9000e-01	3.9630e+00
    6.0000e-01	3.9400e+00
    6.1000e-01	3.9180e+00
    6.2000e-01	3.8980e+00
    6.3000e-01	3.8790e+00
    6.4000e-01	3.8610e+00
    6.5000e-01	3.8440e+00
    6.6000e-01	3.8280e+00
    6.7000e-01	3.8130e+00
    6.8000e-01	3.7980e+00
    6.9000e-01	3.7840e+00
    7.0000e-01	3.7720e+00
    7.1000e-01	3.7590e+00
    7.2000e-01	3.7480e+00
    7.3000e-01	3.7370e+00
    7.4000e-01	3.7270e+00
    7.5000e-01	3.7170e+00
    7.6000e-01	3.7080e+00
    7.7000e-01	3.6990e+00
    7.8000e-01	3.6910e+00
    7.9000e-01	3.6830e+00
    8.0000e-01	3.6750e+00
    8.1000e-01	3.6680e+00
    8.2000e-01	3.6610e+00
    8.3000e-01	3.6540e+00
    8.4000e-01	3.6470e+00
    8.5000e-01	3.6410e+00
    8.6000e-01	3.6350e+00
    8.7000e-01	3.6300e+00
    8.8000e-01	3.6240e+00
    8.9000e-01	3.6190e+00
    9.0000e-01	3.6140e+00
    9.1000e-01	3.6090e+00
    9.2000e-01	3.6040e+00
    9.3000e-01	3.6000e+00
    9.4000e-01	3.5950e+00
    9.5000e-01	3.5910e+00
    9.6000e-01	3.5870e+00
    9.7000e-01	3.5830e+00
    9.8000e-01	3.5790e+00
    9.9000e-01	3.5750e+00
    1.0000e+00	3.5720e+00
    1.0100e+00	3.5680e+00
    1.0200e+00	3.5650e+00
    1.0300e+00	3.5620e+00
    1.0400e+00	3.5590e+00
    1.0500e+00	3.5560e+00
    1.0600e+00	3.5530e+00
    1.0700e+00	3.5500e+00
    1.0800e+00	3.5470e+00
    1.0900e+00	3.5450e+00
    1.1000e+00	3.5420e+00
    1.1100e+00	3.5400e+00
    1.1200e+00	3.5370e+00
    1.1300e+00	3.5350e+00
    1.1400e+00	3.5320e+00
    1.1500e+00	3.5300e+00
    1.1600e+00	3.5280e+00
    1.1700e+00	3.5260e+00
    1.1800e+00	3.5240e+00
    1.1900e+00	3.5220e+00
    1.2000e+00	3.5200e+00
    1.2100e+00	3.5180e+00
    1.2200e+00	3.5170e+00
    1.2300e+00	3.5150e+00
    1.2400e+00	3.5130e+00
    1.2500e+00	3.5110e+00
    1.2600e+00	3.5090e+00
    1.2700e+00	3.5080e+00
    1.2800e+00	3.5060e+00
    1.2900e+00	3.5050e+00
    1.3000e+00	3.5030e+00
    1.3100e+00	3.5020e+00
    1.3200e+00	3.5000e+00
    1.3300e+00	3.4990e+00
    1.3400e+00	3.4970e+00
    1.3500e+00	3.4960e+00
    1.3600e+00	3.4950e+00
    1.3700e+00	3.4940e+00
    1.3800e+00	3.4920e+00
    1.3900e+00	3.4910e+00
    1.4000e+00	3.4900e+00
    1.4100e+00	3.4890e+00
    1.4200e+00	3.4880e+00
    1.4300e+00	3.4870e+00
    1.4400e+00	3.4860e+00
    1.4500e+00	3.4850e+00
    2.0000e+00	3.4850e+00];
wk = [2.5000e-01	3.6650e+00
    2.6000e-01	4.0840e+00
    2.7000e-01	4.6800e+00
    2.8000e-01	5.2870e+00
    2.9000e-01	5.2860e+00
    3.0000e-01	4.2340e+00
    3.1000e-01	3.5980e+00
    3.2000e-01	3.3030e+00
    3.3000e-01	3.1000e+00
    3.4000e-01	2.9770e+00
    3.5000e-01	2.9380e+00
    3.6000e-01	2.9660e+00
    3.7000e-01	2.1710e+00
    3.8000e-01	9.4600e-01
    3.9000e-01	4.4500e-01
    4.0000e-01	2.9600e-01
    4.1000e-01	2.2700e-01
    4.2000e-01	1.7600e-01
    4.3000e-01	1.3800e-01
    4.4000e-01	1.0700e-01
    4.5000e-01	8.6302e-02
    4.6000e-01	7.1381e-02
    4.7000e-01	6.2086e-02
    4.8000e-01	5.5004e-02
    4.9000e-01	4.9131e-02
    5.0000e-01	4.4165e-02
    5.1000e-01	3.9367e-02
    5.2000e-01	3.6415e-02
    5.3000e-01	3.3108e-02
    5.4000e-01	3.0295e-02
    5.5000e-01	2.7968e-02
    5.6000e-01	2.5758e-02
    5.7000e-01	2.4131e-02
    5.8000e-01	2.2524e-02
    5.9000e-01	2.1081e-02
    6.0000e-01	1.9934e-02
    6.1000e-01	1.8446e-02
    6.2000e-01	1.7367e-02
    6.3000e-01	1.6444e-02
    6.4000e-01	1.5432e-02
    6.5000e-01	1.4431e-02
    6.6000e-01	1.3498e-02
    6.7000e-01	1.2743e-02
    6.8000e-01	1.1905e-02
    6.9000e-01	1.1201e-02
    7.0000e-01	1.0528e-02
    7.1000e-01	1.0057e-02
    7.2000e-01	9.6257e-03
    7.3000e-01	8.9461e-03
    7.4000e-01	8.3620e-03
    7.5000e-01	7.8185e-03
    7.6000e-01	7.1970e-03
    7.7000e-01	6.7402e-03
    7.8000e-01	6.3933e-03
    7.9000e-01	5.8340e-03
    8.0000e-01	5.4113e-03
    8.1000e-01	4.9955e-03
    8.2000e-01	4.6134e-03
    8.3000e-01	4.2734e-03
    8.4000e-01	3.9439e-03
    8.5000e-01	3.6120e-03
    8.6000e-01	3.2781e-03
    8.7000e-01	2.9839e-03
    8.8000e-01	2.6821e-03
    8.9000e-01	2.4293e-03
    9.0000e-01	2.1701e-03
    9.1000e-01	1.9625e-03
    9.2000e-01	1.7571e-03
    9.3000e-01	1.5467e-03
    9.4000e-01	1.3689e-03
    9.5000e-01	1.1793e-03
    9.6000e-01	1.0237e-03
    9.7000e-01	8.7225e-04
    9.8000e-01	7.4866e-04
    9.9000e-01	6.2238e-04
    1.0000e+00	5.0930e-04
    1.0100e+00	4.1071e-04
    1.0200e+00	3.2386e-04
    1.0300e+00	2.4753e-04
    1.0400e+00	1.8704e-04
    1.0500e+00	1.3620e-04
    1.0600e+00	9.3631e-05
    1.0700e+00	6.8118e-05
    1.0800e+00	5.3285e-05
    1.0900e+00	4.0768e-05
    1.1000e+00	3.0637e-05
    1.1100e+00	2.3849e-05
    1.1200e+00	1.7825e-05
    1.1300e+00	1.3488e-05
    1.1400e+00	9.0718e-06
    1.1500e+00	6.2230e-06
    1.1600e+00	3.8770e-06
    1.1700e+00	2.0483e-06
    1.1800e+00	6.1036e-07
    1.1900e+00	3.4091e-07
    1.2000e+00	2.1008e-07
    1.2100e+00	1.2518e-07
    1.2200e+00	7.9609e-08
    1.2300e+00	4.6004e-08
    1.2400e+00	2.3682e-08
    1.2500e+00	9.9472e-09
    1.2600e+00	3.6096e-09
    1.2700e+00	2.0213e-09
    1.2800e+00	1.2223e-09
    1.2900e+00	7.2885e-10
    1.3000e+00	4.6553e-10
    1.3100e+00	2.8147e-10
    1.3200e+00	1.6807e-10
    1.3300e+00	8.4670e-11
    1.3400e+00	3.7322e-11
    1.3500e+00	1.8263e-11
    1.3600e+00	1.0281e-11
    1.3700e+00	6.5413e-12
    1.3800e+00	4.1730e-12
    1.3900e+00	2.5441e-12
    1.4000e+00	1.5597e-12
    1.4100e+00	9.5374e-13
    1.4200e+00	5.6500e-13
    1.4300e+00	2.8449e-13
    1.4400e+00	2.0626e-13
    1.4500e+00	0
    2.0000e+00	0];
n = interp1(wn(:,1),wn(:,2),wl*1e6,'linear');
k = interp1(wk(:,1),wk(:,2),wl*1e6,'linear');
alpha = 4*pi*k/wl*1e-2;  %(1/cm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n, k, alpha] = get_Ge_n_k(wl)
wn = [1.879E-1	8.351E-1
    1.881E-1	8.394E-1
    1.884E-1	8.434E-1
    1.887E-1	8.475E-1
    1.890E-1	8.518E-1
    1.893E-1	8.558E-1
    1.896E-1	8.598E-1
    1.899E-1	8.640E-1
    1.902E-1	8.679E-1
    1.905E-1	8.720E-1
    1.907E-1	8.760E-1
    1.910E-1	8.800E-1
    1.913E-1	8.840E-1
    1.916E-1	8.880E-1
    1.919E-1	8.918E-1
    1.922E-1	8.958E-1
    1.925E-1	8.998E-1
    1.928E-1	9.035E-1
    1.931E-1	9.074E-1
    1.934E-1	9.113E-1
    1.937E-1	9.152E-1
    1.940E-1	9.189E-1
    1.943E-1	9.228E-1
    1.946E-1	9.265E-1
    1.949E-1	9.302E-1
    1.953E-1	9.340E-1
    1.956E-1	9.377E-1
    1.959E-1	9.414E-1
    1.962E-1	9.450E-1
    1.965E-1	9.487E-1
    1.968E-1	9.523E-1
    1.971E-1	9.560E-1
    1.974E-1	9.598E-1
    1.977E-1	9.633E-1
    1.981E-1	9.668E-1
    1.984E-1	9.705E-1
    1.987E-1	9.741E-1
    1.990E-1	9.779E-1
    1.993E-1	9.814E-1
    1.997E-1	9.851E-1
    2.000E-1	9.888E-1
    2.003E-1	9.925E-1
    2.006E-1	9.964E-1
    2.009E-1	1.000E+0
    2.013E-1	1.004E+0
    2.016E-1	1.008E+0
    2.019E-1	1.012E+0
    2.023E-1	1.016E+0
    2.026E-1	1.021E+0
    2.029E-1	1.025E+0
    2.033E-1	1.030E+0
    2.036E-1	1.034E+0
    2.039E-1	1.039E+0
    2.043E-1	1.045E+0
    2.046E-1	1.050E+0
    2.049E-1	1.056E+0
    2.053E-1	1.061E+0
    2.056E-1	1.067E+0
    2.060E-1	1.074E+0
    2.063E-1	1.080E+0
    2.066E-1	1.087E+0
    2.070E-1	1.094E+0
    2.073E-1	1.102E+0
    2.077E-1	1.109E+0
    2.080E-1	1.117E+0
    2.084E-1	1.125E+0
    2.087E-1	1.134E+0
    2.091E-1	1.143E+0
    2.094E-1	1.152E+0
    2.098E-1	1.161E+0
    2.101E-1	1.170E+0
    2.105E-1	1.180E+0
    2.109E-1	1.190E+0
    2.112E-1	1.200E+0
    2.116E-1	1.210E+0
    2.119E-1	1.220E+0
    2.123E-1	1.231E+0
    2.127E-1	1.241E+0
    2.130E-1	1.251E+0
    2.134E-1	1.262E+0
    2.138E-1	1.272E+0
    2.141E-1	1.283E+0
    2.145E-1	1.293E+0
    2.149E-1	1.303E+0
    2.153E-1	1.312E+0
    2.156E-1	1.322E+0
    2.160E-1	1.332E+0
    2.164E-1	1.341E+0
    2.168E-1	1.350E+0
    2.171E-1	1.358E+0
    2.175E-1	1.366E+0
    2.179E-1	1.375E+0
    2.183E-1	1.382E+0
    2.187E-1	1.389E+0
    2.191E-1	1.396E+0
    2.194E-1	1.402E+0
    2.198E-1	1.408E+0
    2.202E-1	1.414E+0
    2.206E-1	1.420E+0
    2.210E-1	1.425E+0
    2.214E-1	1.429E+0
    2.218E-1	1.433E+0
    2.222E-1	1.437E+0
    2.226E-1	1.440E+0
    2.230E-1	1.444E+0
    2.234E-1	1.446E+0
    2.238E-1	1.449E+0
    2.242E-1	1.451E+0
    2.246E-1	1.453E+0
    2.250E-1	1.455E+0
    2.254E-1	1.456E+0
    2.258E-1	1.457E+0
    2.262E-1	1.458E+0
    2.267E-1	1.459E+0
    2.271E-1	1.459E+0
    2.275E-1	1.460E+0
    2.279E-1	1.460E+0
    2.283E-1	1.460E+0
    2.288E-1	1.460E+0
    2.292E-1	1.460E+0
    2.296E-1	1.459E+0
    2.300E-1	1.459E+0
    2.305E-1	1.459E+0
    2.309E-1	1.458E+0
    2.313E-1	1.457E+0
    2.317E-1	1.457E+0
    2.322E-1	1.456E+0
    2.326E-1	1.456E+0
    2.331E-1	1.455E+0
    2.335E-1	1.455E+0
    2.339E-1	1.454E+0
    2.344E-1	1.453E+0
    2.348E-1	1.453E+0
    2.353E-1	1.452E+0
    2.357E-1	1.452E+0
    2.362E-1	1.451E+0
    2.366E-1	1.451E+0
    2.371E-1	1.450E+0
    2.375E-1	1.450E+0
    2.380E-1	1.450E+0
    2.384E-1	1.450E+0
    2.389E-1	1.450E+0
    2.394E-1	1.450E+0
    2.398E-1	1.450E+0
    2.403E-1	1.450E+0
    2.407E-1	1.451E+0
    2.412E-1	1.451E+0
    2.417E-1	1.452E+0
    2.422E-1	1.453E+0
    2.426E-1	1.454E+0
    2.431E-1	1.455E+0
    2.436E-1	1.456E+0
    2.441E-1	1.457E+0
    2.445E-1	1.459E+0
    2.450E-1	1.461E+0
    2.455E-1	1.463E+0
    2.460E-1	1.465E+0
    2.465E-1	1.467E+0
    2.470E-1	1.470E+0
    2.475E-1	1.472E+0
    2.480E-1	1.475E+0
    2.485E-1	1.478E+0
    2.490E-1	1.481E+0
    2.495E-1	1.485E+0
    2.500E-1	1.489E+0
    2.505E-1	1.493E+0
    2.510E-1	1.497E+0
    2.515E-1	1.502E+0
    2.520E-1	1.506E+0
    2.525E-1	1.511E+0
    2.530E-1	1.516E+0
    2.535E-1	1.522E+0
    2.541E-1	1.527E+0
    2.546E-1	1.534E+0
    2.551E-1	1.540E+0
    2.556E-1	1.546E+0
    2.562E-1	1.553E+0
    2.567E-1	1.560E+0
    2.572E-1	1.568E+0
    2.578E-1	1.575E+0
    2.583E-1	1.583E+0
    2.588E-1	1.592E+0
    2.594E-1	1.600E+0
    2.599E-1	1.609E+0
    2.605E-1	1.618E+0
    2.610E-1	1.627E+0
    2.616E-1	1.637E+0
    2.621E-1	1.647E+0
    2.627E-1	1.657E+0
    2.632E-1	1.667E+0
    2.638E-1	1.678E+0
    2.644E-1	1.689E+0
    2.649E-1	1.701E+0
    2.655E-1	1.712E+0
    2.661E-1	1.724E+0
    2.666E-1	1.737E+0
    2.672E-1	1.749E+0
    2.678E-1	1.762E+0
    2.684E-1	1.776E+0
    2.689E-1	1.791E+0
    2.695E-1	1.806E+0
    2.701E-1	1.822E+0
    2.707E-1	1.839E+0
    2.713E-1	1.857E+0
    2.719E-1	1.876E+0
    2.725E-1	1.898E+0
    2.731E-1	1.921E+0
    2.737E-1	1.947E+0
    2.743E-1	1.975E+0
    2.749E-1	2.006E+0
    2.755E-1	2.041E+0
    2.761E-1	2.079E+0
    2.768E-1	2.120E+0
    2.774E-1	2.166E+0
    2.780E-1	2.216E+0
    2.786E-1	2.270E+0
    2.792E-1	2.329E+0
    2.799E-1	2.391E+0
    2.805E-1	2.457E+0
    2.811E-1	2.527E+0
    2.818E-1	2.600E+0
    2.824E-1	2.676E+0
    2.831E-1	2.754E+0
    2.837E-1	2.833E+0
    2.844E-1	2.913E+0
    2.850E-1	2.993E+0
    2.857E-1	3.073E+0
    2.863E-1	3.152E+0
    2.870E-1	3.228E+0
    2.877E-1	3.302E+0
    2.883E-1	3.372E+0
    2.890E-1	3.439E+0
    2.897E-1	3.501E+0
    2.904E-1	3.559E+0
    2.910E-1	3.613E+0
    2.917E-1	3.661E+0
    2.924E-1	3.704E+0
    2.931E-1	3.742E+0
    2.938E-1	3.776E+0
    2.945E-1	3.805E+0
    2.952E-1	3.829E+0
    2.959E-1	3.850E+0
    2.966E-1	3.867E+0
    2.973E-1	3.882E+0
    2.980E-1	3.893E+0
    2.988E-1	3.903E+0
    2.995E-1	3.910E+0
    3.002E-1	3.916E+0
    3.009E-1	3.921E+0
    3.017E-1	3.925E+0
    3.024E-1	3.929E+0
    3.031E-1	3.932E+0
    3.039E-1	3.935E+0
    3.046E-1	3.938E+0
    3.054E-1	3.941E+0
    3.061E-1	3.944E+0
    3.069E-1	3.947E+0
    3.077E-1	3.950E+0
    3.084E-1	3.953E+0
    3.092E-1	3.957E+0
    3.100E-1	3.960E+0
    3.107E-1	3.963E+0
    3.115E-1	3.967E+0
    3.123E-1	3.970E+0
    3.131E-1	3.973E+0
    3.139E-1	3.977E+0
    3.147E-1	3.980E+0
    3.155E-1	3.983E+0
    3.163E-1	3.986E+0
    3.171E-1	3.989E+0
    3.179E-1	3.992E+0
    3.187E-1	3.994E+0
    3.195E-1	3.997E+0
    3.204E-1	3.999E+0
    3.212E-1	4.002E+0
    3.220E-1	4.004E+0
    3.229E-1	4.006E+0
    3.237E-1	4.008E+0
    3.246E-1	4.010E+0
    3.254E-1	4.011E+0
    3.263E-1	4.013E+0
    3.271E-1	4.015E+0
    3.280E-1	4.016E+0
    3.289E-1	4.018E+0
    3.297E-1	4.019E+0
    3.306E-1	4.020E+0
    3.315E-1	4.022E+0
    3.324E-1	4.023E+0
    3.333E-1	4.024E+0
    3.342E-1	4.025E+0
    3.351E-1	4.026E+0
    3.360E-1	4.028E+0
    3.369E-1	4.029E+0
    3.378E-1	4.030E+0
    3.388E-1	4.031E+0
    3.397E-1	4.033E+0
    3.406E-1	4.034E+0
    3.416E-1	4.036E+0
    3.425E-1	4.037E+0
    3.434E-1	4.039E+0
    3.444E-1	4.041E+0
    3.454E-1	4.043E+0
    3.463E-1	4.045E+0
    3.473E-1	4.047E+0
    3.483E-1	4.050E+0
    3.493E-1	4.053E+0
    3.502E-1	4.056E+0
    3.512E-1	4.059E+0
    3.522E-1	4.063E+0
    3.532E-1	4.067E+0
    3.542E-1	4.071E+0
    3.553E-1	4.076E+0
    3.563E-1	4.081E+0
    3.573E-1	4.086E+0
    3.583E-1	4.092E+0
    3.594E-1	4.098E+0
    3.604E-1	4.104E+0
    3.615E-1	4.111E+0
    3.625E-1	4.117E+0
    3.636E-1	4.124E+0
    3.647E-1	4.131E+0
    3.657E-1	4.138E+0
    3.668E-1	4.144E+0
    3.679E-1	4.151E+0
    3.690E-1	4.157E+0
    3.701E-1	4.163E+0
    3.712E-1	4.168E+0
    3.723E-1	4.173E+0
    3.734E-1	4.177E+0
    3.746E-1	4.181E+0
    3.757E-1	4.184E+0
    3.769E-1	4.187E+0
    3.780E-1	4.189E+0
    3.792E-1	4.191E+0
    3.803E-1	4.194E+0
    3.815E-1	4.196E+0
    3.827E-1	4.199E+0
    3.839E-1	4.202E+0
    3.850E-1	4.204E+0
    3.862E-1	4.206E+0
    3.875E-1	4.208E+0
    3.887E-1	4.210E+0
    3.899E-1	4.211E+0
    3.911E-1	4.211E+0
    3.924E-1	4.211E+0
    3.936E-1	4.209E+0
    3.949E-1	4.207E+0
    3.961E-1	4.205E+0
    3.974E-1	4.201E+0
    3.987E-1	4.197E+0
    3.999E-1	4.192E+0
    4.012E-1	4.187E+0
    4.025E-1	4.181E+0
    4.039E-1	4.175E+0
    4.052E-1	4.168E+0
    4.065E-1	4.161E+0
    4.078E-1	4.154E+0
    4.092E-1	4.147E+0
    4.105E-1	4.140E+0
    4.119E-1	4.133E+0
    4.133E-1	4.126E+0
    4.147E-1	4.119E+0
    4.161E-1	4.113E+0
    4.175E-1	4.107E+0
    4.189E-1	4.101E+0
    4.203E-1	4.096E+0
    4.217E-1	4.091E+0
    4.232E-1	4.086E+0
    4.246E-1	4.082E+0
    4.261E-1	4.079E+0
    4.275E-1	4.076E+0
    4.290E-1	4.073E+0
    4.305E-1	4.071E+0
    4.320E-1	4.070E+0
    4.335E-1	4.069E+0
    4.350E-1	4.068E+0
    4.366E-1	4.069E+0
    4.381E-1	4.069E+0
    4.397E-1	4.071E+0
    4.412E-1	4.072E+0
    4.428E-1	4.075E+0
    4.444E-1	4.077E+0
    4.460E-1	4.081E+0
    4.476E-1	4.084E+0
    4.492E-1	4.089E+0
    4.509E-1	4.094E+0
    4.525E-1	4.099E+0
    4.542E-1	4.105E+0
    4.558E-1	4.111E+0
    4.575E-1	4.117E+0
    4.592E-1	4.125E+0
    4.609E-1	4.132E+0
    4.626E-1	4.140E+0
    4.644E-1	4.149E+0
    4.661E-1	4.158E+0
    4.679E-1	4.167E+0
    4.696E-1	4.177E+0
    4.714E-1	4.188E+0
    4.732E-1	4.199E+0
    4.750E-1	4.210E+0
    4.769E-1	4.222E+0
    4.787E-1	4.235E+0
    4.806E-1	4.248E+0
    4.824E-1	4.261E+0
    4.843E-1	4.276E+0
    4.862E-1	4.290E+0
    4.881E-1	4.306E+0
    4.901E-1	4.322E+0
    4.920E-1	4.339E+0
    4.940E-1	4.357E+0
    4.959E-1	4.376E+0
    4.979E-1	4.395E+0
    4.999E-1	4.416E+0
    5.020E-1	4.439E+0
    5.040E-1	4.463E+0
    5.061E-1	4.488E+0
    5.081E-1	4.516E+0
    5.102E-1	4.546E+0
    5.123E-1	4.578E+0
    5.145E-1	4.613E+0
    5.166E-1	4.650E+0
    5.188E-1	4.690E+0
    5.209E-1	4.731E+0
    5.231E-1	4.775E+0
    5.254E-1	4.820E+0
    5.276E-1	4.865E+0
    5.298E-1	4.910E+0
    5.321E-1	4.953E+0
    5.344E-1	4.996E+0
    5.367E-1	5.035E+0
    5.391E-1	5.072E+0
    5.414E-1	5.105E+0
    5.438E-1	5.134E+0
    5.462E-1	5.160E+0
    5.486E-1	5.183E+0
    5.510E-1	5.202E+0
    5.535E-1	5.219E+0
    5.560E-1	5.236E+0
    5.585E-1	5.252E+0
    5.610E-1	5.270E+0
    5.636E-1	5.292E+0
    5.661E-1	5.319E+0
    5.687E-1	5.353E+0
    5.714E-1	5.396E+0
    5.740E-1	5.446E+0
    5.767E-1	5.504E+0
    5.794E-1	5.564E+0
    5.821E-1	5.625E+0
    5.848E-1	5.680E+0
    5.876E-1	5.725E+0
    5.904E-1	5.758E+0
    5.932E-1	5.778E+0
    5.961E-1	5.785E+0
    5.990E-1	5.780E+0
    6.019E-1	5.766E+0
    6.048E-1	5.745E+0
    6.078E-1	5.719E+0
    6.108E-1	5.690E+0
    6.138E-1	5.660E+0
    6.168E-1	5.628E+0
    6.199E-1	5.596E+0
    6.230E-1	5.564E+0
    6.262E-1	5.532E+0
    6.294E-1	5.500E+0
    6.326E-1	5.468E+0
    6.358E-1	5.437E+0
    6.391E-1	5.407E+0
    6.424E-1	5.377E+0
    6.458E-1	5.348E+0
    6.491E-1	5.319E+0
    6.525E-1	5.292E+0
    6.560E-1	5.265E+0
    6.595E-1	5.239E+0
    6.630E-1	5.214E+0
    6.666E-1	5.190E+0
    6.702E-1	5.167E+0
    6.738E-1	5.144E+0
    6.775E-1	5.122E+0
    6.812E-1	5.101E+0
    6.850E-1	5.081E+0
    6.888E-1	5.061E+0
    6.926E-1	5.041E+0
    6.965E-1	5.022E+0
    7.005E-1	5.004E+0
    7.045E-1	4.986E+0
    7.085E-1	4.969E+0
    7.126E-1	4.951E+0
    7.167E-1	4.935E+0
    7.208E-1	4.918E+0
    7.251E-1	4.902E+0
    7.293E-1	4.887E+0
    7.336E-1	4.871E+0
    7.380E-1	4.856E+0
    7.424E-1	4.841E+0
    7.469E-1	4.827E+0
    7.514E-1	4.812E+0
    7.560E-1	4.798E+0
    7.606E-1	4.785E+0
    7.653E-1	4.771E+0
    7.701E-1	4.758E+0
    7.749E-1	4.745E+0
    7.798E-1	4.733E+0
    7.847E-1	4.721E+0
    7.897E-1	4.708E+0
    7.948E-1	4.697E+0
    7.999E-1	4.685E+0
    8.051E-1	4.674E+0
    8.104E-1	4.663E+0
    8.157E-1	4.652E+0
    8.211E-1	4.641E+0
    8.266E-1	4.631E+0
    8.321E-1	4.621E+0
    8.377E-1	4.611E+0
    8.434E-1	4.601E+0
    8.492E-1	4.591E+0
    8.551E-1	4.582E+0
    8.610E-1	4.573E+0
    8.670E-1	4.564E+0
    8.731E-1	4.555E+0
    8.793E-1	4.546E+0
    8.856E-1	4.538E+0
    8.920E-1	4.530E+0
    8.984E-1	4.522E+0
    9.050E-1	4.514E+0
    9.116E-1	4.506E+0
    9.184E-1	4.499E+0
    9.253E-1	4.492E+0
    9.322E-1	4.484E+0
    9.393E-1	4.477E+0
    9.464E-1	4.470E+0
    9.537E-1	4.464E+0
    9.611E-1	4.457E+0
    9.686E-1	4.450E+0
    9.763E-1	4.444E+0
    9.840E-1	4.438E+0
    9.919E-1	4.432E+0
    9.999E-1	4.426E+0
    1.008E+0	4.420E+0
    1.016E+0	4.414E+0
    1.025E+0	4.408E+0
    1.033E+0	4.403E+0
    1.042E+0	4.397E+0
    1.051E+0	4.392E+0
    1.060E+0	4.387E+0
    1.069E+0	4.382E+0
    1.078E+0	4.377E+0
    1.088E+0	4.372E+0
    1.097E+0	4.368E+0
    1.107E+0	4.364E+0
    1.117E+0	4.360E+0
    1.127E+0	4.355E+0
    1.137E+0	4.350E+0
    1.148E+0	4.344E+0
    1.159E+0	4.337E+0
    1.170E+0	4.331E+0
    1.181E+0	4.326E+0
    1.192E+0	4.320E+0
    1.204E+0	4.315E+0
    1.216E+0	4.310E+0
    1.228E+0	4.306E+0
    1.240E+0	4.301E+0
    1.252E+0	4.296E+0
    1.265E+0	4.292E+0
    1.278E+0	4.287E+0
    1.292E+0	4.283E+0
    1.305E+0	4.279E+0
    1.319E+0	4.275E+0
    1.333E+0	4.271E+0
    1.348E+0	4.267E+0
    1.362E+0	4.263E+0
    1.378E+0	4.259E+0
    1.393E+0	4.256E+0
    1.409E+0	4.253E+0
    1.425E+0	4.249E+0
    1.442E+0	4.246E+0
    1.459E+0	4.244E+0
    1.476E+0	4.242E+0
    1.494E+0	4.241E+0
    1.512E+0	4.243E+0
    1.531E+0	4.249E+0
    1.550E+0	4.250E+0
    1.569E+0	4.234E+0
    1.590E+0	4.216E+0
    1.610E+0	4.204E+0
    1.631E+0	4.194E+0
    1.653E+0	4.186E+0
    1.675E+0	4.178E+0
    1.698E+0	4.171E+0
    1.722E+0	4.164E+0
    1.746E+0	4.158E+0
    1.771E+0	4.152E+0
    1.797E+0	4.146E+0
    1.823E+0	4.141E+0
    1.851E+0	4.136E+0
    1.879E+0	4.131E+0
    1.907E+0	4.126E+0
    1.937E+0	4.122E+0
    1.968E+0	4.117E+0
    2.000E+0	4.113E+0
    2.033E+0	4.109E+0
    2.066E+0	4.105E+0
    2.101E+0	4.101E+0
    2.138E+0	4.097E+0
    2.175E+0	4.094E+0
    2.214E+0	4.090E+0
    2.254E+0	4.087E+0
    2.296E+0	4.083E+0
    2.339E+0	4.080E+0
    2.384E+0	4.077E+0
    2.431E+0	4.074E+0
    2.480E+0	4.071E+0];
wk = [1.879E-1	2.567E+0
    1.881E-1	2.570E+0
    1.884E-1	2.573E+0
    1.887E-1	2.576E+0
    1.890E-1	2.579E+0
    1.893E-1	2.582E+0
    1.896E-1	2.585E+0
    1.899E-1	2.588E+0
    1.902E-1	2.592E+0
    1.905E-1	2.595E+0
    1.907E-1	2.598E+0
    1.910E-1	2.601E+0
    1.913E-1	2.604E+0
    1.916E-1	2.607E+0
    1.919E-1	2.611E+0
    1.922E-1	2.614E+0
    1.925E-1	2.617E+0
    1.928E-1	2.621E+0
    1.931E-1	2.624E+0
    1.934E-1	2.628E+0
    1.937E-1	2.631E+0
    1.940E-1	2.635E+0
    1.943E-1	2.638E+0
    1.946E-1	2.642E+0
    1.949E-1	2.646E+0
    1.953E-1	2.649E+0
    1.956E-1	2.653E+0
    1.959E-1	2.657E+0
    1.962E-1	2.661E+0
    1.965E-1	2.665E+0
    1.968E-1	2.669E+0
    1.971E-1	2.673E+0
    1.974E-1	2.676E+0
    1.977E-1	2.681E+0
    1.981E-1	2.686E+0
    1.984E-1	2.690E+0
    1.987E-1	2.695E+0
    1.990E-1	2.699E+0
    1.993E-1	2.704E+0
    1.997E-1	2.708E+0
    2.000E-1	2.713E+0
    2.003E-1	2.719E+0
    2.006E-1	2.724E+0
    2.009E-1	2.729E+0
    2.013E-1	2.734E+0
    2.016E-1	2.740E+0
    2.019E-1	2.746E+0
    2.023E-1	2.752E+0
    2.026E-1	2.757E+0
    2.029E-1	2.763E+0
    2.033E-1	2.770E+0
    2.036E-1	2.777E+0
    2.039E-1	2.783E+0
    2.043E-1	2.789E+0
    2.046E-1	2.796E+0
    2.049E-1	2.802E+0
    2.053E-1	2.809E+0
    2.056E-1	2.816E+0
    2.060E-1	2.823E+0
    2.063E-1	2.830E+0
    2.066E-1	2.837E+0
    2.070E-1	2.843E+0
    2.073E-1	2.850E+0
    2.077E-1	2.857E+0
    2.080E-1	2.863E+0
    2.084E-1	2.870E+0
    2.087E-1	2.876E+0
    2.091E-1	2.882E+0
    2.094E-1	2.888E+0
    2.098E-1	2.893E+0
    2.101E-1	2.899E+0
    2.105E-1	2.904E+0
    2.109E-1	2.908E+0
    2.112E-1	2.913E+0
    2.116E-1	2.917E+0
    2.119E-1	2.921E+0
    2.123E-1	2.924E+0
    2.127E-1	2.927E+0
    2.130E-1	2.930E+0
    2.134E-1	2.932E+0
    2.138E-1	2.934E+0
    2.141E-1	2.936E+0
    2.145E-1	2.936E+0
    2.149E-1	2.937E+0
    2.153E-1	2.938E+0
    2.156E-1	2.937E+0
    2.160E-1	2.937E+0
    2.164E-1	2.937E+0
    2.168E-1	2.936E+0
    2.171E-1	2.934E+0
    2.175E-1	2.933E+0
    2.179E-1	2.931E+0
    2.183E-1	2.930E+0
    2.187E-1	2.928E+0
    2.191E-1	2.926E+0
    2.194E-1	2.924E+0
    2.198E-1	2.922E+0
    2.202E-1	2.919E+0
    2.206E-1	2.917E+0
    2.210E-1	2.915E+0
    2.214E-1	2.913E+0
    2.218E-1	2.910E+0
    2.222E-1	2.908E+0
    2.226E-1	2.907E+0
    2.230E-1	2.905E+0
    2.234E-1	2.903E+0
    2.238E-1	2.902E+0
    2.242E-1	2.901E+0
    2.246E-1	2.899E+0
    2.250E-1	2.899E+0
    2.254E-1	2.898E+0
    2.258E-1	2.898E+0
    2.262E-1	2.898E+0
    2.267E-1	2.898E+0
    2.271E-1	2.898E+0
    2.275E-1	2.899E+0
    2.279E-1	2.900E+0
    2.283E-1	2.901E+0
    2.288E-1	2.903E+0
    2.292E-1	2.905E+0
    2.296E-1	2.907E+0
    2.300E-1	2.909E+0
    2.305E-1	2.912E+0
    2.309E-1	2.915E+0
    2.313E-1	2.919E+0
    2.317E-1	2.923E+0
    2.322E-1	2.927E+0
    2.326E-1	2.931E+0
    2.331E-1	2.935E+0
    2.335E-1	2.940E+0
    2.339E-1	2.945E+0
    2.344E-1	2.950E+0
    2.348E-1	2.956E+0
    2.353E-1	2.962E+0
    2.357E-1	2.968E+0
    2.362E-1	2.975E+0
    2.366E-1	2.982E+0
    2.371E-1	2.989E+0
    2.375E-1	2.996E+0
    2.380E-1	3.004E+0
    2.384E-1	3.012E+0
    2.389E-1	3.020E+0
    2.394E-1	3.028E+0
    2.398E-1	3.037E+0
    2.403E-1	3.046E+0
    2.407E-1	3.055E+0
    2.412E-1	3.065E+0
    2.417E-1	3.075E+0
    2.422E-1	3.085E+0
    2.426E-1	3.095E+0
    2.431E-1	3.106E+0
    2.436E-1	3.116E+0
    2.441E-1	3.127E+0
    2.445E-1	3.139E+0
    2.450E-1	3.151E+0
    2.455E-1	3.162E+0
    2.460E-1	3.174E+0
    2.465E-1	3.187E+0
    2.470E-1	3.199E+0
    2.475E-1	3.212E+0
    2.480E-1	3.225E+0
    2.485E-1	3.238E+0
    2.490E-1	3.252E+0
    2.495E-1	3.265E+0
    2.500E-1	3.279E+0
    2.505E-1	3.294E+0
    2.510E-1	3.308E+0
    2.515E-1	3.323E+0
    2.520E-1	3.337E+0
    2.525E-1	3.353E+0
    2.530E-1	3.368E+0
    2.535E-1	3.383E+0
    2.541E-1	3.399E+0
    2.546E-1	3.415E+0
    2.551E-1	3.431E+0
    2.556E-1	3.448E+0
    2.562E-1	3.464E+0
    2.567E-1	3.481E+0
    2.572E-1	3.498E+0
    2.578E-1	3.516E+0
    2.583E-1	3.533E+0
    2.588E-1	3.551E+0
    2.594E-1	3.569E+0
    2.599E-1	3.587E+0
    2.605E-1	3.606E+0
    2.610E-1	3.625E+0
    2.616E-1	3.644E+0
    2.621E-1	3.663E+0
    2.627E-1	3.683E+0
    2.632E-1	3.703E+0
    2.638E-1	3.724E+0
    2.644E-1	3.745E+0
    2.649E-1	3.766E+0
    2.655E-1	3.788E+0
    2.661E-1	3.810E+0
    2.666E-1	3.833E+0
    2.672E-1	3.857E+0
    2.678E-1	3.882E+0
    2.684E-1	3.907E+0
    2.689E-1	3.933E+0
    2.695E-1	3.961E+0
    2.701E-1	3.989E+0
    2.707E-1	4.019E+0
    2.713E-1	4.050E+0
    2.719E-1	4.082E+0
    2.725E-1	4.116E+0
    2.731E-1	4.151E+0
    2.737E-1	4.187E+0
    2.743E-1	4.224E+0
    2.749E-1	4.262E+0
    2.755E-1	4.300E+0
    2.761E-1	4.339E+0
    2.768E-1	4.378E+0
    2.774E-1	4.416E+0
    2.780E-1	4.454E+0
    2.786E-1	4.489E+0
    2.792E-1	4.523E+0
    2.799E-1	4.554E+0
    2.805E-1	4.581E+0
    2.811E-1	4.604E+0
    2.818E-1	4.624E+0
    2.824E-1	4.638E+0
    2.831E-1	4.647E+0
    2.837E-1	4.651E+0
    2.844E-1	4.648E+0
    2.850E-1	4.640E+0
    2.857E-1	4.626E+0
    2.863E-1	4.606E+0
    2.870E-1	4.580E+0
    2.877E-1	4.549E+0
    2.883E-1	4.513E+0
    2.890E-1	4.473E+0
    2.897E-1	4.429E+0
    2.904E-1	4.381E+0
    2.910E-1	4.331E+0
    2.917E-1	4.279E+0
    2.924E-1	4.226E+0
    2.931E-1	4.172E+0
    2.938E-1	4.118E+0
    2.945E-1	4.064E+0
    2.952E-1	4.012E+0
    2.959E-1	3.961E+0
    2.966E-1	3.911E+0
    2.973E-1	3.863E+0
    2.980E-1	3.818E+0
    2.988E-1	3.774E+0
    2.995E-1	3.733E+0
    3.002E-1	3.694E+0
    3.009E-1	3.658E+0
    3.017E-1	3.623E+0
    3.024E-1	3.590E+0
    3.031E-1	3.559E+0
    3.039E-1	3.530E+0
    3.046E-1	3.502E+0
    3.054E-1	3.475E+0
    3.061E-1	3.449E+0
    3.069E-1	3.425E+0
    3.077E-1	3.401E+0
    3.084E-1	3.377E+0
    3.092E-1	3.355E+0
    3.100E-1	3.333E+0
    3.107E-1	3.311E+0
    3.115E-1	3.290E+0
    3.123E-1	3.270E+0
    3.131E-1	3.249E+0
    3.139E-1	3.229E+0
    3.147E-1	3.210E+0
    3.155E-1	3.191E+0
    3.163E-1	3.172E+0
    3.171E-1	3.153E+0
    3.179E-1	3.135E+0
    3.187E-1	3.117E+0
    3.195E-1	3.099E+0
    3.204E-1	3.082E+0
    3.212E-1	3.065E+0
    3.220E-1	3.048E+0
    3.229E-1	3.031E+0
    3.237E-1	3.015E+0
    3.246E-1	2.999E+0
    3.254E-1	2.984E+0
    3.263E-1	2.969E+0
    3.271E-1	2.954E+0
    3.280E-1	2.939E+0
    3.289E-1	2.925E+0
    3.297E-1	2.911E+0
    3.306E-1	2.897E+0
    3.315E-1	2.884E+0
    3.324E-1	2.871E+0
    3.333E-1	2.858E+0
    3.342E-1	2.845E+0
    3.351E-1	2.833E+0
    3.360E-1	2.821E+0
    3.369E-1	2.810E+0
    3.378E-1	2.798E+0
    3.388E-1	2.787E+0
    3.397E-1	2.777E+0
    3.406E-1	2.766E+0
    3.416E-1	2.756E+0
    3.425E-1	2.746E+0
    3.434E-1	2.736E+0
    3.444E-1	2.727E+0
    3.454E-1	2.718E+0
    3.463E-1	2.709E+0
    3.473E-1	2.700E+0
    3.483E-1	2.692E+0
    3.493E-1	2.684E+0
    3.502E-1	2.675E+0
    3.512E-1	2.667E+0
    3.522E-1	2.659E+0
    3.532E-1	2.652E+0
    3.542E-1	2.644E+0
    3.553E-1	2.636E+0
    3.563E-1	2.628E+0
    3.573E-1	2.620E+0
    3.583E-1	2.612E+0
    3.594E-1	2.604E+0
    3.604E-1	2.595E+0
    3.615E-1	2.586E+0
    3.625E-1	2.577E+0
    3.636E-1	2.567E+0
    3.647E-1	2.557E+0
    3.657E-1	2.546E+0
    3.668E-1	2.535E+0
    3.679E-1	2.523E+0
    3.690E-1	2.511E+0
    3.701E-1	2.499E+0
    3.712E-1	2.486E+0
    3.723E-1	2.473E+0
    3.734E-1	2.459E+0
    3.746E-1	2.446E+0
    3.757E-1	2.433E+0
    3.769E-1	2.421E+0
    3.780E-1	2.408E+0
    3.792E-1	2.396E+0
    3.803E-1	2.384E+0
    3.815E-1	2.372E+0
    3.827E-1	2.360E+0
    3.839E-1	2.348E+0
    3.850E-1	2.336E+0
    3.862E-1	2.323E+0
    3.875E-1	2.310E+0
    3.887E-1	2.297E+0
    3.899E-1	2.283E+0
    3.911E-1	2.269E+0
    3.924E-1	2.256E+0
    3.936E-1	2.242E+0
    3.949E-1	2.228E+0
    3.961E-1	2.215E+0
    3.974E-1	2.202E+0
    3.987E-1	2.190E+0
    3.999E-1	2.178E+0
    4.012E-1	2.167E+0
    4.025E-1	2.157E+0
    4.039E-1	2.147E+0
    4.052E-1	2.138E+0
    4.065E-1	2.131E+0
    4.078E-1	2.124E+0
    4.092E-1	2.117E+0
    4.105E-1	2.112E+0
    4.119E-1	2.107E+0
    4.133E-1	2.104E+0
    4.147E-1	2.101E+0
    4.161E-1	2.099E+0
    4.175E-1	2.097E+0
    4.189E-1	2.097E+0
    4.203E-1	2.097E+0
    4.217E-1	2.097E+0
    4.232E-1	2.098E+0
    4.246E-1	2.100E+0
    4.261E-1	2.102E+0
    4.275E-1	2.105E+0
    4.290E-1	2.107E+0
    4.305E-1	2.111E+0
    4.320E-1	2.115E+0
    4.335E-1	2.119E+0
    4.350E-1	2.123E+0
    4.366E-1	2.127E+0
    4.381E-1	2.132E+0
    4.397E-1	2.137E+0
    4.412E-1	2.142E+0
    4.428E-1	2.147E+0
    4.444E-1	2.152E+0
    4.460E-1	2.158E+0
    4.476E-1	2.164E+0
    4.492E-1	2.169E+0
    4.509E-1	2.175E+0
    4.525E-1	2.181E+0
    4.542E-1	2.187E+0
    4.558E-1	2.193E+0
    4.575E-1	2.199E+0
    4.592E-1	2.204E+0
    4.609E-1	2.211E+0
    4.626E-1	2.217E+0
    4.644E-1	2.223E+0
    4.661E-1	2.229E+0
    4.679E-1	2.235E+0
    4.696E-1	2.242E+0
    4.714E-1	2.248E+0
    4.732E-1	2.254E+0
    4.750E-1	2.261E+0
    4.769E-1	2.267E+0
    4.787E-1	2.274E+0
    4.806E-1	2.280E+0
    4.824E-1	2.287E+0
    4.843E-1	2.294E+0
    4.862E-1	2.300E+0
    4.881E-1	2.307E+0
    4.901E-1	2.314E+0
    4.920E-1	2.322E+0
    4.940E-1	2.329E+0
    4.959E-1	2.337E+0
    4.979E-1	2.345E+0
    4.999E-1	2.353E+0
    5.020E-1	2.361E+0
    5.040E-1	2.369E+0
    5.061E-1	2.378E+0
    5.081E-1	2.386E+0
    5.102E-1	2.394E+0
    5.123E-1	2.400E+0
    5.145E-1	2.406E+0
    5.166E-1	2.410E+0
    5.188E-1	2.411E+0
    5.209E-1	2.410E+0
    5.231E-1	2.405E+0
    5.254E-1	2.397E+0
    5.276E-1	2.385E+0
    5.298E-1	2.368E+0
    5.321E-1	2.348E+0
    5.344E-1	2.323E+0
    5.367E-1	2.296E+0
    5.391E-1	2.265E+0
    5.414E-1	2.233E+0
    5.438E-1	2.200E+0
    5.462E-1	2.167E+0
    5.486E-1	2.134E+0
    5.510E-1	2.104E+0
    5.535E-1	2.076E+0
    5.560E-1	2.051E+0
    5.585E-1	2.031E+0
    5.610E-1	2.014E+0
    5.636E-1	2.000E+0
    5.661E-1	1.989E+0
    5.687E-1	1.979E+0
    5.714E-1	1.966E+0
    5.740E-1	1.948E+0
    5.767E-1	1.920E+0
    5.794E-1	1.880E+0
    5.821E-1	1.826E+0
    5.848E-1	1.759E+0
    5.876E-1	1.681E+0
    5.904E-1	1.595E+0
    5.932E-1	1.505E+0
    5.961E-1	1.416E+0
    5.990E-1	1.329E+0
    6.019E-1	1.249E+0
    6.048E-1	1.175E+0
    6.078E-1	1.108E+0
    6.108E-1	1.047E+0
    6.138E-1	9.913E-1
    6.168E-1	9.408E-1
    6.199E-1	8.947E-1
    6.230E-1	8.524E-1
    6.262E-1	8.137E-1
    6.294E-1	7.780E-1
    6.326E-1	7.453E-1
    6.358E-1	7.151E-1
    6.391E-1	6.874E-1
    6.424E-1	6.618E-1
    6.458E-1	6.382E-1
    6.491E-1	6.164E-1
    6.525E-1	5.961E-1
    6.560E-1	5.774E-1
    6.595E-1	5.599E-1
    6.630E-1	5.435E-1
    6.666E-1	5.282E-1
    6.702E-1	5.137E-1
    6.738E-1	5.001E-1
    6.775E-1	4.871E-1
    6.812E-1	4.748E-1
    6.850E-1	4.630E-1
    6.888E-1	4.517E-1
    6.926E-1	4.409E-1
    6.965E-1	4.305E-1
    7.005E-1	4.204E-1
    7.045E-1	4.108E-1
    7.085E-1	4.015E-1
    7.126E-1	3.925E-1
    7.167E-1	3.838E-1
    7.208E-1	3.755E-1
    7.251E-1	3.674E-1
    7.293E-1	3.597E-1
    7.336E-1	3.522E-1
    7.380E-1	3.449E-1
    7.424E-1	3.380E-1
    7.469E-1	3.313E-1
    7.514E-1	3.248E-1
    7.560E-1	3.186E-1
    7.606E-1	3.126E-1
    7.653E-1	3.068E-1
    7.701E-1	3.012E-1
    7.749E-1	2.958E-1
    7.798E-1	2.906E-1
    7.847E-1	2.855E-1
    7.897E-1	2.807E-1
    7.948E-1	2.760E-1
    7.999E-1	2.714E-1
    8.051E-1	2.670E-1
    8.104E-1	2.628E-1
    8.157E-1	2.587E-1
    8.211E-1	2.547E-1
    8.266E-1	2.508E-1
    8.321E-1	2.470E-1
    8.377E-1	2.433E-1
    8.434E-1	2.397E-1
    8.492E-1	2.362E-1
    8.551E-1	2.328E-1
    8.610E-1	2.295E-1
    8.670E-1	2.262E-1
    8.731E-1	2.230E-1
    8.793E-1	2.199E-1
    8.856E-1	2.168E-1
    8.920E-1	2.137E-1
    8.984E-1	2.108E-1
    9.050E-1	2.078E-1
    9.116E-1	2.049E-1
    9.184E-1	2.020E-1
    9.253E-1	1.991E-1
    9.322E-1	1.963E-1
    9.393E-1	1.935E-1
    9.464E-1	1.907E-1
    9.537E-1	1.880E-1
    9.611E-1	1.852E-1
    9.686E-1	1.825E-1
    9.763E-1	1.797E-1
    9.840E-1	1.770E-1
    9.919E-1	1.743E-1
    9.999E-1	1.716E-1
    1.008E+0	1.689E-1
    1.016E+0	1.661E-1
    1.025E+0	1.634E-1
    1.033E+0	1.607E-1
    1.042E+0	1.580E-1
    1.051E+0	1.553E-1
    1.060E+0	1.525E-1
    1.069E+0	1.498E-1
    1.078E+0	1.470E-1
    1.088E+0	1.442E-1
    1.097E+0	1.412E-1
    1.107E+0	1.378E-1
    1.117E+0	1.338E-1
    1.127E+0	1.292E-1
    1.137E+0	1.245E-1
    1.148E+0	1.204E-1
    1.159E+0	1.171E-1
    1.170E+0	1.143E-1
    1.181E+0	1.118E-1
    1.192E+0	1.093E-1
    1.204E+0	1.068E-1
    1.216E+0	1.044E-1
    1.228E+0	1.019E-1
    1.240E+0	9.945E-2
    1.252E+0	9.698E-2
    1.265E+0	9.449E-2
    1.278E+0	9.200E-2
    1.292E+0	8.950E-2
    1.305E+0	8.699E-2
    1.319E+0	8.447E-2
    1.333E+0	8.194E-2
    1.348E+0	7.941E-2
    1.362E+0	7.687E-2
    1.378E+0	7.433E-2
    1.393E+0	7.180E-2
    1.409E+0	6.928E-2
    1.425E+0	6.677E-2
    1.442E+0	6.429E-2
    1.459E+0	6.183E-2
    1.476E+0	5.941E-2
    1.494E+0	5.700E-2
    1.512E+0	5.406E-2
    1.531E+0	4.564E-2
    1.550E+0	2.519E-2
    1.569E+0	7.651E-3
    1.590E+0	2.758E-3
    1.610E+0	1.915E-3
    1.631E+0	1.465E-3
    1.653E+0	1.106E-3
    1.675E+0	8.197E-4
    1.698E+0	5.967E-4
    1.722E+0	4.261E-4
    1.746E+0	2.984E-4
    1.771E+0	2.048E-4
    1.797E+0	1.376E-4
    1.823E+0	9.060E-5
    1.851E+0	5.835E-5
    1.879E+0	3.677E-5
    1.907E+0	2.267E-5
    1.937E+0	1.366E-5
    1.968E+0	8.048E-6
    2.000E+0	4.634E-6
    2.033E+0	2.607E-6
    2.066E+0	1.433E-6
    2.101E+0	7.691E-7
    2.138E+0	4.043E-7
    2.175E+0	2.065E-7
    2.214E+0	0.000E+0
    2.254E+0	0.000E+0
    2.296E+0	0.000E+0
    2.339E+0	0.000E+0
    2.384E+0	0.000E+0
    2.431E+0	0.000E+0
    2.480E+0	0.000E+0];
n = interp1(wn(:,1),wn(:,2),wl*1e6,'linear');
k = interp1(wk(:,1),wk(:,2),wl*1e6,'linear');
alpha = 4*pi*k/wl*1e-2;  %(1/cm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n, k, alpha] = get_GaN_n_k(wl)
wn = [1.3051e-01	1.1488e+00
    1.3074e-01	1.1588e+00
    1.3097e-01	1.1688e+00
    1.3119e-01	1.1789e+00
    1.3142e-01	1.1890e+00
    1.3166e-01	1.1992e+00
    1.3189e-01	1.2095e+00
    1.3212e-01	1.2197e+00
    1.3235e-01	1.2300e+00
    1.3259e-01	1.2403e+00
    1.3282e-01	1.2506e+00
    1.3306e-01	1.2609e+00
    1.3329e-01	1.2711e+00
    1.3353e-01	1.2813e+00
    1.3377e-01	1.2914e+00
    1.3401e-01	1.3015e+00
    1.3425e-01	1.3115e+00
    1.3449e-01	1.3215e+00
    1.3473e-01	1.3313e+00
    1.3497e-01	1.3410e+00
    1.3522e-01	1.3506e+00
    1.3546e-01	1.3600e+00
    1.3571e-01	1.3693e+00
    1.3595e-01	1.3784e+00
    1.3620e-01	1.3874e+00
    1.3645e-01	1.3961e+00
    1.3669e-01	1.4047e+00
    1.3694e-01	1.4130e+00
    1.3720e-01	1.4212e+00
    1.3745e-01	1.4290e+00
    1.3770e-01	1.4367e+00
    1.3795e-01	1.4441e+00
    1.3821e-01	1.4512e+00
    1.3846e-01	1.4581e+00
    1.3872e-01	1.4647e+00
    1.3897e-01	1.4710e+00
    1.3923e-01	1.4771e+00
    1.3949e-01	1.4828e+00
    1.3975e-01	1.4883e+00
    1.4001e-01	1.4935e+00
    1.4027e-01	1.4984e+00
    1.4054e-01	1.5030e+00
    1.4080e-01	1.5073e+00
    1.4107e-01	1.5114e+00
    1.4133e-01	1.5152e+00
    1.4160e-01	1.5187e+00
    1.4187e-01	1.5219e+00
    1.4214e-01	1.5249e+00
    1.4241e-01	1.5277e+00
    1.4268e-01	1.5302e+00
    1.4295e-01	1.5325e+00
    1.4322e-01	1.5345e+00
    1.4350e-01	1.5364e+00
    1.4377e-01	1.5381e+00
    1.4405e-01	1.5397e+00
    1.4432e-01	1.5411e+00
    1.4460e-01	1.5424e+00
    1.4488e-01	1.5435e+00
    1.4516e-01	1.5446e+00
    1.4544e-01	1.5457e+00
    1.4573e-01	1.5467e+00
    1.4601e-01	1.5477e+00
    1.4629e-01	1.5488e+00
    1.4658e-01	1.5499e+00
    1.4687e-01	1.5510e+00
    1.4716e-01	1.5523e+00
    1.4745e-01	1.5538e+00
    1.4774e-01	1.5554e+00
    1.4803e-01	1.5573e+00
    1.4832e-01	1.5595e+00
    1.4861e-01	1.5619e+00
    1.4891e-01	1.5647e+00
    1.4921e-01	1.5679e+00
    1.4950e-01	1.5716e+00
    1.4980e-01	1.5757e+00
    1.5010e-01	1.5804e+00
    1.5040e-01	1.5857e+00
    1.5070e-01	1.5915e+00
    1.5101e-01	1.5981e+00
    1.5131e-01	1.6053e+00
    1.5162e-01	1.6133e+00
    1.5193e-01	1.6220e+00
    1.5223e-01	1.6316e+00
    1.5254e-01	1.6419e+00
    1.5286e-01	1.6530e+00
    1.5317e-01	1.6649e+00
    1.5348e-01	1.6775e+00
    1.5380e-01	1.6908e+00
    1.5411e-01	1.7047e+00
    1.5443e-01	1.7191e+00
    1.5475e-01	1.7339e+00
    1.5507e-01	1.7490e+00
    1.5539e-01	1.7641e+00
    1.5571e-01	1.7793e+00
    1.5604e-01	1.7942e+00
    1.5636e-01	1.8086e+00
    1.5669e-01	1.8226e+00
    1.5702e-01	1.8358e+00
    1.5735e-01	1.8482e+00
    1.5768e-01	1.8597e+00
    1.5801e-01	1.8702e+00
    1.5834e-01	1.8796e+00
    1.5868e-01	1.8880e+00
    1.5901e-01	1.8953e+00
    1.5935e-01	1.9017e+00
    1.5969e-01	1.9072e+00
    1.6003e-01	1.9120e+00
    1.6037e-01	1.9160e+00
    1.6072e-01	1.9195e+00
    1.6106e-01	1.9226e+00
    1.6141e-01	1.9254e+00
    1.6176e-01	1.9280e+00
    1.6211e-01	1.9305e+00
    1.6246e-01	1.9331e+00
    1.6281e-01	1.9358e+00
    1.6317e-01	1.9388e+00
    1.6352e-01	1.9420e+00
    1.6388e-01	1.9457e+00
    1.6424e-01	1.9498e+00
    1.6460e-01	1.9544e+00
    1.6496e-01	1.9596e+00
    1.6532e-01	1.9654e+00
    1.6569e-01	1.9718e+00
    1.6606e-01	1.9788e+00
    1.6642e-01	1.9865e+00
    1.6679e-01	1.9950e+00
    1.6717e-01	2.0041e+00
    1.6754e-01	2.0139e+00
    1.6791e-01	2.0244e+00
    1.6829e-01	2.0357e+00
    1.6867e-01	2.0476e+00
    1.6905e-01	2.0603e+00
    1.6943e-01	2.0737e+00
    1.6982e-01	2.0878e+00
    1.7020e-01	2.1026e+00
    1.7059e-01	2.1180e+00
    1.7098e-01	2.1341e+00
    1.7137e-01	2.1508e+00
    1.7176e-01	2.1681e+00
    1.7215e-01	2.1861e+00
    1.7255e-01	2.2046e+00
    1.7295e-01	2.2236e+00
    1.7335e-01	2.2431e+00
    1.7375e-01	2.2631e+00
    1.7415e-01	2.2836e+00
    1.7456e-01	2.3044e+00
    1.7497e-01	2.3257e+00
    1.7538e-01	2.3472e+00
    1.7579e-01	2.3691e+00
    1.7620e-01	2.3911e+00
    1.7661e-01	2.4134e+00
    1.7703e-01	2.4358e+00
    1.7745e-01	2.4583e+00
    1.7787e-01	2.4809e+00
    1.7829e-01	2.5034e+00
    1.7872e-01	2.5259e+00
    1.7915e-01	2.5483e+00
    1.7958e-01	2.5705e+00
    1.8001e-01	2.5925e+00
    1.8044e-01	2.6142e+00
    1.8087e-01	2.6357e+00
    1.8131e-01	2.6567e+00
    1.8175e-01	2.6774e+00
    1.8219e-01	2.6976e+00
    1.8264e-01	2.7173e+00
    1.8308e-01	2.7365e+00
    1.8353e-01	2.7552e+00
    1.8398e-01	2.7732e+00
    1.8443e-01	2.7905e+00
    1.8489e-01	2.8072e+00
    1.8534e-01	2.8233e+00
    1.8580e-01	2.8386e+00
    1.8627e-01	2.8532e+00
    1.8673e-01	2.8671e+00
    1.8720e-01	2.8802e+00
    1.8766e-01	2.8925e+00
    1.8813e-01	2.9041e+00
    1.8861e-01	2.9150e+00
    1.8908e-01	2.9251e+00
    1.8956e-01	2.9344e+00
    1.9004e-01	2.9430e+00
    1.9052e-01	2.9509e+00
    1.9101e-01	2.9581e+00
    1.9150e-01	2.9646e+00
    1.9199e-01	2.9704e+00
    1.9248e-01	2.9755e+00
    1.9298e-01	2.9800e+00
    1.9347e-01	2.9838e+00
    1.9397e-01	2.9871e+00
    1.9448e-01	2.9898e+00
    1.9498e-01	2.9919e+00
    1.9549e-01	2.9935e+00
    1.9600e-01	2.9946e+00
    1.9652e-01	2.9952e+00
    1.9703e-01	2.9954e+00
    1.9755e-01	2.9951e+00
    1.9807e-01	2.9944e+00
    1.9860e-01	2.9933e+00
    1.9912e-01	2.9919e+00
    1.9966e-01	2.9900e+00
    2.0019e-01	2.9879e+00
    2.0072e-01	2.9855e+00
    2.0126e-01	2.9827e+00
    2.0180e-01	2.9798e+00
    2.0235e-01	2.9765e+00
    2.0290e-01	2.9730e+00
    2.0345e-01	2.9694e+00
    2.0400e-01	2.9655e+00
    2.0456e-01	2.9614e+00
    2.0512e-01	2.9572e+00
    2.0568e-01	2.9528e+00
    2.0624e-01	2.9483e+00
    2.0681e-01	2.9437e+00
    2.0738e-01	2.9389e+00
    2.0796e-01	2.9340e+00
    2.0854e-01	2.9291e+00
    2.0912e-01	2.9241e+00
    2.0970e-01	2.9190e+00
    2.1029e-01	2.9138e+00
    2.1088e-01	2.9086e+00
    2.1148e-01	2.9033e+00
    2.1208e-01	2.8980e+00
    2.1268e-01	2.8927e+00
    2.1328e-01	2.8873e+00
    2.1389e-01	2.8819e+00
    2.1450e-01	2.8766e+00
    2.1512e-01	2.8712e+00
    2.1574e-01	2.8658e+00
    2.1636e-01	2.8604e+00
    2.1699e-01	2.8550e+00
    2.1762e-01	2.8496e+00
    2.1825e-01	2.8443e+00
    2.1889e-01	2.8390e+00
    2.1953e-01	2.8337e+00
    2.2017e-01	2.8284e+00
    2.2082e-01	2.8232e+00
    2.2147e-01	2.8179e+00
    2.2213e-01	2.8128e+00
    2.2279e-01	2.8076e+00
    2.2345e-01	2.8025e+00
    2.2412e-01	2.7975e+00
    2.2479e-01	2.7925e+00
    2.2547e-01	2.7875e+00
    2.2615e-01	2.7826e+00
    2.2683e-01	2.7777e+00
    2.2752e-01	2.7729e+00
    2.2821e-01	2.7682e+00
    2.2891e-01	2.7635e+00
    2.2961e-01	2.7588e+00
    2.3031e-01	2.7542e+00
    2.3102e-01	2.7497e+00
    2.3174e-01	2.7452e+00
    2.3246e-01	2.7408e+00
    2.3318e-01	2.7364e+00
    2.3391e-01	2.7321e+00
    2.3464e-01	2.7278e+00
    2.3537e-01	2.7236e+00
    2.3612e-01	2.7195e+00
    2.3686e-01	2.7154e+00
    2.3761e-01	2.7114e+00
    2.3837e-01	2.7074e+00
    2.3913e-01	2.7035e+00
    2.3989e-01	2.6997e+00
    2.4066e-01	2.6959e+00
    2.4144e-01	2.6922e+00
    2.4222e-01	2.6886e+00
    2.4300e-01	2.6850e+00
    2.4379e-01	2.6814e+00
    2.4459e-01	2.6780e+00
    2.4539e-01	2.6746e+00
    2.4619e-01	2.6712e+00
    2.4700e-01	2.6679e+00
    2.4782e-01	2.6647e+00
    2.4864e-01	2.6616e+00
    2.4947e-01	2.6585e+00
    2.5030e-01	2.6554e+00
    2.5114e-01	2.6524e+00
    2.5198e-01	2.6495e+00
    2.5283e-01	2.6466e+00
    2.5369e-01	2.6438e+00
    2.5455e-01	2.6411e+00
    2.5542e-01	2.6384e+00
    2.5629e-01	2.6358e+00
    2.5717e-01	2.6332e+00
    2.5805e-01	2.6307e+00
    2.5894e-01	2.6283e+00
    2.5984e-01	2.6259e+00
    2.6074e-01	2.6236e+00
    2.6165e-01	2.6213e+00
    2.6257e-01	2.6191e+00
    2.6349e-01	2.6169e+00
    2.6442e-01	2.6148e+00
    2.6536e-01	2.6128e+00
    2.6630e-01	2.6108e+00
    2.6725e-01	2.6089e+00
    2.6820e-01	2.6070e+00
    2.6917e-01	2.6052e+00
    2.7014e-01	2.6034e+00
    2.7111e-01	2.6017e+00
    2.7210e-01	2.6000e+00
    2.7309e-01	2.5984e+00
    2.7409e-01	2.5969e+00
    2.7509e-01	2.5954e+00
    2.7610e-01	2.5940e+00
    2.7712e-01	2.5926e+00
    2.7815e-01	2.5913e+00
    2.7919e-01	2.5900e+00
    2.8023e-01	2.5888e+00
    2.8128e-01	2.5876e+00
    2.8234e-01	2.5865e+00
    2.8341e-01	2.5854e+00
    2.8448e-01	2.5844e+00
    2.8557e-01	2.5834e+00
    2.8666e-01	2.5825e+00
    2.8776e-01	2.5816e+00
    2.8887e-01	2.5808e+00
    2.8998e-01	2.5801e+00
    2.9111e-01	2.5793e+00
    2.9224e-01	2.5787e+00
    2.9339e-01	2.5780e+00
    2.9454e-01	2.5774e+00
    2.9570e-01	2.5769e+00
    2.9687e-01	2.5764e+00
    2.9805e-01	2.5760e+00
    2.9924e-01	2.5755e+00
    3.0044e-01	2.5752e+00
    3.0165e-01	2.5749e+00
    3.0287e-01	2.5746e+00
    3.0409e-01	2.5743e+00
    3.0533e-01	2.5741e+00
    3.0658e-01	2.5739e+00
    3.0784e-01	2.5738e+00
    3.0911e-01	2.5737e+00
    3.1039e-01	2.5736e+00
    3.1168e-01	2.5735e+00
    3.1298e-01	2.5735e+00
    3.1429e-01	2.5735e+00
    3.1561e-01	2.5735e+00
    3.1695e-01	2.5735e+00
    3.1829e-01	2.5735e+00
    3.1965e-01	2.5735e+00
    3.2102e-01	2.5736e+00
    3.2240e-01	2.5736e+00
    3.2379e-01	2.5736e+00
    3.2519e-01	2.5736e+00
    3.2661e-01	2.5736e+00
    3.2804e-01	2.5736e+00
    3.2948e-01	2.5734e+00
    3.3093e-01	2.5733e+00
    3.3240e-01	2.5730e+00
    3.3388e-01	2.5727e+00
    3.3537e-01	2.5723e+00
    3.3688e-01	2.5717e+00
    3.3840e-01	2.5710e+00
    3.3993e-01	2.5701e+00
    3.4148e-01	2.5689e+00
    3.4304e-01	2.5675e+00
    3.4462e-01	2.5657e+00
    3.4621e-01	2.5636e+00
    3.4782e-01	2.5610e+00
    3.4944e-01	2.5578e+00
    3.5107e-01	2.5540e+00
    3.5273e-01	2.5494e+00
    3.5439e-01	2.5440e+00
    3.5608e-01	2.5377e+00
    3.5777e-01	2.5307e+00
    3.5949e-01	2.5236e+00
    3.6122e-01	2.5179e+00
    3.6297e-01	2.5172e+00
    3.6473e-01	2.5273e+00
    3.6652e-01	2.5547e+00
    3.6832e-01	2.5975e+00
    3.7014e-01	2.6402e+00
    3.7197e-01	2.6656e+00
    3.7383e-01	2.6710e+00
    3.7570e-01	2.6634e+00
    3.7759e-01	2.6499e+00
    3.7950e-01	2.6347e+00
    3.8143e-01	2.6195e+00
    3.8338e-01	2.6051e+00
    3.8535e-01	2.5918e+00
    3.8734e-01	2.5795e+00
    3.8935e-01	2.5682e+00
    3.9138e-01	2.5578e+00
    3.9344e-01	2.5481e+00
    3.9551e-01	2.5391e+00
    3.9761e-01	2.5307e+00
    3.9973e-01	2.5228e+00
    4.0187e-01	2.5153e+00
    4.0404e-01	2.5082e+00
    4.0623e-01	2.5015e+00
    4.0844e-01	2.4951e+00
    4.1067e-01	2.4890e+00
    4.1294e-01	2.4831e+00
    4.1522e-01	2.4775e+00
    4.1753e-01	2.4721e+00
    4.1987e-01	2.4668e+00
    4.2224e-01	2.4618e+00
    4.2463e-01	2.4569e+00
    4.2704e-01	2.4522e+00
    4.2949e-01	2.4476e+00
    4.3196e-01	2.4431e+00
    4.3447e-01	2.4388e+00
    4.3700e-01	2.4346e+00
    4.3956e-01	2.4304e+00
    4.4215e-01	2.4264e+00
    4.4477e-01	2.4225e+00
    4.4743e-01	2.4187e+00
    4.5011e-01	2.4149e+00
    4.5283e-01	2.4113e+00
    4.5558e-01	2.4077e+00
    4.5837e-01	2.4041e+00
    4.6119e-01	2.4007e+00
    4.6404e-01	2.3973e+00
    4.6693e-01	2.3940e+00
    4.6985e-01	2.3907e+00
    4.7282e-01	2.3875e+00
    4.7582e-01	2.3844e+00
    4.7886e-01	2.3813e+00
    4.8193e-01	2.3782e+00
    4.8505e-01	2.3752e+00
    4.8821e-01	2.3723e+00
    4.9141e-01	2.3694e+00
    4.9465e-01	2.3665e+00
    4.9793e-01	2.3637e+00
    5.0126e-01	2.3609e+00
    5.0463e-01	2.3582e+00
    5.0805e-01	2.3555e+00
    5.1152e-01	2.3528e+00
    5.1503e-01	2.3502e+00
    5.1859e-01	2.3476e+00
    5.2220e-01	2.3450e+00
    5.2587e-01	2.3425e+00
    5.2958e-01	2.3400e+00
    5.3335e-01	2.3376e+00
    5.3717e-01	2.3351e+00
    5.4104e-01	2.3327e+00
    5.4497e-01	2.3304e+00
    5.4896e-01	2.3280e+00
    5.5301e-01	2.3257e+00
    5.5712e-01	2.3234e+00
    5.6129e-01	2.3211e+00
    5.6552e-01	2.3189e+00
    5.6982e-01	2.3167e+00
    5.7418e-01	2.3145e+00
    5.7861e-01	2.3123e+00
    5.8311e-01	2.3102e+00
    5.8768e-01	2.3081e+00
    5.9232e-01	2.3060e+00
    5.9704e-01	2.3039e+00
    6.0183e-01	2.3019e+00
    6.0670e-01	2.2998e+00
    6.1165e-01	2.2978e+00
    6.1668e-01	2.2958e+00
    6.2179e-01	2.2939e+00
    6.2699e-01	2.2919e+00
    6.3228e-01	2.2900e+00
    6.3765e-01	2.2881e+00
    6.4312e-01	2.2862e+00
    6.4868e-01	2.2843e+00
    6.5434e-01	2.2825e+00
    6.6010e-01	2.2806e+00
    6.6596e-01	2.2788e+00
    6.7193e-01	2.2770e+00
    6.7801e-01	2.2752e+00
    6.8419e-01	2.2735e+00
    6.9049e-01	2.2717e+00
    6.9691e-01	2.2700e+00
    7.0345e-01	2.2682e+00
    7.1011e-01	2.2665e+00
    7.1690e-01	2.2649e+00
    7.2382e-01	2.2632e+00
    7.3087e-01	2.2615e+00
    7.3806e-01	2.2599e+00
    7.4540e-01	2.2582e+00
    7.5288e-01	2.2566e+00
    7.6052e-01	2.2550e+00
    7.6831e-01	2.2534e+00
    7.7626e-01	2.2519e+00
    7.8438e-01	2.2503e+00
    7.9267e-01	2.2488e+00
    8.0114e-01	2.2472e+00
    8.0979e-01	2.2457e+00
    8.1863e-01	2.2442e+00
    8.2767e-01	2.2427e+00
    8.3690e-01	2.2412e+00
    8.4635e-01	2.2397e+00
    8.5601e-01	2.2383e+00
    8.6589e-01	2.2368e+00
    8.7601e-01	2.2354e+00
    8.8636e-01	2.2340e+00
    8.9696e-01	2.2326e+00
    9.0782e-01	2.2312e+00
    9.1895e-01	2.2298e+00
    9.3035e-01	2.2284e+00
    9.4203e-01	2.2270e+00
    9.5402e-01	2.2257e+00
    9.6631e-01	2.2244e+00
    9.7893e-01	2.2230e+00
    9.9187e-01	2.2217e+00];

wk = [1.3051e-01	1.6661e+00
    1.3074e-01	1.6697e+00
    1.3097e-01	1.6732e+00
    1.3119e-01	1.6765e+00
    1.3142e-01	1.6795e+00
    1.3166e-01	1.6824e+00
    1.3189e-01	1.6850e+00
    1.3212e-01	1.6874e+00
    1.3235e-01	1.6895e+00
    1.3259e-01	1.6915e+00
    1.3282e-01	1.6932e+00
    1.3306e-01	1.6946e+00
    1.3329e-01	1.6958e+00
    1.3353e-01	1.6968e+00
    1.3377e-01	1.6976e+00
    1.3401e-01	1.6981e+00
    1.3425e-01	1.6984e+00
    1.3449e-01	1.6984e+00
    1.3473e-01	1.6982e+00
    1.3497e-01	1.6978e+00
    1.3522e-01	1.6971e+00
    1.3546e-01	1.6963e+00
    1.3571e-01	1.6952e+00
    1.3595e-01	1.6939e+00
    1.3620e-01	1.6924e+00
    1.3645e-01	1.6908e+00
    1.3669e-01	1.6890e+00
    1.3694e-01	1.6870e+00
    1.3720e-01	1.6848e+00
    1.3745e-01	1.6825e+00
    1.3770e-01	1.6801e+00
    1.3795e-01	1.6776e+00
    1.3821e-01	1.6750e+00
    1.3846e-01	1.6722e+00
    1.3872e-01	1.6695e+00
    1.3897e-01	1.6666e+00
    1.3923e-01	1.6638e+00
    1.3949e-01	1.6609e+00
    1.3975e-01	1.6580e+00
    1.4001e-01	1.6552e+00
    1.4027e-01	1.6524e+00
    1.4054e-01	1.6497e+00
    1.4080e-01	1.6470e+00
    1.4107e-01	1.6444e+00
    1.4133e-01	1.6420e+00
    1.4160e-01	1.6397e+00
    1.4187e-01	1.6375e+00
    1.4214e-01	1.6355e+00
    1.4241e-01	1.6337e+00
    1.4268e-01	1.6322e+00
    1.4295e-01	1.6308e+00
    1.4322e-01	1.6297e+00
    1.4350e-01	1.6289e+00
    1.4377e-01	1.6284e+00
    1.4405e-01	1.6281e+00
    1.4432e-01	1.6282e+00
    1.4460e-01	1.6286e+00
    1.4488e-01	1.6293e+00
    1.4516e-01	1.6305e+00
    1.4544e-01	1.6319e+00
    1.4573e-01	1.6338e+00
    1.4601e-01	1.6360e+00
    1.4629e-01	1.6387e+00
    1.4658e-01	1.6418e+00
    1.4687e-01	1.6452e+00
    1.4716e-01	1.6491e+00
    1.4745e-01	1.6534e+00
    1.4774e-01	1.6582e+00
    1.4803e-01	1.6633e+00
    1.4832e-01	1.6689e+00
    1.4861e-01	1.6749e+00
    1.4891e-01	1.6813e+00
    1.4921e-01	1.6880e+00
    1.4950e-01	1.6951e+00
    1.4980e-01	1.7025e+00
    1.5010e-01	1.7102e+00
    1.5040e-01	1.7182e+00
    1.5070e-01	1.7264e+00
    1.5101e-01	1.7347e+00
    1.5131e-01	1.7431e+00
    1.5162e-01	1.7516e+00
    1.5193e-01	1.7600e+00
    1.5223e-01	1.7682e+00
    1.5254e-01	1.7762e+00
    1.5286e-01	1.7839e+00
    1.5317e-01	1.7911e+00
    1.5348e-01	1.7978e+00
    1.5380e-01	1.8038e+00
    1.5411e-01	1.8091e+00
    1.5443e-01	1.8135e+00
    1.5475e-01	1.8169e+00
    1.5507e-01	1.8193e+00
    1.5539e-01	1.8207e+00
    1.5571e-01	1.8210e+00
    1.5604e-01	1.8202e+00
    1.5636e-01	1.8183e+00
    1.5669e-01	1.8155e+00
    1.5702e-01	1.8118e+00
    1.5735e-01	1.8075e+00
    1.5768e-01	1.8025e+00
    1.5801e-01	1.7971e+00
    1.5834e-01	1.7915e+00
    1.5868e-01	1.7859e+00
    1.5901e-01	1.7804e+00
    1.5935e-01	1.7751e+00
    1.5969e-01	1.7703e+00
    1.6003e-01	1.7661e+00
    1.6037e-01	1.7624e+00
    1.6072e-01	1.7595e+00
    1.6106e-01	1.7574e+00
    1.6141e-01	1.7561e+00
    1.6176e-01	1.7556e+00
    1.6211e-01	1.7559e+00
    1.6246e-01	1.7571e+00
    1.6281e-01	1.7590e+00
    1.6317e-01	1.7616e+00
    1.6352e-01	1.7650e+00
    1.6388e-01	1.7690e+00
    1.6424e-01	1.7736e+00
    1.6460e-01	1.7788e+00
    1.6496e-01	1.7844e+00
    1.6532e-01	1.7904e+00
    1.6569e-01	1.7969e+00
    1.6606e-01	1.8036e+00
    1.6642e-01	1.8105e+00
    1.6679e-01	1.8177e+00
    1.6717e-01	1.8250e+00
    1.6754e-01	1.8323e+00
    1.6791e-01	1.8397e+00
    1.6829e-01	1.8471e+00
    1.6867e-01	1.8543e+00
    1.6905e-01	1.8615e+00
    1.6943e-01	1.8685e+00
    1.6982e-01	1.8752e+00
    1.7020e-01	1.8817e+00
    1.7059e-01	1.8878e+00
    1.7098e-01	1.8936e+00
    1.7137e-01	1.8990e+00
    1.7176e-01	1.9039e+00
    1.7215e-01	1.9084e+00
    1.7255e-01	1.9123e+00
    1.7295e-01	1.9156e+00
    1.7335e-01	1.9184e+00
    1.7375e-01	1.9205e+00
    1.7415e-01	1.9219e+00
    1.7456e-01	1.9227e+00
    1.7497e-01	1.9227e+00
    1.7538e-01	1.9219e+00
    1.7579e-01	1.9204e+00
    1.7620e-01	1.9181e+00
    1.7661e-01	1.9149e+00
    1.7703e-01	1.9109e+00
    1.7745e-01	1.9061e+00
    1.7787e-01	1.9004e+00
    1.7829e-01	1.8938e+00
    1.7872e-01	1.8863e+00
    1.7915e-01	1.8780e+00
    1.7958e-01	1.8688e+00
    1.8001e-01	1.8588e+00
    1.8044e-01	1.8479e+00
    1.8087e-01	1.8362e+00
    1.8131e-01	1.8237e+00
    1.8175e-01	1.8104e+00
    1.8219e-01	1.7964e+00
    1.8264e-01	1.7816e+00
    1.8308e-01	1.7661e+00
    1.8353e-01	1.7500e+00
    1.8398e-01	1.7332e+00
    1.8443e-01	1.7159e+00
    1.8489e-01	1.6980e+00
    1.8534e-01	1.6797e+00
    1.8580e-01	1.6608e+00
    1.8627e-01	1.6416e+00
    1.8673e-01	1.6220e+00
    1.8720e-01	1.6021e+00
    1.8766e-01	1.5819e+00
    1.8813e-01	1.5615e+00
    1.8861e-01	1.5409e+00
    1.8908e-01	1.5202e+00
    1.8956e-01	1.4994e+00
    1.9004e-01	1.4785e+00
    1.9052e-01	1.4575e+00
    1.9101e-01	1.4366e+00
    1.9150e-01	1.4158e+00
    1.9199e-01	1.3950e+00
    1.9248e-01	1.3743e+00
    1.9298e-01	1.3538e+00
    1.9347e-01	1.3335e+00
    1.9397e-01	1.3133e+00
    1.9448e-01	1.2933e+00
    1.9498e-01	1.2736e+00
    1.9549e-01	1.2541e+00
    1.9600e-01	1.2349e+00
    1.9652e-01	1.2160e+00
    1.9703e-01	1.1974e+00
    1.9755e-01	1.1790e+00
    1.9807e-01	1.1610e+00
    1.9860e-01	1.1433e+00
    1.9912e-01	1.1259e+00
    1.9966e-01	1.1089e+00
    2.0019e-01	1.0921e+00
    2.0072e-01	1.0758e+00
    2.0126e-01	1.0597e+00
    2.0180e-01	1.0440e+00
    2.0235e-01	1.0287e+00
    2.0290e-01	1.0136e+00
    2.0345e-01	9.9896e-01
    2.0400e-01	9.8461e-01
    2.0456e-01	9.7058e-01
    2.0512e-01	9.5689e-01
    2.0568e-01	9.4351e-01
    2.0624e-01	9.3045e-01
    2.0681e-01	9.1771e-01
    2.0738e-01	9.0528e-01
    2.0796e-01	8.9315e-01
    2.0854e-01	8.8131e-01
    2.0912e-01	8.6977e-01
    2.0970e-01	8.5852e-01
    2.1029e-01	8.4754e-01
    2.1088e-01	8.3685e-01
    2.1148e-01	8.2642e-01
    2.1208e-01	8.1625e-01
    2.1268e-01	8.0634e-01
    2.1328e-01	7.9668e-01
    2.1389e-01	7.8727e-01
    2.1450e-01	7.7809e-01
    2.1512e-01	7.6915e-01
    2.1574e-01	7.6043e-01
    2.1636e-01	7.5194e-01
    2.1699e-01	7.4366e-01
    2.1762e-01	7.3559e-01
    2.1825e-01	7.2773e-01
    2.1889e-01	7.2006e-01
    2.1953e-01	7.1259e-01
    2.2017e-01	7.0531e-01
    2.2082e-01	6.9821e-01
    2.2147e-01	6.9129e-01
    2.2213e-01	6.8454e-01
    2.2279e-01	6.7796e-01
    2.2345e-01	6.7154e-01
    2.2412e-01	6.6528e-01
    2.2479e-01	6.5918e-01
    2.2547e-01	6.5323e-01
    2.2615e-01	6.4742e-01
    2.2683e-01	6.4176e-01
    2.2752e-01	6.3623e-01
    2.2821e-01	6.3083e-01
    2.2891e-01	6.2557e-01
    2.2961e-01	6.2043e-01
    2.3031e-01	6.1542e-01
    2.3102e-01	6.1052e-01
    2.3174e-01	6.0574e-01
    2.3246e-01	6.0106e-01
    2.3318e-01	5.9650e-01
    2.3391e-01	5.9205e-01
    2.3464e-01	5.8769e-01
    2.3537e-01	5.8343e-01
    2.3612e-01	5.7927e-01
    2.3686e-01	5.7521e-01
    2.3761e-01	5.7123e-01
    2.3837e-01	5.6734e-01
    2.3913e-01	5.6353e-01
    2.3989e-01	5.5981e-01
    2.4066e-01	5.5616e-01
    2.4144e-01	5.5260e-01
    2.4222e-01	5.4910e-01
    2.4300e-01	5.4568e-01
    2.4379e-01	5.4233e-01
    2.4459e-01	5.3904e-01
    2.4539e-01	5.3582e-01
    2.4619e-01	5.3267e-01
    2.4700e-01	5.2957e-01
    2.4782e-01	5.2653e-01
    2.4864e-01	5.2355e-01
    2.4947e-01	5.2062e-01
    2.5030e-01	5.1775e-01
    2.5114e-01	5.1492e-01
    2.5198e-01	5.1215e-01
    2.5283e-01	5.0942e-01
    2.5369e-01	5.0673e-01
    2.5455e-01	5.0409e-01
    2.5542e-01	5.0149e-01
    2.5629e-01	4.9893e-01
    2.5717e-01	4.9641e-01
    2.5805e-01	4.9392e-01
    2.5894e-01	4.9147e-01
    2.5984e-01	4.8905e-01
    2.6074e-01	4.8666e-01
    2.6165e-01	4.8430e-01
    2.6257e-01	4.8197e-01
    2.6349e-01	4.7967e-01
    2.6442e-01	4.7739e-01
    2.6536e-01	4.7513e-01
    2.6630e-01	4.7289e-01
    2.6725e-01	4.7068e-01
    2.6820e-01	4.6848e-01
    2.6917e-01	4.6630e-01
    2.7014e-01	4.6413e-01
    2.7111e-01	4.6198e-01
    2.7210e-01	4.5984e-01
    2.7309e-01	4.5771e-01
    2.7409e-01	4.5560e-01
    2.7509e-01	4.5348e-01
    2.7610e-01	4.5138e-01
    2.7712e-01	4.4928e-01
    2.7815e-01	4.4718e-01
    2.7919e-01	4.4509e-01
    2.8023e-01	4.4299e-01
    2.8128e-01	4.4090e-01
    2.8234e-01	4.3880e-01
    2.8341e-01	4.3669e-01
    2.8448e-01	4.3458e-01
    2.8557e-01	4.3246e-01
    2.8666e-01	4.3033e-01
    2.8776e-01	4.2819e-01
    2.8887e-01	4.2604e-01
    2.8998e-01	4.2387e-01
    2.9111e-01	4.2169e-01
    2.9224e-01	4.1948e-01
    2.9339e-01	4.1726e-01
    2.9454e-01	4.1502e-01
    2.9570e-01	4.1275e-01
    2.9687e-01	4.1045e-01
    2.9805e-01	4.0813e-01
    2.9924e-01	4.0577e-01
    3.0044e-01	4.0338e-01
    3.0165e-01	4.0096e-01
    3.0287e-01	3.9850e-01
    3.0409e-01	3.9601e-01
    3.0533e-01	3.9347e-01
    3.0658e-01	3.9088e-01
    3.0784e-01	3.8825e-01
    3.0911e-01	3.8557e-01
    3.1039e-01	3.8284e-01
    3.1168e-01	3.8006e-01
    3.1298e-01	3.7721e-01
    3.1429e-01	3.7431e-01
    3.1561e-01	3.7134e-01
    3.1695e-01	3.6831e-01
    3.1829e-01	3.6521e-01
    3.1965e-01	3.6203e-01
    3.2102e-01	3.5878e-01
    3.2240e-01	3.5545e-01
    3.2379e-01	3.5204e-01
    3.2519e-01	3.4854e-01
    3.2661e-01	3.4495e-01
    3.2804e-01	3.4128e-01
    3.2948e-01	3.3750e-01
    3.3093e-01	3.3364e-01
    3.3240e-01	3.2968e-01
    3.3388e-01	3.2562e-01
    3.3537e-01	3.2146e-01
    3.3688e-01	3.1721e-01
    3.3840e-01	3.1288e-01
    3.3993e-01	3.0847e-01
    3.4148e-01	3.0400e-01
    3.4304e-01	2.9949e-01
    3.4462e-01	2.9498e-01
    3.4621e-01	2.9052e-01
    3.4782e-01	2.8618e-01
    3.4944e-01	2.8208e-01
    3.5107e-01	2.7838e-01
    3.5273e-01	2.7533e-01
    3.5439e-01	2.7330e-01
    3.5608e-01	2.7287e-01
    3.5777e-01	2.7487e-01
    3.5949e-01	2.8051e-01
    3.6122e-01	2.9135e-01
    3.6297e-01	3.0862e-01
    3.6473e-01	3.3112e-01
    3.6652e-01	3.5132e-01
    3.6832e-01	3.5496e-01
    3.7014e-01	3.3280e-01
    3.7197e-01	2.9297e-01
    3.7383e-01	2.5123e-01
    3.7570e-01	2.1649e-01
    3.7759e-01	1.9015e-01
    3.7950e-01	1.7070e-01
    3.8143e-01	1.5626e-01
    3.8338e-01	1.4536e-01
    3.8535e-01	1.3694e-01
    3.8734e-01	1.3029e-01
    3.8935e-01	1.2491e-01
    3.9138e-01	1.2049e-01
    3.9344e-01	1.1677e-01
    3.9551e-01	1.1360e-01
    3.9761e-01	1.1085e-01
    3.9973e-01	1.0844e-01
    4.0187e-01	1.0629e-01
    4.0404e-01	1.0436e-01
    4.0623e-01	1.0262e-01
    4.0844e-01	1.0102e-01
    4.1067e-01	9.9546e-02
    4.1294e-01	9.8179e-02
    4.1522e-01	9.6903e-02
    4.1753e-01	9.5705e-02
    4.1987e-01	9.4575e-02
    4.2224e-01	9.3505e-02
    4.2463e-01	9.2486e-02
    4.2704e-01	9.1515e-02
    4.2949e-01	9.0584e-02
    4.3196e-01	8.9690e-02
    4.3447e-01	8.8830e-02
    4.3700e-01	8.8000e-02
    4.3956e-01	8.7197e-02
    4.4215e-01	8.6420e-02
    4.4477e-01	8.5665e-02
    4.4743e-01	8.4932e-02
    4.5011e-01	8.4218e-02
    4.5283e-01	8.3522e-02
    4.5558e-01	8.2843e-02
    4.5837e-01	8.2180e-02
    4.6119e-01	8.1531e-02
    4.6404e-01	8.0897e-02
    4.6693e-01	8.0275e-02
    4.6985e-01	7.9666e-02
    4.7282e-01	7.9069e-02
    4.7582e-01	7.8482e-02
    4.7886e-01	7.7906e-02
    4.8193e-01	7.7340e-02
    4.8505e-01	7.6784e-02
    4.8821e-01	7.6236e-02
    4.9141e-01	7.5698e-02
    4.9465e-01	7.5167e-02
    4.9793e-01	7.4645e-02
    5.0126e-01	7.4130e-02
    5.0463e-01	7.3623e-02
    5.0805e-01	7.3123e-02
    5.1152e-01	7.2630e-02
    5.1503e-01	7.2143e-02
    5.1859e-01	7.1663e-02
    5.2220e-01	7.1189e-02
    5.2587e-01	7.0722e-02
    5.2958e-01	7.0260e-02
    5.3335e-01	6.9804e-02
    5.3717e-01	6.9353e-02
    5.4104e-01	6.8908e-02
    5.4497e-01	6.8468e-02
    5.4896e-01	6.8034e-02
    5.5301e-01	6.7604e-02
    5.5712e-01	6.7179e-02
    5.6129e-01	6.6759e-02
    5.6552e-01	6.6343e-02
    5.6982e-01	6.5932e-02
    5.7418e-01	6.5526e-02
    5.7861e-01	6.5124e-02
    5.8311e-01	6.4726e-02
    5.8768e-01	6.4332e-02
    5.9232e-01	6.3942e-02
    5.9704e-01	6.3557e-02
    6.0183e-01	6.3175e-02
    6.0670e-01	6.2797e-02
    6.1165e-01	6.2423e-02
    6.1668e-01	6.2052e-02
    6.2179e-01	6.1686e-02
    6.2699e-01	6.1322e-02
    6.3228e-01	6.0963e-02
    6.3765e-01	6.0606e-02
    6.4312e-01	6.0253e-02
    6.4868e-01	5.9904e-02
    6.5434e-01	5.9557e-02
    6.6010e-01	5.9214e-02
    6.6596e-01	5.8874e-02
    6.7193e-01	5.8537e-02
    6.7801e-01	5.8203e-02
    6.8419e-01	5.7873e-02
    6.9049e-01	5.7545e-02
    6.9691e-01	5.7220e-02
    7.0345e-01	5.6898e-02
    7.1011e-01	5.6578e-02
    7.1690e-01	5.6262e-02
    7.2382e-01	5.5948e-02
    7.3087e-01	5.5637e-02
    7.3806e-01	5.5329e-02
    7.4540e-01	5.5023e-02
    7.5288e-01	5.4720e-02
    7.6052e-01	5.4420e-02
    7.6831e-01	5.4121e-02
    7.7626e-01	5.3826e-02
    7.8438e-01	5.3533e-02
    7.9267e-01	5.3242e-02
    8.0114e-01	5.2954e-02
    8.0979e-01	5.2668e-02
    8.1863e-01	5.2384e-02
    8.2767e-01	5.2103e-02
    8.3690e-01	5.1824e-02
    8.4635e-01	5.1547e-02
    8.5601e-01	5.1273e-02
    8.6589e-01	5.1000e-02
    8.7601e-01	5.0730e-02
    8.8636e-01	5.0462e-02
    8.9696e-01	5.0196e-02
    9.0782e-01	4.9932e-02
    9.1895e-01	4.9670e-02
    9.3035e-01	4.9410e-02
    9.4203e-01	4.9152e-02
    9.5402e-01	4.8896e-02
    9.6631e-01	4.8642e-02
    9.7893e-01	4.8390e-02
    9.9187e-01	4.8140e-02];
if wl*1e6<0.99
    n = interp1(wn(:,1),wn(:,2),wl*1e6,'linear');
    k = interp1(wk(:,1),wk(:,2),wl*1e6,'linear');
else
    n = wn(end,2);
    k = wk(end,2);
end
alpha = 4*pi*k/wl*1e-2;  %(1/cm)

