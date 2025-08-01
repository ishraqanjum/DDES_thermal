% v3.0
function [Sol] = Obtain_Static_Solution(Design, params,P00)
Sol = [];
Sol.converged = 0;
% for parallel computing, P00 is needed
% for regular computation, "Design and params" are enough
if nargin < 3
    P00 = params.P0_factor;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the default values and formats
set(0,'defaultlinelinewidth',1.8);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultTextFontSize',18);
format short e;
warning('off')
mkdir('Results');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the constants
constant;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bias = params.Bias;
Impact = params.Impact; % impact ionization or not
Thermionic_e = params.Thermionic_e;
Thermionic_h = params.Thermionic_h;
FranzKeldysh = params.FranzKeldysh;
%%% Excitation
wavelength = params.wavelength;         %wavelength of the light
params.P0_factor = P00;
P0_factor = params.P0_factor;
%%%
R_load = params.R_load;
% Set R_load to a small value of 0 entered, to make the code more stable
if R_load == 0
    R_load = 1e-3;
end
Temperature = params.Temperature;            % normlization temperature
% params Dynamic Calculation
dt0 = params.dt;
Totaltime0 = params.Totaltime;
Dmax = params.Dmax;                 % Parameter to control the mesh quality (nm)
% Device params
Device_Diameter = params.Device_Diameter;
Beam_Diameter = params.Beam_Diameter;
% Solver settings (used in params_NR)
In = params.In;                     % exponent power, m in Eq. 23 of Yue Hu's paper
gamma = params.gamma;               %
% Settings to plot and save the results
filenum = params.filenum;
plot_results = params.plot_results;
Solution_Name  = params.Solution_Name;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% STEADY STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pulse_normalizer = 1;
P0=P0_factor;
Parameters_NR % parameter setting

alpha0 = mean(PD_Str1.alpha(PD_Str1.alpha>100));
NX = mesh1.NX;

Biase=Bias;
Eph=const.h*const.c/wavelength;
A=pi*0.25*min([Beam_Diameter, Device_Diameter])^2;
rs = linspace(0,Device_Diameter,1001);
ff = exp(-(rs/Beam_Diameter).^2);
f_ave = sum(ff*(rs(2)-rs(1)))/Device_Diameter/1.001;
clear rs;
clear ff;

G0= f_ave*P0./(A*Eph)*alpha0*NP.NT^2;
% disp(G0);

structure_loss = ones(N,1);
if params.side ==1  %  n-side illumination
    cell_loss = exp(-PD_Str1.alphas_abs.*NX.*[mesh1.dx; 0]);
    structure_loss = flipud(cumprod(flipud(cell_loss)));
elseif params.side ==2  %  p-side illumination
    cell_loss = exp(-PD_Str1.alphas_abs.*NX.*[mesh1.dx; 0]);
    structure_loss = cumprod(cell_loss);
elseif params.side ==3  % WG mode perp illumination
    structure_loss = structure_loss.*[1; PD_Str1.Eb_re./max(PD_Str1.Eb_re)];
end

GL(1:N,1)=G0.*structure_loss; % per volume
djterm = DeviceDiameter(1:end-1).^2./sum(mesh1.dx)*0.25*pi; % pi*r^2/L = pi*R^2/4/L
converged = 0;
initial_guess_method = params.initial_guess_method;
VT=NP.VT;
while converged ==0
    w=zeros(N,1); % w is the potential unit V
    p1=NV.*exp(-DeltaEA(1)./NP.VT);
    n1=NC.*exp(-DeltaED(end)./NP.VT);
    Js=0;
    switch initial_guess_method
        case -1
            disp('initial guess method: -1')
            p(1:N1,1)=(-p1(1:N1)+sqrt(p1(1:N1).^2+4*p1(1:N1).*abs(Tol(1:N1))))/2;
            n(N2+1:N,1)=(n1(N2+1:N)+sqrt(n1(N2+1:N).^2+4*n1(N2+1:N).*abs(Tol(N2+1:N))))/2;
            n(N1+1:N2,1)=ni(N1+1:N2);% initial condition of the electron density for the second region. It depends on the donor doping density.
            p(N1+1:N2,1)=ni(N1+1:N2);% in the equilibrium condition, p*n=ni^2
            p(N2+1:N,1)=ni(N2+1:N,1).^2./n(N2+1:N,1);
            n(1:N1,1)=ni(1:N1,1).^2./p(1:N1,1);
            
            pb1=p(1);
            nb1=ni(1)^2/pb1;
            nbn=n(end);
            pbn=ni(N)^2/nbn;
                        
            % boundary condition for p, n, and w
            pm=p(2:end-1);
            nm=n(2:end-1);
            w0a = -Eg(1)/VT;
            w0b = Eg(end)/VT;
            
            w(1:N1,1) = w0a;
            w(N2+1:N,1) = w0b ;
            w(N1:N2+1,1) = linspace(w0a,w0b,N2-N1+2);
            wbf=w(1);
            wbn=w(end);
            wb1=w(1);
                        
            dBias = Biase/10;
            for Bias0=0:dBias:Biase
                w(1)=wbf-Bias0/VT;
                % boundary condition for p, n, and w
                wb1=w(1);
                if Bias0>0
                    ny1=N1-round(mesh1.Length/100);
                    ny2=N2+round(mesh1.Length/100);
                    if ny1<3
                        ny1 = N1;
                    end
                    if ny2>N
                        ny2 = N2;
                    end
                    w(2:ny1)=w(2:ny1)-dBias/VT;
                    w(ny1:ny2)=w(ny1)+(w(ny2)-w(ny1))./(ny2-ny1).*([ny1:ny2]-ny1);
                else
                    w(2:N-1)=((w(N)-w(1))/(N-1))*(1:N-2)+w(1);
                end
                x0=[pm;nm;w(2:end-1)];% initial guess of the solution
                fd1=fd(mesh1);
                
                ss1 = Semiconductor_Device(NP,Temperature,Length,mesh1,fd1,Tol,Beta,...
                    un,up,tn,tp,ni,w,Vps,Vns,Ep,Diameter,DeviceDiameter,gamma,n1,p1,...
                    Impact,Thermionic_e,Thermionic_h,FranzKeldysh,Ae,Be,Ah,Bh,Eg,CnAu,CpAu,Br,DeltaEA, DeltaED,...
                    Ne_ref, Np_ref, const.eps_r_silicon,Xi,k_L); % Create the Device and initial condition.
                %                
                DD0=Drift_Diff_N(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,0*GL,wavelength);
                opts = optimoptions('fsolve','Display','off','TolFun',1e-8,'MaxIter',200,'Tolx',1e-6...
                    ,'Diagnostics','off','Jacobian','on','DerivativeCheck','off');
                [x,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x0,opts);
                pm=x(1:N-2);
                nm=x(N-1:2*N-4);
                wm=x(2*N-3:3*N-6);
                w=[wb1;wm;wbn];
            end                 
        case 0
            disp('initial guess method: 0')            
            p(1:N1,1)=(-p1(1:N1)+sqrt(p1(1:N1).^2+8*p1(1:N1).*abs(Tol(1:N1))))/4;
            n(N2+1:N,1)=(-n1(N2+1:N)+sqrt(n1(N2+1:N).^2+16*n1(N2+1:N).*abs(Tol(N2+1:N))))/8;
            n(N1+1:N2,1)=ni(N1+1:N2);% initial condition of the electron density for the second region. It depends on the donor doping density.
            p(N1+1:N2,1)=ni(N1+1:N2);% in the equilibrium condition, p*n=ni^2
            p(N2+1:N,1)=ni(N2+1:N,1).^2./n(N2+1:N,1);
            n(1:N1,1)=ni(1:N1,1).^2./p(1:N1,1);
            
            pb1=p(1);
            nb1=ni(1)^2/pb1;
            nbn=n(end);
            pbn=ni(N)^2/nbn;
                        
            % boundary condition for p, n, and w
            pm=p(2:end-1);
            nm=n(2:end-1);
            w0a = -Eg(1)/VT;
            w0b = Eg(end)/VT;
            
            w(1:N1,1) = w0a;
            w(N2+1:N,1) = w0b ;
            w(N1:N2+1,1) = linspace(w0a,w0b,N2-N1+2);
            wbf=w(1);
            wbn=w(end);
            wb1=w(1);
                        
            dBias = Biase/10;
            for Bias0=0:dBias:Biase
                w(1)=wbf-Bias0/VT;
                % boundary condition for p, n, and w
                wb1=w(1);
                if Bias0>0
                    ny1=N1-round(mesh1.Length/100);
                    ny2=N2+round(mesh1.Length/100);
                    if ny1<3
                        ny1 = N1;
                    end
                    if ny2>N
                        ny2 = N2;
                    end
                    w(2:ny1)=w(2:ny1)-dBias/VT;
                    w(ny1:ny2)=w(ny1)+(w(ny2)-w(ny1))./(ny2-ny1).*([ny1:ny2]-ny1);
                else
                    w(2:N-1)=((w(N)-w(1))/(N-1))*(1:N-2)+w(1);
                end
                x0=[pm;nm;w(2:end-1)];% initial guess of the solution
                fd1=fd(mesh1);
                
                ss1 = Semiconductor_Device(NP,Temperature,Length,mesh1,fd1,Tol,Beta,...
                    un,up,tn,tp,ni,w,Vps,Vns,Ep,Diameter,DeviceDiameter,gamma,n1,p1,...
                    Impact,Thermionic_e,Thermionic_h,FranzKeldysh,Ae,Be,Ah,Bh,Eg,CnAu,CpAu,Br,DeltaEA, DeltaED,...
                    Ne_ref, Np_ref, const.eps_r_silicon,Xi,k_L); % Create the Device and initial condition.
                %
                DD0=Drift_Diff_N(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,0*GL,wavelength);
                opts = optimoptions('fsolve','Display','off','TolFun',1e-8,'MaxIter',200,'Tolx',1e-6...
                    ,'Diagnostics','off','Jacobian','on','DerivativeCheck','off');
                [x,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x0,opts);
                pm=x(1:N-2);
                nm=x(N-1:2*N-4);
                wm=x(2*N-3:3*N-6);
                w=[wb1;wm;wbn];
            end   
        case 1
            disp('initial guess method: 1')            
            p(1:N1,1)=(-p1(1:N1)+sqrt(p1(1:N1).^2+16*p1(1:N1:N1).*abs(Tol(1:N1))))/8;
            n(N2+1:N,1)=(-n1(N2+1:N)+sqrt(n1(N2+1:N).^2+8*n1(N2+1:N).*abs(Tol(N2+1:N))))/4;
            n(N1+1:N2,1)=ni(N1+1:N2);% initial condition of the electron density for the second region. It depends on the donor doping density.
            p(N1+1:N2,1)=ni(N1+1:N2);% in the equilibrium condition, p*n=ni^2
            p(N2+1:N,1)=ni(N2+1:N,1).^2./n(N2+1:N,1);
            n(1:N1,1)=ni(1:N1,1).^2./p(1:N1,1);
            
            pb1=p(1);
            nb1=ni(1)^2/pb1;
            nbn=n(end);
            pbn=ni(N)^2/nbn;            
            
            % boundary condition for p, n, and w
            pm=p(2:end-1);
            nm=n(2:end-1);
            
            w0a = -Eg(1)/NP.VT/2-Xi(1)/NP.VT;
            w0b = -Xi(end)/NP.VT;
            w(1:N1,1) = w0a;
            w(N2+1:N,1) = w0b ;
            w(N1:N2+1,1) = linspace(w0a,w0b,N2-N1+2);
            wbf=w(1);
            wbn=w(end);
            wb1=w(1);
            
            dBias = Biase/10;
            for Bias0=0:dBias:Biase
                w(1)=wbf-Bias0/VT;
                % boundary condition for p, n, and w
                wb1=w(1);
                if Bias0>0
                    ny1=N1-round(mesh1.Length/100);
                    ny2=N2+round(mesh1.Length/100);
                    if ny1<3
                        ny1 = N1;
                    end
                    if ny2>N
                        ny2 = N2;
                    end
                    w(2:ny1)=w(2:ny1)-dBias/VT;
                    w(ny1:ny2)=w(ny1)+(w(ny2)-w(ny1))./(ny2-ny1).*([ny1:ny2]-ny1);
                else
                    w(2:N-1)=((w(N)-w(1))/(N-1))*(1:N-2)+w(1);
                end
                x0=[pm;nm;w(2:end-1)];% initial guess of the solution
                fd1=fd(mesh1);
                
                ss1 = Semiconductor_Device(NP,Temperature,Length,mesh1,fd1,Tol,Beta,...
                    un,up,tn,tp,ni,w,Vps,Vns,Ep,Diameter,DeviceDiameter,gamma,n1,p1,...
                    Impact,Thermionic_e,Thermionic_h,FranzKeldysh,Ae,Be,Ah,Bh,Eg,CnAu,CpAu,Br,DeltaEA, DeltaED,...
                    Ne_ref, Np_ref, const.eps_r_silicon,Xi,k_L); % Create the Device and initial condition.
                %
                DD0=Drift_Diff_N(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,0*GL,wavelength);
                opts = optimoptions('fsolve','Display','off','TolFun',1e-8,'MaxIter',200,'Tolx',1e-6...
                    ,'Diagnostics','off','Jacobian','on','DerivativeCheck','off');
                [x,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x0,opts);
                pm=x(1:N-2);
                nm=x(N-1:2*N-4);
                wm=x(2*N-3:3*N-6);
                w=[wb1;wm;wbn];
            end            
        case 2
            % disp('initial guess method: 2')
            p(1:N1,1)=(-p1(1:N1)+sqrt(p1(1:N1).^2+8*p1(1:N1).*abs(Tol(1:N1))))/4;
            n(N2+1:N,1)=(-n1(N2+1:N)+sqrt(n1(N2+1:N).^2+16*n1(N2+1:N).*abs(Tol(N2+1:N))))/8;
            n(N1+1:N2,1)=ni(N1+1:N2);% initial condition of the electron density for the second region. It depends on the donor doping density.
            p(N1+1:N2,1)=ni(N1+1:N2);% in the equilibrium condition, p*n=ni^2
            p(N2+1:N,1)=ni(N2+1:N,1).^2./n(N2+1:N,1);
            n(1:N1,1)=ni(1:N1,1).^2./p(1:N1,1);
            
            pb1=p(1);
            nb1=ni(1)^2/pb1;
            nbn=n(end);
            pbn=ni(N)^2/nbn;            
            
            % boundary condition for p, n, and w
            pm=p(2:end-1);
            nm=n(2:end-1);

            w0a = -(PD_Str1.delta_Fermis(end)-PD_Str1.delta_Fermis(1))/VT;
            w0b = 0;
            w(1:N1,1) = w0a;
            w(N2+1:N,1) = w0b ;
            w(N1:N2+1,1) = linspace(w0a,w0b,N2-N1+2);
            wbf=w(1);
            wbn=w(end);
            wb1=w(1);
                        
            dBias = Biase/10;
            for Bias0=0:dBias:Biase
                w(1)=wbf-Bias0/VT;
                % boundary condition for p, n, and w
                wb1=w(1);
                if Bias0>0
                    ny1=N1-round(mesh1.Length/100);
                    ny2=N2+round(mesh1.Length/100);
                    if ny1<3
                        ny1 = N1;
                    end
                    if ny2>N
                        ny2 = N2;
                    end
                    w(2:ny1)=w(2:ny1)-dBias/VT;
                    w(ny1:ny2)=w(ny1)+(w(ny2)-w(ny1))./(ny2-ny1).*([ny1:ny2]-ny1);
                else
                    w(2:N-1)=((w(N)-w(1))/(N-1))*(1:N-2)+w(1);
                end
                x0=[pm;nm;w(2:end-1)];% initial guess of the solution
                fd1=fd(mesh1);
                
                ss1 = Semiconductor_Device(NP,Temperature,Length,mesh1,fd1,Tol,Beta,...
                    un,up,tn,tp,ni,w,Vps,Vns,Ep,Diameter,DeviceDiameter,gamma,n1,p1,...
                    Impact,Thermionic_e,Thermionic_h,FranzKeldysh,Ae,Be,Ah,Bh,Eg,CnAu,CpAu,Br,DeltaEA, DeltaED,...
                    Ne_ref, Np_ref, const.eps_r_silicon,Xi,k_L); % Create the Device and initial condition.
                %
                DD0=Drift_Diff_N(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,0*GL,wavelength);
                opts = optimoptions('fsolve','Display','off','TolFun',1e-8,'MaxIter',200,'Tolx',1e-6...
                    ,'Diagnostics','off','Jacobian','on','DerivativeCheck','off');
                [x,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x0,opts);
                pm=x(1:N-2);
                nm=x(N-1:2*N-4);
                wm=x(2*N-3:3*N-6);
                w=[wb1;wm;wbn];
            end
        case 3
            disp('initial guess method: 3')
            
            p(1:N1,1)=(-p1(1:N1)+sqrt(p1(1:N1).^2+16*p1(1:N1:N1).*abs(Tol(1:N1))))/8;
            n(N2+1:N,1)=(-n1(N2+1:N)+sqrt(n1(N2+1:N).^2+8*n1(N2+1:N).*abs(Tol(N2+1:N))))/4;
            n(N1+1:N2,1)=ni(N1+1:N2);% initial condition of the electron density for the second region. It depends on the donor doping density.
            p(N1+1:N2,1)=ni(N1+1:N2);% in the equilibrium condition, p*n=ni^2
            p(N2+1:N,1)=ni(N2+1:N,1).^2./n(N2+1:N,1);
            n(1:N1,1)=ni(1:N1,1).^2./p(1:N1,1);
            
            pb1=p(1);
            nb1=ni(1)^2/pb1;
            nbn=n(end);
            pbn=ni(N)^2/nbn;
                        
            % boundary condition for p, n, and w
            pm=p(2:end-1);
            nm=n(2:end-1);
            
            w0a = -Eg(1)/NP.VT/2-Xi(1)/NP.VT;
            w0b = -Xi(end)/NP.VT;
            w(1:N1,1) = w0a;
            w(N2+1:N,1) = w0b ;
            w(N1:N2+1,1) = linspace(w0a,w0b,N2-N1+2);
            wbf=w(1);
            wbn=w(end);
            wb1=w(1);
            
            dBias = Biase/10;
            for Bias0=0:dBias:Biase
                w(1)=wbf-Bias0/VT;
                % boundary condition for p, n, and w
                wb1=w(1);
                if Bias0>0
                    ny1=N1-round(mesh1.Length/100);
                    ny2=N2+round(mesh1.Length/100);
                    if ny1<3
                        ny1 = N1;
                    end
                    if ny2>N
                        ny2 = N2;
                    end
                    w(2:ny1)=w(2:ny1)-dBias/VT;
                    w(ny1:ny2)=w(ny1)+(w(ny2)-w(ny1))./(ny2-ny1).*([ny1:ny2]-ny1);
                else
                    w(2:N-1)=((w(N)-w(1))/(N-1))*(1:N-2)+w(1);
                end
                x0=[pm;nm;w(2:end-1)];% initial guess of the solution
                fd1=fd(mesh1);
                
                ss1 = Semiconductor_Device(NP,Temperature,Length,mesh1,fd1,Tol,Beta,...
                    un,up,tn,tp,ni,w,Vps,Vns,Ep,Diameter,DeviceDiameter,gamma,n1,p1,...
                    Impact,Thermionic_e,Thermionic_h,FranzKeldysh,Ae,Be,Ah,Bh,Eg,CnAu,CpAu,Br,DeltaEA, DeltaED,...
                    Ne_ref, Np_ref, const.eps_r_silicon,Xi,k_L); % Create the Device and initial condition.
                %
                DD0=Drift_Diff_N(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,0*GL,wavelength);
                opts = optimoptions('fsolve','Display','off','TolFun',1e-8,'MaxIter',200,'Tolx',1e-6...
                    ,'Diagnostics','off','Jacobian','on','DerivativeCheck','off');
                [x,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x0,opts);
                pm=x(1:N-2);
                nm=x(N-1:2*N-4);
                wm=x(2*N-3:3*N-6);
                w=[wb1;wm;wbn];
            end
        case 4
            disp('initial guess method: 4')
            pb1=(-p1(1)+sqrt(p1(1)^2+16*p1(1)*abs(Tol(1))))/8;
            nb1=ni(1)^2/pb1;
            nbn=(-n1(end)+sqrt(n1(end)^2+8*n1(end)*abs(Tol(end))))/4;
            pbn=ni(N)^2/nbn;
            %
            p(1:N1,1)=pb1;  % initial condition of the hole density for the first region. It depends on the acceptor doping density.
            n(1:N1,1)=ni(1:N1,1).^2./p(1:N1,1); % in the equilibrium condition, p*n=ni^2.
            %
            n(N1+1:N2,1)=ni(N1+1:N2);% initial condition of the electron density for the second region. It depends on the donor doping density.
            p(N1+1:N2,1)=ni(N1+1:N2);% in the equilibrium condition, p*n=ni^2
            %
            n(N2+1:N,1)=nbn;% initial condition of the electron density for the third region. It depends on the donor doping density.
            p(N2+1:N,1)=ni(N2+1:N,1).^2./n(N2+1:N,1);% in the equilibrium condition, p*n=ni^2
            % initial and boundary condition for the potential
            w=zeros(N,1); % w is the potential unit V
            w(1:N1) = -cumsum(log(-ni(1:N1)./Tol(1:N1))*NP.VT)/N;
            w(N1+1:N)= cumsum(log(abs(Tol(N1+1:N))./ni(N1+1:N))*NP.VT)/N;
            
            % boundary condition for p, n, and w
            w=w-Xi-Eg/2;
            w=w/VT;
            wbf=w(1);
            wbn=w(end);
            %
            pm=p(2:end-1);
            nm=n(2:end-1);
            Js=0;
            dBias = Biase/10;
            for Bias0=0:dBias:Biase
                w(1)=wbf-(Bias0)/VT;
                % boundary condition for p, n, and w
                wb1=w(1);
                if Bias0>0
                    ny1=N1-round(mesh1.Length/100);
                    ny2=N2+round(mesh1.Length/100);
                    if ny1<3
                        ny1 = N1;
                    end
                    if ny2>N
                        ny2 = N2;
                    end
                    
                    w(2:ny1)=w(2:ny1)-dBias/VT;
                    w(ny1:ny2)=w(ny1)+(w(ny2)-w(ny1))./(ny2-ny1).*([ny1:ny2]-ny1);
                else
                    w(2:N-1)=((w(N)-w(1))/(N-1))*(1:N-2)+w(1);
                end
                x0=[pm;nm;w(2:end-1)];% initial guess of the solution
                fd1=fd(mesh1);
                %
                ss1 = Semiconductor_Device(NP,Temperature,Length,mesh1,fd1,Tol,Beta,...
                    un,up,tn,tp,ni,w,Vps,Vns,Ep,Diameter,DeviceDiameter,gamma,n1,p1,...
                    Impact,Thermionic_e,Thermionic_h,FranzKeldysh,Ae,Be,Ah,Bh,Eg,CnAu,CpAu,Br,DeltaEA, DeltaED,...
                    Ne_ref, Np_ref, const.eps_r_silicon,Xi,k_L); % Create the Device and initial condition.
                %
                DD0=Drift_Diff_N(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,0*GL,wavelength);
                opts = optimoptions('fsolve','Display','off','TolFun',1e-8,'MaxIter',100,'Tolx',1e-6,'MaxIter',100 ...
                    ,'Diagnostics','off','Jacobian','on','DerivativeCheck','off');
                [x,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x0,opts);
                pm=x(1:N-2);
                nm=x(N-1:2*N-4);
                wm=x(2*N-3:3*N-6);
                w=[wb1;wm;wbn];
            end
            
            wb10=w(1);
            w(1)=wbf-Biase/VT;
            w(1)=wbf;
            dw1=wb10-w(1);
            ny1=N1-round(mesh1.Length/100);
            ny2=N2+round(mesh1.Length/100);
            if ny1<3
                ny1 = N1;
            end
            if ny2>N
                ny2 = N2;
            end
            
            w(2:ny1)=w(2:ny1)-dw1;
            w(ny1:ny2)=w(ny1)+(w(ny2)-w(ny1))./(ny2-ny1).*((ny1:ny2)-ny1);
            % boundary condition for p, n, and w
            wb1=w(1);
            w(2:N-1)=((w(N)-w(1))/(N-1)^2)*((1:N-2).^2)+w(1);
        case 5
            disp('initial guess method: 5')            
            pb1=(-p1(1)+sqrt(p1(1)^2+16*p1(1)*abs(Tol(1))))/8;
            nb1=ni(1)^2/pb1;
            nbn=(-n1(end)+sqrt(n1(end)^2+8*n1(end)*abs(Tol(end))))/4;
            pbn=ni(N)^2/nbn;
            %
            p(1:N1,1)=pb1;  % initial condition of the hole density for the first region. It depends on the acceptor doping density.
            n(1:N1,1)=ni(1:N1,1).^2./p(1:N1,1); % in the equilibrium condition, p*n=ni^2.
            %
            n(N1+1:N2,1)=ni(N1+1:N2);% initial condition of the electron density for the second region. It depends on the donor doping density.
            p(N1+1:N2,1)=ni(N1+1:N2);% in the equilibrium condition, p*n=ni^2
            %
            n(N2+1:N,1)=nbn;% initial condition of the electron density for the third region. It depends on the donor doping density.
            p(N2+1:N,1)=ni(N2+1:N,1).^2./n(N2+1:N,1);% in the equilibrium condition, p*n=ni^2
            % initial and boundary condition for the potential
            w=zeros(N,1); % w is the potential unit V
            w(1:N1) = -cumsum(log(-ni(1:N1)./Tol(1:N1))*NP.VT)/N;
            w(N1+1:N)= cumsum(log(abs(Tol(N1+1:N))./ni(N1+1:N))*NP.VT)/N;
            
            % boundary condition for p, n, and w
            w=w-Xi-Eg/2;
            w=w/VT;
            wbf=w(1);
            wbn=w(end);
            %
            pm=p(2:end-1);
            nm=n(2:end-1);
            Js=0;
            dBias = 4*Biase/10;
            for Bias0=0:dBias:4*Biase
                w(1)=wbf-(Bias0)/VT;
                % boundary condition for p, n, and w
                wb1=w(1);
                if Bias0>0
                    ny1=N1-round(mesh1.Length/100);
                    ny2=N2+round(mesh1.Length/100);
                    if ny1<3
                        ny1 = N1;
                    end
                    if ny2>N
                        ny2 = N2;
                    end
                    
                    w(2:ny1)=w(2:ny1)-dBias/VT;
                    w(ny1:ny2)=w(ny1)+(w(ny2)-w(ny1))./(ny2-ny1).*([ny1:ny2]-ny1);
                else
                    w(2:N-1)=((w(N)-w(1))/(N-1))*(1:N-2)+w(1);
                end
                x0=[pm;nm;w(2:end-1)];% initial guess of the solution
                fd1=fd(mesh1);
                %
                ss1 = Semiconductor_Device(NP,Temperature,Length,mesh1,fd1,Tol,Beta,...
                    un,up,tn,tp,ni,w,Vps,Vns,Ep,Diameter,DeviceDiameter,gamma,n1,p1,...
                    Impact,Thermionic_e,Thermionic_h,FranzKeldysh,Ae,Be,Ah,Bh,Eg,CnAu,CpAu,Br,DeltaEA, DeltaED,...
                    Ne_ref, Np_ref, const.eps_r_silicon,Xi,k_L); % Create the Device and initial condition.
                %
                DD0=Drift_Diff_N(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,0*GL,wavelength);
                opts = optimoptions('fsolve','Display','off','TolFun',1e-8,'MaxIter',100,'Tolx',1e-6,'MaxIter',100 ...
                    ,'Diagnostics','off','Jacobian','on','DerivativeCheck','off');
                [x,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x0,opts);
                pm=x(1:N-2);
                nm=x(N-1:2*N-4);
                wm=x(2*N-3:3*N-6);
                w=[wb1;wm;wbn];
            end
            for Bias0=linspace(4*Biase, Biase,20)
                w(1)=wbf-(Bias0)/VT;
                % boundary condition for p, n, and w
                wb1=w(1);
                if Bias0>0
                    ny1=N1-round(mesh1.Length/100);
                    ny2=N2+round(mesh1.Length/100);
                    if ny1<3
                        ny1 = N1;
                    end
                    if ny2>N
                        ny2 = N2;
                    end
                    
                    w(2:ny1)=w(2:ny1)-dBias/VT;
                    w(ny1:ny2)=w(ny1)+(w(ny2)-w(ny1))./(ny2-ny1).*([ny1:ny2]-ny1);
                else
                    w(2:N-1)=((w(N)-w(1))/(N-1))*(1:N-2)+w(1);
                end
                x0=[pm;nm;w(2:end-1)];% initial guess of the solution
                fd1=fd(mesh1);
                %
                ss1 = Semiconductor_Device(NP,Temperature,Length,mesh1,fd1,Tol,Beta,...
                    un,up,tn,tp,ni,w,Vps,Vns,Ep,Diameter,DeviceDiameter,gamma,n1,p1,...
                    Impact,Thermionic_e,Thermionic_h,FranzKeldysh,Ae,Be,Ah,Bh,Eg,CnAu,CpAu,Br,DeltaEA, DeltaED,...
                    Ne_ref, Np_ref, const.eps_r_silicon,Xi,k_L); % Create the Device and initial condition.
                %
                DD0=Drift_Diff_N(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,0*GL,wavelength);
                opts = optimoptions('fsolve','Display','off','TolFun',1e-8,'MaxIter',100,'Tolx',1e-6,'MaxIter',100 ...
                    ,'Diagnostics','off','Jacobian','on','DerivativeCheck','off');
                [x,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x0,opts);
                pm=x(1:N-2);
                nm=x(N-1:2*N-4);
                wm=x(2*N-3:3*N-6);
                w=[wb1;wm;wbn];
            end            
            wb10=w(1);
            w(1)=wbf-Biase/VT;
            w(1)=wbf;
            dw1=wb10-w(1);
            ny1=N1-round(mesh1.Length/100);
            ny2=N2+round(mesh1.Length/100);
            if ny1<3
                ny1 = N1;
            end
            if ny2>N
                ny2 = N2;
            end
            
            w(2:ny1)=w(2:ny1)-dw1;
            w(ny1:ny2)=w(ny1)+(w(ny2)-w(ny1))./(ny2-ny1).*((ny1:ny2)-ny1);
            % boundary condition for p, n, and w
            wb1=w(1);
            w(2:N-1)=((w(N)-w(1))/(N-1)^2)*((1:N-2).^2)+w(1);
    end
    x0=[pm;nm;w(2:end-1)];% initial guess of the solution
    fd1=fd(mesh1);
    %%%%%
    ss1 = Semiconductor_Device(NP,Temperature,Length,mesh1,fd1,Tol,Beta,...
        un,up,tn,tp,ni,w,Vps,Vns,Ep,Diameter,DeviceDiameter,gamma,n1,p1,...
        Impact,Thermionic_e,Thermionic_h,FranzKeldysh,Ae,Be,Ah,Bh,Eg,CnAu,CpAu,Br,DeltaEA, DeltaED,...
        Ne_ref, Np_ref, const.eps_silicon,Xi,k_L); % Create the Device and initial condition.
    
    
    werror=1;
    jjj=1;
    while werror>1e-4 && jjj<20
        jjj=jjj+1;
        DD0=Drift_Diff_N(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,GL,wavelength); % last parameter is generation rate
        opts = optimoptions('fsolve','Display','iter-detailed','TolFun',1e-8,'MaxIter',100,'Tolx',1e-6,'MaxIter',100 ...
            ,'Diagnostics','off','Jacobian','on','DerivativeCheck','off');
        [x,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x0,opts);
        pm=x(1:N-2);
        nm=x(N-1:2*N-4);
        wm=x(2*N-3:3*N-6);
        [~,~,~,Jp,Jn]=DD0.Cal_Current(pm,nm,wm,GL);
        J0=NP.J0;
        Js_ini1=trapz(mesh1.Lx_half,(Jp+Jn)*J0.*djterm);
        
        wb10=wb1;
        if jjj>10
            wb1=wbf-(Biase+Js_ini1*R_load)/VT;
        else
            wb1=wbf-(Biase+(jjj/10)*Js_ini1*R_load)/VT;
        end
        werror=abs(wb1-wb10)./abs(wb1);
    end
%     w_c = [wb1;wm;wbn];
%     tm = 1+abs(fd1.diffmp(w_c))/2.458504453260142e+06;
%     tm = smooth(tm, 50, 'loess');
% %     ff = 1e7*max(E_g)/1.000000018710562;
% %     tm = abs(E_g(1:end-1)+1)/ff;
%     tm=1*tm(1:end-1);
%     x0=[pm;nm;wm;tm];
    fd1=fd(mesh1);
    ss1 = Semiconductor_Device(NP,Temperature,Length,mesh1,fd1,Tol,Beta,...
        un,up,tn,tp,ni,w,Vps,Vns,Ep,Diameter,DeviceDiameter,gamma,n1,p1,...
        Impact,Thermionic_e,Thermionic_h,FranzKeldysh,Ae,Be,Ah,Bh,Eg,CnAu,CpAu,Br,DeltaEA, DeltaED,...
        Ne_ref, Np_ref, const.eps_r_silicon,Xi,k_L); % Create the Device and initial condition.
%     werror=1;
%     jjj=1;
%     while werror>1e-4 && jjj<20
%         jjj=jjj+1;
%         DD0=Drift_Diff_N_T(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,GL,wavelength); % last parameter is generation rate
%         opts = optimoptions('fsolve', 'Algorithm', 'trust-region-reflective','Display','iter-detailed','TolFun',1e-8,'MaxIter',200,'Tolx',1e-9,'MaxIter',100 ...
%             ,'Diagnostics','on','Jacobian','off','DerivativeCheck','off');
%         [x,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x0,opts);
%         pm=x(1:N-2);
%         nm=x(N-1:2*N-4);
%         wm=x(2*N-3:3*N-6);
%         tm=x(3*N-5:4*N-8);
%         [~,~,~,~,Jp,Jn]=DD0.Cal_Current(pm,nm,wm,tm,GL);
%         J0=NP.J0;
%         Js_ini1=trapz(mesh1.Lx_half,(Jp+Jn)*J0.*djterm);
%         wb10=wb1;
%         if jjj>10
%             wb1=wbf-(Biase+Js_ini1*R_load)/VT;
%         else
%             wb1=wbf-(Biase+(jjj/10)*Js_ini1*R_load)/VT;
%         end
%         werror=abs(wb1-wb10)./abs(wb1);
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w(1)=wbf-(Bias+0.98*Js_ini1*R_load)/VT;
    tm = ones(3409, 1);
    w_c = [wb1;wm;wbn];
    E_temp = abs(fd1.diffmp(w_c));
    H = E_temp*mean(Jn+Jp);
%     tm = 1+abs(fd1.diffmp(w_c))/2.458504453260142e+06;
%     tm = smooth(tm, 50, 'loess');
%     ff = 1e7*max(E_g)/1.000000018710562;
%     tm = abs(E_g(1:end-1)+1)/ff;
%     tm=8e-10*tm(1:end-1);
%     tm=8e-10*tm(1:end-1);
%     H=1*tm(1:end-1);
%     tm=heat_solver(H);
%     tm=heat_solver_neumann(H);
% %     tm=heat_solver_newton(H);
% %     tm = heat_solver_newton_neumann(H);
% %     dm=100*(tm(1:end-1)-1);
% %     tm=tm(1:end-1) + dm;
% %     tm(end+1)=tm(end);
%     dm=1*(tm-1);
%     tm= 1 + dm;
    tm(:)=1;
%     tm(3409:1)=1;
%     figure; plot(tm);
%     tm=tm';
%     aaa=1:3409;
%     figure; plot(aaa*0.948,tm*300);
    x0=[pm;nm;wm;tm];
    fd1=fd(mesh1);
    ss1 = Semiconductor_Device(NP,Temperature,Length,mesh1,fd1,Tol,Beta,...
        un,up,tn,tp,ni,w,Vps,Vns,Ep,Diameter,DeviceDiameter,gamma,n1,p1,...
        Impact,Thermionic_e,Thermionic_h,FranzKeldysh,Ae,Be,Ah,Bh,Eg,CnAu,CpAu,Br,DeltaEA, DeltaED,...
        Ne_ref, Np_ref, const.eps_r_silicon,Xi,k_L); % Create the Device and initial condition.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     opts = optimoptions('fsolve','Algorithm','trust-region','ScaleProblem','jacobian','Display','iter-detailed','FunctionTolerance',1e-35,'StepTolerance',1e-5,'MaxIter',0 ...
%             ,'Diagnostics','on','Jacobian','off','CheckGradients',false);

%     opts = optimoptions('fsolve', ...
%     'Display', 'iter-detailed', ...
%     'FunctionTolerance', 1e-30, ...
%     'StepTolerance',1e-95, ...
%     'MaxIter', 100, ...
%     'TolX', 1e-6, ...
%     'Diagnostics', 'off', ...
%     'Jacobian','off', ...
%     'SpecifyObjectiveGradient', false, ... % MATLAB computes Jacobian automatically
%     'DerivativeCheck', 'off');  % No need for Jacobian checking

opts = optimoptions('fsolve', ...
    'Display', 'iter-detailed', ...
    'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10, ...
    'MaxIter', 100, ...
    'TolX', 1e-8, ...
    'Jacobian', 'on',...
    'CheckGradients', false);


%     DD0=Drift_Diff_N_T_Neumann(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,GL*1e-3,wavelength);
    DD0=Drift_Diff_N_T(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,GL*1e-3,wavelength);
    [y,M]=DD0.Newton_FD(x0);
%     full_M = full(M);
%     figure; plot(tm);
%     [x0,fval,exitflag,output]=fsolve(@DD0.Newton_FD,x0,opts);
    [x0,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x0,opts);
    pm=x0(1:N-2);
    nm=x0(N-1:2*N-4);
    wm=x0(2*N-3:3*N-6);
    tm=x0(3*N-5:4*N-8);
    x0=[pm;nm;wm;tm];
%     figure; plot(tm);
    if mean(tm)<1
        tm = 2-tm;
        [x0,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x0,opts);
    end
%     tm=x0(3*N-5:4*N-8);
    L=NX*mesh1.Lx;
    x1 = L*1e4;
    figure; plot(x1(2:end-1)*1000,300*tm);
%       figure; plot(tm*300);
    xlim([0 3229]);
    xlabel('Length (nm)');
    ylabel('Temperature (K)');
    keyboard;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    opts = optimoptions('fsolve', ...
    'Display', 'iter-detailed', ...
    'FunctionTolerance', 1e-30, ...
    'StepTolerance',1e-95, ...
    'MaxIter', 0, ...
    'TolX', 1e-6, ...
    'Diagnostics', 'on', ...
    'SpecifyObjectiveGradient', false, ... % MATLAB computes Jacobian automatically
    'DerivativeCheck', 'off');  % No need for Jacobian checking

%     opts = optimoptions('fsolve','Display','off','TolFun',1e-8,'MaxIter',100,'Tolx',1e-6,'MaxIter',100 ...
%             ,'Diagnostics','on','Jacobian','off','DerivativeCheck','off');
    DD0=Drift_Diff_N_T(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,GL,wavelength);
    [x0,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FDT,x0,opts);
    pm=x0(1:N-2);
    nm=x0(N-1:2*N-4);
    wm=x0(2*N-3:3*N-6);
    tm=x0(3*N-5:4*N-8);
    
%     tm(1:3409)=1;
%     x0=[pm;nm;wm;tm'];
%     opts = optimoptions('fsolve', ...
%     'FunctionTolerance', 1e-30, ... % Much smaller than 1e-15 to ensure strict convergence
%     'StepTolerance', 1e-30, ... % Ensure small step size tolerance
%     'OptimalityTolerance', 1e-30, ... % Ensure solver keeps refining the solution
%     'MaxIterations', 1000, ... % Allow more iterations if necessary
%     'MaxFunctionEvaluations', 100, ... % Allow more function evaluations
%     'Display', 'iter-detailed'); % Show iteration progress
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DD0=Drift_Diff_N_R(pb1,pbn,nb1,nbn,wbf,wbn,fd1,mesh1,PD_Str1,NP,ss1,GL,Bias,R_load,DeviceDiameter,wavelength);
    x1=[x;wb1];
    [x1,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x1,opts);
    pm=x1(1:N-2);
    nm=x1(N-1:2*N-4);
    wm=x1(2*N-3:3*N-6);
    wb1=x1(end);

 if max(abs(imag(wm)))>0 || min(pm)<0 || min(nm)<0
     disp([length(find(pm<0)), length(find(nm<0))])
 end
    
    w(1)=wbf-(Bias+0.98*Js_ini1*R_load)/VT;
    x0=[pm;nm;wm];
    fd1=fd(mesh1);
    ss1 = Semiconductor_Device(NP,Temperature,Length,mesh1,fd1,Tol,Beta,...
        un,up,tn,tp,ni,w,Vps,Vns,Ep,Diameter,DeviceDiameter,gamma,n1,p1,...
        Impact,Thermionic_e,Thermionic_h,FranzKeldysh,Ae,Be,Ah,Bh,Eg,CnAu,CpAu,Br,DeltaEA, DeltaED,...
        Ne_ref, Np_ref, const.eps_r_silicon,Xi,k_L); % Create the Device and initial condition.
    
    DD0=Drift_Diff_N_B(pb1,pbn,nb1,nbn,wb1,wbn,fd1,mesh1,PD_Str1,NP,ss1,GL,wavelength);
    [x0,fval,exitflag,output,jacobian]=fsolve(@DD0.Newton_FD,x0,opts);
    pm=x0(1:N-2);
    nm=x0(N-1:2*N-4);
    wm=x0(2*N-3:3*N-6);

 if max(abs(imag(wm)))>0 || min(pm)<0 || min(nm)<0
     disp([length(find(pm<0)), length(find(nm<0))])
     % keyboard
 end
%     [~,~,~,Jp,Jn,ave_alpha]=DD0.Cal_Current(pm,nm,wm,GL);
    [~,~,~,Jp,Jn,~]=DD0.Cal_Current(pm,nm,wm,GL);

    J0=NP.J0;
    NX=NP.NX;
    Ni=NP.Ni;
    p=[pb1;pm;pbn];
    n=[nb1;nm;nbn];
    w=VT*[wb1;wm;wbn];
    L=NX*mesh1.Lx;
    
    DD2 = Drift_Diff_N_C_B_R(pb1,pbn,nb1,nbn,wbf,wbn,fd1,mesh1,ss1,NP,GL,PD_Str1,...
        TL,Bias,R_load,Device_Diameter,wavelength);
    x1=[x0;wb1];
    opts = optimoptions('fsolve','Display','off','TolFun',1e-8,'MaxIter',100,'Tolx',1e-6,'MaxIter',100 ...
        ,'Diagnostics','off','Jacobian','on','DerivativeCheck','off');
    [x1,fval,exitflag,output,jacobian]=fsolve(@DD2.Newton_FD,x1,opts);
    
    pm=x1(1:N-2);
    nm=x1(N-1:2*N-4);
    wm=x1(2*N-3:3*N-6);
    
    wb1=x1(end);
    w=VT*[wb1;wm;wbn];
    p=[pb1;pm;pbn];
    n=[nb1;nm;nbn];
    [~,~,~,~,Jp,Jn,Jpd,Jpdu,Doping_ion,ave_alpha,G,R,Vn_a,Vp_a,NAP,NDP]=DD2.Cal_Currentv2(pm,nm,wm,wb1,GL);
    Js_n = trapz(mesh1.Lx_half, Jn*J0.*djterm);
    Js_p = trapz(mesh1.Lx_half, Jp*J0.*djterm);
    Js_C = Js_n + Js_p;
    %
    p_c=Ni*[pb1;pm;pbn];
    n_c=Ni*[nb1;nm;nbn];
    w_c=VT*[wb1;wm;wbn];

    if max(abs(imag(wm)))>0 || min(pm)<0 || min(nm)<0
       disp([length(find(pm<0)), length(find(nm<0))])
        initial_guess_method = initial_guess_method +1;
        
        if initial_guess_method ==6
            disp('No Steady State Solution!');    
            converged =0;
            Sol.converged = 0;
            w0a = -Bias;
            w0b = 0;
            w_c(1:N1,1) = w0a;
            w_c(N2+1:N,1) = w0b ;
            w_c(N1:N2+1,1) = linspace(w0a,w0b,N2-N1+2);

            xx = mesh1.Lx*mesh1.NX*1e4;
            f1 = -w_c+PD_Str1.Eg/2-PD_Str1.delta_Fermis/2;
            f2 = -w_c-PD_Str1.Eg/2-PD_Str1.delta_Fermis/2;
            f1 = f1-f2(1);
            f2 = f2-f2(1);
            Sol.BandC = f1;
            Sol.BandV = f2;
            Sol.BandX =xx;                       
            break
        end
    else
        converged = 1;

    end
    % % keyboard
end
% %%%%%%%%%%%%%
E_c=-(fd1.diffmp(w_c))/NX;
%%%%%%%%%%%
Sol.pn_mean = mean(p_c(mesh1.Node(interfaces(2)):end));
Sol.L=L;
Sol.p=p_c;
Sol.n=n_c;
Sol.w=w_c;
Sol.pm=pm;
Sol.nm=nm;
Sol.wm=wm;
Sol.Js_C=Js_C;
Sol.Js_n=Js_n;
Sol.Js_p=Js_p;
Sol.Jp=Jp;
Sol.Jn=Jn;
Sol.E=E_c;

Sol.djterm = djterm;
Sol.J0 = NP.J0;

Sol.PD_Str = PD_Str1;
Sol.params = params;
Sol.Design = Design;
Sol.Boundary.pb1=pb1;
Sol.Boundary.pbn=pbn;
Sol.Boundary.nb1=nb1;
Sol.Boundary.nbn=nbn;
Sol.Boundary.wb1=wb1;
Sol.Boundary.wbn=wbn;
Sol.Boundary.wbf=wbf;

Sol.Parameters.NP=NP;
Sol.Parameters.fd1=fd1;
Sol.Parameters.mesh1=mesh1;
Sol.Parameters.ss1=ss1;
Sol.GL=GL;
Sol.TL=TL;
Sol.Parameters.R_load=R_load;
Sol.Parameters.DeltaEA=DeltaEA;
Sol.Parameters.DeltaED=DeltaED;
Sol.NAP=NAP;
Sol.NDP=NDP;
Sol.Stats.ave_alpha=ave_alpha;
Sol.converged = converged;
%%%%%%%%%%%

xx = mesh1.Lx*mesh1.NX*1e4;
f1 = -w_c+PD_Str1.Eg/2-PD_Str1.delta_Fermis/2;
f2 = -w_c-PD_Str1.Eg/2-PD_Str1.delta_Fermis/2;

f1 = f1-f2(1);
f2 = f2-f2(1);
Sol.BandC = f1;
Sol.BandV = f2;
Sol.BandX =xx;

if plot_results == 1 && converged ==1

    Emax_plot = max([50, 10*(round(max(max(abs(E_c)/1e3)/10)))+10]);
    low_lim = min([-20, -round(max(max(abs(E_c)/1e3)/10))]);

    xs = [0 cumsum(mesh1.zis)]*1e-3;
    Face_Colors = [0.9 0.6 0.8; 1 0.5 0; 0.3 0.6 0.2; 0.6 0.4 0.2; 0.75 0.67 0.89; rand(9,3)];
    Face_Colors2 = [227, 242, 253; 237, 231, 246; 255, 235, 238]/255;

    ys2 = [0 Emax_plot];
    ys = [low_lim  Emax_plot-low_lim];
    ys3 = [low_lim abs(low_lim)];
    %
    x1 = L*1e4;
    x2 = L(1:end-1)*1e4;
    %
    SurfaceArea = DeviceDiameter(1:end-1).^2*0.25*pi;
    figure(filenum); clf;
    subplot(231);
    plot(x2,1e3*Jn.*J0.*SurfaceArea,x2,1e3*Jp.*J0.*SurfaceArea,'r');
    legend('\it{I_n}','\it{I_p}','Location','East'); grid on;
    ylabel('Currents (mA)')
    xlabel('Position (\mu{m})');
    xlim([0 max(x2)])

    subplot(232);
    semilogy(x1(2:end-1),n_c(2:end-1),x1(2:end-1),p_c(2:end-1),'linewidth',3); grid on;
    xlabel('Position (\mu{m})');
    ylabel('Conc. (cm^{-3})')
    xlim([0 max(x1)])
    legend('Electron','Hole','Location','East')

    subplot(233);
    plot(x1,log10(abs(PD_Str1.Tol)),'linewidth',3); grid on;
    ylabel('Doping (cm^{-3})')
    xlabel('Position (\mu{m})');
    xlim([0 max(x1)])

    subplot(234);
    for ii = 1:PD_Str1.Num_Layer
        rectangle('Position',[xs(ii),ys3(1),xs(ii+1)-xs(ii),ys3(2)],'FaceColor',Face_Colors(PD_Str1.what_sc(ii),:),'EdgeColor','w','LineWidth',1); hold on
        rectangle('Position',[xs(ii),ys2(1),xs(ii+1)-xs(ii),ys2(2)],'FaceColor',Face_Colors2(PD_Str1.Type_Layer(ii)+2,:),'EdgeColor',Face_Colors2(PD_Str1.Type_Layer(ii)+2,:)); hold on
    end
    plot(x2,abs(E_c)/1e3,'linewidth',3,'color','k'); grid on;
    xlabel('Position (\mu{m})');
    ylabel('|{\bf{E}}| (KV/cm)')
    xlim([0 max(x2)]); ylim([low_lim Emax_plot]);

    subplot(235);
    plot(x1,w_c-w_c(end),'linewidth',3); grid on;
    ylabel('Potential (V)')
    xlabel('Position (\mu{m})');
    xlim([0 max(x1)])

    subplot(236);
    plot(xx,f1,xx,f2); grid on;
    xlim([0 xx(end)]);
    xlabel('Position (\mu{m})');
    ylabel('Energy Band (eV)');
    
    if params.print_plots == 1
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.04, 0.04, 0.6, 0.6]);
        oldscreenunits = get(gcf,'Units');
        oldpaperunits = get(gcf,'PaperUnits');
        oldpaperpos = get(gcf,'PaperPosition');
        set(gcf,'Units','pixels');
        scrpos = get(gcf,'Position');
        newpos = scrpos/100;
        set(gcf,'PaperUnits','inches','PaperPosition',newpos)
        print('-dpng', [Solution_Name '_static'], '-r300');
        drawnow
        set(gcf,'Units',oldscreenunits,'PaperUnits',oldpaperunits,'PaperPosition',oldpaperpos)
%         !mv *.png ./figures/
    end
end

function Temp = heat_solver(H)
% Solves the 1D heat conduction equation with given heat source H
% H: source term (Nx x 1 vector)
% T: temperature distribution (Nx x 1 vector)

% Define parameters
Nx = length(H);  % Number of grid points
L = 86.4;    % Length of domain (arbitrary units)
dx = L / (Nx - 1); % Grid spacing

% Thermal conductivity (assumed constant for simplicity)
kappa_L = 2.1; 

% Construct finite difference matrix for 1D Laplacian
main_diag = -2 * ones(Nx, 1);
off_diag = ones(Nx-1, 1);
A = (1/dx^2) * (diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1));

% Apply boundary conditions (Dirichlet: T = 1 at both ends)
A(1, :) = 0; A(1, 1) = 1;
A(end, :) = 0; A(end, end) = 1;
H(1) = 1;  % Corresponding Dirichlet values
H(end) = 1;

% Solve system Ax = H for temperature T
Temp = A \ H;

% Plot results
% figure;
% plot(linspace(0, L, Nx), T, 'b', 'LineWidth', 1.5);
% xlabel('Position');
% ylabel('Temperature');
% title('Numerical Solution of Heat Equation');
% grid on;
% end

function Temp = heat_solver_neumann(H)
% Solves 1D steady-state heat conduction equation:
%   -d²T/dx² = H(x)/kappa, with:
%   Dirichlet BC (T = 1) at x = 0 (left boundary)
%   Neumann BC (dT/dx = 0) at x = L (right boundary)
%
% Input:
%   H    : Heat source term vector (Nx x 1)
% Output:
%   Temp : Temperature distribution (Nx x 1)

% Parameters
Nx = length(H);       % Number of grid points
L = 86.4;             % Length of the domain (units)
dx = L / (Nx - 1);    % Grid spacing
kappa = 2.1;          % Thermal conductivity (not used for matrix here)

% Construct finite-difference matrix (1D Laplacian)
main_diag = -2 * ones(Nx, 1);
off_diag = ones(Nx - 1, 1);
A = (1 / dx^2) * (diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1));

% Apply Dirichlet BC at left end (x = 0): T(1) = 1
A(1, :) = 0;
A(1, 1) = 1;
H(1) = 1;

% Apply Neumann BC at right end (x = L): dT/dx = 0 -> T(end) - T(end-1) = 0
A(end, :) = 0;
A(end, end-1) = -1;
A(end, end) = 1;
H(end) = 0;

% Solve linear system A * Temp = H
Temp = A \ H;



function T = heat_solver_newton(H)
% Solves -d²T/dx² = H(x)/kappa using Newton's method
% Inputs:
%   H : heat generation vector (Nx x 1)
% Output:
%   T : temperature profile (Nx x 1)

% --- Parameters ---
Nx = length(H);             % Number of grid points
L = 86.4;                   % Length of domain
dx = L / (Nx - 1);          % Grid spacing
kappa_L = 2.1;              % Thermal conductivity (constant)

% --- Initial guess for temperature ---
T = ones(Nx, 1);            % Start with uniform guess

% --- Newton iteration settings ---
max_iter = 100;
tol = 1e-6;

% --- Newton Loop ---
for iter = 1:max_iter

    % Step 1: Construct residual vector F
    F = zeros(Nx, 1);

    % Interior nodes: apply central difference and subtract RHS H/kappa
    for i = 2:Nx-1
        Ti = T(i);
        % Discretized Laplacian: (T_{i+1} - 2T_i + T_{i-1}) / dx²
        % Residual: -Laplacian - H/kappa
        F(i) = (-T(i+1) + 2*T(i) - T(i-1)) / dx^2 - H(i) / kappa_L;
    end

    % Dirichlet BC at left: T(1) = 1
    F(1) = T(1) - 1;

    % Dirichlet BC at right: T(end) = 1
    F(end) = T(end) - 1;

    % Step 2: Construct Jacobian matrix J = dF/dT
    J = zeros(Nx, Nx);

    for i = 2:Nx-1
        % Partial derivatives of residual F(i) w.r.t. neighboring T values
        J(i, i-1) = -1 / dx^2;
        J(i, i)   =  2 / dx^2;
        J(i, i+1) = -1 / dx^2;
    end

    % Jacobian rows for Dirichlet BCs
    J(1, :) = 0; J(1, 1) = 1;
    J(end, :) = 0; J(end, end) = 1;

    % Step 3: Newton update
    dT = -J \ F;
    T = T + dT;

    % Step 4: Convergence check
    if norm(dT, inf) < tol
        fprintf('Converged in %d iterations\n', iter);
        return
    end
end

error('Newton method did not converge within %d iterations.', max_iter);


function T = heat_solver_newton_neumann(H)
% Solves -d²T/dx² = H(x)/kappa using Newton's method
% BCs: Neumann at x=0 (∂T/∂x = 0), Dirichlet at x=L (T = 1)

% --- Parameters ---
Nx = length(H);             % Number of grid points
L = 86.4;                   % Domain length
dx = L / (Nx - 1);          % Grid spacing
kappa_L = 2.1;              % Constant thermal conductivity

% --- Initial guess ---
T = ones(Nx, 1);            % Initial temperature guess

% --- Newton iteration settings ---
max_iter = 100;
tol = 1e-15;

% --- Newton loop ---
for iter = 1:max_iter

    % Step 1: Residual vector F
    F = zeros(Nx, 1);

    % Neumann BC at left (∂T/∂x = 0) ⇒ T(2) - T(1) = 0
    F(1) = T(2) - T(1);

    % Interior nodes: apply discretized PDE
    for i = 2:Nx-1
        F(i) = (-T(i+1) + 2*T(i) - T(i-1)) / dx^2 - H(i) / kappa_L;
    end

    % Dirichlet BC at right (T = 1)
    F(end) = T(end) - 1;

    % Step 2: Jacobian matrix J = dF/dT
    J = zeros(Nx, Nx);

    % Neumann BC row: F(1) = T(2) - T(1)
    J(1,1) = -1;
    J(1,2) = 1;

    % Interior nodes
    for i = 2:Nx-1
        J(i, i-1) = -1 / dx^2;
        J(i, i)   =  2 / dx^2;
        J(i, i+1) = -1 / dx^2;
    end

    % Dirichlet BC row
    J(end, :) = 0;
    J(end, end) = 1;

    % Step 3: Newton update
    dT = -J \ F;
    T = T + dT;

    % Step 4: Check convergence
    if norm(dT, inf) < tol
        fprintf('Converged in %d iterations\n', iter);
        return
    end
end

error('Newton did not converge in %d iterations', max_iter);



