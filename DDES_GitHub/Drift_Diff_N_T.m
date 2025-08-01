classdef Drift_Diff_N_T
    %DRIFT_DIFF Summary of this class goes here
    properties        
        pb1% boundary condition for hole
        pbn% boundary condition(end) for hole
        nb1% boundary condition for electron
        nbn% boundary condition(end) for electron
        wb1% boundary condition for potential
        wbn% boundary condition(end) for potential
        fd
        wavelength
        NP
        mesh
        PD_Str
        Semiconductor_Device
        GL % generation rate        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function self=Drift_Diff_N_T(pb1,pbn,nb1,nbn,wb1,wbn,fd,mesh,PD_Str,NP,Semiconductor_Device,GL,wavelength)
            self.pb1=pb1;
            self.wavelength = wavelength;
            self.pbn=pbn;
            self.nb1=nb1;
            self.nbn=nbn;
            self.wb1=wb1;
            self.wbn=wbn;
            self.fd = fd;
            self.mesh=mesh;
            self.PD_Str=PD_Str;
            self.NP=NP;
            self.Semiconductor_Device = Semiconductor_Device;
            self.GL=GL;
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Fp,Fn,Fw,Ft,J_Hole,J_Elec,alpha]=Cal_Current(self,pm,nm,phim,tm,GL)
            % calculate the function for drift&diffusion equation and
            % possion equation
            % the following are the parameters for the device.
            p=[self.pb1;pm;self.pbn];
            n=[self.nb1;nm;self.nbn];
            phi=[self.wb1;phim;self.wbn];
            T=[1;tm;1];
            N1=self.mesh.N1;
            N=self.mesh.N;
%             kB = self.NP.KB;
            kB = 1;
%             q = self.NP.q;
            q = 1;
%             k_L=2.1;
            k_L = self.Semiconductor_Device.k_L;
            tau_n = self.Semiconductor_Device.Lifetime_n;
            ni = self.Semiconductor_Device.Intrinsic_Density;
            Doping = self.Semiconductor_Device.Doping_Profile;
            n1=self.Semiconductor_Device.n1;
            p1=self.Semiconductor_Device.p1;
            % incomplete ionization model
            Doping_ion=zeros(N,1);
            Doping_ion(1:N1)=Doping(1:N1)./(1+4*p(1:N1)./p1(1:N1));
            Doping_ion(N1+1:end)=Doping(N1+1:end)./(1+2*n(N1+1:end)./n1(N1+1:end));
            Electron_Mobility = self.Semiconductor_Device.Electron_Mobility;
            Hole_Mobility = self.Semiconductor_Device.Hole_Mobility;
            %%%%%%%%%%%%            
            E = -self.fd.diffmp(phi);% electrical field
            %%%%%%%%%%%%            
            Electron_Velocity = E.*Electron_Mobility(2:end); % Electron Velocity
            Hole_Velocity = E.*Hole_Mobility(2:end); % Hole Velocity
            %%%%%%%%%%%%            
            Recombination=(n.*p)./(tau_n.*(n+p+2*ni)); % Recombination-ni.*ni            
            VT=self.Semiconductor_Device.NP.VT;
            NX=self.Semiconductor_Device.NP.NX;
%             Generation =0*GL;
            E_main = -self.fd.diff3p(phi);
            E_main_all=[E_main(1);E_main;E_main(end)];
            VT=self.Semiconductor_Device.NP.VT;
            NX=self.Semiconductor_Device.NP.NX;
            alpha0=Absorp_FK(abs(E_main_all)*VT/NX*100,self.wavelength,self.PD_Str,self.Semiconductor_Device.FranzKeldysh);                        
            absorbers = ones(length(alpha0),1);
            absorbers(alpha0==0) = 0;
            cell_loss = exp(-alpha0.*NX.*[self.mesh.dx; 0]);
            structure_loss = flipud(cumprod(flipud(cell_loss)));
            Generation = GL.*structure_loss.*absorbers;  

            E_main = -self.fd.diff3p(phi); % the electric field at the main point
            G_imp=0;                        
            J_Hole = Hole_Velocity.*self.fd.averagemp(p)...% Hole current,there should be a q times, but it cancelled in the Fp.
                - Hole_Mobility(2:end).*self.fd.diffmp(p) - kB*Hole_Mobility(2:end).*self.fd.averagemp(p).*self.fd.diffmp(T)/q;
            %%%%%%%%%%%%%%%%%%%%            
            if numel(J_Hole) == 3409
                disp('J_Hole has a size of 3409. Entering debug mode...');
                keyboard; % Pauses execution for interactive debugging
%             else
%                 disp('J_Hole does not have a size of 3409. Continuing execution...');
            end
            %%%%%%%%%%%%%%%%%%%%
            J_Elec_drift = Electron_Velocity.*self.fd.averagemp(n);
            J_Elec_diff = Electron_Mobility(2:end).*self.fd.diffmp(n);
            J_Elec_th = kB*Electron_Mobility(2:end).*self.fd.averagemp(n).*self.fd.diffmp(T)/q;
            J_Elec = J_Elec_drift + J_Elec_diff+ J_Elec_th; % electron current
            Fp = -self.fd.diffmidp(J_Hole)-Recombination(2:end-1)+G_imp+Generation(2:end-1); % drift % diffusion current equation for hole
            Fn = self.fd.diffmidp(J_Elec)-Recombination(2:end-1)+G_imp+Generation(2:end-1);% drift % diffusion current equation for electron
            dx=self.fd.mesh.dx;
            dx1=self.fd.mesh.dx1;
            Fw = (Doping_ion(2:end-1)+p(2:end-1)-n(2:end-1))... % ni(2:end-1)is norm
                +(1./((dx(1:end-1)).*(dx1(1:end)))).*phi(1:end-2)...
                -((1./dx(1:end-1)+1./dx(2:end))./dx1(1:end)).*phi(2:end-1)...
                +(1./((dx(2:end)).*(dx1(1:end)))).*phi(3:end);% possion equation 
%             Ft = k_L*((1./((dx(1:end-1)).*(dx1(1:end)))).*T(1:end-2)...         % 4.2e10
%                 -((1./dx(1:end-1)+1./dx(2:end))./dx1(1:end)).*T(2:end-1)...
%                 +(1./((dx(2:end)).*(dx1(1:end)))).*T(3:end)) + ((J_Elec(2:end)+J_Hole(2:end)).*E(2:end)+(J_Elec(1:end-1)+J_Hole(1:end-1)).*E(1:end-1))/2;

            Ft = k_L(2:end-1).*((1./((dx(1:end-1)).*(dx1(1:end)))).*T(1:end-2)...         % 4.2e10
                -((1./dx(1:end-1)+1./dx(2:end))./dx1(1:end)).*T(2:end-1)...
                +(1./((dx(2:end)).*(dx1(1:end)))).*T(3:end)) + abs((E(2:end).*(J_Elec(2:end)+J_Hole(2:end)) + E(1:end-1).*(J_Elec(1:end-1)+J_Hole(1:end-1))))/2;
%               Ft = k_L*((1./((dx(1:end-1)).*(dx1(1:end)))).*T(1:end-2)...         % 4.2e10
%                 -((1./dx(1:end-1)+1./dx(2:end))./dx1(1:end)).*T(2:end-1)...
%                 +(1./((dx(2:end)).*(dx1(1:end)))).*T(3:end)) + ((E(2:end).*(J_Elec(2:end)) + E(1:end-1).*(J_Elec(1:end-1))))/2;
%             figure; plot(abs((E(2:end).*(J_Elec(2:end)+J_Hole(2:end)) + E(1:end-1).*(J_Elec(1:end-1)+J_Hole(1:end-1))))/2);
%               figure; plot(abs((E(2:end).*(J_Elec(2:end)+J_Hole(2:end)) + E(1:end-1).*(J_Elec(1:end-1)+J_Hole(1:end-1))))/2);
%             hold on;
%             Ft = k_L*((1./((dx(1:end-1)).*(dx1(1:end)))).*T(1:end-2)...         % 4.2e10
%                 -((1./dx(1:end-1)+1./dx(2:end))./dx1(1:end)).*T(2:end-1)...
%                 +(1./((dx(2:end)).*(dx1(1:end)))).*T(3:end));
          
%             Ft = k_L*((1./((dx(1:end-1)).*(dx1(1:end)))).*T(1:end-2)...         % 4.2e10
%                 -((1./dx(1:end-1)+1./dx(2:end))./dx1(1:end)).*T(2:end-1)...
%                 +(1./((dx(2:end)).*(dx1(1:end)))).*T(3:end)) + ((E(2:end).*(J_Elec(2:end)) + E(1:end-1).*(J_Elec(1:end-1))))/2;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [JacobM]=ConsJacobian(self,pm,nm,phim,tm)
            p=[self.pb1;pm;self.pbn];
            n=[self.nb1;nm;self.nbn];
            phi=[self.wb1;phim;self.wbn];
            T=[1;tm;1];
            tau_n = self.Semiconductor_Device.Lifetime_n;
            tau_p = self.Semiconductor_Device.Lifetime_p;
            ni = self.Semiconductor_Device.Intrinsic_Density;
            Doping = self.Semiconductor_Device.Doping_Profile;
            n1=self.Semiconductor_Device.n1;
            p1=self.Semiconductor_Device.p1;
            N=self.mesh.N;
            N1=self.mesh.N1;
%             kB = self.NP.KB;
            kB = 1;
%             q = self.NP.q;
            q = 1;
%             k_L = 2.1;
            k_L = self.Semiconductor_Device.k_L;
            % incomplete ionization model
            Doping_ion=zeros(N,1);
            Doping_ion(1:N1)=Doping(1:N1)./(1+4*p(1:N1)./p1(1:N1));
            Doping_ion(N1+1:end)=Doping(N1+1:end)./(1+2*n(N1+1:end)./n1(N1+1:end));
            % Doing_ion  derive to the p and n
            Doping_ion_p=zeros(N,1);
            Doping_ion_n=zeros(N,1);
            Doping_ion_p(1:N1)=-Doping_ion(1:N1)./(p1(1:N1)/4+p(1:N1));
            Doping_ion_n(N1+1:end)=-Doping_ion(N1+1:end)./(n1(N1+1:end)/2+n(N1+1:end));
            
            E = -self.fd.diffmp(phi);  
            Electron_Mobility = self.Semiconductor_Device.Electron_Mobility;
            Electron_Velocity = E.*Electron_Mobility(2:end); % Electron Velocity
            Hole_Mobility = self.Semiconductor_Device.Hole_Mobility;
            Hole_Velocity = E.*Hole_Mobility(2:end); % Hole Velocity
            E = -self.fd.diffmp(phi);   
            
            J_Hole = Hole_Velocity.*self.fd.averagemp(p)...% Hole current,there should be a q times, but it cancelled in the Fp.
                - Hole_Mobility(2:end).*self.fd.diffmp(p) - kB*Hole_Mobility(2:end).*self.fd.averagemp(p).*self.fd.diffmp(T)/q;
            J_Elec_drift = Electron_Velocity.*self.fd.averagemp(n);
            J_Elec_diff = Electron_Mobility(2:end).*self.fd.diffmp(n);
            J_Elec_th = kB*Electron_Mobility(2:end).*self.fd.averagemp(n).*self.fd.diffmp(T)/q;
            J_Elec = J_Elec_drift + J_Elec_diff+ J_Elec_th; % electron current
            
            Recombination=(n.*p)./(tau_n.*(n+p+2*ni));
            DUDp=(n-tau_n.*Recombination)./(tau_n.*(p+ni)+tau_p.*(n+ni));
            % first derivative  Recombination(i) to Hole_Density(i)
            DUDn=(p-tau_p.*Recombination)./(tau_n.*(p+ni)+tau_p.*(n+ni));
            % this vn and vp should be on the main point;
            DG_impDn=0;
            DG_impDp=0;
            DG_impDwm=zeros(N-2,1);
            DG_impDwp=zeros(N-2,1);
                            
            %CreateJacobM()
            Electron_Velocity=E.*Electron_Mobility(2:end); % Electron Velocity
            Hole_Velocity=E.*Hole_Mobility(2:end); % Hole Velocity
            DJpDp = Hole_Velocity/2+Hole_Mobility(2:end).*self.fd.differ + kB*Hole_Mobility(2:end).*self.fd.diffmp(T)/2/q;% first derivative hole current J_Hole(i) to Hole_Density(i)
            DJpDp2 = Hole_Velocity/2-Hole_Mobility(2:end).*self.fd.differ + kB*Hole_Mobility(2:end).*self.fd.diffmp(T)/2/q;% first derivative J_Hole(i) to Hole_Density(i+1)            
            DJpDw = Hole_Mobility(2:end).*self.fd.averagemp(p).*self.fd.differ; % first derivative J_Hole(i) to Potential(i)
            DJpDw2 = -DJpDw; %  % first derivative J_Hole(i) to Potential(i+1)
            
            % J_Hole first derivative to T  
            DJpDT = -kB*Hole_Mobility(2:end).*self.fd.averagemp(p).*self.fd.differ/q;
            DJpDT2 = kB*Hole_Mobility(2:end).*self.fd.averagemp(p).*self.fd.differ/q;
            %
            DJnDn = Electron_Velocity/2-Electron_Mobility(2:end).*self.fd.differ + kB*Electron_Mobility(2:end).*self.fd.diffmp(T)/2/q; % % first derivative electron current J_Elec(i) to Elec_Density(i)
            DJnDn2 = Electron_Velocity/2+Electron_Mobility(2:end).*self.fd.differ + kB*Electron_Mobility(2:end).*self.fd.diffmp(T)/2/q; % first derivative current J_Elec(i) to Elec_Density(i+1)
            DJnDw = Electron_Mobility(2:end).*(self.fd.averagemp(n)).*self.fd.differ;  % first derivative  J_Elec(i) to potential(i)
            DJnDw2 = -DJnDw; % first derivative  J_Elec(i) to potential(i+1)
            
            % J_electron first derivative to T  
            DJnDT = -kB*Electron_Mobility(2:end).*self.fd.averagemp(n).*self.fd.differ/q;
            DJnDT2 = kB*Electron_Mobility(2:end).*self.fd.averagemp(n).*self.fd.differ/q;
            
            %calculate the current of hole and electrons
            %so the Jacob matrix***DFp/Dp,DFp/Dn,DFp/Dw************
            %**********************DFn/Dp,DFn/Dn,DFn/Dw************
            %**********************DFw/Dp,DFw/Dn,DFw/Dw************
%             dx=self.fd.differ;
            dx=self.fd.mesh.dx;
            dx1=self.mesh.dx1;
            Num=length(dx1)+2;
            
            % Fp derivative to p , n and phi
            a=[DJpDp(2:end-1)./(dx1(2:end));0];
            e=(-DJpDp(2:end)+DJpDp2(1:end-1))./(dx1(1:end))-DUDp(2:end-1)+DG_impDp;%(2:N-1)
            b=[0;-DJpDp2(2:end-1)./(dx1(1:end-1))];
            M1=spdiags([a e b],-1:1,Num-2,Num-2); % Fp first derivative to p(hole denstiy)
            M2=spdiags(-DUDn(2:end-1)+DG_impDn,0,Num-2,Num-2); % Fp first derivative to n (electron denstiy)
            a=[DJpDw(2:end-1)./(dx1(2:end))+DG_impDwm(2:end);0];
            e=(-DJpDw(2:end)+DJpDw2(1:end-1))./(dx1(1:end));
            b=[0;-DJpDw2(2:end-1)./(dx1(1:end-1))+DG_impDwp(1:end-1)];
            M3=spdiags([a e b],-1:1,Num-2,Num-2);% Fp first derivative to potential
            a = [-DJpDT(2:end-1)./(dx1(2:end));0];
            e = (DJpDT(2:end)-DJpDT2(1:end-1))./(dx1(1:end));
            b=[0;DJpDT2(2:end-1)./(dx1(1:end-1))];
            M4=spdiags([a e b],-1:1,Num-2,Num-2); % Fp first derivative to T(temperature)
            MFp=[M1,M2,M3,M4]; % Fp first derivative to p ,n and potential
            
            % Fn derive to p n w
            M1=spdiags(-DUDp(2:end-1)+DG_impDp,0,Num-2,Num-2); % Fn first derivative to p(hole denstiy)
            a=[-DJnDn(2:end-1)./(dx1(2:end));0];
            e=(DJnDn(2:end)-DJnDn2(1:end-1))./(dx1(1:end))-DUDn(2:end-1)+DG_impDn;
            b=[0;DJnDn2(2:end-1)./(dx1(1:end-1))];
            M2=spdiags([a e b],-1:1,Num-2,Num-2); % Fn first derivative to n(electron denstiy)
            a=[-DJnDw(2:end-1)./(dx1(2:end))+DG_impDwm(2:end);0]; % diag-1,row :i, col; i-1
            e=(DJnDw(2:end)-DJnDw2(1:end-1))./(dx1(1:end)); % diag
            b=[0;DJnDw2(2:end-1)./(dx1(1:end-1))+DG_impDwp(1:end-1)];
            M3=spdiags([a e b],-1:1,Num-2,Num-2); % Fn first derivative to potential
            a = [-DJnDT(2:end-1)./(dx1(2:end));0];
            e = (DJnDT(2:end)-DJnDT2(1:end-1))./(dx1(1:end));
            b=[0;DJnDT2(2:end-1)./(dx1(1:end-1))];
            M4=spdiags([a e b],-1:1,Num-2,Num-2); % Fn first derivative to T(temperature)
            MFn=[M1,M2,M3,M4];
            
            % Fw first derivative to p n w T
            DFphiDp = self.fd.CreateZOM+spdiags(Doping_ion_p(2:end-1),0,Num-2,Num-2);% Fphi first derivative to p
            DFphiDn = -self.fd.CreateZOM+spdiags(Doping_ion_n(2:end-1),0,Num-2,Num-2); % Fphi first derivative to n
            DFphiDphi= self.fd.CreateSOM; % Fp first derivative to phi
            DFphiDT = zeros(Num-2,Num-2);
            MFw=[DFphiDp,DFphiDn,DFphiDphi,DFphiDT];
            
%             Ft derivative to p n w t

%             DHDp = DJpDp(2:end);
% %             DHDp2 = DJpDp2*self.fd.diffmp(w) - (DG_impDp2 + DG_imp2Dp2)/2 + (DUDp(2:end+1) + Dre_auDp(2:end+1) + Dre_radDp(2:end+1))/2;
% %             DHDp_1 = DG_impDp_1/2;
%             
%             DH_1Dp = DJpDp2(1:end-1);
% %             DH_1Dp2 = - (DG_impDp2)/2;
%             DH_1Dp_1 = DJpDp2(2:end-1);
            
%             DHDp3 = -DG_imp2Dp3/2;
%             DH_1DP_2 = -DG_1impDp_2/2;
            
            e = sign((E(2:end).*(J_Elec(2:end)+J_Hole(2:end)) + E(1:end-1).*(J_Elec(1:end-1)+J_Hole(1:end-1)))).*(DJpDp(2:end).*E(2:end) + DJpDp2(1:end-1).*E(1:end-1))/2;
            a = [sign((E(3:end).*(J_Elec(3:end)+J_Hole(3:end)) + E(2:end-1).*(J_Elec(2:end-1)+J_Hole(2:end-1)))).*((DJpDp(2:end-1)).*E(2:end-1)/2);0];
%             aa = DHDp3;
            b = [0;sign((E(2:end-1).*(J_Elec(2:end-1)+J_Hole(2:end-1)) + E(1:end-2).*(J_Elec(1:end-2)+J_Hole(1:end-2)))).*(DJpDp2(2:end-1).*E(2:end-1)/2)];
%             bb = DH_1DP_2/2;
            M1=spdiags([a e b],-1:1,Num-2,Num-2); 
            
            
%             DHDp = DJpDp(2:end);
% %             DHDp2 = DJpDp2*self.fd.diffmp(w) - (DG_impDp2 + DG_imp2Dp2)/2 + (DUDp(2:end+1) + Dre_auDp(2:end+1) + Dre_radDp(2:end+1))/2;
% %             DHDp_1 = DG_impDp_1/2;
%             
%             DH_1Dp = DJpDp2(1:end-1);
% %             DH_1Dp2 = - (DG_impDp2)/2;
%             DH_1Dp_1 = DJpDp2(2:end-1);
%             
% %             DHDp3 = -DG_imp2Dp3/2;
% %             DH_1DP_2 = -DG_1impDp_2/2;
            
            e = sign((E(2:end).*(J_Elec(2:end)+J_Hole(2:end)) + E(1:end-1).*(J_Elec(1:end-1)+J_Hole(1:end-1)))).*((DJnDn(2:end).*E(2:end) + DJnDn2(1:end-1).*E(1:end-1))/2);
            a = [sign((E(3:end).*(J_Elec(3:end)+J_Hole(3:end)) + E(2:end-1).*(J_Elec(2:end-1)+J_Hole(2:end-1)))).*((DJnDn(2:end-1)).*E(2:end-1)/2);0];
%             aa = DHDp3;
            b = [0;sign((E(2:end-1).*(J_Elec(2:end-1)+J_Hole(2:end-1)) + E(1:end-2).*(J_Elec(1:end-2)+J_Hole(1:end-2)))).*(DJnDn2(2:end-1).*E(2:end-1)/2)];
%             bb = DH_1DP_2/2;
            M2=spdiags([a e b],-1:1,Num-2,Num-2);    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             DHDw = (J_Hole + J_Elec).*self.fd.differ;
%             DHDw_1 = DG_impDwm/2;
%             DHDw2 = DG_impDwp/2 + diffmp(w)*(DJpDw2 + DJnDw2) + (Jp + Jn)*self.fd.differ;
%             DHDw3 = DG_impDw3/2;
%             
%             DH_1Dw = (DG_1impDw)/2 + diffmp(w)*(DJpDw2 + DJnDw2) + (Jp + Jn)*self.fd.differ;
%             DH_1Dw_2 = (DG_1impDwmm)/2;
%             DH_1Dw2 = DG_impDwp/2;
%             DH_1Dw_1 = DG_impDwm/2 + diffmp(w)*(DJpDw + DJnDw) - (Jp + Jn)*self.fd.differ;
            denom = dx;
            % Full Joule heating without absolute value
%             e = (0.5*(J_Hole(2:end)./denom(2:end)-J_Hole(1:end-1)./denom(1:end-1)) + 0.5*(J_Elec(2:end)./denom(2:end)-J_Elec(1:end-1)./denom(1:end-1))...        
%                 + 0.5*E(2:end).*(DJpDw(2:end)+DJnDw(2:end)) + 0.5*E(1:end-1).*(DJpDw2(1:end-1)+DJnDw2(1:end-1)));
%             a = [(0.5*J_Hole(2:end-1)./denom(2:end-1) + 0.5*J_Elec(2:end-1)./denom(2:end-1)... 
%                 + 0.5*E(2:end-1).*(DJpDw(2:end-1)+DJnDw(2:end-1)));0];
%             b = [0;(-0.5*J_Hole(2:end-1)./denom(2:end-1) - 0.5*J_Elec(2:end-1)./denom(2:end-1)... 
%                 + 0.5*E(2:end-1).*(DJpDw2(2:end-1)+DJnDw2(2:end-1)))]; 
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            e = sign((E(2:end).*(J_Elec(2:end)+J_Hole(2:end)) + E(1:end-1).*(J_Elec(1:end-1)+J_Hole(1:end-1)))).*...
                (0.5*(J_Hole(2:end)./denom(2:end)-J_Hole(1:end-1)./denom(1:end-1)) + 0.5*(J_Elec(2:end)./denom(2:end)-J_Elec(1:end-1)./denom(1:end-1))...        
                + 0.5*E(2:end).*(DJpDw(2:end)+DJnDw(2:end)) + 0.5*E(1:end-1).*(DJpDw2(1:end-1)+DJnDw2(1:end-1)));
            
            a = [sign((E(3:end).*(J_Elec(3:end)+J_Hole(3:end)) + E(2:end-1).*(J_Elec(2:end-1)+J_Hole(2:end-1)))).*...
                (0.5*J_Hole(2:end-1)./denom(2:end-1) + 0.5*J_Elec(2:end-1)./denom(2:end-1)... 
                + 0.5*E(2:end-1).*(DJpDw(2:end-1)+DJnDw(2:end-1)));0];
            
            b = [0;sign((E(2:end-1).*(J_Elec(2:end-1)+J_Hole(2:end-1)) + E(1:end-2).*(J_Elec(1:end-2)+J_Hole(1:end-2)))).*(-0.5*J_Hole(2:end-1)./denom(2:end-1) - 0.5*J_Elec(2:end-1)./denom(2:end-1)... 
                + 0.5*E(2:end-1).*(DJpDw2(2:end-1)+DJnDw2(2:end-1)))]; 
            
            % Full Joule heating with absolute value
%             e = sign((E(2:end).*(J_Elec(2:end)+J_Hole(2:end)) + E(1:end-1).*(J_Elec(1:end-1)+J_Hole(1:end-1)))).*...
%                 (0.5*(J_Hole(2:end)./denom(2:end)-J_Hole(1:end-1)./denom(1:end-1)) + 0.5*(J_Elec(2:end)./denom(2:end)-J_Elec(1:end-1)./denom(1:end-1))...        
%                 + 0.5*E(2:end).*(DJpDw(2:end)+DJnDw(2:end)) + 0.5*E(1:end-1).*(DJpDw2(1:end-1)+DJnDw2(1:end-1)));
%             a = [sign((E(3:end).*(J_Elec(3:end)+J_Hole(3:end)) + E(2:end-1).*(J_Elec(2:end-1)+J_Hole(2:end-1))))...
%                 .*(0.5*J_Hole(2:end-1)./denom(2:end-1) + 0.5*J_Elec(2:end-1)./denom(2:end-1) + 0.5*E(2:end-1).*(DJpDw(2:end-1)+DJnDw(2:end-1)));0]; 
%             b = [0;sign(((E(2:end-1).*(J_Elec(2:end-1)+J_Hole(2:end-1)) + E(1:end-2).*(J_Elec(1:end-2)+J_Hole(1:end-2)))))...
%                 .*(-0.5*J_Hole(2:end-1)./denom(2:end-1) - 0.5*J_Elec(2:end-1)./denom(2:end-1) + 0.5*E(2:end-1).*(DJpDw2(2:end-1)+DJnDw2(2:end-1)))];
            %%%%%%%%%%%%%%%%%%%%%%%%%
%             e = 0.5*(1./denom(2:end).*(J_Elec(2:end))-1./denom(1:end-1).*(J_Elec(1:end-1)))... 
%                 + 0.5*E(2:end).*(DJnDw(2:end)) + 0.5*E(1:end-1).*(DJnDw2(1:end-1));
%             a = [0.5*1./denom(2:end-1).*(J_Elec(2:end-1)) + 0.5*E(2:end-1).*(DJnDw(2:end-1));0];
%             b = [0;- 0.5*1./denom(2:end-1).*(J_Elec(2:end-1)) + 0.5*E(2:end-1).*(DJnDw2(2:end-1))]; 
            %%%%%%%%%%%%%%%%%%%%%%%%%
%             e = 0.5*(1./denom(2:end).*(J_Elec(2:end))-1./denom(1:end-1).*(J_Elec(1:end-1)))... 
%                 + 0.5*E(2:end).*(DJnDw(2:end)) + 0.5*E(1:end-1).*(DJnDw2(1:end-1));
% %             e = 0.5*E(2:end).*(DJnDw(2:end)) + 0.5*E(1:end-1).*(DJnDw2(1:end-1));
%             a = [0.5*1./denom(2:end-1).*(J_Elec(2:end-1)) + 0.5*E(2:end-1).*(DJnDw(2:end-1));0];
%             b = [0;- 0.5*1./denom(2:end-1).*(J_Elec(2:end-1)) + 0.5*E(2:end-1).*(DJnDw2(2:end-1))];
            %%%%%%%%%%%%%%%%%%%%%%%%%
%             e = (- DHDw(2:end) + DHDw(1:end-1))/2;
%             a = [-DHDw(2:end-1)/2;0];
% %             aa = DHDw3;
%             b = [0;DHDw(2:end-1)/2];
%             bb = DH_1Dw_2/2;
            M3=spdiags([a e b],-1:1,Num-2,Num-2); 
        %%%%%%%%%%%%%%%%%%%%%%%%%
%             e = (k_L./(dx1(1:end))).*(-1./dx(2:end) - 1./dx(1:end-1));
%             a = [(k_L./dx1(2:end)).*(1./dx(2:end-1));0];
%             b = [0;k_L./(dx1(1:end-1)).*(1./dx(2:end-1))];
%             denom = self.fd.differ;
            e = (k_L(2:end-1)./(dx1(1:end))).*(-1./dx(2:end) - 1./dx(1:end-1)) + sign((E(2:end).*(J_Elec(2:end)+J_Hole(2:end)) + E(1:end-1).*(J_Elec(1:end-1)+J_Hole(1:end-1)))).*(0.5.*E(2:end).*(DJnDT(2:end)+DJpDT2(2:end)) + 0.5.*E(1:end-1).*(DJnDT2(1:end-1)+DJpDT(1:end-1)));
            a = [(k_L(3:end-1)./dx1(2:end)).*(1./dx(2:end-1)) + sign((E(3:end).*(J_Elec(3:end)+J_Hole(3:end)) + E(2:end-1).*(J_Elec(2:end-1)+J_Hole(2:end-1)))).*(0.5.*E(2:end-1).*(DJnDT(2:end-1) + DJpDT2(2:end-1)));0];
            b = [0;k_L(2:end-2)./(dx1(1:end-1)).*(1./dx(2:end-1)) + sign((E(2:end-1).*(J_Elec(2:end-1)+J_Hole(2:end-1)) + E(1:end-2).*(J_Elec(1:end-2)+J_Hole(1:end-2)))).*(0.5.*E(2:end-1).*(DJnDT2(2:end-1) + DJpDT(2:end-1)))];
            %%%%%%%%%%%%%%%%%%%%%%%%%%
%             e = (k_L./(dx1(1:end))).*(-1./dx(2:end) - 1./dx(1:end-1)) + 0.5.*E(2:end).*(DJnDT(2:end)) + 0.5.*E(1:end-1).*(DJnDT2(1:end-1));
%             a = [(k_L./dx1(2:end)).*(1./dx(2:end-1)) + 0.5.*E(2:end-1).*(DJnDT(2:end-1));0];
%             b = [0;k_L./(dx1(1:end-1)).*(1./dx(2:end-1)) + 0.5.*E(2:end-1).*(DJnDT2(2:end-1))];
            
            M4=spdiags([a e b],-1:1,Num-2,Num-2); 
            MFt=[M1,M2,M3,M4];
            JacobM = [MFp;MFn;MFw;MFt];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [y,M]=Newton_FD(self,x)
            N=length(self.mesh.dx)+1;
            pm = x(1:N-2);
            nm = x(N-1:2*N-4);
            wm = x(2*N-3:3*N-6);
            tm = x(3*N - 5 : 4*N - 8);
            
            [Fp,Fn,Fw,Ft]=self.Cal_Current(pm,nm,wm,tm,self.GL);
            y = [Fp;Fn;Fw;Ft];
            M = self.ConsJacobian(pm,nm,wm,tm);
        end
         function [y]=Newton_FDT(self,x)
            N=length(self.mesh.dx)+1;
            pm = x(1:N-2);
            nm = x(N-1:2*N-4);
            wm = x(2*N-3:3*N-6);
            tm = x(3*N - 5 : 4*N - 8);
            
            [Fp,Fn,Fw,Ft]=self.Cal_Current(pm,nm,wm,tm,self.GL);
            y = [Fp;Fn;Fw;Ft];
%             M = self.ConsJacobian(pm,nm,wm,tm);
        end

    end    
end

