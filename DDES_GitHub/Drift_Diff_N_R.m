classdef Drift_Diff_N_R
    %DRIFT_DIFF Summary of this class goes here
    %   Detailed explanation goes here
    %8/16/2011,add generation term
    % add potential at the p side as a variable 3/22/2015
    properties
        wavelength
        pb1% boundary condition for hole
        pbn% boundary condition(end) for hole
        nb1% boundary condition for electron
        nbn% boundary condition(end) for electron
        wbf% boundary condition for potential at zero bias
        wbn% boundary condition(end) for potential
        fd
        NP
        mesh
        PD_Str
        Semiconductor_Device
        GL % generation rate
        Bias
        R_Load
        DeviceDiameter        
    end
    
    methods
        function self=Drift_Diff_N_R(pb1,pbn,nb1,nbn,wbf,wbn,fd,mesh,PD_Str,NP,Semiconductor_Device,GL,Bias,R_Load,DeviceDiameter,wavelength)
            self.pb1=pb1;
            self.wavelength= wavelength;
            self.pbn=pbn;
            self.nb1=nb1;
            self.nbn=nbn;
            self.wbf=wbf;
            self.wbn=wbn;
            self.fd = fd;
            self.mesh=mesh;
            self.PD_Str=PD_Str;
            self.NP=NP;
            self.Semiconductor_Device = Semiconductor_Device;
            self.GL=GL;
            self.Bias=Bias;
            self.R_Load=R_Load;
            self.DeviceDiameter=DeviceDiameter;
            
        end
        
        function [Fp,Fn,Fw,Fb,J_Hole,J_Elec,alpha]=Cal_Current(self,pm,nm,phim,wb1,GL)
            p=[self.pb1;pm;self.pbn];
            n=[self.nb1;nm;self.nbn];
            phi=[wb1;phim;self.wbn];
            N1=self.mesh.N1;
            N=self.mesh.N;
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
            %
            Electron_Velocity = E.*Electron_Mobility(2:end); % Electron Velocity
            Hole_Velocity = E.*Hole_Mobility(2:end); % Hole Velocity
            Recombination=(n.*p)./(tau_n.*(n+p+2*ni)); % Recombination-ni.*ni
            
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
            G_imp=0;            
            J_Hole = Hole_Velocity.*self.fd.averagemp(p)...% Hole current,there should be a q times, but it cancelled in the Fp.
                - Hole_Mobility(2:end).*self.fd.diffmp(p);
            J_Elec_drift = Electron_Velocity.*self.fd.averagemp(n);
            J_Elec_diff = Electron_Mobility(2:end).*self.fd.diffmp(n);
            J_Elec = J_Elec_drift + J_Elec_diff; % electron current
            J_T=J_Hole+J_Elec;
            Fp = -self.fd.diffmidp(J_Hole)-Recombination(2:end-1)+G_imp+Generation(2:end-1); % drift % diffusion current equation for hole
            Fn = self.fd.diffmidp(J_Elec)-Recombination(2:end-1)+G_imp+Generation(2:end-1);% drift % diffusion current equation for electron
            dx=self.fd.mesh.dx;
            dx1=self.fd.mesh.dx1;
            Fw = (Doping_ion(2:end-1)+p(2:end-1)-n(2:end-1))... % ni(2:end-1)is norm
                +(1./((dx(1:end-1)).*(dx1(1:end)))).*phi(1:end-2)...
                -((1./dx(1:end-1)+1./dx(2:end))./dx1(1:end)).*phi(2:end-1)...
                +(1./((dx(2:end)).*(dx1(1:end)))).*phi(3:end);% possion equation
            
            Jave=mean(J_T(1))*self.NP.J0*pi*(self.DeviceDiameter(1))^2*0.25;
            Fb = self.wbf-wb1-(self.Bias+(Jave)*self.R_Load)/self.NP.VT;
        end
        function [JacobM]=ConsJacobian(self,pm,nm,phim,wb1)
            p=[self.pb1;pm;self.pbn];
            n=[self.nb1;nm;self.nbn];
            phi=[wb1;phim;self.wbn];
            
            tau_n = self.Semiconductor_Device.Lifetime_n;
            tau_p = self.Semiconductor_Device.Lifetime_p;
            ni = self.Semiconductor_Device.Intrinsic_Density;
            Doping = self.Semiconductor_Device.Doping_Profile;
            n1=self.Semiconductor_Device.n1;
            p1=self.Semiconductor_Device.p1;
            N=self.mesh.N;
            N1=self.mesh.N1;
            % incomplete ionization model
            Doping_ion=zeros(N,1);
            Doping_ion(1:N1)=Doping(1:N1)./(1+4*p(1:N1)./p1(1:N1));
            Doping_ion(N1+1:end)=Doping(N1+1:end)./(1+2*n(N1+1:end)./n1(N1+1:end));
            % Doing_ion  derive to the p and n
            Doping_ion_p=zeros(N,1);
            Doping_ion_n=zeros(N,1);
            Doping_ion_p(1:N1)=-Doping_ion(1:N1)./(p1(1:N1)/4+p(1:N1));
            Doping_ion_n(N1+1:end)=-Doping_ion(N1+1:end)./(n1(N1+1:end)/2+n(N1+1:end));
            
            Electron_Mobility = self.Semiconductor_Device.Electron_Mobility;
            Hole_Mobility = self.Semiconductor_Device.Hole_Mobility;
            E = -self.fd.diffmp(phi);
%            
            Recombination=(n.*p)./(tau_n.*(n+p+2*ni));
            DUDp=(n-tau_n.*Recombination)./...
                (tau_n.*(p+ni)...
                +tau_p.*(n+ni));
            % first derivative  Recombination(i) to Hole_Density(i)
            DUDn=(p-tau_p.*Recombination)./...
                (tau_n.*(p+ni)+tau_p.*(n+ni));
            %             DUDn=n*0;
            %           % the impact ionization generation
            % G_imp=alph_An*n*abs(vn)+alpha_Ap*p*abs(vp)
            % this vn and vp should be on the main point;

                DG_impDn=0;
                DG_impDp=0;
                DG_impDwm=zeros(N-2,1);
                DG_impDwp=zeros(N-2,1);
                
            Electron_Velocity=E.*Electron_Mobility(2:end); % Electron Velocity
            Hole_Velocity=E.*Hole_Mobility(2:end); % Hole Velocity
            % J_Hole and J_ Electron first derivative to p ,n and phi
            DJpDp = Hole_Velocity/2+Hole_Mobility(2:end).*self.fd.differ;% first derivative hole current J_Hole(i) to Hole_Density(i)
            DJpDp2 = Hole_Velocity/2-Hole_Mobility(2:end).*self.fd.differ;% first derivative J_Hole(i) to Hole_Density(i+1)
            
            DJpDw = Hole_Mobility(2:end).*self.fd.averagemp(p).*self.fd.differ; % first derivative J_Hole(i) to Potential(i)
            DJpDw2 = -DJpDw; %  % first derivative J_Hole(i) to Potential(i+1)            
            %
            DJnDn = Electron_Velocity/2-Electron_Mobility(2:end).*self.fd.differ; % % first derivative electron current J_Elec(i) to Elec_Density(i)
            DJnDn2 = Electron_Velocity/2+Electron_Mobility(2:end).*self.fd.differ; % first derivative current J_Elec(i) to Elec_Density(i+1)
            DJnDw = Electron_Mobility(2:end).*(self.fd.averagemp(n)).*self.fd.differ;  % first derivative  J_Elec(i) to potential(i)
            DJnDw2 = -DJnDw; % first derivative  J_Elec(i) to potential(i+1)
            %calculate the current of hole and electrons
            %so the Jacob matrix***DFp/Dp,DFp/Dn,DFp/Dw************
            %**********************DFn/Dp,DFn/Dn,DFn/Dw************
            %**********************DFw/Dp,DFw/Dn,DFw/Dw************
            dx1=self.mesh.dx1;
            dx=self.mesh.dx;
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
            M4=zeros(Num-2,1);
            M4(1)=DJpDw(1)./dx1(1);
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
            M4=zeros(Num-2,1);
            M4(1)=-DJnDw(1)./dx1(1);
            MFn=[M1,M2,M3,M4];
            
            % Fw derivative
            DFphiDp = self.fd.CreateZOM+spdiags(Doping_ion_p(2:end-1),0,Num-2,Num-2);% Fphi first derivative to p
            DFphiDn = -self.fd.CreateZOM+spdiags(Doping_ion_n(2:end-1),0,Num-2,Num-2); % Fphi first derivative to n
            DFphiDphi= self.fd.CreateSOM; % Fp first derivative to phi
            M4=zeros(Num-2,1);
            M4(1)=1./(dx1(1)*dx(1));
            MFb=zeros(1,3*(Num-2)+1);
            MFb(end)=-1-(DJpDw(1)+DJnDw(1))*self.NP.J0*pi*(self.DeviceDiameter(1))^2*0.25*self.R_Load/self.NP.VT;
            MFb(1)=-DJpDp2(1)*self.NP.J0*pi*(self.DeviceDiameter(1))^2*0.25*self.R_Load/self.NP.VT;
            MFb(N-1)=-DJnDn2(1)*self.NP.J0*pi*(self.DeviceDiameter(1))^2*0.25*self.R_Load/self.NP.VT;
            MFb(2*N-3)=-(DJpDw2(1)+DJnDw2(1))*self.NP.J0*pi*(self.DeviceDiameter(1))^2*0.25*self.R_Load/self.NP.VT;
            JacobM = [MFp;MFn;DFphiDp,DFphiDn,DFphiDphi,M4;MFb];
        end
        function [y,M]=Newton_FD(self,x)
            N=length(self.mesh.dx)+1;
            pm = x(1:N-2);
            nm = x(N-1:2*N-4);
            wm = x(2*N-3:3*N-6);
            wb1 = x(end);
            
            [Fp,Fn,Fw,Fb]=self.Cal_Current(pm,nm,wm,wb1,self.GL);
            y = [Fp;Fn;Fw;Fb];
            M = self.ConsJacobian(pm,nm,wm,wb1);
        end
    end
    
end

