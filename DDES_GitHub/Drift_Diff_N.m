classdef Drift_Diff_N
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
        function self=Drift_Diff_N(pb1,pbn,nb1,nbn,wb1,wbn,fd,mesh,PD_Str,NP,Semiconductor_Device,GL,wavelength)
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
        function [Fp,Fn,Fw,J_Hole,J_Elec,alpha]=Cal_Current(self,pm,nm,phim,GL)
            % calculate the function for drift&diffusion equation and
            % possion equation
            % the following are the parameters for the device.
            p=[self.pb1;pm;self.pbn];
            n=[self.nb1;nm;self.nbn];
            phi=[self.wb1;phim;self.wbn];
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
            %%%%%%%%%%%%            
            Electron_Velocity = E.*Electron_Mobility(2:end); % Electron Velocity
            Hole_Velocity = E.*Hole_Mobility(2:end); % Hole Velocity
            %%%%%%%%%%%%            
            Recombination=(n.*p)./(tau_n.*(n+p+2*ni)); % Recombination-ni.*ni            
            VT=self.Semiconductor_Device.NP.VT;
            NX=self.Semiconductor_Device.NP.NX;
            Generation =0*GL;

            E_main = -self.fd.diff3p(phi); % the electric field at the main point
            G_imp=0;                        
            J_Hole = Hole_Velocity.*self.fd.averagemp(p)...% Hole current,there should be a q times, but it cancelled in the Fp.
                - Hole_Mobility(2:end).*self.fd.diffmp(p);
            J_Elec_drift = Electron_Velocity.*self.fd.averagemp(n);
            J_Elec_diff = Electron_Mobility(2:end).*self.fd.diffmp(n);
            J_Elec = J_Elec_drift + J_Elec_diff; % electron current
            Fp = -self.fd.diffmidp(J_Hole)-Recombination(2:end-1)+G_imp+Generation(2:end-1); % drift % diffusion current equation for hole
            Fn = self.fd.diffmidp(J_Elec)-Recombination(2:end-1)+G_imp+Generation(2:end-1);% drift % diffusion current equation for electron
            dx=self.fd.mesh.dx;
            dx1=self.fd.mesh.dx1;
            Fw = (Doping_ion(2:end-1)+p(2:end-1)-n(2:end-1))... % ni(2:end-1)is norm
                +(1./((dx(1:end-1)).*(dx1(1:end)))).*phi(1:end-2)...
                -((1./dx(1:end-1)+1./dx(2:end))./dx1(1:end)).*phi(2:end-1)...
                +(1./((dx(2:end)).*(dx1(1:end)))).*phi(3:end);% possion equation 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [JacobM]=ConsJacobian(self,pm,nm,phim)
            p=[self.pb1;pm;self.pbn];
            n=[self.nb1;nm;self.nbn];
            phi=[self.wb1;phim;self.wbn];            
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
            MFp=[M1,M2,M3]; % Fp first derivative to p ,n and potential
            
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
            MFn=[M1,M2,M3];
            
            % Fn first derivative to phi
            DFphiDp = self.fd.CreateZOM+spdiags(Doping_ion_p(2:end-1),0,Num-2,Num-2);% Fphi first derivative to p
            DFphiDn = -self.fd.CreateZOM+spdiags(Doping_ion_n(2:end-1),0,Num-2,Num-2); % Fphi first derivative to n
            DFphiDphi= self.fd.CreateSOM; % Fp first derivative to phi                        
            JacobM = [MFp;MFn;DFphiDp,DFphiDn,DFphiDphi];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [y,M]=Newton_FD(self,x)
            N=length(self.mesh.dx)+1;
            pm = x(1:N-2);
            nm = x(N-1:2*N-4);
            wm = x(2*N-3:3*N-6);
            
            [Fp,Fn,Fw]=self.Cal_Current(pm,nm,wm,self.GL);
            y = [Fp;Fn;Fw];
            M = self.ConsJacobian(pm,nm,wm);
        end

    end    
end

