classdef Drift_Diff_N_B
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
    
    methods
        function self=Drift_Diff_N_B(pb1,pbn,nb1,nbn,wb1,wbn,fd,mesh,PD_Str,NP,Semiconductor_Device,GL,wavelength)
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
        
        function [Fp,Fn,Fw,J_Hole,J_Elec,alpha]=Cal_Current(self,pm,nm,phim,GL)
            % calculate the function for drift&diffusion equation and
            % possion equation
            % the following are the parameters for the device.
            q = 1.602e-19;
            me = 9.1e-31;
            hbar = 1.05457e-34;
            TL = self.PD_Str.TL;
            NC = self.PD_Str.NC*self.PD_Str.Ni;
            NV = self.PD_Str.NV*self.PD_Str.Ni;
            Xi = self.PD_Str.Xi;
            Eg = self.PD_Str.Eg;
            mstar_n = self.PD_Str.mstar_n;
            mstar_p = self.PD_Str.mstar_p;
            VT=self.Semiconductor_Device.NP.VT;
            NX=self.NP.NX;
            KB = self.NP.KB;
            Nmu = self.NP.Nmu;
            
            p=[self.pb1;pm;self.pbn];
            n=[self.nb1;nm;self.nbn];
            phi=[self.wb1;phim;self.wbn];
            N1=self.mesh.N1;
            N=self.mesh.N;
            In = self.PD_Str.In;
            tau_n = self.Semiconductor_Device.Lifetime_n;
            ni = self.Semiconductor_Device.Intrinsic_Density;
            Doping = self.Semiconductor_Device.Doping_Profile;
            n1 = self.Semiconductor_Device.n1;
            p1 = self.Semiconductor_Device.p1;
            % incomplete ionization model
            Doping_ion=zeros(N,1);
            Doping_ion(1:N1)=Doping(1:N1)./(1+4*p(1:N1)./p1(1:N1));
            Doping_ion(N1+1:end)=Doping(N1+1:end)./(1+2*n(N1+1:end)./n1(N1+1:end));
            %T_Vol = self.Semiconductor_Device.T_Vol;
            Electron_Mobility = self.Semiconductor_Device.Electron_Mobility;
            Hole_Mobility = self.Semiconductor_Device.Hole_Mobility;
            un=Electron_Mobility;
            up=Hole_Mobility;
            %%%%%%%%%%%%
            E = -self.fd.diffmp(phi);% electrical field
            %
            Electron_Velocity = E.*Electron_Mobility(2:end); % Electron Velocity
            Hole_Velocity = E.*Hole_Mobility(2:end); % Hole Velocity
            %
            Recombination=(n.*p)./(tau_n.*(n+p+2*ni)); % Recombination-ni.*ni
            %
            alpha=Absorp_FK(abs([E(1);E])*VT/NX*100,self.wavelength,self.PD_Str,self.Semiconductor_Device.FranzKeldysh);
            absorbers = ones(length(alpha),1);
            absorbers(alpha==0) = 0;
            cell_loss = exp(-alpha.*NX.*[self.mesh.dx; 0]);
            structure_loss = flipud(cumprod(flipud(cell_loss)));
            Generation = GL.*structure_loss.*absorbers;
            E_main = -self.fd.diff3p(phi); % the electric field at the main point
            if self.Semiconductor_Device.Impact_Ion
                Ae = self.Semiconductor_Device.Ae.*self.Semiconductor_Device.NP.NX;
                Be = self.Semiconductor_Device.Be.*self.Semiconductor_Device.NP.NX/self.Semiconductor_Device.NP.VT;
                Ah = self.Semiconductor_Device.Ah.*self.Semiconductor_Device.NP.NX;
                Bh = self.Semiconductor_Device.Bh.*self.Semiconductor_Device.NP.NX/self.Semiconductor_Device.NP.VT;
                                
                Vn_main = E_main.*(un(2:end-1)); % Electron Velocity at the main point
                Vp_main = E_main.*(up(2:end-1)); % Hole Velocity at the main point
                Dnz_main = un(2:end-1);
                Dpz_main = up(2:end-1);
                alpha_An = Ae(2:end-1).*exp(-(Be(2:end-1)./abs(E_main)).^In);
                alpha_Ap = Ah(2:end-1).*exp(-Bh(2:end-1)./abs(E_main));
                re = 4e-6./self.NP.NX;
                i1 = self.mesh.N1;
                i2 = self.mesh.N2;
                for ii = i1:i2
                    f_ce = 2/(sqrt(pi)*re)*exp(-(self.mesh.Lx(2:ii)-self.mesh.Lx(ii)).^2./re^2);
                    a_temp=f_ce.*E_main(1:ii-1).*self.mesh.dx(1:ii-1);
                    E_eff=sum(a_temp);
                    alpha_An(ii) = Ae(ii).*exp(-(Be(ii)./abs(E_eff)).^In);
                    alpha_Ap(ii) = Ah(ii).*exp(-(Bh(ii)./abs(E_eff)).^In);
                end
                G_imp = alpha_An.*abs(n(2:end-1).*abs(Vn_main)+Dnz_main.*self.fd.diff3p(n(1:end)));
                G_imp = alpha_Ap.*abs(p(2:end-1).*abs(Vp_main)+Dpz_main.*self.fd.diff3p(p(1:end)))+G_imp;
                G_imp(i2+1:end)=0;
                G_imp(1:i1)=0;
            else
                G_imp=0;
            end
            J_Hole_drift = Hole_Velocity.*self.fd.averagemp(p);
            J_Hole_diff  = - Hole_Mobility(2:end).*self.fd.diffmp(p);
            J_Elec_drift = Electron_Velocity.*self.fd.averagemp(n);
            J_Elec_diff = Electron_Mobility(2:end).*self.fd.diffmp(n);
            
            hetjunctions = find(diff(self.PD_Str.what_sc)~=0);
            nn = max([abs(self.PD_Str.Tol(self.mesh.Node(hetjunctions))) abs(self.PD_Str.Tol(self.mesh.Node(hetjunctions)+1))].').';
            epn = min([(self.PD_Str.Epsilons(self.mesh.Node(hetjunctions))) (self.PD_Str.Epsilons(self.mesh.Node(hetjunctions)+1))].').';
            mes = min([(self.PD_Str.mes(self.mesh.Node(hetjunctions))) (self.PD_Str.mes(self.mesh.Node(hetjunctions)+1))].').';
            mps = min([(self.PD_Str.mes(self.mesh.Node(hetjunctions))) (self.PD_Str.mes(self.mesh.Node(hetjunctions)+1))].').';
            E00n = hbar/2*sqrt(nn./epn./mes);
            E00p = hbar/2*sqrt(nn./epn./mps);
            hetjunctions1 = hetjunctions(KB*self.NP.T>E00n);
            hetjunctions2 = hetjunctions(KB*self.NP.T>E00p);
            if self.Semiconductor_Device.Thermionic_e ==1 &&  isempty(hetjunctions1) ==0
                for i56 = 1:length(hetjunctions1)
                    yy = hetjunctions1(i56);
                    he=self.mesh.Node(yy);
                    
                    me_min = min([mstar_n(yy) mstar_n(yy+1)])*me;
                    Astar=1e-4*me_min*KB^2/(2*pi^2*hbar^3);
                    vn1=Astar*TL(he)^2/NC(he)*NX/VT/Nmu;
                    vn2=Astar*TL(he+1)^2/NC(he+1)*NX/VT/Nmu;
                    Ec1_Ec2 = Xi(he)-Xi(he+1);
                    m_factor = exp(-abs(Ec1_Ec2)/8.617333262e-5/self.PD_Str.TL(he));
                    if Ec1_Ec2 <0
                        J_Elec_diff(he)=J_Elec_diff(he)-vn1*n(he)+vn2*n(he+1)*m_factor;
                    else
                        J_Elec_diff(he)=J_Elec_diff(he)+vn1*n(he)*m_factor-vn2*n(he+1);
                    end
                end
            end
            if self.Semiconductor_Device.Thermionic_h ==1 && isempty(hetjunctions2) ==0
                for i56 = 1:length(hetjunctions2)
                    yy = hetjunctions2(i56);
                    he=self.mesh.Node(yy);
                    mh_min = min([mstar_p(yy) mstar_p(yy+1)])*me;
                    Astar=1e-4*mh_min*KB^2/(2*pi^2*hbar^3);
                    vp1=Astar*TL(he)^2/NV(he)*NX/VT/Nmu;
                    vp2=Astar*TL(he+1)^2/NV(he+1)*NX/VT/Nmu;
                    
                    Ev1_Ev2 = Eg(he)+Xi(he)-Eg(he+1)-Xi(he+1);
                    m_factor = exp(-abs(Ev1_Ev2)/8.617333262e-5/self.PD_Str.TL(he));
                    if Ev1_Ev2 <0
                        J_Hole_diff(he)= J_Hole_diff(he)+vp1*p(he)*m_factor-vp2*p(he+1);
                    else
                        J_Hole_diff(he)= J_Hole_diff(he)-vp1*p(he)+vp2*p(he+1)*m_factor;
                    end
                end
            end
            J_Hole = J_Hole_drift + J_Hole_diff;
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [JacobM]=ConsJacobian(self,pm,nm,phim)
            p=[self.pb1;pm;self.pbn];
            n=[self.nb1;nm;self.nbn];
            phi=[self.wb1;phim;self.wbn];
            hbar=1.05457e-34;
            q = 1.602e-19;
            me = 9.1e-31;
            TL = self.PD_Str.TL;
            NC = self.PD_Str.NC*self.PD_Str.Ni;
            NV = self.PD_Str.NV*self.PD_Str.Ni;
            Xi = self.PD_Str.Xi;
            Eg = self.PD_Str.Eg;
            mstar_n = self.PD_Str.mstar_n;
            mstar_p = self.PD_Str.mstar_p;
            VT=self.Semiconductor_Device.NP.VT;
            NX=self.NP.NX;
            KB = self.NP.KB;
            Nmu = self.NP.Nmu;
            
            tau_n = self.Semiconductor_Device.Lifetime_n;
            tau_p = self.Semiconductor_Device.Lifetime_p;
            ni = self.Semiconductor_Device.Intrinsic_Density;
            Doping = self.Semiconductor_Device.Doping_Profile;
            n1=self.Semiconductor_Device.n1;
            p1=self.Semiconductor_Device.p1;
            N=self.mesh.N;
            N1=self.mesh.N1;
            Vps=self.Semiconductor_Device.P_Vsaturate;
            Vns=self.Semiconductor_Device.N_Vsaturate;
            beta=self.Semiconductor_Device.Beta;
            Ep=self.Semiconductor_Device.Ep;
            gamma=self.Semiconductor_Device.gamma;
            In = self.PD_Str.In;
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
            un=Electron_Mobility;
            up=Hole_Mobility;
            E = -self.fd.diffmp(phi);
            %
            Recombination=(n.*p)./(tau_n.*(n+p+2*ni));
            DUDp=(n-tau_n.*Recombination)./...
                (tau_n.*(p+ni)...
                +tau_p.*(n+ni));
            % first derivative  Recombination(i) to Hole_Density(i)
            DUDn=(p-tau_p.*Recombination)./...
                (tau_n.*(p+ni)+tau_p.*(n+ni));
            %           % the impact ionization generation
            if  self.Semiconductor_Device.Impact_Ion
                Ae = self.Semiconductor_Device.Ae*self.Semiconductor_Device.NP.NX;
                Be = self.Semiconductor_Device.Be*self.Semiconductor_Device.NP.NX/self.Semiconductor_Device.NP.VT;
                Ah = self.Semiconductor_Device.Ah*self.Semiconductor_Device.NP.NX;
                Bh = self.Semiconductor_Device.Bh*self.Semiconductor_Device.NP.NX/self.Semiconductor_Device.NP.VT;
                
                
                E_main = -self.fd.diff3p(phi); % the electric field at the main point
                Vn_main = E_main.*(un(2:end-1)); % Electron Velocity at the main point
                Vp_main = E_main.*(up(2:end-1)); % Hole Velocity at the main point
                Dnz_main = un(2:end-1);
                Dpz_main = (Vps(2:end-1).*up(2:end-1)./((Vps(2:end-1).^gamma+up(2:end-1).^gamma.*(abs(E_main).^gamma)).^(1/gamma))); % Hole Velocity at the main point
                alpha_An = Ae(2:end-1).*exp(-(Be(2:end-1)./abs(E_main)).^In);
                alpha_Ap = Ah(2:end-1).*exp(-Bh(2:end-1)./abs(E_main));
                % the derive of G_imp to n,p and w
                DVnmainDE = (2*Vns(2:end-1).*beta(2:end-1).*abs(E_main) - un(2:end-1).*beta(2:end-1).*E_main.^2+un(2:end-1))./(1+beta(2:end-1).*E_main.^2).^2;
                DDnzDE =  un(2:end-1).*E_main.*(abs(Ep(1:end-1))-abs(E_main))./abs(Ep(1:end-1)).^3./(1-2*E_main.^2./Ep(1:end-1).^2+4/3*E_main.^2.*abs(E_main)./abs(Ep(1:end-1)).^3).^1.25;
                %
                DVpmainDE = up(2:end-1).*Vps(2:end-1).^(1+gamma).*(up(2:end-1).^gamma.*abs(E_main).^gamma+Vps(2:end-1)).^(-1/gamma-1);
                DDpzDE = -up(2:end-1).^(1+gamma).*Vps(2:end-1).*abs(E_main).^gamma./E_main.*(up(2:end-1).^gamma.*abs(E_main).^gamma+Vps(2:end-1).^gamma).^(-1/gamma-1);
                
                theta = n(2:end-1).*Vn_main+Dnz_main.*self.fd.diff3p(n(1:end));
                DG_impDn = alpha_An.*(Vn_main+Dnz_main.*self.fd.diff2p(n(1:end))).*sign(theta);
                theta = p(2:end-1).*Vp_main-Dpz_main.*self.fd.diff3p(p(1:end));
                DG_impDp = alpha_Ap.*(Vp_main-Dnz_main.*self.fd.diff2p(p(1:end))).*sign(theta);
                
                theta = n(2:end-1).*Vn_main+Dnz_main.*self.fd.diff3p(n(1:end));
                theta2 = p(2:end-1).*Vp_main-Dpz_main.*self.fd.diff3p(p(1:end));
                DG_impDE = alpha_An.*(n(2:end-1).*DVnmainDE+DDnzDE.*self.fd.diff3p(n(1:end))).*sign(theta)+...
                    alpha_Ap.*(p(2:end-1).*DVpmainDE-DDpzDE.*self.fd.diff3p(p(1:end))).*sign(theta2);
                
                DG_impDwm = DG_impDE./(self.mesh.dx(1:end-1)+self.mesh.dx(2:end));% DG_imp to w_i-1
                DG_impDwp = -DG_impDwm;
                DG_impDn(self.mesh.N2+1:end)=0;
                DG_impDp(self.mesh.N2+1:end)=0;
                DG_impDwm(self.mesh.N2+1:end)=0;
                DG_impDwp(self.mesh.N2+1:end)=0;
                DG_impDn(1:self.mesh.N1)=0;
                DG_impDp(1:self.mesh.N1)=0;
                DG_impDwm(1:self.mesh.N1)=0;
                DG_impDwp(1:self.mesh.N1)=0;
            else
                DG_impDn=0;
                DG_impDp=0;
                DG_impDwm=zeros(N-2,1);
                DG_impDwp=zeros(N-2,1);
            end
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
            
            hetjunctions = find(diff(self.PD_Str.what_sc)~=0);
            nn = max([abs(self.PD_Str.Tol(self.mesh.Node(hetjunctions))) abs(self.PD_Str.Tol(self.mesh.Node(hetjunctions)+1))].').';
            epn = min([(self.PD_Str.Epsilons(self.mesh.Node(hetjunctions))) (self.PD_Str.Epsilons(self.mesh.Node(hetjunctions)+1))].').';
            mes = min([(self.PD_Str.mes(self.mesh.Node(hetjunctions))) (self.PD_Str.mes(self.mesh.Node(hetjunctions)+1))].').';
            mps = min([(self.PD_Str.mes(self.mesh.Node(hetjunctions))) (self.PD_Str.mes(self.mesh.Node(hetjunctions)+1))].').';
            E00n = hbar/2*sqrt(nn./epn./mes);
            E00p = hbar/2*sqrt(nn./epn./mps);
            hetjunctions1 = hetjunctions(KB*self.NP.T>E00n);
            hetjunctions2 = hetjunctions(KB*self.NP.T>E00p);
            if self.Semiconductor_Device.Thermionic_e ==1 && isempty(hetjunctions1) == 0
                for i56 = 1:length(hetjunctions1)
                    yy = hetjunctions1(i56);
                    he=self.mesh.Node(yy);
                    me_min = min([mstar_n(yy) mstar_n(yy+1)])*me;
                    Astar=1e-4*me_min*KB^2/(2*pi^2*hbar^3);
                    vn1=Astar*TL(he)^2/NC(he)*NX/VT/Nmu;
                    vn2=Astar*TL(he+1)^2/NC(he+1)*NX/VT/Nmu;
                    Ec1_Ec2 = Xi(he)-Xi(he+1);
                    if Ec1_Ec2 <0
                        m_factor = exp(Ec1_Ec2/8.617333262e-5/TL(he+1));
                        DJnDn(he)=-vn1;
                        DJnDn2(he)=vn2*m_factor;
                    else
                        m_factor = exp(-Ec1_Ec2/8.617333262e-5/TL(he));
                        DJnDn(he)=vn1*m_factor;
                        DJnDn2(he)=-vn2;
                    end
                end
            end
            if self.Semiconductor_Device.Thermionic_h ==1 && isempty(hetjunctions2) == 0
                for i56 = 1:length(hetjunctions2)
                    yy = hetjunctions2(i56);
                    he=self.mesh.Node(yy);
                    mh_min = min([mstar_p(yy) mstar_p(yy+1)])*me;
                    Astar=1e-4*mh_min*KB^2/(2*pi^2*hbar^3);
                    vp1=Astar*TL(he)^2/NV(he)*NX/VT/Nmu;
                    vp2=Astar*TL(he+1)^2/NV(he+1)*NX/VT/Nmu;
                    Ev1_Ev2 = Eg(he)+Xi(he)-Eg(he+1)-Xi(he+1);
                    m_factor = exp(-abs(Ev1_Ev2)/8.617333262e-5/TL(he));
                    if Ev1_Ev2 <0
                        DJpDp(he)= DJpDp(he) + vp1*m_factor;
                        DJpDp2(he)= DJpDp2(he)-vp2;
                    else
                        DJpDp(he)= DJpDp(he)-vp1;
                        DJpDp2(he)= DJpDp(he)+vp2*m_factor;
                    end
                end
            end
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

