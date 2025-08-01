classdef Drift_Diff_N_C_B_R
    properties
        pb1% boundary condition for hole
        wavelength
        pbn% boundary condition(end) for hole
        nb1% boundary condition for electron
        nbn% boundary condition(end) for electron
        wbf% boundary condition for potential without bias
        wbn% boundary condition(end) for potential
        fd
        mesh
        Semiconductor_Device
        NP% normolized parameter class
        GL% light internsity generation
        BF% boundary conditionflag 1 for true using boundary, 0 for no
        % the coeificient for the diffusion current
        PD_Str % PD structure parameters
        TL % temperature
        Tg % temperature 2
        Bias
        R_Load
        Diameter
        dir00
        pulse % input pulse properties
        E_pre % last step e field
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function self=Drift_Diff_N_C_B_R(pb1,pbn,nb1,nbn,wbf,wbn,fd,mesh,...
                Semiconductor_Device,NP,GL,PD_Str1,TL,Bias,...
                R_Load,Diameter,wavelength)
            self.pb1=pb1;
            self.wavelength = wavelength;
            self.pbn=pbn;
            self.nb1=nb1;
            self.nbn=nbn;
            self.wbf=wbf;
            self.wbn=wbn;
            self.fd = fd;
            self.mesh=mesh;
            self.Semiconductor_Device = Semiconductor_Device;
            self.NP=NP;
            self.GL=GL;
            self.PD_Str=PD_Str1;
            self.TL = TL;
            self.Bias = Bias;
            self.R_Load = R_Load;
            self.Diameter = Diameter;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Fp,Fn,Fw,Fb,J_Hole,J_Elec,J_Elec_drift,J_Elec_diff,Doping_ion,ave_alpha,Generation,Recombination,Electron_Velocity,Hole_Velocity...
                ,NAp,NDp]=Cal_Currentv2(self,pm,nm,phim,wb1,GL)
            N1=self.mesh.N1;
            N2=self.mesh.N2;
            N=self.mesh.N;
            hbar = 1.05457e-34;
            q = 1.602e-19;
            me = 9.1e-31;
            TL = self.TL;
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
            phi=[wb1;phim;self.wbn];
            Vps=self.Semiconductor_Device.P_Vsaturate;
            Vns=self.Semiconductor_Device.N_Vsaturate;
            beta=self.Semiconductor_Device.Beta;
            Ep=self.Semiconductor_Device.Ep;
            gamma=self.Semiconductor_Device.gamma;            
            In=self.PD_Str.In;
            CnAu=self.Semiconductor_Device.CnAu;
            CpAu=self.Semiconductor_Device.CpAu;
            Br=self.Semiconductor_Device.Br;
            Eb_re=self.PD_Str.Eb_re;
            
            tau_n = self.Semiconductor_Device.Lifetime_n;
            tau_p = self.Semiconductor_Device.Lifetime_p;
            ni = self.Semiconductor_Device.Intrinsic_Density;
            Doping = self.Semiconductor_Device.Doping_Profile;
            n1=self.Semiconductor_Device.n1; % used in the incomplete ionization model
            p1=self.Semiconductor_Device.p1; % used in the incomplete ionization model

            % incomplete ionization model
            Doping_ion=zeros(N,1);
            NAp=zeros(N,1);
            NDp=zeros(N,1);
            Doping_ion(self.PD_Str.nt==2)=Doping(self.PD_Str.nt==2)./(1+4*p(self.PD_Str.nt==2)./p1(self.PD_Str.nt==2));
            Doping_ion(self.PD_Str.nt==1) = Doping(self.PD_Str.nt==1)./(1+2*n(self.PD_Str.nt==1)./n1(self.PD_Str.nt==1));
            NAp(self.PD_Str.nt==2)=abs(Doping_ion(self.PD_Str.nt==2));
            NDp(self.PD_Str.nt==1)=abs(Doping_ion(self.PD_Str.nt==1));
            
            un = self.Semiconductor_Device.Electron_Mobility;
            up = self.Semiconductor_Device.Hole_Mobility;
            %%%%%%%%%%%%
            hetjunctions = find(diff(self.PD_Str.what_sc)~=0);

            Ec = phi*VT+Xi;
            Ev = phi*VT+Xi+Eg;            
            if self.Semiconductor_Device.Thermionic_e ==1 && isempty(hetjunctions) == 0
                for i56 = 1:length(hetjunctions)
                    yy = hetjunctions(i56);
                    he=self.mesh.Node(yy);
                    Ec1_Ec2 = Ec(he)-Ec(he+1);
                    m_factor = exp(-abs(Ec1_Ec2)/8.617333262e-5/TL(he+1));                    
                    if abs(Ec1_Ec2) >0.01
                        if Ec1_Ec2 >0
                            un(he)=un(he)*m_factor;
                            beta(he) = 0;
                            Vns(he)= Vns(he)*m_factor;
                        end
                    end
                end
            end
            if self.Semiconductor_Device.Thermionic_h ==1 && isempty(hetjunctions) == 0
                for i56 = 1:length(hetjunctions)
                    yy = hetjunctions(i56);
                    he=self.mesh.Node(yy);
                    Ev1_Ev2 = Ev(he)-Ev(he+1);
                    m_factor = exp(-abs(Ev1_Ev2)/8.617333262e-5/TL(he));                    
                    if abs(Ev1_Ev2)>0.01
                        if Ev1_Ev2 >0
                            up(he+1)= up(he+1)*m_factor;
                            Vps(he+1)= Vps(he+1)*1e100;
                        end
                    end
                end
            end                
            
            E = -self.fd.diffmp(phi);% electrical field
            Electron_Velocity = E.*((un(2:end)+Vns(2:end).*beta(2:end).*abs(E))./(1+beta(2:end).*(E.^2)));
            Hole_Velocity = E.*(Vps(2:end).*up(2:end)./((Vps(2:end).^gamma+up(2:end).^gamma.*(abs(E).^gamma)).^(1/gamma))); % Hole Velocity
            Electron_Diffusion = TL(2:end)./300.*(un(2:end)./((1-2*((E)./Ep).^2+4/3*(abs(E)./Ep).^3).^(1/4)));
            Hole_Diffusion = TL(2:end)./300.*(up(2:end).*Vps(2:end)./((Vps(2:end).^gamma+up(2:end).^gamma.*(abs(E).^gamma)).^(1/gamma)));
            E_main = -self.fd.diff3p(phi); % the electric field at the main point
            
            if self.Semiconductor_Device.Impact_Ion
                Ae = self.Semiconductor_Device.Ae*NX; 
                Be = self.Semiconductor_Device.Be*NX/VT; 
                Ah = self.Semiconductor_Device.Ah*NX; 
                Bh = self.Semiconductor_Device.Bh*NX/VT; 
                %
                Vn_main = E_main.*((un(2:end-1)+Vns(2:end-1).*beta(2:end-1).*abs(E_main))./(1+beta(2:end-1).*(E_main.^2))); % Electron Velocity at the main point
                Vp_main = E_main.*(Vps(2:end-1).*up(2:end-1)./((Vps(2:end-1).^gamma+up(2:end-1).^gamma.*(abs(E_main).^gamma)).^(1/gamma))); % Hole Velocity at the main point
                Dnz_main = (un(2:end-1)./((1-2*(abs(E_main)./Ep(1:end-1)).^2+4/3*(abs(E_main)./Ep(1:end-1)).^3).^(1/4)));
                Dpz_main = Vps(2:end-1).*up(2:end-1)./((Vps(2:end-1).^gamma+up(2:end-1).^gamma.*(abs(E_main).^gamma)).^(1/gamma));
                alpha_An = Ae(2:end-1).*exp(-(Be(2:end-1)./abs(E_main)).^In);
                alpha_Ap = Ah(2:end-1).*exp(-(Bh(2:end-1)./abs(E_main)).^In);
                G_imp = alpha_An.*abs(n(2:end-1).*Vn_main+Dnz_main.*self.fd.diff3p(n(1:end)));
                G_imp = alpha_Ap.*abs(p(2:end-1).*Vp_main-Dpz_main.*self.fd.diff3p(p(1:end)))+G_imp;
            else
                G_imp=0;
            end
            
            alpha0=Absorp_FK(abs([E(1);E])*VT/NX*100,self.wavelength,self.PD_Str,self.Semiconductor_Device.FranzKeldysh);
            ave_alpha = mean(alpha0(alpha0>0));
            absorbers = ones(length(alpha0),1);
            absorbers(alpha0<10) = 0;
            Generation = GL.*absorbers;

            Recombination=(n.*p-ni.^2)./(tau_p.*(n+ni)+tau_n.*(p+ni));
            Re_Aug=(CnAu.*n+CpAu.*p).*(n.*p-ni.^2);     % Auger recombination
            Re_Rad=Br.*(n.*p-ni.^2);    % Radiative recombination
            %
         
            J_Hole_drift =     Hole_Velocity.*self.fd.averagemp(p);% drift current
            J_Hole_diff =    -Hole_Diffusion.*self.fd.diffmp(p);% diffusion current
            J_Elec_drift=      Electron_Velocity.*self.fd.averagemp(n);
            J_Elec_diff=    Electron_Diffusion.*self.fd.diffmp(n);
            
            J_Hole = J_Hole_drift+J_Hole_diff;% Hole current,there should be a q times, but it cancelled in the Fp.
            J_Elec = J_Elec_drift+J_Elec_diff;% electron current
            J_T=J_Hole+J_Elec;
            Fp = -self.fd.diffmidp(J_Hole)-Recombination(2:end-1)-Re_Aug(2:end-1)-Re_Rad(2:end-1)+Generation(2:end-1)+G_imp; % drift % diffusion current equation for hole
            Fn = self.fd.diffmidp(J_Elec)-Recombination(2:end-1)-Re_Aug(2:end-1)-Re_Rad(2:end-1)+Generation(2:end-1)+G_imp;% drift % diffusion current equation for electron
            dx=self.fd.mesh.dx;
            dx1=self.fd.mesh.dx1;
            %            
            Fw = Eb_re(1:end-1).*(Doping_ion(2:end-1)+p(2:end-1)-n(2:end-1))... % ni(2:end-1)is norm
                +(1./((dx(1:end-1)).*(dx1(1:end)))).*phi(1:end-2)...
                -((1./dx(1:end-1)+1./dx(2:end))./dx1(1:end)).*phi(2:end-1)...
                +(1./((dx(2:end)).*(dx1(1:end)))).*phi(3:end);% possion equation
            %
            Jave=mean(J_T(1))*self.NP.J0*pi*(self.Diameter(1))^2*0.25;
            Fb = self.wbf-wb1-(self.Bias+(Jave)*self.R_Load)/self.NP.VT;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [y,M]=Newton_FD(self,x)
            N=length(self.mesh.dx)+1;
            pm = x(1:N-2);
            nm = x(N-1:2*N-4);
            wm = x(2*N-3:3*N-6);
            wb1 = x(end);
            [Fp,Fn,Fw,Fb]=self.Cal_Currentv2(pm,nm,wm,wb1,self.GL);
            y = [Fp;Fn;Fw;Fb];
            M = self.ConsJacobian(pm,nm,wm,wb1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [JacobM,MNADt]=ConsJacobian(self,pm,nm,phim,wb1)
            % Nmu=self.Semiconductor_Device.NP.Nmu; % normlized mobolity;
            q = 1.602e-19;
            hbar = 1.05457e-34;
            me = 9.1e-31;
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
            D0=self.Semiconductor_Device.NP.D0; % normlized diffusion
            p=[self.pb1;pm;self.pbn];
            n=[self.nb1;nm;self.nbn];
            phi=[wb1;phim;self.wbn];
            Vps=self.Semiconductor_Device.P_Vsaturate;
            Vns=self.Semiconductor_Device.N_Vsaturate;
            beta=self.Semiconductor_Device.Beta;
            Ep=self.Semiconductor_Device.Ep;
            gamma=self.Semiconductor_Device.gamma;
            tau_n = self.Semiconductor_Device.Lifetime_n;
            tau_p = self.Semiconductor_Device.Lifetime_p;
            ni = self.Semiconductor_Device.Intrinsic_Density;
            Doping = self.Semiconductor_Device.Doping_Profile;
            CnAu=self.Semiconductor_Device.CnAu;
            CpAu=self.Semiconductor_Device.CpAu;
            Br=self.Semiconductor_Device.Br;
            up = self.Semiconductor_Device.Hole_Mobility;
            un = self.Semiconductor_Device.Electron_Mobility;            
            N=self.mesh.N;
            N1=self.mesh.N1;
            % incomplete ionization model
            n1=self.Semiconductor_Device.n1; % used in the incomplete ionization model
            p1=self.Semiconductor_Device.p1; % used in the incomplete ionization model
            Doping_ion=zeros(N,1);
            % Doing_ion  derive to the p and n
            Doping_ion_p=zeros(N,1);
            Doping_ion_n=zeros(N,1);
            pdoped = find(sign(Doping)==-1);
            ndoped = find(sign(Doping)==1);
            TL = self.PD_Str.TL;
            hetjunctions = find(diff(self.PD_Str.what_sc)~=0);
            
            Ec = phi*VT+Xi;
            Ev = phi*VT+Xi+Eg;

            
            if self.Semiconductor_Device.Thermionic_e ==1 && isempty(hetjunctions) == 0
                for i56 = 1:length(hetjunctions)
                    yy = hetjunctions(i56);
                    he=self.mesh.Node(yy);
                    Ec1_Ec2 = Ec(he)-Ec(he+1);
                    m_factor = exp(-abs(Ec1_Ec2)/8.617333262e-5/TL(he+1));                    
                    if abs(Ec1_Ec2) >0.01
                        if Ec1_Ec2 >0
                            un(he)=un(he)*m_factor;
                            beta(he) = 0;
                            Vns(he)=Vns(he)*m_factor;
                        end
                    end
                end
            end
            if self.Semiconductor_Device.Thermionic_h ==1 && isempty(hetjunctions) == 0
                for i56 = 1:length(hetjunctions)
                    yy = hetjunctions(i56);
                    he=self.mesh.Node(yy);
                    Ev1_Ev2 = Ev(he)-Ev(he+1);
                    m_factor = exp(-abs(Ev1_Ev2)/8.617333262e-5/TL(he));
                    
                    if abs(Ev1_Ev2)>0.01
                        if Ev1_Ev2 >0
                            up(he+1)= up(he+1)*m_factor;
                            Vps(he+1)= Vps(he+1)*1e100;
                        end
                    end
                end
            end            

            Doping_ion(self.PD_Str.nt==2)=Doping(self.PD_Str.nt==2)./(1+4*p(self.PD_Str.nt==2)./p1(self.PD_Str.nt==2));
            Doping_ion(self.PD_Str.nt==1) = Doping(self.PD_Str.nt==1)./(1+2*n(self.PD_Str.nt==1)./n1(self.PD_Str.nt==1));
            Doping_ion_p(pdoped)=-4*Doping_ion(pdoped)./p1(pdoped)./(p1(pdoped)+4*p(pdoped)).^2;
            Doping_ion_n(ndoped)=-2*Doping_ion(ndoped)./n1(ndoped)./(n1(ndoped)+2*n(ndoped)).^2;

            tempa=[-Doping_ion_p(2:N-1);Doping_ion_n(2:N-1);zeros(N-2+1,1)];
            MNADt=spdiags(tempa,0,3*(N-2)+1,3*(N-2)+1);
            E = -self.fd.diffmp(phi);
            In = self.PD_Str.In;
            
            DUDn=(p+ni).*(p.*tau_n+tau_p.*ni)./(tau_p.*(n+ni)+tau_n.*(p+ni)).^2;
            DUDp=(n+ni).*(n.*tau_p+tau_n.*ni)./(tau_p.*(n+ni)+tau_n.*(p+ni)).^2;

            % Auger recombination
            Dre_auDp=(CpAu.*(n.*p-ni.^2)+(CnAu.*n+CnAu.*p).*n);
            Dre_auDn=(CnAu.*(n.*p-ni.^2)+(CnAu.*n+CnAu.*p).*p);
            % Radiative recombination
            Dre_radDp=Br.*n;
            Dre_radDn=Br.*p;
            Dre_auDp=Dre_auDp+Dre_radDp;
            Dre_auDn=Dre_auDn+Dre_radDn;
            %
            %CreateJacobM()
            % the impact ionization generation
            % G_imp=alph_An*n*abs(vn)+alpha_Ap*p*abs(vp)
            % this vn and vp should be on the main point;
            if self.Semiconductor_Device.Impact_Ion
                Ae = self.Semiconductor_Device.Ae*NX;
                Be = self.Semiconductor_Device.Be*NX/VT;
                Ah = self.Semiconductor_Device.Ah*NX;
                Bh = self.Semiconductor_Device.Bh*NX/VT;
                
                E_main = -self.fd.diff3p(phi); % the electric field at the main point
                Vn_main = E_main.*((un(2:end-1)+Vns(2:end-1).*beta(2:end-1).*abs(E_main))./(1+beta(2:end-1).*(E_main.^2))); % Electron Velocity at the main point
                Vp_main = E_main.*(Vps(2:end-1).*up(2:end-1)./((Vps(2:end-1).^gamma+up(2:end-1).^gamma.*(abs(E_main).^gamma)).^(1/gamma))); % Hole Velocity at the main point
                Dnz_main = (un(2:end-1)./((1-2*(abs(E_main)./Ep(1:end-1)).^2+4/3*(abs(E_main)./Ep(1:end-1)).^3).^(1/4)));
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
            else
                DG_impDn=0;
                DG_impDp=0;
                DG_impDwm=zeros(N-2,1);
                DG_impDwp=zeros(N-2,1);                
            end

            Electron_Velocity = E.*((un(2:end)+Vns(2:end).*beta(2:end).*abs(E))./(1+beta(2:end).*(E.^2)));
            Hole_Velocity = E.*(Vps(2:end).*up(2:end)./((Vps(2:end).^gamma+up(2:end).^gamma.*(abs(E).^gamma)).^(1/gamma))); % Hole Velocity
            Electron_Diffusion = self.TL(2:end)./300.*(un(2:end)./((1-2*((E)./Ep).^2+4/3*(abs(E)./Ep).^3).^(1/4)));
            Hole_Diffusion = self.TL(2:end)./300.*(up(2:end).*Vps(2:end)./((Vps(2:end).^gamma+up(2:end).^gamma.*(abs(E).^gamma)).^(1/gamma)));
            DVnDE = (1-beta(2:end).*(E.^2)).*(un(2:end)+Vns(2:end).*beta(2:end).*abs(E))./((1+beta(2:end).*(E.^2)).^2)+...
                E.*(Vns(2:end).*beta(2:end).*sign(E))./(1+beta(2:end).*(E.^2));
            DVpDE = up(2:end).*Vps(2:end)./((Vps(2:end).^gamma+up(2:end).^gamma.*(abs(E).^gamma)).^(1/gamma))...
                -up(2:end).^(gamma+1).*Vps(2:end).*(abs(E).^gamma)./((Vps(2:end).^gamma+up(2:end).^gamma.*(abs(E).^gamma)).^((gamma+1)/gamma));
            DDnDE = self.TL(2:end)./300.*un(2:end).*((E)./Ep.^2-sign(E).*(E.^2)./Ep.^3)./((1-2*(E./Ep).^2+4/3*(abs(E)./Ep).^3).^(5/4));
            DDpDE = -self.TL(2:end)./300.*up(2:end).^(gamma+1).*Vps(2:end).*(abs(E).^(gamma-1)).*sign(E)./((Vps(2:end).^(gamma)+up(2:end).^(gamma).*(abs(E).^(gamma))).^((gamma+1)/gamma));
            %%
            % J_Hole and J_ Electron first derivative to p ,n and phi            
            DJpDp = Hole_Velocity/2+Hole_Diffusion.*self.fd.differ;% first derivative hole current J_Hole(i) to Hole_Density(i)
            DJpDp2 = Hole_Velocity/2-Hole_Diffusion.*self.fd.differ;% first derivative J_Hole(i) to Hole_Density(i+1)
            
            DJpDw = DVpDE.*self.fd.averagemp(p).*self.fd.differ-DDpDE.*self.fd.diffmp(p).*self.fd.differ; % first derivative J_Hole(i) to Potential(i)
            DJpDw2 = -DJpDw; %  % first derivative J_Hole(i) to Potential(i+1)
            % boundary condition
            
            DJnDn = Electron_Velocity/2-Electron_Diffusion.*self.fd.differ; % % first derivative electron current J_Elec(i) to Elec_Density(i)
            DJnDn2 = Electron_Velocity/2+Electron_Diffusion.*self.fd.differ; % first derivative current J_Elec(i) to Elec_Density(i+1)
            DJnDw = DVnDE.*(self.fd.averagemp(n)).*self.fd.differ+DDnDE.*self.fd.diffmp(n).*self.fd.differ;  % first derivative  J_Elec(i) to potential(i)
            DJnDw2 = -DJnDw; % first derivative  J_Elec(i) to potential(i+1)
            
            %%%

            %**********************DFw/Dp,DFw/Dn,DFw/Dw************            
            dx1=self.mesh.dx1;
            Num=length(dx1)+2;
            % Fp derivative to p , n and phi
            a=[DJpDp(2:end-1)./(dx1(2:end));0];
            e=(-DJpDp(2:end)+DJpDp2(1:end-1))./(dx1(1:end))-DUDp(2:end-1)-Dre_auDp(2:end-1)+DG_impDp;%(2:N-1)
            b=[0;-DJpDp2(2:end-1)./(dx1(1:end-1))];
            M1=spdiags([a e b],-1:1,Num-2,Num-2); % Fp first derivative to p(hole denstiy)
            M2=spdiags(-DUDn(2:end-1)-Dre_auDn(2:end-1)+DG_impDn,0,Num-2,Num-2); % Fp first derivative to n (electron denstiy)
            a=[DJpDw(2:end-1)./(dx1(2:end))+DG_impDwm(2:end);0];
            e=(-DJpDw(2:end)+DJpDw2(1:end-1))./(dx1(1:end));
            b=[0;-DJpDw2(2:end-1)./(dx1(1:end-1))+DG_impDwp(1:end-1)];
            M3=spdiags([a e b],-1:1,Num-2,Num-2);% Fp first derivative to potential
            M4=zeros(Num-2,1);
            M4(1)=DJpDw(1)./dx1(1);
            MFp=[M1,M2,M3,M4];% Fp first derivative to p ,n and potential
            
            % Fn derive to p n w
            M1=spdiags(-DUDp(2:end-1)-Dre_auDp(2:end-1)+DG_impDp,0,Num-2,Num-2); % Fn first derivative to p(hole denstiy)
            a=[-DJnDn(2:end-1)./(dx1(2:end));0];
            e=(DJnDn(2:end)-DJnDn2(1:end-1))./(dx1(1:end))-DUDn(2:end-1)-Dre_auDn(2:end-1)+DG_impDn;
            b=[0;DJnDn2(2:end-1)./(dx1(1:end-1))];
            M2=spdiags([a e b],-1:1,Num-2,Num-2); % Fn first derivative to n(electron denstiy)
            a=[-DJnDw(2:end-1)./(dx1(2:end))+DG_impDwm(2:end);0]; % diag-1,row :i, col; i-1
            e=(DJnDw(2:end)-DJnDw2(1:end-1))./(dx1(1:end)); % diag
            b=[0;DJnDw2(2:end-1)./(dx1(1:end-1))+DG_impDwp(1:end-1)];
            M3=spdiags([a e b],-1:1,Num-2,Num-2); % Fn first derivative to potential
            M4=zeros(Num-2,1);
            M4(1)=-DJnDw(1)./dx1(1);
            MFn=[M1,M2,M3,M4];
            % Phi first derivative to p,n and phi
            DFphiDp = self.fd.CreateZOM+spdiags(Doping_ion_p(2:end-1),0,Num-2,Num-2);% Fphi first derivative to p
            permit_matrix = spdiags(self.PD_Str.Eb_re,0,Num-2,Num-2);
            DFphiDn = -self.fd.CreateZOM-spdiags(Doping_ion_n(2:end-1),0,Num-2,Num-2); % Fphi first derivative to n            
            DFphiDn = permit_matrix*DFphiDn;
            DFphiDp = permit_matrix*DFphiDp;
            DFphiDphi= self.fd.CreateSOM; % Fp first derivative to phi
            dx=self.fd.mesh.dx;
            dx1=self.fd.mesh.dx1;
            
            M4=zeros(Num-2,1);
            M4(1)=1./(dx1(1)*dx(1));
            MFb=zeros(1,3*(Num-2)+1);
            divterm = pi*self.Diameter(1)^2*0.25*self.R_Load/self.NP.VT;
            MFb(end)=-1-(DJpDw(1)+DJnDw(1))*self.NP.J0*divterm;
            MFb(1)=-DJpDp2(1)*self.NP.J0*divterm;
            MFb(N-1)=-DJnDn2(1)*self.NP.J0*divterm;
            MFb(2*N-3)=-(DJpDw2(1)+DJnDw2(1))*self.NP.J0*divterm;
            JacobM = [MFp;MFn;DFphiDp,DFphiDn,DFphiDphi,M4;MFb];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [y]=Newton_FD_y(self,t,x)
            N=length(self.mesh.dx)+1;
            pm = x(1:N-2);
            nm = x(N-1:2*N-4);
            wm = x(2*N-3:3*N-6);
            wb1=x(end);
            G0=self.GL;
            [Fp,Fn,Fw,Fb]=self.Cal_Currentv2(pm,nm,wm,wb1,G0);
            y = [Fp;Fn;Fw;Fb];
            %M = self.ConsJacobian(pm,nm,wm);
        end
        function [M]=Newton_FD_J(self,t,x)
            N=length(self.mesh.dx)+1;
            pm = x(1:N-2);
            nm = x(N-1:2*N-4);
            wm = x(2*N-3:3*N-6);
            wb1 = x(end);
            M = self.ConsJacobian(pm,nm,wm,wb1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [y,M]=Newton_FD_Time(self,x)
            Deltat=self.deltat;
            N=length(self.mesh.dx)+1;
            pm = x(1:N-2);
            nm = x(N-1:2*N-4);
            wm = x(2*N-3:3*N-6);
            x0=x(3*N-5:end);
            pm0 = x0(1:N-2);
            nm0 = x0(N-1:2*N-4);
            Ft=([pm-pm0;nm-nm0;wm-wm])/Deltat;
            
            [Fp,Fn,Fw]=self.Cal_Currentv2(pm,nm,wm);
            y = [Fp;Fn;Fw]-Ft;
            
            M = self.ConsJacobian(pm,nm,wm)+self.Mt;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
