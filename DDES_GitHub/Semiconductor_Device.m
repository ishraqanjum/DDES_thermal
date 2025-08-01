classdef Semiconductor_Device < handle
    %SEMICONDUCTOR_DEVICE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
         % total nodes
        Length% the device length unit cm
        mesh % the class mesh
        %dx % grid profile array(N-1), the distance between Nodes, unit cm
        fd % Finite difference object        
        Doping_Profile% doping profile of the device ,
        %it should be a array. unit cm-3, acceptor is minus, donor is positive        
        Electron_Mobility % the electron mobility unit,cm^2/(V*s)
        Hole_Mobility % the hole mobility unit,cm^2/(V*s)
        Beta % fitting parameter for the electron velocity        
        Lifetime_n % recombination lifetime for electron        
        Lifetime_p % recombination lifetime for hole     
        Temperature % Temperature unit K
        T_Vol % electron voltage. unit V, T_Vol=KB*T/q; q is the charge        
        Hole_Density% hole density profile array of N, unit cm-3
        Elec_Density% electron density profile array of N unit cm-3
        Potential% Potential Voltage unitV
        Elec_Field % Electrical Field unit V/cm
        Ep  % calculated saturted electrical field parameter
        Intrinsic_Density % intrinsic density
        P_Vsaturate
        N_Vsaturate
        Diameter  % Diameter of the Beam
        DeviceDiameter % diameter to calculate the current
        NP % normlized parameter class
        gamma % the parameter in the hole velocity, to control the increase of the velocity as a function of E
        n1 % the ionization density
        p1 % the ionization denstiy
	    Impact_Ion % 1 for the impact ionization model
        Thermionic_e % 1 for thermionic field emission for electrons
        Thermionic_h % 1 for thermionic field emission for holes
        FranzKeldysh % 1 for FK effect
        Ae % impact ionization parameter unit cm-1
        Be % impact ionization parameter unit V/cm 
        Ah % impact ionization parameter unit cm-1
        Bh % impact ionization parameter unit V/cm
        k_L
        Eg1 % EG used in the Franz-keldysh
        CnAu % Auger Recombination para
        CpAu % Auger Recombination para
        Br % radiative recombination para
        DeltaEA 
        DeltaED
        Ne_ref 
        Np_ref
        Xi
        epsilon_reference
    end
    properties (Constant)
        KB=1.38e-23;% Boltzmann constant,unit,J/K
        q=1.6e-19;
    end
    methods
        function obj=Semiconductor_Device(NP,Temperature,Length,mesh,fd,Doping_Profile,Beta,Electron_Mobility,Hole_Mobility,...
                Lifetime_n,Lifetime_p,Intrinsic_Density,Potential,P_Vsaturate,N_Vsaturate,Ep,Beam_Dia,Device_Dia,gamma,n1,p1,...
                Impact_Ion,Thermionic_e,Thermionic_h,FranzKeldysh,Ae,Be,Ah,Bh,Eg1,CnAu,CpAu,Br,DeltaEA, DeltaED, Ne_ref, Np_ref,epsilon_reference,Xi,k_L)
            % the Doping_Profile should be the array. And the length of
            % Doping Profile should be equal to the "Part"
            obj.epsilon_reference = epsilon_reference;
            obj.Xi= Xi;
            obj.DeltaEA = DeltaEA;
            obj.DeltaED = DeltaED;            
            obj.Ne_ref = Ne_ref;
            obj.Np_ref = Np_ref;
            obj.k_L=k_L;
            obj.NP=NP;
            N=length(mesh.dx)+1;
            obj.Temperature=Temperature;
            obj.Length=Length;
            obj.mesh=mesh;
            obj.fd=fd;
            obj.Doping_Profile=Doping_Profile;
            obj.T_Vol=obj.KB*Temperature/obj.q;
            obj.Electron_Mobility=Electron_Mobility;
            obj.Hole_Mobility=Hole_Mobility; 
            obj.Beta=Beta; % the fitting parameter
            obj.Lifetime_n=Lifetime_n;
            obj.Lifetime_p=Lifetime_p;
            obj.Intrinsic_Density=Intrinsic_Density;
            obj.Hole_Density=abs(min(Doping_Profile,zeros(N,1)));% initial the Hole density to the acceptor doping
            obj.Elec_Density=abs(max(Doping_Profile,zeros(N,1)));% initial the Hole density to the donor doping
            obj.Potential=Potential;
            obj.Elec_Field=-fd.diffmp(obj.Potential);
            obj.P_Vsaturate=P_Vsaturate;
            obj.N_Vsaturate=N_Vsaturate;
            obj.Ep=Ep;
            obj.Diameter=Beam_Dia;
            obj.DeviceDiameter=Device_Dia;
            obj.gamma=gamma;
            obj.n1=n1;
            obj.p1=p1;
            obj.Impact_Ion=Impact_Ion;
            obj.FranzKeldysh=FranzKeldysh;
            obj.Thermionic_e = Thermionic_e;
            obj.Thermionic_h = Thermionic_h;
            obj.Ae = Ae;
            obj.Be = Be;
            obj.Ah= Ah;
            obj.Bh= Bh;
            obj.Eg1= Eg1;
            obj.CnAu=CnAu;
            obj.CpAu=CpAu;
            obj.Br=Br;
        end
        

        function Update(self,p,n,phi)
            % update the carrier density and potential
            self.Hole_Density=p;
            self.Elec_Density=n;
            self.Potential=phi;
        end
        function Plotn(self,x)
            % x is the variable you want to plot
            % the normal coordinalator
            plot(self.mesh.Lx,x);
        end
        function Plotl(self,x)
            % x is the variable you want to plot
            % plot log function
            plot(self.mesh.Lx,log10(x));
        end
        function [p,n,phi]=Get(self)
            % to get the carrier and potential function 
            p=self.Hole_Density;
            n=self.Elec_Density;
            phi=self.Potential;
        end
    end    
end

