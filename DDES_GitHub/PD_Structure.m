classdef PD_Structure < handle
    properties
        Num_Layer % number of layers
        what_sc % the type of the semiconductor for layer (n, p, Q doesn't matter)
        Type_sc %  the type of the semiconductor for each node
        Layer_Lx % the position of each layer
        Type_Layer % the type of each layer 1for p,-1for n, 0 for i;
        Layer_Doping % Doping information
        % Doping_type % doping type constant or graded
        Node % mesh including node and distance
        Lx
        NC % conduction band density  in numerical
        NV % valance band density  in numerical
        ni % numerical intrisic density
        Ni % for normalization
        tn % electron life time numerical
        tp % hole life time numerical
        Tol % doping density of layer in numerical
        mstar_n % effective mass of electron for each layer
        mstar_p % effective mass of hole for each layer
        mes % effective mass of electron for each node
        mhs % effective mass of hole for each node
        Ne_ref   % reference electron concentration
        Np_ref   % reference hole concentration
        TL % temperature
        Tr % temperature for the alpha absorption
        Vns % electron saturated velocity
        Vps % hole saturated velocity
        un % electron mobility as Temperature dependent at low field
        up % hole mobility as Temperature dependent at low field
        %Doping % numerical doping information
        Beta
        CnAu % AUger recombination
        CpAu % AUger recombination
        Br % radiative recombiantion
        Eg % energy gap
        Eb_re % relative permitivity to InGaAs
        Epsilons % epsilon_relative*8.854 for each layer
        nt  % 2 for p doped, 1 for negative doped
        Ep
        DeltaEA
        DeltaED
        Ae
        Be
        Ah
        Bh
        k_L
        Xi
        eta_n
        eta_p
        In % the exponent power, m in Eq. 23 Yue Hu 10.1109/JLT.2014.2315740
        alpha
        alphas_abs
        delta_Fermis
    end
    %
    methods
        function self=PD_Structure(Ni,Num_Layer, Layer_Lx, ...
                Type_Layer,Layer_Doping,Node,Lx,what_sc,In,...
                mes, mhs, Egs, CnAus, CpAus, Brs, NCs, NVs, Aes, Bes, ...
                Ahs, Bhs, betas, Eps, DeltaEAs, DeltaEDs, Ne_refs,  ...
                Np_refs, Vnsats, Vpsats, Ebslons, uns, ups, tor1s, tor2s, ...
                n_is, Xis, alphas, alphas_abs,eta_ns, eta_ps, TL,Tr,...
                epsilon_reference,Qvalues,delta_Fermi,k_L)
            
            self.Ni = Ni;
            self.Num_Layer = Num_Layer;
            self.Layer_Lx = Layer_Lx;
            self.Type_Layer = Type_Layer;
            self.Layer_Doping = Layer_Doping;
            self.Node = Node;
            self.Lx = Lx;            
            self.what_sc = what_sc;
            self.In = In;            
            self.Type_sc = ones(Num_Layer,1);

            self.alpha = ones(Num_Layer,1);
            self.alphas_abs= ones(Num_Layer,1);
            self.eta_n = zeros(Num_Layer,1);
            self.eta_p = zeros(Num_Layer,1);
            self.mstar_n = zeros(Num_Layer,1);
            self.mstar_p = zeros(Num_Layer,1);
            self.Eg = zeros(Node(end),1);
            self.delta_Fermis = zeros(Node(end),1);
            self.mhs = zeros(Node(end),1);
            self.mes = zeros(Node(end),1);
            self.CnAu = zeros(Node(end),1);
            self.CpAu = zeros(Node(end),1);
            self.Br = zeros(Node(end),1);
            self.NC = zeros(Node(end),1);
            self.NV = zeros(Node(end),1);
            self.Ae = zeros(Node(end),1);
            self.Be = zeros(Node(end),1);
            self.Ah = zeros(Node(end),1);
            self.Bh = zeros(Node(end),1);
            self.k_L = zeros(Node(end),1);
            self.Beta = zeros(Node(end),1);
            self.Ep = ones(Node(end)-1,1);
            self.DeltaEA = zeros(Num_Layer,1);
            self.DeltaED = zeros(Num_Layer,1);
            self.Ne_ref= zeros(Node(end),1);
            self.Np_ref = zeros(Node(end),1);
            self.Vns = zeros(Node(end),1);
            self.Vps = zeros(Node(end),1);
            self.Epsilons = ones(Node(end)-1,1);
            self.un = zeros(Node(end),1);
            self.up = zeros(Node(end),1);
            self.tn = zeros(Node(end),1);
            self.tp = zeros(Node(end),1);
            self.ni = zeros(Node(end),1);
            self.Xi = zeros(Node(end),1);            
            self.TL = TL;
            self.Tr = Tr;
            self.Eb_re= ones(Node(end)-1,1);    % relative permittivity normalized by eps_InGaAs
            self.nt = zeros(Node(end),1);            
            self.Tol = zeros(Node(end),1);
            for ii=1:Num_Layer
                first_node=1;
                first_node_shifted=1;
                last_node_shifted = Node(ii)-1;
                
                if ii>1; first_node = Node(ii-1)+1; first_node_shifted = Node(ii-1); end
                self.DeltaEA(ii)=DeltaEAs(ii);
                self.DeltaED(ii)=DeltaEDs(ii);               
                self.eta_n(ii) = eta_ns(ii); 
                self.eta_p(ii) = eta_ps(ii); 
                self.mstar_n(ii) = mes(ii);
                self.mstar_p(ii) = mhs(ii);
                self.mes(first_node:Node(ii)) = mes(ii);
                self.mhs(first_node:Node(ii)) = mhs(ii);
                self.Type_sc(first_node:Node(ii)) = what_sc(ii);
                self.alpha(first_node:Node(ii)) = alphas(ii);
                self.alphas_abs(first_node:Node(ii)) =alphas_abs(ii);
                self.Eg(first_node:Node(ii))=Egs(ii);
                self.delta_Fermis(first_node:Node(ii))=delta_Fermi(ii);
                self.CnAu(first_node:Node(ii))=CnAus(ii);
                self.CpAu(first_node:Node(ii))=CpAus(ii);
                self.Br(first_node:Node(ii))=Brs(ii);
                self.NC(first_node:Node(ii))=NCs(ii);
                self.NV(first_node:Node(ii))=NVs(ii);
                self.Ae(first_node:Node(ii))=Aes(ii);
                self.Be(first_node:Node(ii))=Bes(ii);
                self.Ah(first_node:Node(ii))=Ahs(ii);
                self.Bh(first_node:Node(ii))=Bhs(ii);
                self.k_L(first_node:Node(ii))=k_L(ii);
                self.Beta(first_node:Node(ii))=betas(ii);
                self.Ep(first_node_shifted:Node(ii)-1)=Eps(ii);
                self.Ne_ref(first_node:Node(ii))=Ne_refs(ii);
                self.Np_ref(first_node:Node(ii))=Np_refs(ii);
                self.Vns(first_node:Node(ii))=Vnsats(ii);
                self.Vps(first_node:Node(ii))=Vpsats(ii);
                self.Epsilons(first_node_shifted:Node(ii)-1)=Ebslons(ii);
                if Layer_Doping(ii)<-1  % for undoped region, use 1 as doping! 
                    self.nt(first_node:Node(ii))=2;
                elseif Layer_Doping(ii)>1  % for undoped region, use 1 as doping! 
                    self.nt(first_node:Node(ii))=1;
                end
                self.un(first_node:Node(ii))=uns(ii)/(1+(abs(Layer_Doping(ii))/Ne_refs(ii))^eta_ns(ii));
                self.up(first_node:Node(ii))=ups(ii)/(1+(abs(Layer_Doping(ii))/Np_refs(ii))^eta_ps(ii));
                self.ni(first_node:Node(ii))=n_is(ii);
                % self.Eb_re(first_node_shifted:last_node_shifted)=Ebslons(ii)/epsilon_reference;
                % 2023/4/27, we are dividing with permittivity, not
                % multipliying! since in other parts we do Eb_re*..., here
                % we need to do eb_reference/eb_real
                self.Eb_re(first_node_shifted:last_node_shifted)=epsilon_reference/Ebslons(ii);
                self.Tol(first_node:Node(ii))=Layer_Doping(ii);
                if Qvalues(ii)> 10 || Qvalues(ii) <-10
                    self.Tol(first_node:Node(ii))=linspace(Layer_Doping(ii),Qvalues(ii),length(first_node:Node(ii)));
                end
%                 if Qvalues(ii)> 10 || Qvalues(ii) <-10                    
%                     self.Tol(first_node:Node(ii))=logspace(log10(abs(Layer_Doping(ii))),log10(abs(Qvalues(ii))),length(first_node:Node(ii)));
%                     self.Tol(first_node:Node(ii)) = sign(Layer_Doping(ii))*self.Tol(first_node:Node(ii));
%                 end
                self.Xi(first_node:Node(ii)) = Xis(ii);    
% tor2 = 1e-10;    % in region-i
% tor1 = 1e-8;   % o/wise                
                if Type_Layer(ii)==0
                    self.Beta(first_node:Node(ii))=0.75*betas(ii);
                    self.tn(first_node:Node(ii)) = tor1s(ii);
                    self.tp(first_node:Node(ii)) = tor1s(ii)*10;
                else
                    self.tn(first_node:Node(ii)) = tor2s(ii);
                    self.tp(first_node:Node(ii)) = tor2s(ii)*10;
                end
            end
            
        end
    end
end
