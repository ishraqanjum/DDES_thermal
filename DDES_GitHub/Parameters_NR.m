%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KB=1.38e-23;% Boltzmann constant,unit,J/K
q=1.6e-19;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tg=Temperature+0*abs(Bias)/8*50; % this is the temperature in the intrinsic region 300-460 15ma, 300-480, 20ma, 300-390, 10ma
Tr=Temperature+0*abs(Bias)/8*30; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Type_Layer = [Design{:,1}];
thicknesses = [Design{:,3}];
dopings = [Design{:,4}];
Qvalues = [Design{:,5}];
Layer_Lx = [0 cumsum(thicknesses)];
Length = Layer_Lx(end);
Layer_Doping = dopings;
neg_ns = find(Type_Layer==1);
asd = find(Type_Layer==1 & Qvalues>10); % find linearly graded p-type layers
Qvalues(asd) = -Qvalues(asd);
Layer_Doping(neg_ns) = -abs(dopings(neg_ns));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NP=Norm_parameter(KB,q,const.eps_silicon,Ni,Temperature,Diff_Coef);
%%
Num_Layer = length(Type_Layer);
what_sc = ones(1,Num_Layer);
% find the unique Q values and assign
% them to different numbers
Qlevels = unique([Design{:,5}]);
Qlevels(Qlevels<0.5) = [];
Qlevels(Qlevels>2) = [];
for iiii = 1:Num_Layer
    material = Design{iiii,2};
    switch material
        case 'InP'
            what_sc(iiii) = 2;
        case 'InGaAsP'
           if Design{iiii,5} == Qlevels(1)
                what_sc(iiii) = 3;
           else
               asd = find(Qlevels == Design{iiii,5});
                what_sc(iiii) = 4+asd;
           end
        case 'Si'
       what_sc(iiii) = 4;
        case 'Ge'
       what_sc(iiii) = 5;
   end
end

mes = zeros(1,Num_Layer);
eta_ns = zeros(1,Num_Layer);
eta_ps = zeros(1,Num_Layer);
mhs = mes;
Egs = mes;
CnAus = mes;
CpAus = mes; 
Brs = mes;
NCs = mes;
NVs = mes;
Aes = mes;
Bes = mes;
Ahs = mes;
Bhs = mes;
k_L = mes;
betas = mes;
Eps = mes;
DeltaEAs = mes;
DeltaEDs = mes;
Ne_refs = mes;
Np_refs = mes;
Vnsats = mes;
Vpsats = mes;
Ebslons = mes;
uns = mes;
ups = mes;
tor1s = mes;
tor2s = mes;
n_is = mes;
Xis = mes;
alphas = mes;
alphas_abs = mes;
delta_Fermi = mes;
% keyboard
for io = 1:Num_Layer
    [mes(io), mhs(io), Egs(io), CnAus(io), CpAus(io), Brs(io), NCs(io),...
        NVs(io), Aes(io), Bes(io),Ahs(io), Bhs(io), betas(io), ...
     Eps(io), DeltaEAs(io), DeltaEDs(io), Ne_refs(io), Np_refs(io), ...
    Vnsats(io), Vpsats(io), Ebslons(io), uns(io), ups(io), tor1s(io),...
    tor2s(io), n_is(io), Xis(io), alphas(io), eta_ns(io), eta_ps(io), ...
    alphas_abs(io), resistivity(io), delta_Fermi(io), k_L(io)] = get_SC_parameters(Design{io,2},Type_Layer(io), dopings(io), Qvalues(io), Ni, wavelength,NP,Design{io,6});
    % disp([io, alphas(io), alphas_abs(io)])
end

% Jin_Li design has a p-doped thin layer, so I had to change
% using "interfaces = find(diff(Type_Layer)==-1);"
interface1 = find(diff(cumprod(Type_Layer))==-1);
interface2 = Num_Layer-find(diff(cumprod(fliplr(-Type_Layer)))==-1);


interfaces = [interface1 interface2];
%%%mesh
mesh1=Variablemesh(Dmax,Length,NP.NX,thicknesses,interfaces);
N=mesh1.N;% The total nodes 
TL=Tg*ones(N,1); % temperature distribution 
TLv=Tr*ones(N,1);
%%%%%%%%%
% keyboard
PD_Str1 = PD_Structure(Ni,Num_Layer, Layer_Lx, ...
    Type_Layer,Layer_Doping,mesh1.Node,mesh1.Lx,what_sc,In,...
    mes, mhs, Egs, CnAus, CpAus, Brs, NCs, NVs, Aes, Bes, ...
    Ahs, Bhs, betas, Eps,DeltaEAs, DeltaEDs, Ne_refs, Np_refs, ...
    Vnsats, Vpsats, Ebslons, uns, ups, tor1s, tor2s, n_is, Xis,...
    alphas, alphas_abs, eta_ns, eta_ps,TL, Tr,const.eps_silicon,Qvalues,delta_Fermi,k_L);

R_um = Device_Diameter*1e4;
Area_um2 = pi*R_um^2/4;
thicknesses_um = thicknesses(Type_Layer==0)*1e-3;
capacitances = 1e-4*Ebslons(Type_Layer==0)*Area_um2./thicknesses_um; % 1e-4* because eps in F/cm
                                                                     % needs to be converted to F/um
total_capacitance_fF = 1/(sum(1./capacitances))*1e15;
total_resistance = sum(resistivity(Type_Layer==0).*thicknesses_um);
PD_parallel_load = R_load*total_resistance/(total_resistance+R_load);
%     disp(['Surface Area: ' num2str(Area_um2) ' um^2']);
    disp(['Capacitance: ' num2str(total_capacitance_fF) ' (fF)'])
    disp(['Resistance: ' num2str(total_resistance) ' (Ohm)'])


% % If you want to plot band diagrams, uncomment the next line
% xx = mesh1.Lx*mesh1.NX;
% xxx = mesh1.Lx(mesh1.Node)*mesh1.NX;
% f3 = -PD_Str1.Xi-PD_Str1.Eg/2+PD_Str1.delta_Fermis;
% f4 = f3-max(f3);
% f1 = -PD_Str1.Xi-f4;
% f2 = -PD_Str1.Xi-PD_Str1.Eg-f4;
% figure(77); 
% plot(xx,f1,xx,f2)
% xlim([0 mesh1.Lx(end)*mesh1.NX]);
% xlabel('z (m)');
% ylabel('Band Diagram (eV)');
% % % text(1.5e-4, -5.6,'E_v');
% % % text(1.5e-4, -4.5,'E_c');
%  keyboard
In = PD_Str1.In;
Xi=PD_Str1.Xi;
me=PD_Str1.mstar_n;
mh=PD_Str1.mstar_p;
Ae=PD_Str1.Ae;
Be=PD_Str1.Be;
Ah=PD_Str1.Ah;
Bh=PD_Str1.Bh;
k_L=PD_Str1.k_L;
DeltaEA=PD_Str1.DeltaEA;
DeltaED=PD_Str1.DeltaED;
Ne_ref=PD_Str1.Ne_ref;
Np_ref=PD_Str1.Np_ref;
Ebslon=PD_Str1.Epsilons;
NC=PD_Str1.NC;
NV=PD_Str1.NV;
ni=PD_Str1.ni;
tn=PD_Str1.tn;
tp=PD_Str1.tp;
Tol=PD_Str1.Tol/NP.Ni;
Vns=PD_Str1.Vns;
Vps=PD_Str1.Vps;
un=PD_Str1.un;
up=PD_Str1.up;
Beta=PD_Str1.Beta;
CnAu=PD_Str1.CnAu;
CpAu=PD_Str1.CpAu;
Br=PD_Str1.Br;
Eg=PD_Str1.Eg;
Xi=PD_Str1.Xi;
Ep=PD_Str1.Ep;
Diameter(1:N,1)= Beam_Diameter; % Beam Diameter unit cm
DeviceDiameter(1:N,1)= Device_Diameter; 
N1=mesh1.Node(interfaces(1)); % N1 is the p-i interface
N2=mesh1.Node(interfaces(end)); % N2 is the i-n interface;