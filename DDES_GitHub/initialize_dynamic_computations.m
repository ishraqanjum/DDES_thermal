% Set the default values and formats
set(0,'defaultlinelinewidth',1.8);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultTextFontSize',18);
format short e;
warning('off')
%%%% initialization %%%
Sol_Dyn.phasenoise = 0;
Sol_Dyn.phasenoise_comb = 0;
Sol_Dyn.P1notShifted = 0;
Sol_Dyn.P1bare = 0;
Sol_Dyn.Ne = 0;
Sol_Dyn.Np = 0;
Sol_Dyn.Responsivity = 0;
Sol_Dyn.Jvt = 0;
Sol_Dyn.Jdis = 0;
Sol_Dyn.J_hole = 0;
Sol_Dyn.J_elec = 0;
Sol_Dyn.Js = 0;
Sol_Dyn.Jsn = 0;
Sol_Dyn.Jsp = 0;
% outputs
Sol_Dyn.Ave_I = 0;
Sol_Dyn.Qeff = 0;
Sol_Dyn.RFpower = 0;
Sol_Dyn.f_0 = 0;
Sol_Dyn.p1 = 0;
Sol_Dyn.completed = 0;

if min(Sol.p)<0 || min(Sol.n)<0 || max(abs(imag(Sol.p)))>0 || max(abs(imag(Sol.n)))
    Sol_Dyn.completed = 0;
    return
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load the constants
    constant;
    % load parameters for the Dynamic solution
    Totaltime0 = parameters.Totaltime;
    dt0 = parameters.dt;
    Bias = parameters.Bias;
    Impact = parameters.Impact; % impact ionization or not
    %%% Excitation
    wavelength = parameters.wavelength;         %wavelength of the light
    R_load = parameters.R_load;
    % Set R_load to a small value of 0 entered, to make the code more stable
    if R_load == 0
        R_load = 1e-3;
    end
    Temperature = parameters.Temperature;            % normlization temperature

    % Device parameters
    Device_Diameter = parameters.Device_Diameter;
    Beam_Diameter = parameters.Beam_Diameter;
    Dmax = parameters.Dmax;                 % Parameter to control the mesh quality (nm)
    % Solver settings
    In = parameters.In;                     % exponent power, m in Eq. 23 of Yue Hu's paper
    gamma = parameters.gamma;               %
    % Settings to plot and save the results
    filenum = parameters.filenum;
    plot_results = parameters.plot_results;
    Solution_Name  = parameters.Solution_Name;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load Steady-State Solution variables
    pm = Sol.pm;
    nm = Sol.nm;
    wm = Sol.wm;
    Js = Sol.Js_C;
    Jn = Sol.Js_n;
    Jp = Sol.Js_p;
    TL = Sol.TL;

    E = Sol.E;
    PD_Str1 = Sol.PD_Str;
    pb1 = Sol.Boundary.pb1;
    pbn = Sol.Boundary.pbn;
    nb1 = Sol.Boundary.nb1;
    nbn = Sol.Boundary.nbn;
    wb1 = Sol.Boundary.wb1;
    wbn = Sol.Boundary.wbn;
    wbf = Sol.Boundary.wbf;

    %keyboard
    NP = Sol.Parameters.NP;
    fd1 = Sol.Parameters.fd1;
    mesh1 = Sol.Parameters.mesh1;
    ss1 = Sol.Parameters.ss1;     
    NAP = Sol.NAP;
    NDP = Sol.NDP;

    % 
    alpha0 = mean(PD_Str1.alpha(PD_Str1.alpha>100));
    xxgh = 1;
    Sol_Dyn.Mfactor = 1;
    Sol_Dyn.xx = xxgh;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % if this solution will be used for pulsed excitation, we need
    % to take care of normalization
    pulse_normalizer = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P0=parameters.P0_factor;   %  note that pulse_normalizer is 1, if it is a cont. excitation
    Parameters_NR % parameter setting
    NX = mesh1.NX;
    %
    Eph=const.h*const.c/wavelength;
    D_effective = min([Beam_Diameter, Device_Diameter]); % Note that here the unit is CM
    A=pi*0.25*D_effective^2;

    rs = linspace(0,Device_Diameter,1001);
    ff = exp(-(rs/Beam_Diameter).^2);
    f_ave = sum(ff*(rs(2)-rs(1)))/Device_Diameter/1.001;
    clear rs;
    clear ff;

    P0_normalizer = 1/(A*Eph)*alpha0*NP.NT^2;
    G0 = f_ave*P0.*P0_normalizer;
    Sol_Dyn.G0= G0;
    Sol_Dyn.f_ave = f_ave;

    structure_loss = ones(N,1);
    if parameters.side ==1  %  n-side illumination
        cell_loss = exp(-PD_Str1.alphas_abs.*NX.*[mesh1.dx; 0]);
        structure_loss = flipud(cumprod(flipud(cell_loss)));
    elseif parameters.side ==2  %  p-side illumination
        cell_loss = exp(-PD_Str1.alphas_abs.*NX.*[mesh1.dx; 0]);
        structure_loss = cumprod(cell_loss);
    elseif parameters.side ==3  % WG mode perp illumination
        structure_loss = structure_loss.*[1; PD_Str1.Eb_re./max(PD_Str1.Eb_re)];
    end

    GL(1:N,1)=G0.*structure_loss; % per volume
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% DYNAMIC PART %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Jsn=Jn;
    Jsp=Jp;
    pm0=pm;
    nm0=nm;
    E0=E;
    NAP0 = NAP;
    NDP0 = NDP;
    
    NT=NP.NT;               % normalized time
    J0=NP.J0;               % normalized current

    Steps=round(Totaltime0/dt0);
    time_rec = linspace(0,Totaltime0,Steps+1);
    Sol_Dyn.time_rec = time_rec;

    dt0s = diff(time_rec);
    Steptime1s=dt0s/NT;
    DeviceDiameter=ss1.DeviceDiameter;
    Ebslon=PD_Str1.Epsilons;

    L1=(length(pm));
    J_elec=zeros(1,Steps);
    J_hole=zeros(1,Steps);
    Jdis=zeros(1,Steps);
    djterm = DeviceDiameter(1:end-1).^2./sum(mesh1.dx)*0.25*pi; % pi*r^2/L = pi*R^2/4/L
end