clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Design = {
   1, 'InGaAs', 50, 2e19, 0, 'Zn'; 
1,'InP',  500, 1.5e18, 0, 'Zn';
1,'InGaAsP', 10, 2e18, 1.10, 'Zn';
1,'InGaAsP', 10, 2e18, 1.40, 'Zn';
1,'InGaAs', 50, 2e18, 2e17, 'Zn';
1,'InGaAs', 40, 1e17, 0, 'Zn';
1,'InGaAs', 10, 4e17, 0, 'Zn';
0, 'InGaAs', 10, 1, 0, 'Si';
0, 'InGaAsP', 10, 1, 1.50, 'Si'; 
0, 'InGaAsP', 10, 1, 1.15, 'Si';  
0, 'InP', 10, 5e17 , 0, 'Si';  
0, 'InP', 40, 5e16 , 0, 'Si';  
0, 'InP',  100, 1e16, 0, 'Si';  
-1,'InGaAsP', 350, 1.5e18, 1.25, 'Si';  
-1,'InP', 1000, 1e19, 0, 'Si'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Device parameters
params.Beam_Diameter = 9.57e-4;        % (cm) r0 in Eq. 9 of Yue's paper
params.Device_Diameter = 9.57e-4;      % (cm) D0 in Eq. 10 of Yue's paper%%%%%%%%%%%%%%%%%%%%%%%%%%
params.Bias= 1.2;                                  % reverse bias voltage (V)
params.R_load=50;                                % load resistance (ohms)
params.Temperature=300;                     % normalization temperature (F)

% Physical Effects to be included
params.Impact=1;                    % 1 for impact ionization or 0 for no Imp. Ion.
params.Thermionic_e=1;              % 1 for thermionic field emission for electrons
params.Thermionic_h=1;              % 1 for thermionic field emission for holes (CAUSED NAN)
params.FranzKeldysh = 1;            % 1 to include FK effect
params.Bleaching = 0;               % 1 include bleaching; 0: not include
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Excitation
params.wavelength = 1.55e-6;        % wavelength of the light (m)
params.P0_factor = 1e-5;               % Excitation strength Determines Generation

params.mod_freq = 1e9;
params.mod_depth = 0.04;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver settings
params.initial_guess_method = 0;
% Simulation Settings
params.Dmax=0.1;                           %  Mesh quality (unit nm)

% phase noise parameters and pulse settings                                        
params.comb_numbers = 1:100;            % phase noise comb parameters
params.dt = 1/params.mod_freq/max(params.comb_numbers)/10; 
params.Totaltime= 1000*params.dt;    
params.pulse_pos = params.Totaltime*0.4;                % center of the pulse (s)
params.pulse_width = 10*params.dt/2.634;  % full width at half-maximum pulse duration for sech input
params.harmonic_order = 5;              % n in Eq. (7) of Ehsan's OPEX paper
params.t_c = params.Totaltime/5;      % t_c in Eq. (9) of Ehsan's OPEX paper

% Broadband Parameters
params.bb_mode = 1;                                 % broadband mode ==> 1: windows, 2: sinc
params.window_type = 'BlackmanHarris';     % required if mode == 1
params.window_max_frequency = 50e9;

% Settings to calculate, plot and save
params.filenum = 10;                   % figure number
params.plot_results = 1;                % 1: plot results. 0: do not plot results
params.print_plots = 1;                 % set this to 1 if you want to print the main figure
params.show_timestepping = 1;           % set this to 1 if you want to see the progress
params.Solution_Name = 'z_Grzeslo';

% DO NOT CHANGE THESE PARAMETERS (NOT IMPLEMENTED YET)
params.side = 1;           % 1: n-side, 2:p-side, 3:wg mode perp, 4: wg mode parallel
params.In = 1;              % exponent power, m in Eq. 23 of Yue Hu's paper
params.gamma=1;        % exponent power, gamma in Eqs. 4 and 5 of Yue Hu's paper
% Yue Hu's paper: 10.1109/JLT.2014.2315740
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sol] = Obtain_Static_Solution(Design, params);   
if Sol.converged ==1       
    [Sol_Dyn2] = Obtain_Broadband_Solution(Design, params,Sol);    

    [Sol_Dyn1] = Obtain_Phase_Noise(Design, params,Sol);
    % let's change the modulation depth seriously, which requires much
    % finer Delta_T
    params.mod_depth = 0.78;                
    params.dt = 1/max([params.mod_freq, 1e5])/1000;   % for lin integration, delta_t
    params.Totaltime= 10000*params.dt;                         % total simulation time
    [Sol_Dyn3] = Obtain_Modulated_Solution(Design, params,Sol);        
end