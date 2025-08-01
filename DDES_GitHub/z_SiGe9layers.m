clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Type SemiCond thickness Doping Q-value or Doping for graded layers
Design = {1, 'Ge',  50,  8.0e18, 4e18, 'B';
          1, 'Ge',  150,  2.0e18, 6e17, 'B';
          1, 'Ge',  175,  6.0e17, 7e16, 'B';
          1, 'Ge',  75,  7.0e16, 8e16, 'B';
          1, 'Ge',  100,  8.0e16, 6e14, 'B';       
          1, 'Ge',  30,  6e14, 7e15, 'B';          
          1, 'Ge',  20,  7e15, 6e14, 'B';          
          0, 'Si',  300,  1e14,   1e14, 'As';
         -1, 'Si', 1000,  1e18,    0, 'As'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Device parameters
params.Beam_Diameter = 15e-4;        % (cm) r0 in Eq. 9 of Yue's paper
params.Device_Diameter = 15e-4;      % (cm) D0 in Eq. 10 of Yue's paper%%%%%%%%%%%%%%%%%%%%%%%%%%
params.wavelength = 1.55e-6;        % wavelength of the light (m)
params.Bias= 1;                     % Bias voltage (V)
params.R_load=50;                   % load resistance (ohms)
params.Temperature=300;             % normalization temperature (F)

% Physical Effects to be included
params.Impact=1;                    % 1 for impact ionization or 0 for no Imp. Ion.
params.Thermionic_e=1;              % 1 for thermionic field emission for electrons
params.Thermionic_h=1;              % 1 for thermionic field emission for holes (CAUSED NAN)
params.FranzKeldysh = 1;            % 1 to include FK effect
params.Bleaching = 0;               % 1 include bleaching; 0: not include
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Excitation
params.wavelength = 1.55e-6;        % wavelength of the light (m)
params.P0_factor = 1e-4;               % Excitation strength Determines Generation

params.mod_freq = 1e9;
params.mod_depth = 0.04;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver settings
params.initial_guess_method = -1;
% Simulation Settings
params.Dmax=1;                           %  Mesh quality (unit nm)

% phase noise parameters and pulse settings                                        
params.comb_numbers = 1:100;            % phase noise comb parameters
params.dt = 1/params.mod_freq/max(params.comb_numbers)/20;   % for lin integration, delta_t
params.Totaltime= 1000*params.dt;                         % total simulation time 
params.pulse_pos = params.Totaltime*0.4;                % center of the pulse (s)
params.pulse_width = 1e-12/2.634;  % full width at half-maximum pulse duration for sech input
params.harmonic_order = 5;              % n in Eq. (7) of Ehsan's OPEX paper
params.t_c = params.Totaltime/5;      % t_c in Eq. (9) of Ehsan's OPEX paper

% Broadband Parameters
params.bb_mode = 1;                                 % broadband mode ==> 1: windows, 2: sinc
params.window_type = 'BlackmanHarris';     % required if mode == 1
params.window_max_frequency = 50e9;

% Settings to calculate, plot and save
params.filenum = 304;                   % figure number
params.plot_results = 1;                % 1: plot results. 0: do not plot results
params.print_plots = 1;                 % set this to 1 if you want to print the main figure
params.show_timestepping = 1;           % set this to 1 if you want to see the progress
params.Solution_Name = 'z_SiGe';

% DO NOT CHANGE THESE PARAMETERS (NOT IMPLEMENTED YET)
params.side = 1;           % 1: n-side, 2:p-side, 3:wg mode perp, 4: wg mode parallel
params.In = 1;              % exponent power, m in Eq. 23 of Yue Hu's paper
params.gamma=1;        % exponent power, gamma in Eqs. 4 and 5 of Yue Hu's paper
% Yue Hu's paper: 10.1109/JLT.2014.2315740
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sol] = Obtain_Static_Solution(Design, params);   
if Sol.converged ==1    
    [Sol_Dyn1] = Obtain_Phase_Noise(Design, params,Sol);
    [Sol_Dyn2] = Obtain_Broadband_Solution(Design, params,Sol);
    
    % let's change the modulation depth and time stepping for seriously
    % modulated CW excitation
    params.mod_depth = 0.78;                
    params.dt = 1/params.mod_freq/100; 
    params.Totaltime= 1000*params.dt;    
    [Sol_Dyn3] = Obtain_Modulated_Solution(Design, params,Sol);        
end