clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Type SemiCond thickness Doping Q-value or Doping for graded layers
Design = {1, 'InGaAs',  50,  2.0e19,    0, 'Zn';  % 1
    1, 'InP',    100,  1.5e18,    0, 'Zn';  % 2
    1, 'InGaAsP', 15,  2.0e18,  1.1, 'Zn';  % 3
    1, 'InGaAsP', 15,  2.0e18,  1.4, 'Zn';  % 4
    1, 'InGaAs', 100,  2.0e18,    0, 'Zn';  % 5
    1, 'InGaAs', 150,  1.2e18,    0, 'Zn';  % 6
    1, 'InGaAs', 200,  8.0e17,    0, 'Zn';  % 7
    1, 'InGaAs', 250,  5.0e17,    0, 'Zn';  % 8
    0, 'InGaAs', 150,  1.0e16,    0, 'Si';  % 8    
    0, 'InGaAsP', 15,  1e16,  1.4, 'Si';  % 9
    0, 'InGaAsP', 15,  1e16,  1.1, 'Si';  % 10
    0, 'InP',     50,  1.4e17,    0, 'Si';  % 11
    0, 'InP',    900,  1.0e16,    0, 'Si';  % 12
    -1, 'InP',    100,  1.0e18,    0, 'Si';  % 13
    -1, 'InP',    900,  1.0e19,    0, 'Si';  % 14
    -1, 'InGaAs', 20,  1.0e19,    0, 'Si';  % 15    
    -1, 'InP',    200,  1.0e19,    0, 'Si'}; % 16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Device parameters
params.Beam_Diameter = 30e-4;        % (cm) r0 in Eq. 9 of Yue's paper
params.Device_Diameter = 28e-4;      % (cm) D0 in Eq. 10 of Yue's paper
%%%%%%%%%%%%%%%%%%%%%%%%%%
params.Bias= 5;                            % Reverse Bias voltage (V)
params.R_load=50;                       % load resistance (ohms)
params.Temperature=300;            % normalization temperature (F)

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

params.mod_freq = 0;                % not used in broadband calculations
params.mod_depth = 0;              % not used in broadband calculations

% Broadband Parameters
params.bb_mode = 1;                                 % broadband mode ==> 1: windows, 2: sinc
params.window_type = 'BlackmanHarris';     % required if mode == 1
params.window_max_frequency = 50e9;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver settings
params.initial_guess_method = 1;
% Simulation Settings
params.Dmax=1;                                                          %  Mesh quality (unit nm)
params.dt = 1/params.window_max_frequency/200;      % for lin integration, delta_t
params.Totaltime= 2000*params.dt;                         % total simulation time 

% Settings to calculate, plot and save
params.filenum = 10;                   % figure number
params.plot_results = 1;                % 1: plot results. 0: do not plot results
params.print_plots = 1;                 % set this to 1 if you want to print the main figure
params.show_timestepping = 1;           % set this to 1 if you want to see the progress
params.Solution_Name = 'z_mutc4broadband';


% DO NOT CHANGE THESE PARAMETERS (NOT IMPLEMENTED YET)
params.side = 1;           % 1: n-side, 2:p-side, 3:wg mode perp, 4: wg mode parallel
params.In = 1;              % exponent power, m in Eq. 23 of Yue Hu's paper
params.gamma=1;        % exponent power, gamma in Eqs. 4 and 5 of Yue Hu's paper
% Yue Hu's paper: 10.1109/JLT.2014.2315740
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sol] = Obtain_Static_Solution(Design, params);   
if Sol.converged ==1
    [Sol_Dyn] = Obtain_Broadband_Solution(Design, params,Sol);
end
