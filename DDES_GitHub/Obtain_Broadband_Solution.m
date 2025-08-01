function Sol_Dyn = Obtain_Broadband_Solution(Design, parameters, Sol)
% This code calculates the broadband response by adding a broadband
% modulation to incident ligth, i.e. InputLight = G0*(1+broadband_pulse).
% Pulse is centered at T/4
% This code shouldn't be used for single frequency modulation or phase
% noise calculations
%%%%%%%%%% default settings %%%%%%%
disp('Calculations started: CW excitation modulated with a broadband window');
initialize_dynamic_computations;
BroadbandMod = broadband_pulse(parameters.window_max_frequency,time_rec,Totaltime0*0.25,parameters.window_type);
BroadbandMod = BroadbandMod-min(BroadbandMod);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For sensitivity calculation use smaller dt for accurate integration
dt = dt0/10;
tt = 0:dt:Totaltime0;
Wbhs = broadband_pulse(parameters.window_max_frequency,time_rec,Totaltime0*0.25,parameters.window_type);
Wbhs = Wbhs-min(Wbhs);
Gint = G0.*(1+Wbhs);
% Np = f_ave*(sum(Gint.*dt)-Gint(1)*dt(1)/2)/P0_normalizer;
% note that G0 already includes f_ave!
Np = (sum(Gint.*dt)-Gint(1)*dt(1)/2)/NP.VT;
Sol_Dyn.Np = Np;
Np = (sum(Gint.*dt)-Gint(1)*dt(1)/2)/NP.VT;
Sol_Dyn.Np = Np;
% clear the variables that won't be used
clear tt;
clear Gint;
clear Wbhs;
% %%%%%%%%%%%%%%%%%%%%%%
mod_function =1+BroadbandMod;       % used in dynamic_loop
dynamic_loop;
%%%% PLOT RESULTS $$$$$$$$$$
x = time_rec(1:length(Jdis))*1e12;
Jvt= J_elec + J_hole + Jdis;
h_norm = Jvt-Js;
c_f = 1/sum(dt0s(1:length(h_norm)).*h_norm);

Sol_Dyn.Jvt = Jvt;
Sol_Dyn.Jdis = Jdis;
Sol_Dyn.J_hole = J_hole;
Sol_Dyn.J_elec = J_elec;
Sol_Dyn.Js = Js;
Sol_Dyn.Jsn = Jsn;
Sol_Dyn.Jsp = Jsp;
Sol_Dyn.c_f = c_f;

if R_load<0
    disp('A convergent issue has been experienced during DTE. Calculations are halted!')
    Sol_Dyn.completed = 0;    
    return
end

fi = 0:parameters.window_max_frequency/100:parameters.window_max_frequency;
[f_0, P1, P1notShifted, P1bare] = get_broadband_RF_Output(parameters,fi,Sol_Dyn, Sol);
Sol_Dyn.P1notShifted = P1notShifted;
Sol_Dyn.P1bare = P1bare;
Sol_Dyn.f_0 = f_0;
Sol_Dyn.p1 = P1;

f3dB = interp1(P1notShifted,f_0, P1notShifted(1)-3,'linear')/1e9;
Sol_Dyn.f3dB = f3dB;
disp(['Bandwidth (GHz):' num2str(f3dB)]);
Sol_Dyn.completed = 1;
if parameters.plot_results==1

    figure(parameters.filenum+1);
    subplot(121);
    plot(x,Jvt,'k',x,J_elec,'r', x,J_hole,'b',x,Jdis,'g');
    legend('total','elec.','hole','disp.','Location','Best');
    grid on;
    xlabel('{\it{t}} (ps)');
    ylabel('Currents Densities');
    subplot(122);
    plot(f_0/1e9,P1notShifted, f3dB, P1notShifted(1)-3,'ro');
    grid on;
    xlabel('{\it{f}} (GHz)')
    ylabel('{\it{P}}_{out} (dB)');

    %disp('684')
    if parameters.print_plots ==1
        figure(parameters.filenum+1);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.04, 0.04, 0.9, 0.5]);
        oldscreenunits = get(gcf,'Units');
        oldpaperunits = get(gcf,'PaperUnits');
        oldpaperpos = get(gcf,'PaperPosition');
        set(gcf,'Units','pixels');
        scrpos = get(gcf,'Position');
        newpos = scrpos/100;
        set(gcf,'PaperUnits','inches','PaperPosition',newpos)
        print('-dpng', [Solution_Name '_broadband'], '-r300');
        drawnow
        set(gcf,'Units',oldscreenunits,'PaperUnits',oldpaperunits,'PaperPosition',oldpaperpos)
        !mv *.png ./figures/
    end
end

% SAVE THE RESULTS
save(['./Results/' parameters.Solution_Name '_broadband'],'Sol_Dyn','Sol')

%%%%%%%%%% BROADBAND PULSE %%%%%%%%%%
