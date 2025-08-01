function Sol_Dyn = Obtain_Modulated_Solution(Design, parameters, Sol)
disp('Calculations started: CW excitation with a single freq modulation');
initialize_dynamic_computations;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For sensitivity calculation use smaller dt for accurate integration
dt = dt0/10;
tt = 0:dt:Totaltime0;
Gint = G0*(1+parameters.mod_depth*sin(2*pi*parameters.mod_freq*tt));
Np = (sum(Gint.*dt)-Gint(1)*dt(1)/2)/NP.VT;
Sol_Dyn.Np = Np;
%keyboard
clear tt;
clear Gint;
% %%%%%%%%%%%%%%%%%%%%%%
mod_function  = 1+parameters.mod_depth*sin(2*pi*parameters.mod_freq*time_rec);       % used in dynamic_loop
dynamic_loop;
%%%% PLOT THE IMPULSE RESPONSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    keyboard

    return

    if parameters.plot_results==1
        figure(parameters.filenum+1000);
        plot(x,Jvt,'k',x,J_elec,'r', x,J_hole,'b',x,Jdis,'g')
        legend('J_{total}','J_{elec}','J_{hole}','J_{dis}');
        xlabel('Time (ps)');
        ylabel('Current Densities');
        title('DID NOT CONVERGE')
    end
else

    Jmean = mean(Jvt(round(length(Jvt)/2):end));
    Jmax = max(abs(Jvt(round(length(Jvt)/2):end)-Jmean));
    Sol_Dyn.RFpower = 10*log10(Jmax^2*parameters.R_load/2e-3);    % 2 comes from RMS

    Ave_I = mean(abs(Jvt.*dt0s));
    Sol_Dyn.Ave_I = Ave_I;
end
%%%%%%% RESPONSIVITY CALCULATION %%%%
% q/(h*c) = 8.0647e+05;
% h*c/q = 1.24e-6;
% Resp = Qeff*8.0647e+05*lambda
Ne = abs(sum((J_elec).*dt0s));
Nh = abs(sum((J_hole).*dt0s));

Sol_Dyn.Nh = Nh;
Sol_Dyn.Ne = Ne;
Sol_Dyn.Qeff = (Ne+Nh)/Np;
Sol_Dyn.Responsivity = Sol_Dyn.Qeff*8.0647e+05*wavelength;

if parameters.show_timestepping == 1
    disp(['Responsivity (A/W) = ' num2str(Sol_Dyn.Responsivity)])
    disp(['Quantum Efficiency = ' num2str(Sol_Dyn.Qeff)])
end
%  keyboard

Sol_Dyn.completed = 1;
if parameters.plot_results==1

    figure(parameters.filenum+20); clf;
    plot(x,Jvt,'k',x,J_elec,'r', x,J_hole,'b',x,Jdis,'g');
    legend('total','elec.','hole','disp.','Location','Best');
    grid on;
    xlabel('{\it{t}} (ps)');
    ylabel('Currents Densities');

    %disp('684')
    if parameters.print_plots ==1
        figure(parameters.filenum+20);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.04, 0.04, 0.5, 0.5]);
        oldscreenunits = get(gcf,'Units');
        oldpaperunits = get(gcf,'PaperUnits');
        oldpaperpos = get(gcf,'PaperPosition');
        set(gcf,'Units','pixels');
        scrpos = get(gcf,'Position');
        newpos = scrpos/100;
        set(gcf,'PaperUnits','inches','PaperPosition',newpos)
        print('-dpng', [Solution_Name '_modulated'], '-r300');
        drawnow
        set(gcf,'Units',oldscreenunits,'PaperUnits',oldpaperunits,'PaperPosition',oldpaperpos)
        !mv *.png ./figures/
    end
end

% SAVE THE RESULTS
save(['./Results/' parameters.Solution_Name '_modulated'],'Sol_Dyn','Sol')
% keyboard
end
