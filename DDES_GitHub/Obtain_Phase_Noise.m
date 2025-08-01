function Sol_Dyn = Obtain_Phase_Noise(Design, parameters, Sol)
disp('Calculations started: CW excitation with a small perturbation');
initialize_dynamic_computations;
f_r = parameters.mod_freq;
T_R = 1/f_r;
t_c = parameters.t_c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For sensitivity calculation use smaller dt for accurate integration
dt = dt0/10;
tt = 0:dt:Totaltime0;
kx = (tt-parameters.pulse_pos)/parameters.pulse_width;
Gint = G0.*(1+parameters.mod_depth*sech(kx));
Np = (sum(Gint.*dt)-Gint(1)*dt(1)/2)/NP.VT;
Sol_Dyn.Np = Np;
%keyboard
clear tt;
clear Gint;
clear kx
% %%%%%%%%%%%%%%%%%%%%%%
mod_function = 1+parameters.mod_depth*sech((time_rec-parameters.pulse_pos)/parameters.pulse_width);
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
    %%%%%%%% PHASE NOISE CALCULATION %%%%%%%%
    t_n = time_rec(time_rec<T_R);
    dt_n = diff(t_n);
    h_t = zeros(1,length(t_n));
    h_t(1:min([ti0, length(t_n)])) = c_f*h_norm(1:min([ti0, length(t_n)]));
    abs_error = 1;
    n = parameters.harmonic_order;
    %
    while abs(abs_error)>=1e-9
        S_n = sin((2*n*pi./T_R).*(t_n-t_c));
        C_n = cos((2*n*pi./T_R).*(t_n-t_c));
        F_n1 = h_t.*S_n;
        F_n2 = -(2*n*pi./T_R).*h_t.*C_n;
        F_n = sum(dt_n.*(F_n1(1:end-1)+F_n1(2:end))/2);
        dF_n = sum(dt_n.*(F_n2(1:end-1)+F_n2(2:end))/2);
        t_c = t_c - F_n/dF_n;
        abs_error=abs(F_n);
    end
    argu = 2*pi*n*(t_n-t_c)./T_R;
    numerator = sin(argu).^2.*h_t;
    G = sum(dt_n.*(numerator(1:end-1)+numerator(2:end))/2);
    denumerator = cos(argu).*h_t;
    G1 = (sum(dt_n.*(denumerator(1:end-1)+denumerator(2:end))/2))^2;

    Ave_I = sum(abs(Jvt.*dt0s))/T_R;
    phasenoise= 10*log10(const.eV*G/2/Ave_I/G1);
    Sol_Dyn.phasenoise = phasenoise;
    Sol_Dyn.Ave_I = Ave_I;

    phasenoise_comb = zeros(1,length(parameters.comb_numbers));
    for n = parameters.comb_numbers
        abs_error = 1;
        while abs(abs_error)>=1e-9
            S_n = sin((2*n*pi./T_R).*(t_n-t_c));
            C_n = cos((2*n*pi./T_R).*(t_n-t_c));
            F_n1 = h_t.*S_n;
            F_n2 = -(2*n*pi./T_R).*h_t.*C_n;
            F_n = sum(dt_n.*(F_n1(1:end-1)+F_n1(2:end))/2);
            dF_n = sum(dt_n.*(F_n2(1:end-1)+F_n2(2:end))/2);
            t_c = t_c - F_n/dF_n;
            abs_error=abs(F_n);
        end
        argu = 2*pi*n*(t_n-t_c)./T_R;
        numerator = sin(argu).^2.*h_t;
        G = sum(dt_n.*(numerator(1:end-1)+numerator(2:end))/2);
        denumerator = cos(argu).*h_t;
        G1 = (sum(dt_n.*(denumerator(1:end-1)+denumerator(2:end))/2))^2;

        Ave_I_calc = sum(abs(Jvt.*dt0s))/T_R;
        phasenoise_comb(n) = 10*log10(const.eV*G/2/Ave_I_calc/G1);
    end
    if length(phasenoise_comb(imag(phasenoise_comb)~=0))>0
        disp('Imaginary phase noise! Calculations are halted!')
    figure(parameters.filenum+10); clf;
    plot(x,-(Jvt-Js),'k',x,-(J_elec-Jsn),'r', x,-(J_hole-Jsp),'b',x,-Jdis,'g');
    % plot(x,-(Jvt),'k',x,-(J_elec),'r', x,-(J_hole),'b',x,-Jdis,'g');
    legend('total','elec.','hole','disp.','Location','Best');
    grid on;
    xlabel('{\it{t}} (ps)');
    ylabel('Currents Perturbations');
        
        keyboard
        Sol_Dyn.completed = 0;
        return
    else
        Sol_Dyn.phasenoise_comb = phasenoise_comb;
    end
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
    disp(['Phase noise (dBc/Hz):' num2str(phasenoise)])
end
%  keyboard

Sol_Dyn.completed = 1;
if parameters.plot_results==1

    figure(parameters.filenum+10);
    subplot(121);
    plot(x,-(Jvt-Js),'k',x,-(J_elec-Jsn),'r', x,-(J_hole-Jsp),'b',x,-Jdis,'g');
    legend('total','elec.','hole','disp.','Location','Best');
    grid on;
    xlabel('{\it{t}} (ps)');
    ylabel('Currents Perturbations');

    xx = parameters.comb_numbers;
    yy = phasenoise_comb;
    n = parameters.harmonic_order;
    xp = 1.05*xx(n)*parameters.mod_freq/1e9;
    yp = yy(n);
    xp2 = mean(xx*parameters.mod_freq/1e9);
    yp2 = min(yy) + (max(yy)-min(yy))/20;

    subplot(122);
    plot(xx*parameters.mod_freq/1e9,yy,'-o')
    grid on;
    xlabel('Comb-Line Frequency (GHz)');
    text(xp,yp, ['\rightarrow ' num2str(round(10*yp)/10)] )
    text(xp2,yp2,['{\it{Q}}_{eff} = ', num2str(round(100*Sol_Dyn.Qeff)/100)]);
    grid on;
    ylabel('{\it{PN}} (dBc/Hz)')

    %disp('684')
    if parameters.print_plots ==1
        figure(parameters.filenum+10);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.04, 0.04, 0.9, 0.5]);
        oldscreenunits = get(gcf,'Units');
        oldpaperunits = get(gcf,'PaperUnits');
        oldpaperpos = get(gcf,'PaperPosition');
        set(gcf,'Units','pixels');
        scrpos = get(gcf,'Position');
        newpos = scrpos/100;
        set(gcf,'PaperUnits','inches','PaperPosition',newpos)
        print('-dpng', [Solution_Name '_single_freq'], '-r300');
        drawnow
        set(gcf,'Units',oldscreenunits,'PaperUnits',oldpaperunits,'PaperPosition',oldpaperpos)
        !mv *.png ./figures/
    end
end

% SAVE THE RESULTS
save(['./Results/' parameters.Solution_Name '_single_freq'],'Sol_Dyn','Sol')
end