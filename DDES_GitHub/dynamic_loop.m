ti = 0;
while ti<Steps
    if parameters.show_timestepping==1 && mod(ti,100)==0; disp([ti, Steps]); end
    ti = ti +1;
    ti0 = ti;
    a=(-1/Steptime1s(ti)).*[ones(2*L1,1);zeros(L1+1,1)];
    Mt1=spdiags(a,0,3*L1+1,3*L1+1);
    Mt=Mt1;
    NTSteptime = NT*Steptime1s(ti);

    %%%%%%%%%%%%%%%%%%% CW or PULSED %%%%%%%%%%%%%%%%%%
    err1 = 1;
    err2_p = 1;
    Gin = GL.*mod_function(ti);

    DD2 = Drift_Diff_N_C_B_R(pb1,pbn,nb1,nbn,wbf,wbn,fd1,mesh1,ss1,NP,Gin,PD_Str1,...
        TL,Bias,R_load,DeviceDiameter,wavelength);
    for i=1:20
        if err1 < 1e-14 && err2_p<1e-11
            break
        end        
        [Fp,Fn,Fw,Fb,J1,J2,~,~,~,~,~,~,~,~,NAP,NDP]=DD2.Cal_Currentv2(pm,nm,wm,wb1,Gin);
        Js11=trapz(mesh1.Lx_half,J1*J0.*djterm);
        Js12=trapz(mesh1.Lx_half,J2*J0.*djterm);
        Fy = -[Fp;Fn;Fw;Fb];
        Ft=([pm-pm0;nm-nm0;wm-wm;0])/Steptime1s(ti);

        Fdt=([NAP(2:end-1)-NAP0(2:end-1);NDP(2:end-1)-NDP0(2:end-1);wm-wm;0])/Steptime1s(ti);
        Fy=Fy+Ft-Fdt;
        [M,M1] = DD2.ConsJacobian(pm,nm,wm,wb1); % Calculate Jacobian
        M=M+Mt+M1./Steptime1s(ti);
        Deltay = M\Fy; % Calculate Newton step
        if any(isnan(Deltay))
            disp('isnan: Loop1')
            keyboard
            Sol_Dyn.completed = 0;
            Jdis = Jdis*0;
            J_hole = J_hole*0;
            J_elec = J_elec*0;
            ti= Steps+101;
            i = 100;
            err1 = 1;
            err2_p = 1;
            Js1 = 1e10;
            R_load = -1e10;
            break
        end

        % Update p,n and w
        pm = pm+Deltay(1:L1);
        nm = nm+Deltay(L1+1:2*L1);
        wm = wm+Deltay(2*L1+1:3*L1);
        wb1= wb1+Deltay(end);

        err1=norm(Deltay)/norm([pm;nm;wm]);
        err2_p=norm(Fp-(pm-pm0-NAP(2:end-1)+NAP0(2:end-1))/Steptime1s(ti));

        w=NP.VT*[wb1;wm;DD2.wbn];
        E=-(fd1.diffmp(w))/NP.NX;
        Jdiss = 0;
        if max(abs(E-E0))/max(abs(E0))>NTSteptime/max(djterm)
            Jdiss=trapz(mesh1.Lx_half,Ebslon.*(E-E0)/NTSteptime.*djterm);
        end

        Js1=Js11+Js12+Jdiss;
    end
    
    if R_load<0
        keyboard
        Js=Js1;
        break;
    end
    Js=Js1;
    Jsn=Js12;
    Jsp=Js11;

    J_elec(ti)=trapz(mesh1.Lx_half,J2*J0.*djterm);%...
    J_hole(ti)=trapz(mesh1.Lx_half,J1*J0.*djterm);%...
    Jdis(ti)=trapz(mesh1.Lx_half,Ebslon.*(E-E0)/NTSteptime.*djterm);
    E0=E;
    pm0=pm;
    nm0=nm;

    NAP0=NAP;
    NDP0=NDP;    
end