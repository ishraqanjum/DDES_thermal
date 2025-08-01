classdef Norm_parameter < handle
    % Parameter normalization
    properties
        KB;
        q;
        Ebslon;
        Ni;
        T;
        Diff_coef;
        VT % temperature
        NX
        D0
        NT
        Nmu
        J0
    end
    properties (Constant)
    end
    
    methods
        function obj=Norm_parameter(KB,q,Ebslon,Ni,T,Diff_coef)
            obj.KB=KB;
            obj.q=q;
            obj.Ebslon=Ebslon;
            obj.Diff_coef = Diff_coef;
            obj.Ni=Ni;
            obj.T=T;
            
            obj.VT=KB*T/q;                                  % normalized voltage or potential
            obj.NX=sqrt(Ebslon/q*obj.VT/Ni);       % normalized position x
            obj.D0=Diff_coef*obj.VT;                    %normalized diffusion coeificient
            obj.NT=obj.NX^2/obj.D0;                   % normalized time t
            obj.Nmu=obj.D0/obj.VT;                     % normalized mobility;
            obj.J0=q*obj.D0*Ni/obj.NX;            
        end
    end
end

