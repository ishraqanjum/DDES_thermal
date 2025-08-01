%%%%%%%%%% RF Output Power %%%%%%%%%%
function [f_0, P1, P1notShifted, P1bare] = get_broadband_RF_Output(params,fi,Sol_Dyn, Sol)
R = params.R_load;     % ohm
ti = Sol_Dyn.time_rec;
p0 = broadband_pulse(params.window_max_frequency,ti,params.Totaltime*0.25, params.window_type);
[fpulse_0, P_pulse_0] = getAmpSpec(p0, Sol);
[f_0, P_0] = getPowerSpec(Sol_Dyn, Sol, R, fi);
normalizer = interp1(fpulse_0,P_pulse_0,f_0,'spline');
P_0= P_0./normalizer.^2;
P1bare=10*log10(P_0/1e-3);
P1 = 10*log10(P_0/1e-3)-10*log10(P_0(1)/1e-3);
P1notShifted = P1+10*log10(Sol_Dyn.Js^2*Sol.params.R_load/2e-3);

% figure(23);
% subplot(221); plot(f_0/1e9, P1bare)
% subplot(222);plot(f_0/1e9, P1); xlim([0 25])
% subplot(223);plot(f_0/1e9, normalizer); xlim([0 25])
% subplot(224);plot(f_0/1e9, P1notShifted)

%%%%%%%%%%  %%%%%%%%%%
function [f, W] = getAmpSpec(p, Sol)
% compute amplitude spectrum of source pulse and corresponding frequency vector
n = 2^nextpow2(length(p));
W = abs(fftshift(fft(p,n)));
W = W./max(W);
fn = 0.5/Sol.params.dt;
df = 2.*fn/n;
f = -fn:df:fn-df;
W = W(n/2+1:end);
f = f(n/2+1:end);

%%%%%%%%%%  %%%%%%%%%%
function [f1, P1] = getPowerSpec(Sol_Dyn, Sol, R,fi)
% this code requires J to be sampled linearly
% however logarithmic integration breaks this rule
% this is why we need to an interpolation
% J = Sol_Dyn.DJvt-Sol_Dyn.DJs;
% normalize broadband current so that it will have same amplitude as the
% single frequency
% XR = max(abs(Sol_Dyn.Jvt-Sol_Dyn.Js));
% J = (Sol_Dyn.Jvt-Sol_Dyn.Js)/XR*abs(Sol_Dyn.Js);
J = (Sol_Dyn.Jvt-Sol_Dyn.Js);
L = length(J);
tt = Sol_Dyn.time_rec;
n_periods = max(tt)*Sol.params.window_max_frequency;

T = tt(2)-tt(1);
Fs = 1/T;

f1 = Fs*(0:(L/2))/L;
Y = abs(fft(J,L));
P1 = 2*R*(Y(1:L/2+1)*n_periods/L).^2;
P1i = interp1(f1, P1, fi,'spline');
f1 = fi;
P1 = P1i;
