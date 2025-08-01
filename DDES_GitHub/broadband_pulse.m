function p = broadband_pulse(fr,t,pulse_center,window_type)
% inputs
% a: window coefficients
% fr: maximum frequency of interst of interest
% t: time-steps, from 0:dt:Tmax
% pulse_center: where we woould like to have the center of pulse
% output
% p: normalized pulse amplitude, i.e. maximum will be equal to 1.
T = 1.14/fr;  % let's choose something sligtly larger than 1/fr
%
switch window_type
    case 'Nuttall3der' % continuous third derivative
        a = [0.338946 0.481973 0.161054  0.018027];
    case 'Nuttall1der' % continuous first derivative
        a = [0.35875 0.487396 0.144232 0.012604];
    case 'BlackmanHarris' % Minimum Four-Term
        a = [0.35322222 0.48829 0.145 0.010222222];
    case 'FlatTop'
        a = [0.21557895 0.41663158 0.277263158 0.083578947 0.006947368];
end

window = zeros(size(t));
for n=0:length(a)-1
    window = window + a(n+1)*cos(2*n*pi*(t-pulse_center)./T);
end
window(abs(t-pulse_center)>=T/2) = 0;
p = window(:)';
p = p./max(abs(p));


