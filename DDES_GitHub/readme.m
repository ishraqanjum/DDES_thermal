% 2023/05/25
% This version assumes the photodetector is reverse biased
% and excited from the n-side.
% If you want to excite the structure from other sides,
% then you will need to make a lot of changes in the Drift_Diff codes
% because the current implementation assumes electron move
% from p side to n side and holes move in the opposite direction
%
%  2023/05/11
% the bug in InGaAsP Q calculation has been fixed
% 
%  2023/05/03
% since when the material is lossy, the E field decays as it propagates
% when we set alpha to 0 for Eg>Einc, the previous implementations had this
% small bug. now it is fixed.
% it should be fixed in other versions as well!
%
% v4.6 what's new
% 1. si and ge material properties are added
% 2. windowing + non-uniform time stepping requires signals to be
% interpolated over uniformly sampled t, to calculate the FFT properly
% ..
% Nth iteration assumes incidence power is P0, 
% then the dynamic solver solves two cases simulatenously
% first one assumes a constant excitation with a perturbation
% second one assumes a constant excitation + broandband signal
