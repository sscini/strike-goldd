% STRIKE-GOLDD Model Creation Script for Cell Calcium Dynamics
% Based on the provided Python model (a two-pool calcium oscillation system).
% This script defines the structural model for identifiability analysis.

clear;

% 2 states (State variables)
% Z: Cytosolic Calcium (IC)
% Y: ER/Pool Calcium (IC)
syms Z Y
x = [Z; Y];

% 1 output (Measured variable)
% Assuming Cytosolic Calcium (Z) is the measured output
h = Z;

% Inputs
u = []; % Known inputs (Assuming no experimental stimuli)
w = []; % Unknown inputs (Assuming no unknown disturbances)

% Unknown parameters (The constants and rates that are typically estimated)
% Hill coefficients included for full structural analysis
syms beta mu0 mu1 vm2 vm3 k2 kr ka k kf n m p_hill
p = [beta; mu0; mu1; vm2; vm3; k2; kr; ka; k; kf; n; m; p_hill];

% --- Intermediate terms (Auxiliary functions) ---
% vin = mu_0 + mu_1 * beta
vin = mu0 + mu1 * beta;

% v2 (Calcium uptake into Y, depends on Z)
% v2 = (vm2 * (Z ^ n)) / (k2 ^ n + (Z) ^ n)
v2 = (vm2 * (Z^2)) / (k2^2 + (Z)^2);

% v3 (Calcium release from Y, depends on Z and Y)
% v3 = ((vm3 * (Z ^ p_hill)) / (ka ^ p_hill + (Z) ^ p_hill)) * ((Y ^ m) / (kr ^ m + Y ^ m))
v3 = ((vm3 * (Z)^p_hill) / (ka^p_hill + (Z)^p_hill)) * ((Y^m) / (kr^m + Y^m));


% --- Dynamic equations (f) ---
% dZ = vin - v2 + v3 + kf * Y - k * Z
% dY = v2 - v3 - kf * Y
f = [vin - v2 + v3 + kf * Y - k * Z;
     v2 - v3 - kf * Y];

% --- Initial Conditions ---
% % Initial conditions that are UNKNOWN parameters (ics)
% % The initial concentration of the non-measured pool (Y) is typically unknown.
% syms Y0 
% ics = [Y0]; % Y(0) is unknown

% Initial conditions that are KNOWN (known_ics)
% The initial concentration of the measured state (Z) is usually considered known (or set to 0).
syms Z_initial Y_initial
% known_ics is a vector of size 'size(x)' where the i-th element is a known value for x(i), or 0 if x(i) is unknown (defined in ics).
known_ics = [Z_initial; Y_initial]; % Z(0) = Z_initial (a known value), Y(0) is parameterized by Y_initial.

% Save the model structure to the file 'calcium_model.mat'
save('twopoolCa_SSC_withhill','x','p','h','u','w','f','known_ics');