%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------- Main code that runs SABEMMT for a given propeller -----------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script : main.m
%   Project: sabemmt
%   Author: Miguel Frade
%   Affiliation at time of publication: Universidad Politecnica de Madrid
%   Date: December 2025
%   License: Creative Commons Attribution-NonCommercial 4.0 (CC BY-NC 4.0)
%
%   Description:
%      Driver script to setup and execute the SABEMMT (Structures And Blade 
%      Element Modified Momentum Theory) solver for a specific propeller 
%      geometry. This script handles:
%         1. Geometry Definition (Approx. APC 11x5.5 Sport Propeller)
%         2. Airfoil Data Loading & Pre-processing
%         3. Simulation Environment Setup
%         4. Solver Execution
%         5. Visualization of Results
% 
%   Dependencies: 
%      - runSABEMMT.m
%      - getAirfoilDataE214_MultiRe_cleaned.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
addpath('functions\');


%% ========================================================================
%  SECTION 0: LOGGING SETUP (SAVE EVERYTHING THAT APPEARS IN THE COMMAND
%             WINDOW AS .TXT)
%  ========================================================================
% 1. Create directory if needed
if ~exist('results', 'dir')
    mkdir('results');
end
% 2. Define log file path
logFile = fullfile('results', 'results.txt');
% 3. Delete existing log so we start fresh (otherwise diary appends)
if exist(logFile, 'file')
    delete(logFile);
end
% 4. Start recording to file
diary(logFile);


%% ========================================================================
%  1. DEFINE SIMULATION PARAMETERS
%  ========================================================================
%  ------------- Name that will appear on top of the plots ----------------
prop_const.name = 'APC 11x5.5 Propeller';

%  ------------------------- Operation Point ------------------------------
v = 15;        % Flight speed [m/s]
rpm = 7000;    % Rotational speed of the propeller [RPM]

fprintf('\n==================== FLIGHT CONDITION =====================\n');
fprintf(' Flight speed and RPM: V=%.1f m/s, RPM=%.0f\n', v, rpm);

%  ----------------- Environment (Sea Level Standard) ---------------------
environment.rho = 1.225;          % Air density [kg/m^3]
environment.nu  = 1.46e-5;        % Kinematic viscosity [m^2/s]
environment.sound_speed = 340.29; % Speed of sound [m/s]

%  ----------- Propeller Geometry (Approx. APC 11x5.5 Sport) --------------
R = (11 * 0.0254) / 2;  % Radius [m] (11 inches diameter)
Nb = 2;                 % Number of Blades

%  --- Define Radial Stations of the aerodynamic blade (Mesh) ---
N_stations = 200;
x0 = 0.1; % Start of the aerodynamic blade (root). 
          % At smaller x, no aerodynamics, but the geometry will be
          % automatically filled for 3D plot.
mesh = linspace(x0, 1, N_stations); % Normalized radius r/R (hub to tip)
dx = gradient(mesh); % SABEMMT also works for vector dx (general mesh)

% --- Define where the structural reinforcement finishes ---
% The sections closer to the root will have structural reinforcement,
% so we won't worry about them.
prop_const.x_finish_reinforcement = 0.18; % This is a realistic value

%  --- Chord Distribution c(r) [m] ---
% Simple quadratic approximation for a sport prop
c_root = 0.022; 
c_max  = 0.034; 
c_tip  = 0.012;
% Shape function peaking around r/R = 0.4
c = spline([x0, 0.4, 1.0], [c_root, c_max, c_tip], mesh);

%  --- Twist Distribution theta(r) [deg] - Geometric Pitch ---
% For an 11x5.5, P ~ 5.5 inches = 0.1397 m
Pitch_geom = 5.5 * 0.0254; 

% Since: P = 2*pi*r(x=0.75)*tan(theta_charact) = 2*pi*0.75*R*tan(theta_c) ,
% we can calculate theta at the characteristic section (x = 0.75) as:
theta_deg_charact = atan(Pitch_geom / (2 * pi * 0.75 * R)) * (180/pi);

% But we will assume that pitch is roughly constant throughout the blade,
% so that we can directly calculate the hole pitch distribution using:
% P = 2*pi*r * tan(gamma) = 2*pi*mesh*R * tan(theta)
theta_deg = atan(Pitch_geom ./ (2 * pi * (mesh * R))) * (180/pi);

%  --- Material Properties (e.g., Glass Fiber Reinforced Nylon) ---
rho_mat = 1600; % [kg/m^3]


%% ========================================================================
%  2. AIRFOIL DATA LOADING
%  ========================================================================
% Load the Eppler E214C-PT data
fprintf('\n===========================================================\n');
fprintf('Loading airfoil data...\n');
prop_const.airfoil_data = getAirfoilDataMultiReE214;
% prop_const.airfoil_data = getAirfoilDataMultiReNACA4412;
% prop_const.airfoil_data = getAirfoilDataMultiReClarkY;


%% ========================================================================
%  3. PRE-CALCULATE AERODYNAMICS AND INERTIAS TO MAKE SABEMMT FASTER
%  ========================================================================
% ------------------ Aerodynamics (Load & Interpolate) --------------------
airfoil_data = prop_const.airfoil_data; 
% 'linear' ensures smooth transitions between alphas and Reynolds numbers
Fc = griddedInterpolant({airfoil_data.alpha_deg*(pi/180), airfoil_data.Re}, ...
                        airfoil_data.cl_matrix, 'linear', 'linear');                      
Fd = griddedInterpolant({airfoil_data.alpha_deg*(pi/180), airfoil_data.Re}, ...
                        airfoil_data.cd_matrix, 'linear', 'linear');

% --------- Linear Aerodynamics (for SABEMMT loop initialization) ---------
% Finds the linear lift slope (Cla) and zero-lift angle (alpha_0)
nom_idx = find(airfoil_data.Re >= 0.6 * max(airfoil_data.Re), 1);
alpha_lin = airfoil_data.alpha_deg * (pi/180);
cl_lin = airfoil_data.cl_matrix(:, nom_idx);
lin_region = alpha_lin >= deg2rad(-1) & alpha_lin <= deg2rad(6);
p_poly = polyfit(alpha_lin(lin_region), cl_lin(lin_region), 1);
cla_data = p_poly(1); 
cl0_data = p_poly(2);

% ---------- Structural Properties (Real Airfoil Geometry) ----------------
% Computes exact A, Ix, Iy, Ixy, and centroid from the coordinate data
% instead of using the equivalent rectangle approximation.
fprintf('\n===========================================================\n');
fprintf('Computing airfoil geometric properties...\n');
geo = getAirfoilInertias(airfoil_data.x_u, airfoil_data.y_u, ...
                             airfoil_data.x_l, airfoil_data.y_l);


%% ========================================================================
%  4. PACK VALUES INTO prop_const STRUCTURE FOR THE SOLVER
%  ========================================================================
% Airfoil aerodynamic characteristics
prop_const.airfoil.Fc = Fc; 
prop_const.airfoil.Fd = Fd;
prop_const.airfoil.cla_data = cla_data; % To init the loop in SABEMMT
prop_const.airfoil.cl0_data = cl0_data;
prop_const.airfoil.alpha_min = min(alpha_lin); 
prop_const.airfoil.alpha_max = max(alpha_lin);
prop_const.airfoil.cd_max_val = max(airfoil_data.cd_matrix(:));

% Airfoil geometry and inertias
prop_const.structural = geo;

% Add mesh, Nb and material to prop_const
prop_const.rho_mat = rho_mat;
prop_const.Nb = Nb;
prop_const.mesh = mesh;
prop_const.dx = dx;


%% ========================================================================
%  5. RUN SOLVER (SABEMMT)
%  ========================================================================
fprintf('\n===========================================================\n');
fprintf('Running SABEMMT Solver for V=%.1f m/s, RPM=%.0f...\n', v, rpm);

try
    result = runSABEMMT(v, rpm, Nb, R, c, theta_deg, prop_const, environment);

    % Only consider the max stress in the region outside the structural
    % reinforcement:
    % Filter: Only consider sections outside the reinforcement zone
    valid_region_mask = result.mesh >= prop_const.x_finish_reinforcement;
    % Create a temporary array for finding the max
    search_stress = result.sigma_total_max;
    % Set stress in the excluded region to -Infinity so they are never picked
    % (unless the entire blade is excluded, which would be an error in inputs)
    search_stress(~valid_region_mask) = -inf;
    % Find the max stress and its index within the valid region
    [max_stress_val, idx_crit] = max(search_stress);
    % If the max stress section is x_finish_reinforcement within a dx 
    % margin (which it probably is), let's make the indication of the most
    % critical section match exactly x_finish_reinforcement:
    % First, we find the corresponding dx (it can be a vector):
    if isscalar(result.dx)
        dx_eff = result.dx;
    else
        dx_eff = result.dx(idx_crit);
    end 
    if result.mesh(idx_crit) > prop_const.x_finish_reinforcement - dx_eff && ...
       result.mesh(idx_crit) < prop_const.x_finish_reinforcement + dx_eff
        x_critical = prop_const.x_finish_reinforcement;
    else
        x_critical = result.mesh(idx_crit);
    end
    
    % Display Summary
    fprintf('\n========================= RESULTS =========================\n');
    fprintf('Thrust (T):         %.4f N\n', result.T);
    fprintf('Torque (Q):         %.4f Nm\n', result.Q);
    fprintf('Power (P):          %.4f W\n', result.P);
    fprintf('Efficiency (eta_p): %.2f %%\n', result.eta_p * 100);
    fprintf('Max Stress:         %.2f MPa\n', max_stress_val / 1e6);
    fprintf('(Max Stress section is: r/R =%.2f )\n', x_critical);
    fprintf('Tip Mach:           %.3f\n', result.M_tip);
    
catch ME
    fprintf('Error during execution: %s\n', ME.message);
    return;
end

% STOP LOGGING
diary off;
fprintf('Full command window output saved to %s\n', logFile);



%% ========================================================================
%  6. PLOTTING
%  ========================================================================
%  The function plotGeometryAndPerformancesSingle saves all plots as .jpg
fprintf('\n===================== PLOT GENERATION =====================\n');
reply = input('\n  > Generate plots for this solution? (y/n): ','s');
if strcmpi(reply,'y')
    fprintf('  > Generating plots...\n');
    plotGeometryAndPerformancesSingle(v, rpm, Nb, R, c, theta_deg, prop_const, environment);
else
    fprintf(' > Plot generation skipped\n');
end


