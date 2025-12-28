%% ---------- Main code that runs SABEMMT for a given propeller -----------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Script : main.m
%  Project: sabemmt
%  Date: 12/27/2025
%  Author: Miguel Frade
%  License: BSD 3-Clause
%
%  Description:
%     Driver script to setup and execute the SABEMMT (Structures And Blade 
%     Element Modified Momentum Theory) solver for a specific propeller 
%     geometry. This script handles:
%        1. Geometry Definition (Approx. APC 11x5.5 Sport Propeller)
%        2. Airfoil Data Loading & Pre-processing
%        3. Simulation Environment Setup
%        4. Solver Execution
%        5. Visualization of Results
%
%  Dependencies: 
%     - runSABEMMT.m
%     - getAirfoilDataE214_MultiRe_cleaned.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% ========================================================================
%  1. DEFINE SIMULATION PARAMETERS
%  ========================================================================
%  Operation Point
v_inf = 15;        % Freestream velocity [m/s]
rpm   = 6000;      % Rotational speed [RPM]

%  Environment (Sea Level Standard)
environment.rho = 1.225;          % Air density [kg/m^3]
environment.nu  = 1.46e-5;        % Kinematic viscosity [m^2/s]
environment.sound_speed = 340.29; % Speed of sound [m/s]

%  Propeller Geometry (Approx. APC 11x5.5 Sport)
R = (11 * 0.0254) / 2;  % Radius [m] (11 inches diameter)
Nb = 2;                 % Number of Blades

% Define Radial Stations (Mesh)
N_stations = 200;
mesh = linspace(0.15, 1, N_stations); % Normalized radius r/R (hub to tip)
dx = mesh(2) - mesh(1);

% Chord Distribution c(r) [m] - Tapered Planform
% Simple quadratic approximation for a sport prop
c_root = 0.022; 
c_max  = 0.034; 
c_tip  = 0.012;
% Shape function peaking around r/R = 0.4
c = spline([0.15, 0.4, 1.0], [c_root, c_max, c_tip], mesh);

% Twist Distribution theta(r) [deg] - Geometric Pitch
% Pitch (P) is roughly constant for fixed pitch props: P = 2*pi*r * tan(theta)
% For an 11x5.5, P ~ 5.5 inches = 0.1397 m
Pitch_geom = 5.5 * 0.0254; 
theta_deg = atan(Pitch_geom ./ (2 * pi * (mesh * R))) * (180/pi);

%  Material Properties (e.g., Glass Fiber Reinforced Nylon)
rho_mat = 1600; % [kg/m^3]


%% ========================================================================
%  2. AIRFOIL DATA LOADING
%  ========================================================================
% Load the Eppler E214C-PT data
fprintf('Loading airfoil data...\n');
prop.airfoil_data = getAirfoilDataE214_MultiRe_cleaned();



%% ========================================================================
%  3. PRE-CALCULATE AERODYNAMICS AND INERTIAS TO MAKE SABEMMT FASTER
%  ========================================================================
% Aerodynamics (Load & Interpolate)
airfoil_data = prop.airfoil_data; 
% 'linear' ensures smooth transitions between Reynolds numbers
Fc = griddedInterpolant({airfoil_data.alpha_deg*(pi/180), airfoil_data.Re}, ...
                        airfoil_data.cl_matrix, 'linear', 'linear');                      
Fd = griddedInterpolant({airfoil_data.alpha_deg*(pi/180), airfoil_data.Re}, ...
                        airfoil_data.cd_matrix, 'linear', 'linear');

% Linear Aerodynamics (for SABEMMT initialization)
% Finds the linear lift slope (Cla) and zero-lift angle (alpha_0)
nom_idx = find(airfoil_data.Re >= 0.6 * max(airfoil_data.Re), 1);
alpha_lin = airfoil_data.alpha_deg * (pi/180);
cl_lin = airfoil_data.cl_matrix(:, nom_idx);
lin_region = alpha_lin >= deg2rad(-1) & alpha_lin <= deg2rad(6);
p_poly = polyfit(alpha_lin(lin_region), cl_lin(lin_region), 1);
cla_data = p_poly(1); 
cl0_data = p_poly(2);

% Structural Properties (Normalized per chord length)
% Assumes a pseudo-shell or solid section based on area
c_norm_struct = 1; 
x_poly = [airfoil_data.x_u, fliplr(airfoil_data.x_l)];
y_poly = [airfoil_data.y_u, fliplr(airfoil_data.y_l)];
A_norm = polyarea(x_poly, y_poly);

% Simplified inertia approximation (Thin airfoil theory / flat plate equivalent)
t_norm = A_norm / c_norm_struct; 
I_y_norm = t_norm / 12;         % Flapwise inertia factor
I_x_norm = t_norm^3 / 12;       % Chordwise (Lead-lag) inertia factor

% Pack values into prop_const structure for the solver
prop_const.airfoil.Fc = Fc; 
prop_const.airfoil.Fd = Fd;
prop_const.airfoil.cla_data = cla_data; 
prop_const.airfoil.cl0_data = cl0_data;
prop_const.airfoil.alpha_min = min(alpha_lin); 
prop_const.airfoil.alpha_max = max(alpha_lin);
prop_const.airfoil.cd_max_val = max(airfoil_data.cd_matrix(:));
prop_const.airfoil.A_norm = A_norm; 
prop_const.airfoil.t_norm = t_norm;
prop_const.airfoil.I_x_norm = I_x_norm; 
prop_const.airfoil.I_y_norm = I_y_norm;

% Add mesh and material to prop_const
prop_const.rho_mat = rho_mat;
prop_const.mesh = mesh;
prop_const.dx = dx;


%% ========================================================================
%  4. RUN SOLVER (SABEMMT)
%  ========================================================================
fprintf('Running BEMT Solver for V=%.1f m/s, RPM=%.0f...\n', v_inf, rpm);

try
    result = runSABEMMT(v_inf, rpm, Nb, R, c, theta_deg, prop_const, environment);
    
    % Display Summary
    fprintf('\n--- RESULTS ---\n');
    fprintf('Thrust (T):       %.4f N\n', result.T);
    fprintf('Torque (Q):       %.4f Nm\n', result.Q);
    fprintf('Power (P):        %.4f W\n', result.P);
    fprintf('Efficiency (eta): %.2f %%\n', result.eta_p * 100);
    fprintf('Max Stress:       %.2f MPa\n', max(abs(result.sigma_total_max)) / 1e6);
    fprintf('Tip Mach:         %.3f\n', result.M_tip);
    
catch ME
    fprintf('Error during execution: %s\n', ME.message);
    return;
end



%% ========================================================================
%  5. PLOTTING
%  ========================================================================
figure('Name', 'SABEMMT Propeller Analysis', 'Color', 'w', 'Position', [100 100 1000 600]);

% Plot 1: Distributed Loads
subplot(2,2,1);
plot(result.r, result.dTdx, 'b-', 'LineWidth', 1.5); hold on;
plot(result.r, result.dFtdx, 'r--', 'LineWidth', 1.5);
xlabel('Radius $r$ [m]', 'Interpreter', 'latex'); 
ylabel('Force Distribution [N/m]', 'Interpreter', 'latex');
legend({'Thrust Dist. ($dT/dx$)', 'Tangential Force Dist. ($dF_t/dx$)'}, ...
       'Interpreter', 'latex', 'Location', 'northwest');
title('Aerodynamic Load Distribution', 'Interpreter', 'latex');
grid on;

% Plot 2: Angle of Attack
subplot(2,2,2);
plot(result.r, rad2deg(result.alpha), 'k-', 'LineWidth', 1.5); hold on;
yline(rad2deg(prop_const.airfoil.alpha_max), 'r--');
yline(rad2deg(prop_const.airfoil.alpha_min), 'r--');
xlabel('Radius $r$ [m]', 'Interpreter', 'latex'); 
ylabel('Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
title('Angle of Attack Distribution', 'Interpreter', 'latex');
legend({'$\alpha$', 'Stall Limits'}, 'Interpreter', 'latex', 'Location', 'northwest');
grid on;

% Plot 3: Structural Stress
subplot(2,2,3);
plot(result.r, result.sigma_total_max / 1e6, 'm-', 'LineWidth', 1.5);
xlabel('Radius $r$ [m]', 'Interpreter', 'latex'); 
ylabel('Max Total Stress $\sigma_{max}$ [MPa]', 'Interpreter', 'latex');
title('Structural Stress Est.', 'Interpreter', 'latex');
grid on;

% Plot 4: Geometry Visualization
subplot(2,2,4);
yyaxis left
plot(result.r, c * 1000, 'b-');
ylabel('Chord $c$ [mm]', 'Interpreter', 'latex');
yyaxis right
plot(result.r, theta_deg, 'r-');
ylabel('Twist $\theta$ [deg]', 'Interpreter', 'latex');
xlabel('Radius $r$ [m]', 'Interpreter', 'latex');
title('Blade Geometry', 'Interpreter', 'latex');
grid on;


%% ========================================================================
%  6. SAVE RESULTS AND FIGURE
%  ========================================================================
% Ensure the results folder exists
resultsDir = fullfile(pwd, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

% Save MATLAB data
save(fullfile(resultsDir, 'SABEMMT_simulation_results.mat'), ...
     'result', 'v_inf', 'rpm', 'Nb', 'R', 'c', 'theta_deg', ...
     'prop_const', 'environment');

% Save figure BEFORE showing it
fig = gcf;
set(fig, 'PaperPositionMode', 'auto');

print(fig, ...
      fullfile(resultsDir, 'SABEMMT_simulation_plot.jpg'), ...
      '-djpeg', '-r300');

% Show figure AFTER it is saved
drawnow;