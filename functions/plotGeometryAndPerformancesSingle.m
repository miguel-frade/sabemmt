%% ------------------- Plot Geometry and Performances ---------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: plotGeometryAndPerformancesSingle
%   Project: sabemmt
%   Author: Miguel Frade
%   Affiliation at time of publication: Universidad Politecnica de Madrid
%   Date: December 2025
%   License: Creative Commons Attribution-NonCommercial 4.0 (CC BY-NC 4.0)
%
%   Description:
%       A comprehensive visualization suite for propeller analysis and 
%       design.
%       Generates and saves the following plots to the 'results/' folder:
%
%       1. Propeller geometry in 3D (with Hub and Transition Shank)
%       2. Propeller Geometry Analysis (Chord/Pitch dist, Planform, 3D Mesh)
%       3. Airfoil Geometry and Aerodynamics (Inflow geometry & Polars)
%       4. Propeller Aerodynamics & Forces (Alpha/Lambda, dT/dx, dFt/dx)
%       5. Stress Distributions (Centrifugal, Bending, and Total Max Stress)
%       6. Airfoil Stress Distribution (Detailed stress map at 3 sections)
%       7. Performance Clouds (Efficiency, Ct, Cp vs Advance Ratio)
%
%   Inputs:
%       v, rpm, Nb, R, c_dist, theta_deg, prop, env
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotGeometryAndPerformancesSingle(v, rpm, Nb, R, c_dist, theta_deg, prop, env)

    % Run SABEMMT to get the outputs neccesary for plots
    outputs = runSABEMMT(v, rpm, Nb, R, c_dist, theta_deg, prop, env);

    % Construct the 3 lines for the title of pictures
    line1 = sprintf('\\underline{\\textbf{%s}}', prop.name);
    line2 = sprintf('Operating condition: $v=%.1f~m/s$, $RPM=%.0f$', v, outputs.rpm);
    line3 = sprintf('Performances: $T=%.2f~N$, $\\eta_p=%.2f$', outputs.T, outputs.eta_p);
    
    % Store combined lines
    titles = {line1; line2; line3};


    % ================= [TOP BLOCK: START] =================
    % Don't show the figures until they are all saved.
    % Tell MATLAB to always restore visibility, even if an error occurs.

    % Store current default settings to restore them later
    defaultVis = get(0, 'DefaultFigureVisible');
    
    % Force all subsequent figures to be invisible initially
    set(0, 'DefaultFigureVisible', 'off');
    
    % Safety mechanism: If code crashes, restore visibility automatically
    cleanupVis = onCleanup(@() set(0, 'DefaultFigureVisible', defaultVis));
    % ================== [TOP BLOCK: END] ==================



    % =======================================================================
    % 1. Entire propeller Geometry in 3D (Hub matches Blade Color)
    % =======================================================================
    figure('Name', '1. Propeller geometry in 3D', 'Color', 'w', ...
           'Renderer', 'opengl'); 
    hold on; grid on; axis equal;
    
    % 1. Setup Geometry Inputs
    mesh = outputs.mesh;
    x0 = outputs.mesh(1);
    R = outputs.R;

    x_hub = 0.07; % I want the hub to be always drawed up to 7% of the radius
    % The aerodynamic sections can start at greater x (x0 approx. 0.2), but
    % a smooth transition will be drawn from the hub to the aerodynamic blade.

    x_u = prop.airfoil_data.x_u;
    x_l = prop.airfoil_data.x_l;
    y_u = prop.airfoil_data.y_u;
    y_l = prop.airfoil_data.y_l;
    x_poly = [x_u, fliplr(x_l)]; 
    y_poly = [y_u, fliplr(y_l)];
    x_axis_norm = 0.25; % Stacking Axis
    
    % Pre-allocate Grids
    num_points_foil = length(x_poly);
    
    % 2. Generate Aerodynamic Blade Mesh (The original blade)
    r_phys_blade = R .* mesh; 
    c_blade = outputs.c;
    theta_rad_blade = outputs.theta;
    
    num_blade_sections = length(r_phys_blade);
    X_blade = zeros(num_points_foil, num_blade_sections);
    Y_blade = zeros(num_points_foil, num_blade_sections);
    Z_blade = zeros(num_points_foil, num_blade_sections);
    for i = 1:num_blade_sections
        radius = r_phys_blade(i);
        chord  = c_blade(i);
        twist  = -theta_rad_blade(i); 
        
        % Scale and Center Airfoil
        xx_local = (x_poly - x_axis_norm) * chord;
        yy_local = y_poly * chord;
        
        % Apply Twist Rotation
        x_rot = xx_local * cos(twist) - yy_local * sin(twist);
        z_rot = xx_local * sin(twist) + yy_local * cos(twist);
        
        % Map to Global Coordinates
        X_blade(:, i) = x_rot;
        Y_blade(:, i) = radius;
        Z_blade(:, i) = z_rot;
    end

    % ---------------------------------------------------------
    % 2B. Generate Transition Shank (Hub to Blade Root)
    % ---------------------------------------------------------
    X_trans = []; Y_trans = []; Z_trans = [];
    
    if x0 > x_hub
        % Settings for the transition
        num_trans_sections = 15; 
        r_trans = linspace(0.91*x_hub*R, x0*R, num_trans_sections);
        % I use 0.91*x_hub*R to make sure there is no gap between hub and
        % transition shank. Don't use less than 0.91, because the blade
        % "enters" the hub and it is not aesthetic.
        
        % Get geometry of the blade root (Target Shape)
        c_root = c_blade(1);
        twist_root = -theta_rad_blade(1);
        
        % 1. Define Root Airfoil Local Coordinates (Target)
        xx_root_local = (x_poly - x_axis_norm) * c_root;
        yy_root_local = y_poly * c_root;
        
        % 2. Define Hub Circle Local Coordinates (Source)
        % We project the airfoil points onto a circle to ensure 
        % point-to-point correspondence (topology matching).
        
        % Calculate thickness of root to sizing the shank cylinder
        root_thickness = max(y_poly) - min(y_poly);
        shank_radius = (root_thickness * c_root) * 1.8; % <--- Radius of intersection (can be modified)
        
        % Calculate angles of the airfoil points relative to centroid
        angles = atan2(yy_root_local, xx_root_local);
        
        % Create circle points based on those angles
        xx_circ_local = shank_radius .* cos(angles);
        yy_circ_local = shank_radius .* sin(angles);
        
        % Initialize grids
        X_trans = zeros(num_points_foil, num_trans_sections);
        Y_trans = zeros(num_points_foil, num_trans_sections);
        Z_trans = zeros(num_points_foil, num_trans_sections);
        
        for i = 1:num_trans_sections
            % Interpolation factor (0 = Hub, 1 = Blade Root)
            % Using a sinusoidal easing for visual smoothness
            s_lin = (i-1) / (num_trans_sections-1);
            s = sin(s_lin * pi/2); % Erase-out
            
            radius = r_trans(i);
            
            % Interpolate Shapes (Morph Circle -> Airfoil)
            xx_current = xx_circ_local * (1-s) + xx_root_local * s;
            yy_current = yy_circ_local * (1-s) + yy_root_local * s;
            
            % Apply Twist
            % We keep the twist constant (equal to root twist) for the 
            % shank so it enters the hub cleanly.
            twist = twist_root; 
            
            x_rot = xx_current * cos(twist) - yy_current * sin(twist);
            z_rot = xx_current * sin(twist) + yy_current * cos(twist);
            
            X_trans(:, i) = x_rot;
            Y_trans(:, i) = radius;
            Z_trans(:, i) = z_rot;
        end
        
        % Remove the last section of transition to avoid duplicate points 
        % at the join with the main blade
        X_trans = X_trans(:, 1:end-1);
        Y_trans = Y_trans(:, 1:end-1);
        Z_trans = Z_trans(:, 1:end-1);
    end
    
    % ---------------------------------------------------------
    % 2C. Combine Grids
    % ---------------------------------------------------------
    X_grid = [X_trans, X_blade];
    Y_grid = [Y_trans, Y_blade];
    Z_grid = [Z_trans, Z_blade];
    
    num_sections = size(X_grid, 2);

    % --- Aesthetic Settings (UPDATED) ---
    upper_color = [0.15 0.35 0.75];   % Blue (upper surface)
    lower_color = [0.80 0.55 0.20];   % Bronze (lower surface)
    
    % Set Hub Color to match Upper Blade Color
    prop_color  = upper_color;         
    
    spec_strength = 0.6;            
    diff_strength = 0.8;            
    
    % --- 3. Draw Blades using PATCH ---
    faces = [];
    vertices = [];
    colors = [];
    
    nv = 0;
    
    for i = 1:num_sections-1
        for j = 1:num_points_foil-1
            % Quad vertices (current blade)
            v1 = [X_grid(j,i),   Y_grid(j,i),   Z_grid(j,i)];
            v2 = [X_grid(j+1,i), Y_grid(j+1,i), Z_grid(j+1,i)];
            v3 = [X_grid(j+1,i+1), Y_grid(j+1,i+1), Z_grid(j+1,i+1)];
            v4 = [X_grid(j,i+1), Y_grid(j,i+1), Z_grid(j,i+1)];
    
            vertices = [vertices; v1; v2; v3; v4];
            faces = [faces; nv+1 nv+2 nv+3 nv+4];
    
            % Upper vs lower by airfoil index
            if j <= length(x_u)
                colors = [colors; upper_color];
            else
                colors = [colors; lower_color];
            end
    
            nv = nv + 4;
        end
    end
    
    patch('Faces',faces,'Vertices',vertices, ...
          'FaceVertexCData',colors, ...
          'FaceColor','flat', ...
          'CDataMapping','direct', ...
          'EdgeColor','none', ...
          'FaceLighting','gouraud', ...
          'BackFaceLighting','reverselit');
    
    % --- Second blade (mirrored) ---
    vertices2 = vertices;
    vertices2(:,1:2) = -vertices2(:,1:2);
    
    patch('Faces',faces,'Vertices',vertices2, ...
          'FaceVertexCData',colors, ...
          'FaceColor','flat', ...
          'CDataMapping','direct', ...
          'EdgeColor','none', ...
          'FaceLighting','gouraud', ...
          'BackFaceLighting','reverselit');
    
    
    % 4. Precision Hub Generation
    % Step A: Find the EXACT vertical extent of the blade root
    % NOTE: We use the transition root (first col of Grid), not the aero root
    root_z_points = Z_grid(:, 1); 
    z_max = max(root_z_points);
    z_min = min(root_z_points);
    
    % Step B: Define Hub Dimensions based on these limits
    hub_height = (z_max - z_min) * 1.05; 
    hub_center_z = (z_max + z_min) / 2;
    
    hub_radius_outer = (R * x_hub); 
    shaft_radius_inner = (R * 0.03); 
    
    % Step C: Draw Outer Hub (Uses prop_color, now Blue)
    [Hx, Hy, Hz] = cylinder(hub_radius_outer, 100); 
    Hz = (Hz - 0.5) * hub_height + hub_center_z; 
    surf(Hx, Hy, Hz, 'FaceColor', prop_color, 'EdgeColor', 'none', ...
         'SpecularStrength', spec_strength, 'DiffuseStrength', diff_strength, ...
         'FaceLighting', 'gouraud');
    
    % Step D: Draw Inner Hole
    [Sx, Sy, Sz] = cylinder(shaft_radius_inner, 50);
    Sz = (Sz - 0.5) * hub_height * 1.02 + hub_center_z; 
    surf(Sx, Sy, Sz, 'FaceColor', 'k', 'EdgeColor', 'none', ...
         'DiffuseStrength', 0, 'SpecularStrength', 0);

    % --- Step E: Draw Hub Covers (Annulus Caps) ---
    % Create a polar grid for the top/bottom caps
    r_cap = linspace(shaft_radius_inner, hub_radius_outer, 2); 
    th_cap = linspace(0, 2*pi, 80);
    [R_cap, TH_cap] = meshgrid(r_cap, th_cap);
    
    X_cap = R_cap .* cos(TH_cap);
    Y_cap = R_cap .* sin(TH_cap);
    
    % Top Cap (Blue)
    Z_cap_top = ones(size(X_cap)) * (hub_center_z + hub_height/2);
    surf(X_cap, Y_cap, Z_cap_top, 'FaceColor', prop_color, 'EdgeColor', 'none', ...
         'SpecularStrength', spec_strength, 'DiffuseStrength', diff_strength, ...
         'FaceLighting', 'gouraud');
         
    % Bottom Cap (Blue)
    Z_cap_bot = ones(size(X_cap)) * (hub_center_z - hub_height/2);
    surf(X_cap, Y_cap, Z_cap_bot, 'FaceColor', prop_color, 'EdgeColor', 'none', ...
         'SpecularStrength', spec_strength, 'DiffuseStrength', diff_strength, ...
         'FaceLighting', 'gouraud');
    
    % 5. Illumination
    % shading flat; 
    lighting gouraud;
    delete(findall(gcf,'Type','light'));
    
    light('Position', [R, R, 2*R], 'Style', 'local', 'Color', [1 1 1]);
    light('Position', [-R, -R, -R], 'Style', 'local', 'Color', [0.6 0.6 0.6]);
    light('Position', [R, -R, 0], 'Style', 'local', 'Color', [0.7 0.7 0.7]);
    
    % 6. Plot Formatting (for LaTeX)
    xlabel('Tangential coordinate, $x$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Spanwise coordinate, $y$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    zlabel('Thrust axis, $z$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    
    mainTitle = sprintf('\\textbf{Propeller geometry in 3D \\textemdash\\ %s} ($R=%.2f$ cm)', ...
                        prop.name, R*100);
    t = title(mainTitle, 'Interpreter', 'latex', 'FontSize', 15);
    t.Units = 'normalized';
    t.Position(2) = t.Position(2) + 0.08;   % increase vertical offset
    % Don't make the offset smaller than 0.05 or the title will overlap the figure.

    
    box on; 
    grid on;
    axis on; 
    set(gca, 'GridAlpha', 0.3); 
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    view(-45, 20); 
    
    limit_span = R * 1.1;
    limit_axial = max(hub_height, max(abs(Z_grid(:)))) * 1.5;
    xlim([-limit_span/2 limit_span/2]); 
    ylim([-limit_span limit_span]);     
    zlim([-limit_axial limit_axial]); 
    
    camlight('headlight');






    % =======================================================================
    % 2. COMPLETE GEOMETRY ANALYSIS (3 Subplots in 1 Figure)
    % =======================================================================
    figure('Name', '2. Propeller Geometry Analysis', 'Position', [100 50 1200 900], 'Color', 'w');
    
    % --- Prepare Data Common to All Plots ---
    % 1. Extract data
    x_u = prop.airfoil_data.x_u;  % u --> upper (extrados)
    y_u = prop.airfoil_data.y_u;
    x_l = prop.airfoil_data.x_l;
    y_l = prop.airfoil_data.y_l;

    R = outputs.R;
    mesh_vec = outputs.mesh; % Non-dimensional (0 to 1)
    r_vec = mesh_vec * R;       % Dimensional radius [m]
    c = outputs.c;
    theta_deg = (180/pi) * outputs.theta;
    n = numel(mesh_vec);
    
    % ---------- SUBPLOT 1 (Top): Chord and Pitch Distribution ------------
    subplot(2, 2, [1, 2]); % Spans the entire top row
    
    % --- Limits Logic (Preserved) ---
    c_ylim = [0, max(c) * 1.1];
    theta_min = min(theta_deg);
    theta_max = max(theta_deg);
    buffer = 0.05 * (theta_max - theta_min);
    theta_ylim = [theta_min - buffer, theta_max + buffer];
    
    % --- Left Axis: Chord ---
    yyaxis left
    plot(mesh_vec, c, 'b-', 'LineWidth', 2, 'DisplayName', 'Chord $c$');
    ylabel('Chord, $c$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    ylim(c_ylim);
    set(gca, 'YColor', 'b'); % Ensure axis text is blue
    
    % --- Right Axis: Pitch ---
    yyaxis right
    plot(mesh_vec, theta_deg, 'r--', 'LineWidth', 2, 'DisplayName', 'Pitch $\theta$');
    ylabel('Pitch, $\theta$ [deg]', 'Interpreter', 'latex', 'FontSize', 12);
    ylim(theta_ylim);
    set(gca, 'YColor', 'r'); % Ensure axis text is red
    
    % --- Common Axis Settings ---
    xlim([mesh_vec(1), 1]);
    xlabel('Radial non-dimensional coordinate, $x = r/R$', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{A. Radial Distribution of Chord and Pitch}', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    
    % --------- SUBPLOT 2 (Bottom Left): 2D Planform (Upper View) ---------
    subplot(2, 2, 3);
    hold on; grid on; axis equal; box on;
    
    % Calculate Leading/Trailing Edges based on 1/4 chord alignment
    X_LE = -0.25 * c;
    X_TE =  0.75 * c;
    
    % Plot Lines
    plot(r_vec, X_LE, 'b-', 'LineWidth', 2);
    plot(r_vec, X_TE, 'r-', 'LineWidth', 2);
    plot(r_vec, zeros(size(r_vec)), 'k--', 'LineWidth', 1); % Quarter chord line
    
    % Fill the blade shape
    fill([r_vec, fliplr(r_vec)], [X_LE, fliplr(X_TE)], [0.9 0.9 0.9], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
    % Formatting
    xlabel('Radius, $r$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Chordwise Position [m]', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{B. 2D Planform View}', 'Interpreter', 'latex', 'FontSize', 14);
    
    % Custom Legend (Manual construction for clarity)
    % We create dummy entries for the legend to avoid cluttering the plot
    h1 = plot(nan, nan, 'b-', 'LineWidth', 2);
    h2 = plot(nan, nan, 'r-', 'LineWidth', 2);
    h3 = plot(nan, nan, 'k--', 'LineWidth', 1);
    legend([h1 h2 h3], {'Leading Edge', 'Trailing Edge', '1/4 Chord'}, ...
        'Interpreter', 'latex', 'Location', 'SouthWest', 'FontSize', 10);
    
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    
    % ----------- SUBPLOT 3 (Bottom Right): 3D Blade Geometry -------------
    subplot(2, 2, 4);
    
    % --- Data Preparation for 3D Surface ---
    % Make airfoil data local and normalized
    x_profile = [x_u, fliplr(x_l)];
    y_profile = [y_u, fliplr(y_l)];
    x_profile = x_profile - min(x_profile);
    x_profile = x_profile ./ max(x_profile);
    Np = length(x_profile);
    
    % Initialize Matrices
    X_surf = zeros(Np, n);
    Y_surf = zeros(Np, n);
    Z_surf = zeros(Np, n);
    
    % Build 3D coordinates
    for j = 1:n
        c_j = c(j);
        th = deg2rad(theta_deg(j));
        
        % Scale Airfoil
        xs = c_j * (x_profile - 0.25);
        ys = c_j * y_profile;
        
        % Rotate for Pitch (Twist)
        % Note: X is Radial, Y is Chordwise (Tangential), Z is Thickness/Vertical
        % Adjusting rotation to match standard "propeller lying flat" or "vertical" view
        
        % Standard rotation logic:
        xr = xs*cos(-th) - ys*sin(-th);
        yr = xs*sin(-th) + ys*cos(-th);
        
        X_surf(:, j) = r_vec(j); % Radial Axis
        Y_surf(:, j) = xr;       % Chord/Tangential Axis
        Z_surf(:, j) = yr;       % Thickness/Twist Axis
    end
    
    % Plot Surface
    surf(X_surf, Y_surf, Z_surf, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 1]);
    
    % Formatting
    axis equal; grid on;
    xlabel('Radius, $r$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Chordwise, $y$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    zlabel('Thickness, $z$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{C. 3D Blade Geometry}', 'Interpreter', 'latex', 'FontSize', 14);
    
    % Lighting for 3D effect
    camlight; lighting gouraud;
    view(3); % Standard isometric view
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

    % Figure Super Title
    mainTitle = sprintf('\\textbf{Propeller Geometry Analysis \\textemdash\\ %s} ($R=%.2f$ cm)', ...
                        prop.name, R*100);
    sgtitle(mainTitle, 'Interpreter', 'latex', 'FontSize', 15);





    % =====================================================================
    % 3. Combined figure of Airfoil Geometry and Aerodynamics
    % =====================================================================
    figure('Name', '3. Airfoil Geometry and Aerodynamics', 'Position', [100 50 1200 800], 'Color', 'w');
    
    % -----------------------------------------------------------------------
    % SUBPLOT 1 (Left Half): Blade-element inflow geometry
    % -----------------------------------------------------------------------
    % ax1 = axes('Position', [0.05, 0.15, 0.47, 0.70]);
    ax1 = axes('Position', [0.05, 0.10, 0.40, 0.75]); 
    % hold on; axis equal; grid on;
    hold on; axis equal tight; grid on;
    
    % -- 1. Select Radial Station (75% Span) --
    [~, idx] = min(abs(outputs.r - 0.75 * outputs.R));
    c_local = outputs.c(idx);
    theta_local = outputs.theta(idx); 
    phi_local = outputs.phi(idx);
    alpha_local = outputs.alpha(idx);
    
    % -- 2. Process Airfoil Geometry --
    x_foil = [prop.airfoil_data.x_u(:); flipud(prop.airfoil_data.x_l(:))];
    y_foil = [prop.airfoil_data.y_u(:); flipud(prop.airfoil_data.y_l(:))];
    
    % Rotate Airfoil by -theta_local 
    theta_rot = -theta_local; 
    R_mat = [cos(theta_rot), -sin(theta_rot); sin(theta_rot), cos(theta_rot)];
    coords_rot = R_mat * [x_foil'*c_local; y_foil'*c_local]; 
    x_rot = coords_rot(1, :);
    y_rot = coords_rot(2, :);
    
    % -- 3. Calculate Zero Lift Line (ZLL) --
    Fc = griddedInterpolant({prop.airfoil_data.alpha_deg*(pi/180), prop.airfoil_data.Re}, ...
                            prop.airfoil_data.cl_matrix, 'linear', 'linear');
    Re_loc = outputs.Re_local(idx);
    find_alpha0 = @(a) Fc(a, Re_loc);
    alpha_0 = fzero(find_alpha0, 0); 
    theta_zll = theta_local - alpha_0; % Note that alpha_0 < 0, thats why we subtract it (to actually add it)
                                       % theta_zll is a clockwise angle.
                                       % Later, we will rotate by -theta_zll
    
    % -- 4. Draw Reference Lines --
    
    % A. Rotor Plane
    line([-c_local*1.5, c_local*2.0], [0, 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.8);
    text(c_local*1.9, 0, 'Rotor Plane', 'FontSize', 9, 'VerticalAlignment', 'bottom', ...
         'HorizontalAlignment', 'right', 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
    
    % B. Chord Line (Prolonged BOTH ways)
    % Forward Extension
    len_chord_fwd = c_local * 2.2;
    line([0, len_chord_fwd*cos(theta_rot)], [0, len_chord_fwd*sin(theta_rot)], ...
        'Color', [0.4 0.4 0.4], 'LineStyle', '-.', 'LineWidth', 1.0);
    % Backward Extension (for Alpha visualization on the left)
    len_chord_back = c_local * 1.5;
    line([0, -0.9*len_chord_back*cos(theta_rot)], [0, -0.9*len_chord_back*sin(theta_rot)], ...
        'Color', [0.4 0.4 0.4], 'LineStyle', '-.', 'LineWidth', 1.0);
    
    % -- C. Zero Lift Line (ZLL) --
    % We define two lengths: one extending forward (upstream/left) and one backward
    len_zll_fwd = c_local * 1.5; % Extension to the front
    len_zll_aft = c_local * 1.8; % Extension to the back
    
    % Calculate coordinates for the endpoints
    % Note: We use -len_zll_fwd for the starting point to extend to the left
    x_zll = [-len_zll_fwd, len_zll_aft] * cos(-theta_zll);
    y_zll = [-len_zll_fwd, len_zll_aft] * sin(-theta_zll);
    
    plot(x_zll, y_zll, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1.2);
    
    % Label: Placed on the forward extension (upstream)
    % 'VerticalAlignment', 'bottom' puts it right above the line
    % 'HorizontalAlignment', 'right' anchors it nicely at the tip
    text(-len_zll_fwd * cos(-theta_zll), -len_zll_fwd * sin(-theta_zll), ...
        'ZLL ', ... % Added a space for padding
        'FontSize', 8, 'Color', [0.5 0.5 0.5], 'Interpreter', 'latex', ...
        'VerticalAlignment', 'bottom', ... 
        'HorizontalAlignment', 'right'); % Aligns text to end at the tip
    
    % D. Resultant Velocity Direction Line (W)
    phi_draw = -phi_local; 
    len_wind_fwd = c_local * 1.8;
    len_wind_back = c_local * 1.5;
    
    % Draw line extending both ways
    line([-0.9*len_wind_back*cos(phi_draw), 0.95*len_wind_fwd*cos(phi_draw)], ...
         [-0.9*len_wind_back*sin(phi_draw), 0.95*len_wind_fwd*sin(phi_draw)], ...
         'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.2);
    
    % Label: "Resultant velocity direction" below the line
    x_lbl = 0.6*len_wind_fwd * cos(phi_draw);
    y_lbl = 1.05*len_wind_fwd * sin(phi_draw);
    text(x_lbl, y_lbl, ' Resultant velocity direction', 'Color', [0, 0.4470, 0.7410], ...
         'FontSize', 10, 'Interpreter', 'latex', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
    % -- 5. Draw Airfoil --
    fill(x_rot, y_rot, [0.2 0.2 0.2], 'FaceAlpha', 0.1, 'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 1.5);
    
    % -- 6. Velocity Triangle --
    Vx = outputs.V_x(idx);
    Vy = outputs.V_y(idx); 
    W_mag = sqrt(Vx^2 + Vy^2);
    scale_fac = (c_local*0.65) / W_mag;         
    
    vec_origin = [-c_local*0.5, -c_local*0.6]; 
    
    % Vectors
    % Omega*r
    quiver(vec_origin(1), vec_origin(2), Vy*scale_fac, 0, ...
        'Color', 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.3, 'AutoScale', 'off'); 
    % V + vi (Down)
    quiver(vec_origin(1)+Vy*scale_fac, vec_origin(2), 0, -Vx*scale_fac, ...
        'Color', [0.85, 0.32, 0.09], 'LineWidth', 1.5, 'MaxHeadSize', 0.3, 'AutoScale', 'off'); 
    % u_R (Resultant)
    quiver(vec_origin(1), vec_origin(2), Vy*scale_fac, -Vx*scale_fac, ...
        'Color', [0, 0.4470, 0.7410], 'LineWidth', 2, 'MaxHeadSize', 0.3, 'AutoScale', 'off'); 
    
    % Labels
    text(vec_origin(1) + Vy*scale_fac/2, vec_origin(2), '$\Omega r$', ...
        'VerticalAlignment', 'bottom', 'Interpreter', 'latex', 'FontSize', 10);
    
    % Change label to (v + vi)
    text(vec_origin(1) + Vy*scale_fac, vec_origin(2) - Vx*scale_fac/2, ' $(v+v_i)$', ...
        'Color', [0.85, 0.32, 0.09], 'HorizontalAlignment', 'left', 'Interpreter', 'latex', 'FontSize', 10);
    
    % Label u_R for resultant
    text(vec_origin(1) + Vy*scale_fac*0.4, vec_origin(2) - Vx*scale_fac*0.6, '$u_R$', ...
        'Color', [0, 0.4470, 0.7410], 'Interpreter', 'latex', 'FontSize', 11, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right');
    
    % Phi Angle in Triangle
    r_tri_arc = Vy*scale_fac * 0.45;
    t_tri_phi = linspace(0, phi_draw, 20);
    plot(vec_origin(1) + r_tri_arc*cos(t_tri_phi), vec_origin(2) + r_tri_arc*sin(t_tri_phi), ...
         'Color', [0, 0.4470, 0.7410], 'LineWidth', 0.8);
    text(vec_origin(1) + r_tri_arc*1.15, vec_origin(2) + r_tri_arc*0.5*sin(phi_draw), '$\phi$', ...
         'Color', [0, 0.4470, 0.7410], 'Interpreter', 'latex', 'FontSize', 9, 'HorizontalAlignment', 'left');
    
    
    % -- 7. Detailed Geometry Angles --
    r_theta = c_local * 0.70;
    r_phi   = c_local * 1.10;
    
    % A. Theta (Pitch) - Right side
    t_theta = linspace(0, theta_rot, 30);
    plot(r_theta * cos(t_theta), r_theta * sin(t_theta), 'k-', 'LineWidth', 0.8);
    text(r_theta * cos(theta_rot/2), r_theta * sin(theta_rot/2), ' $\theta$', ...
        'Color', 'k', 'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'FontSize', 10);
    
    % B. Phi (Inflow) - Right side
    t_phi = linspace(0, phi_draw, 30);
    plot(r_phi * cos(t_phi), r_phi * sin(t_phi), '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 0.8);
    text(r_phi * cos(phi_draw/2), r_phi * sin(phi_draw/2), ' $\phi$', ...
        'Color', [0, 0.4470, 0.7410], 'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'FontSize', 10);
    
    % C. Alpha (AoA) - LEFT SIDE (Negative X)
    % Calculate angles on the left (rotated 180 deg / pi radians)
    angle_chord_left = theta_rot + pi;
    angle_wind_left = phi_draw + pi;
    r_alpha_left = c_local * 1.2;
    
    t_alpha = linspace(angle_wind_left, angle_chord_left, 40);
    plot(r_alpha_left * cos(t_alpha), r_alpha_left * sin(t_alpha), 'r-', 'LineWidth', 1.2);
    
    % Label Alpha on the left
    mid_alpha_left = (angle_wind_left + angle_chord_left) / 2;
    text(1.2*r_alpha_left * cos(mid_alpha_left), 1.2*r_alpha_left*sin(mid_alpha_left), '$\alpha$', ...
        'Color', 'r', 'Interpreter', 'latex', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold');
    
    % -- 8. Formatting Subplot 1 --
    title(sprintf('\\textbf{A. Inflow Geometry at characteristic section} ($r/R = %.2f$)', outputs.r(idx)/outputs.R), 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('Rotor Plane $x$ [m]', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('Axial Direction $y$ [m]', 'Interpreter', 'latex', 'FontSize', 11);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
    
    % Adjust limits to include left side visualization
    all_x = [x_rot, vec_origin(1), x_lbl, -len_chord_back, len_chord_fwd];
    all_y = [y_rot, vec_origin(2)-Vx*scale_fac, y_lbl];
    xlim([min(all_x)-0.02, max(all_x)+0.02]);
    ylim([min(all_y)-0.02, max(all_y)+0.02]);
    
    
    % -----------------------------------------------------------------------
    % SUBPLOT 2 & 3 (Unchanged)
    % -----------------------------------------------------------------------
    % ax2 = axes('Position', [0.53, 0.15, 0.20, 0.65]); 
    ax2 = axes('Position', [0.53, 0.15, 0.20, 0.6]); 
    hold on; grid on;
    airfoil_data = prop.airfoil_data;
    colors = jet(length(airfoil_data.Re));
    legend_list = {};
    for i = 1:length(airfoil_data.Re)
        plot(airfoil_data.alpha_deg, airfoil_data.cl_matrix(:,i), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
        legend_list{end+1} = sprintf('$Re = %.1e$', airfoil_data.Re(i));
    end
    xlabel('$\alpha$ [deg]', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$c_l$', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{B. Lift Curves}', 'Interpreter', 'latex', 'FontSize', 12);
    legend(legend_list, 'Interpreter', 'latex', 'Location', 'SouthEast', 'FontSize', 7);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
    
    % ax3 = axes('Position', [0.78, 0.15, 0.20, 0.65]);
    ax3 = axes('Position', [0.78, 0.15, 0.20, 0.6]);
    hold on; grid on;
    for i = 1:length(airfoil_data.Re)
        plot(airfoil_data.alpha_deg, airfoil_data.cd_matrix(:,i), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    xlabel('$\alpha$ [deg]', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$c_d$', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{C. Drag Polars}', 'Interpreter', 'latex', 'FontSize', 12);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
    
    % =======================================================================
    % GLOBAL TITLE
    % =======================================================================
    if ~exist('titles', 'var'), titles = '3. Airfoil Geometry and Aerodynamics'; end
    hTitle = annotation('textbox', [0.5, 0.88, 0, 0], 'String', titles, ...
        'Interpreter', 'latex', 'FontSize', 13, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'EdgeColor', 'black', 'LineWidth', 1, ...
        'BackgroundColor', 'white', 'FitBoxToText', 'on');
    drawnow; 
    currentPos = get(hTitle, 'Position'); 
    set(hTitle, 'Position', [0.5 - currentPos(3)/2, currentPos(2), currentPos(3), currentPos(4)]);






    % =====================================================================
    % 4. Combined Figure of Aerodynamic Angles and Forces
    % =====================================================================
    
    % 1. Create the main wide figure
    figure('Name', '4. Propeller Aerodynamics and Forces', 'Position', [100 100 1600 550], 'Color', 'w');
    
    % =====================================================================
    % SUBPLOT 1: Aerodynamics (Lambda and Alpha)
    % =====================================================================
    subplot(1, 3, 1);
    
    % --- ADJUST POSITION HERE ---
    pos = get(gca, 'Position');
    pos(4) = 0.60;   % Reduce height
    pos(2) = 0.15;   % Bottom margin
    set(gca, 'Position', pos);
    % ----------------------------
    
    hold on;
    
    lambda = outputs.lambda;
    alpha  = outputs.alpha;
    
    Re_local = outputs.Re_local;
    Fc = prop.airfoil.Fc;
    
    % 1. Extract the original alpha grid points
    grid_alphas = Fc.GridVectors{1};
    num_alphas  = length(grid_alphas);
    
    % 2. Evaluate Fc for every alpha against every target_Re
    % cl_surface is (num_alphas x numel(Re_local))
    cl_surface = Fc({grid_alphas, Re_local});
    
    % 3. Discrete indices for max and min Cl for each Re
    [clmax_disc, idx_clmax] = max(cl_surface, [], 1);
    [clmin_disc, idx_clmin] = min(cl_surface, [], 1);
    
    % 4. Refine using Parabolic Fit (Sub-grid interpolation) for BOTH extrema
    alpha_stall_clmax = zeros(size(Re_local));
    alpha_stall_clmin = zeros(size(Re_local));
    
    % Optional: refined cl values at refined alpha (from parabola vertex)
    clmax_refined = clmax_disc;  % initialize with discrete as fallback
    clmin_refined = clmin_disc;
    
    for k = 1:numel(Re_local)
    
        % ---- CLMAX refinement ----
        idx = idx_clmax(k);
        if idx > 1 && idx < num_alphas
            x_sub = grid_alphas(idx-1:idx+1);
            y_sub = cl_surface(idx-1:idx+1, k);
    
            p = polyfit(x_sub, y_sub, 2);  % y = a x^2 + b x + c
            a = p(1); b = p(2); c = p(3);
    
            if abs(a) > eps  % avoid divide-by-zero / near-linear
                alpha_refined = -b/(2*a);
    
                % sanity: keep within one grid step of the discrete peak
                if abs(alpha_refined - grid_alphas(idx)) > (grid_alphas(idx) - grid_alphas(idx-1))
                    alpha_stall_clmax(k) = grid_alphas(idx);
                else
                    alpha_stall_clmax(k) = alpha_refined;
                    clmax_refined(k) = polyval(p, alpha_refined);
                end
            else
                alpha_stall_clmax(k) = grid_alphas(idx);
            end
        else
            alpha_stall_clmax(k) = grid_alphas(idx);
        end
    
        % ---- CLMIN refinement ----
        idx = idx_clmin(k);
        if idx > 1 && idx < num_alphas
            x_sub = grid_alphas(idx-1:idx+1);
            y_sub = cl_surface(idx-1:idx+1, k);
    
            p = polyfit(x_sub, y_sub, 2);
            a = p(1); b = p(2); c = p(3);
    
            if abs(a) > eps
                alpha_refined = -b/(2*a);
    
                % sanity: keep within one grid step of the discrete minimum
                if abs(alpha_refined - grid_alphas(idx)) > (grid_alphas(idx) - grid_alphas(idx-1))
                    alpha_stall_clmin(k) = grid_alphas(idx);
                else
                    alpha_stall_clmin(k) = alpha_refined;
                    clmin_refined(k) = polyval(p, alpha_refined);
                end
            else
                alpha_stall_clmin(k) = grid_alphas(idx);
            end
        else
            alpha_stall_clmin(k) = grid_alphas(idx);
        end
    
    end
    
    % 5. Store results
    outputs.alpha_stall_clmax_local = alpha_stall_clmax;
    outputs.alpha_stall_clmin_local = alpha_stall_clmin;
    
    % Optional
    outputs.clmax_local = clmax_refined;
    outputs.clmin_local = clmin_refined;

    % Y-limits (include alpha constraints)
    global_min = min([lambda, alpha, alpha_stall_clmin, alpha_stall_clmax]);
    global_max = max([lambda, alpha, alpha_stall_clmin, alpha_stall_clmax]);
    % % If I want to use the alpha data limits:
    % global_min = min([lambda, alpha, prop.airfoil.alpha_min, prop.airfoil.alpha_max]);
    % global_max = max([lambda, alpha, prop.airfoil.alpha_min, prop.airfoil.alpha_max]);
    padding = 0.1 * (global_max - global_min);
    if padding == 0
        padding = 0.1;
    end
    Aero_ylim = [global_min - padding, global_max + padding];
    
    % Left axis in radians
    yyaxis left
    f1 = plot(outputs.mesh, lambda, 'b', ...
        'LineWidth', 2, 'DisplayName', '$\lambda$');
    f2 = plot(outputs.mesh, alpha, 'r-', ...
        'LineWidth', 2, 'DisplayName', '$\alpha$');
    f3 = plot(outputs.mesh, outputs.alpha_stall_clmax_local, 'r:', ...
        'LineWidth', 2, 'DisplayName', '$\alpha_{stall,max}$');
    f4 = plot(outputs.mesh, outputs.alpha_stall_clmin_local, 'r:', ...
        'LineWidth', 2, 'DisplayName', '$\alpha_{stall,min}$');
    % % If I want to use the alpha data limits:
    % yline(prop.airfoil.alpha_max, 'r:', ...
    %     'LineWidth', 1.5, 'DisplayName', '$\alpha_{stall,max}$');
    % yline(prop.airfoil.alpha_min, 'r:', ...
    %     'LineWidth', 1.5, 'DisplayName', '$\alpha_{stall,min}$');
    ylim(Aero_ylim);
    ylabel('Value [rad]', 'Interpreter', 'latex', 'FontSize', 12);
    
    % Right axis in degrees (and red color)
    yyaxis right
    ylim(rad2deg(Aero_ylim));
    ylabel('Value [deg]', 'Interpreter', 'latex', 'FontSize', 12);
    ax = gca;
    ax.YAxis(2).Color = 'r';
    
    % Labels and legend
    xlabel('Radial coordinate, $x$', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{A. }$\mathbf{\lambda}$\textbf{ and }$\mathbf{\alpha}$\textbf{ vs }$\mathbf{x}$', ...
        'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend([f1, f2, f3], { ...
        '$\lambda$', ...
        '$\alpha$', ...
        sprintf([' $\\alpha$ stall limits considering\n', ...
                 ' $Re$ number effect\n']) ...
        }, ...
        'Location', 'best', ...
        'Interpreter', 'latex', ...
        'FontSize', 7);

    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    hold off;

    
    % =======================================================================
    % SUBPLOT 2: Thrust Distribution
    % =======================================================================
    subplot(1, 3, 2);
    
    % --- ADJUST POSITION HERE ---
    pos = get(gca, 'Position');
    pos(4) = 0.60; 
    pos(2) = 0.15; 
    set(gca, 'Position', pos);
    % -----------------------------
    
    dTdx = outputs.dTdx;
    dTdx_ylim = [min(dTdx) * 1.1, max(dTdx) * 1.1];
    
    plot(outputs.mesh, outputs.dTdx, 'k', 'LineWidth', 2);
    xlabel('Radial coordinate, $x$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\frac{dT}{dx}$ [N/m]', 'Interpreter', 'latex', 'FontSize', 12);
    ylim(dTdx_ylim);
    grid on;
    title('\textbf{B. Thrust Distribution} ($dT/dx$)', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    % =======================================================================
    % SUBPLOT 3: Tangential Force Distribution
    % =======================================================================
    subplot(1, 3, 3);
    
    % --- ADJUST POSITION HERE ---
    pos = get(gca, 'Position');
    pos(4) = 0.60; 
    pos(2) = 0.15; 
    set(gca, 'Position', pos);
    % -----------------------------
    
    dFtdx = outputs.dFtdx;
    dFtdx_ylim = [min(dFtdx) * 1.1, max(dFtdx) * 1.1];
    
    plot(outputs.mesh, outputs.dFtdx, 'k', 'LineWidth', 2);
    xlabel('Radial coordinate, $x$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\frac{dF_t}{dx}$ [N/m]', 'Interpreter', 'latex', 'FontSize', 12);
    ylim(dFtdx_ylim);
    grid on;
    title('\textbf{C. Tangential Force} ($dF_t/dx$)', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    % =======================================================================
    % GLOBAL TITLE (Robust Centering)
    % ======================================================================= 
    % 1. Create the annotation with 'FitBoxToText' ON
    % We start at a dummy position [0.5, 0.85] with 0 width.
    % MATLAB will automatically expand the width/height to fit the text.
    hTitle = annotation('textbox', [0.5, 0.88, 0, 0], ...
        'String', titles, ...
        'Interpreter', 'latex', ...
        'FontSize', 13, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'EdgeColor', 'black', ...
        'LineWidth', 1, ...
        'BackgroundColor', 'white', ...
        'FitBoxToText', 'on');
    
    % 2. Force MATLAB to render the text so we can get its true size
    drawnow; 
    
    % 3. Recalculate Position to Center it
    % Get the fitted size [x, y, w, h]
    currentPos = get(hTitle, 'Position'); 
    newWidth = currentPos(3);
    
    % Calculate new X to center it: (0.5 - half_width)
    newX = 0.5 - (newWidth / 2);
    
    % Update the position (Keep the original Y, use calculated X)
    set(hTitle, 'Position', [newX, currentPos(2), newWidth, currentPos(4)]);





    % =======================================================================
    % 5. Combined Figure of Stress Distributions
    % =======================================================================
    
    % 1. Create the main wide figure
    figure('Name', '5. Stress Distributions', 'Position', [100 100 1600 550], 'Color', 'w');

    % --- Shared Definitions for the Reinforcement Line ---
    x_reinf_end = prop.x_finish_reinforcement;
    % Use a Cell Array for multiline LaTeX support
    reinfLabelStr = {'\textbf{Structural reinforcement must be used in this region}', 'Focus only on stresses to the right'};
    % Style: Red dashed line
    reinfLineStyle = { 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5 };
    
    % =======================================================================
    % SUBPLOT 1: Centrifugal Stress
    % =======================================================================
    subplot(1, 3, 1);
    
    % --- ADJUST POSITION (Consistent layout) ---
    pos = get(gca, 'Position');
    pos(4) = 0.60;   % Reduce height
    pos(2) = 0.15;   % Bottom margin
    set(gca, 'Position', pos);
    % -------------------------------------------
    
    % Data Processing
    sigma_cf_MPa = outputs.sigma_cf / 1e6;
    sigma_cf_max = max(sigma_cf_MPa);
    sigma_cf_ylim = [0, sigma_cf_max * 1.1]; % Start from 0 for stress
    
    plot(outputs.mesh, sigma_cf_MPa, 'b-', 'LineWidth', 2);
    hold on
    xline(x_reinf_end, ...
          'Label', reinfLabelStr, ...
          reinfLineStyle{:}, ...
          'Interpreter', 'latex', ...
          'LabelVerticalAlignment', 'bottom', ...
          'LabelHorizontalAlignment', 'left', ...
          'FontSize', 9, ...
          'HandleVisibility', 'off');
    hold off
    xlabel('Radial coordinate, $x$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\sigma_{\mathrm{cf}}$ [MPa]', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{A. Centrifugal Stress}', 'Interpreter', 'latex', 'FontSize', 14);
    ylim(sigma_cf_ylim);
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    
    % =======================================================================
    % SUBPLOT 2: Maximum Bending Stresses
    % =======================================================================
    subplot(1, 3, 2);
    
    % --- ADJUST POSITION ---
    pos = get(gca, 'Position');
    pos(4) = 0.60; 
    pos(2) = 0.15; 
    set(gca, 'Position', pos);
    % -----------------------
    
    % Data Processing
    sigma_bend_Ix = outputs.sigma_bend_max_Ix / 1e6;
    sigma_bend_Iy = outputs.sigma_bend_max_Iy / 1e6;
    
    % Robust limit calculation (using (:) to handle vectors safely)
    all_bend = [sigma_bend_Ix(:); sigma_bend_Iy(:)];
    bend_max = max(all_bend);
    bend_min = min(all_bend);
    padding = (bend_max - bend_min) * 0.1;
    if padding == 0, padding = 0.1; end
    sigma_bend_ylim = [bend_min - padding, bend_max + padding];
    hold on
    plot(outputs.mesh, sigma_bend_Ix, 'b', 'LineWidth', 2, 'DisplayName', '$\sigma_{bend, I_x}$');
    plot(outputs.mesh, sigma_bend_Iy, 'r--', 'LineWidth', 2, 'DisplayName', '$\sigma_{bend, I_y}$');
    xline(x_reinf_end, ...
          'Label', reinfLabelStr, ...
          reinfLineStyle{:}, ...
          'Interpreter', 'latex', ...
          'LabelVerticalAlignment', 'bottom', ...
          'LabelHorizontalAlignment', 'left', ...
          'FontSize', 9, ...
          'HandleVisibility', 'off');
    hold off
    xlabel('Radial coordinate, $x$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\sigma_{\mathrm{bend, max}}$ [MPa]', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{B. Bending Stress at critical points}', 'Interpreter', 'latex', 'FontSize', 14);
    ylim(sigma_bend_ylim);
    grid on;
    legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 11);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    
    % =======================================================================
    % SUBPLOT 3: Maximum Total Stress
    % =======================================================================
    subplot(1, 3, 3);
    
    % --- ADJUST POSITION ---
    pos = get(gca, 'Position');
    pos(4) = 0.60; 
    pos(2) = 0.15; 
    set(gca, 'Position', pos);
    % -----------------------
    
    % Data Processing
    sigma_total_MPa = outputs.sigma_total_max / 1e6;
    sigma_total_max = max(sigma_total_MPa);
    sigma_total_ylim = [0, sigma_total_max * 1.1];
    
    plot(outputs.mesh, sigma_total_MPa, 'k', 'LineWidth', 2); % Black for total
    hold on
    xline(x_reinf_end, ...
          'Label', reinfLabelStr, ...
          reinfLineStyle{:}, ...
          'Interpreter', 'latex', ...
          'LabelVerticalAlignment', 'bottom', ...
          'LabelHorizontalAlignment', 'left', ...
          'FontSize', 9, ...
          'HandleVisibility', 'off');
    hold off
    xlabel('Radial coordinate, $x$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\sigma_{\mathrm{total, max}}$ [MPa]', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{C. Total Stress at critical points}', 'Interpreter', 'latex', 'FontSize', 14);
    ylim(sigma_total_ylim);
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    % =======================================================================
    % GLOBAL TITLE (Robust Centering)
    % =======================================================================
    % 1. Create the annotation with 'FitBoxToText' ON
    hTitle = annotation('textbox', [0.5, 0.88, 0, 0], ...
        'String', titles, ...
        'Interpreter', 'latex', ...
        'FontSize', 13, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'EdgeColor', 'black', ...
        'LineWidth', 1, ...
        'BackgroundColor', 'white', ...
        'FitBoxToText', 'on');
    
    % 2. Force MATLAB to render the text
    drawnow; 
    
    % 3. Recalculate Position to Center it
    currentPos = get(hTitle, 'Position'); 
    newWidth = currentPos(3);
    newX = 0.5 - (newWidth / 2);
    set(hTitle, 'Position', [newX, currentPos(2), newWidth, currentPos(4)]);

    % =====================================================================
    % NORTHWEST NOTE (Top-Left Corner)
    % =====================================================================
    note_str = '* \textbf{Note:} The high stresses observed near the root are theoretical. In real propellers, these stresses are significantly lower due to the structural reinforcement of the root region (approx. up to $20\%$ of the radius).';
    
    % Position: [x y w h] -> Top Left (Northwest)
    % y=0.82 places it below the title (0.88) but above the subplots (0.75)
    annotation('textbox', [0.03, 0.82, 0.3, 0.1], ...
        'String', note_str, ...
        'Interpreter', 'latex', ...
        'FontSize', 10, ...
        'Color', [0.3 0.3 0.3], ... 
        'HorizontalAlignment', 'left', ... % Left aligned text
        'VerticalAlignment', 'top', ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none');





    %% ====================================================================
    %  6. PLOT THE AIRFOIL STRESS DISTRIBUTION (at 3 sections side-by-side)
    %  ====================================================================

    mesh = outputs.mesh;
    dx = outputs.dx;

    % --- Indices for the 3 sections ---
    % A. Identify Critical Section (CONSIDERING THAT IT HAS TO BE AT X > X_FINISH_REINFORCEMENT)
    % Filter: Only consider sections outside the reinforcement zone
    valid_region_mask = outputs.mesh >= prop.x_finish_reinforcement;
    % Create a temporary array for finding the max
    search_stress = outputs.sigma_total_max;
    % Set stress in the excluded region to -Infinity so they are never picked
    % (unless the entire blade is excluded, which would be an error in inputs)
    search_stress(~valid_region_mask) = -inf;
    % Find the max stress and its index within the valid region
    [max_stress_val, idx_crit] = max(search_stress);
    % If the max stress section is x_finish_reinforcement within a dx 
    % margin (which it probably is), let's make the indication of the most
    % critical section match exactly x_finish_reinforcement:
    % First, we find the corresponding dx (it can be a vector):
    if isscalar(outputs.dx)
        dx_eff = outputs.dx;
    else
        dx_eff = outputs.dx(idx_crit);
    end
    if mesh(idx_crit) > prop.x_finish_reinforcement - dx_eff && ...
       mesh(idx_crit) < prop.x_finish_reinforcement + dx_eff
        x_critical = prop.x_finish_reinforcement;
    else
        x_critical = mesh(idx_crit);
    end

    % B. Mid-span section (x=0.40)
    [~, idx_mid]  = min(abs(mesh - 0.40));   % x = 0.4
    % B. Characteristic section (x=0.75)
    [~, idx_075]  = min(abs(mesh - 0.75));   % x = 0.75

    idx_list   = [idx_crit, idx_mid, idx_075];
    
    title_list = { ...
        { ...
          sprintf('\\textbf{A. Most critical section (r/R=%.3f)}', x_critical), ...
          sprintf('(Only aerodynamic zone is considered: $x \\ge %.3f$)', prop.x_finish_reinforcement) ...
        }, ...
        sprintf('\\textbf{B. Mid-span section} (r/R=%.2f)', mesh(idx_mid)), ...
        sprintf('\\textbf{C. Characteristic section} (r/R=%.2f)', mesh(idx_075)) ...
    };

    % --- Common color scale across the three subplots (MPa) ---
    sig_all = [];
    for k = 1:numel(idx_list)
        sig_all = [sig_all; outputs.sigma_total_pts(:, idx_list(k)) / 1e6]; %#ok<AGROW>
    end
    cax = [min(sig_all), max(sig_all)];
    % Optional: make symmetric about zero for cleaner comparison
    % m = max(abs(cax)); cax = [-m, m];
    
    % --- Figure & layout ---
    figure('Name','6. Airfoil stress distribution (3 sections)','Color','w','Renderer','opengl');
    tlo = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    colormap(turbo);
    
    ax_arr = gobjects(1,3);
    for k = 1:3
        ax_arr(k) = nexttile(tlo,k);
        plotAirfoilStressSection(ax_arr(k), idx_list(k), title_list{k}, cax, outputs, prop);
    end
    
    % Shared colorbar (same scale for all axes)
    cb = colorbar(ax_arr(end));        % attach to the 3rd axes
    cb.Location = 'eastoutside';
    cb.TickLabelInterpreter = 'latex';
    cb.Label.String = 'Normal stress, $\sigma$ [MPa]';
    cb.Label.Interpreter = 'latex';
    cb.FontSize = 11;

    % =======================================================================
    % GLOBAL TITLE (Robust Centering)
    % =======================================================================
    % 1. Create the annotation with 'FitBoxToText' ON
    hTitle = annotation('textbox', [0.5, 0.88, 0, 0], ...
        'String', titles, ...
        'Interpreter', 'latex', ...
        'FontSize', 13, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'EdgeColor', 'black', ...
        'LineWidth', 1, ...
        'BackgroundColor', 'white', ...
        'FitBoxToText', 'on');
    
    % 2. Force MATLAB to render the text
    drawnow; 
    
    % 3. Recalculate Position to Center it
    currentPos = get(hTitle, 'Position'); 
    newWidth = currentPos(3);
    newX = 0.5 - (newWidth / 2);
    set(hTitle, 'Position', [newX, currentPos(2), newWidth, currentPos(4)]);
    






    % =====================================================================
    % 7. Performance plots of the optimal propeller: eta_p(J), ct(J), cp(J)
    % =====================================================================

    % When the performance plots of a propeller are obtained experimentally
    % in a wind tunel, the Measurement Procedure (The "Sweep") is:
    % Since the coefficients depend on the Advance Ratio J = V/nD, you must
    % vary J from 0 (static) to a value where thrust typically drops to 
    % zero. This is usually done by holding the RPM constant and changing 
    % the wind speed.
    % 1) Set Fixed RPM: The propeller is spun up to a constant rotational 
    %    speed (n). This is often chosen to match a specific Reynolds 
    %    number or the operational RPM of the real aircraft. 
    % 2) Static Test (J=0): The wind tunnel fan is off (V=0). Thrust and 
    %    Torque are measured. This gives the "static thrust" point on the 
    %    far left of the graph.
    % 3) Velocity Sweep: The wind tunnel speed (V) is gradually increased 
    %    in steps while the propeller controller maintains the fixed RPM.
    % 4) Data Acquisition: At each velocity step, the computer records:
    %    Freestream Velocity (V), Rotational Speed (n), Thrust (T), 
    %    Torque (Q), Air Density (rho) (calculated from P and Temp)
    % 5) Data Processing (Calculating the Coefficients): Once the raw data 
    %    (T, Q, V, n) is collected for the entire sweep, the 
    %    non-dimensional coefficients are calculated for each data point 
    %    using the standard formulas
    %
    % In this code, since the "experiments" are being run by a computer, we
    % don't care about reducing the number of experiments. We don't need to
    % fix RPM and only vary v.
    % We will vary both RPM and v, and that will provide a cloud of points
    % that really represent the behaviour of the propeller at different
    % Reynolds numbers (not just a pre-defined mission profile of the
    % propeller with fixed RPM and varying v).

    fprintf('-------------------------------------------------------------------------\n');
    fprintf('Running BEMT for many v, rpm to generate PERFORMANCES PLOTS...\n'); 
    fprintf('(The code will try some unreal v, rpm. It is normal to see SABEMMT warning).\n');
    fprintf('-------------------------------------------------------------------------\n');

    % pause(5); % Wait 5 seconds to let the user see that Warning is normal

    R = outputs.R;
    c = outputs.c;
    theta_deg = (180/pi)*outputs.theta;

    D = 2 * R;

    % --- 1. Define the velocity and rpm sweep ---
    v_vec = linspace(0, v, 30);
    rpm_vec = linspace(0.25*rpm, rpm, 30);

    n_rot_vec = rpm_vec./60;  % rev/s

    J_min_tested = 0;
    J_max_tested = max(v_vec) / ( min(n_rot_vec)*D );

    num_points = length(v_vec);
    num_rpms = length(rpm_vec);

    % Pre-allocate result vectors
    J_vec     = zeros(1, num_points*num_rpms);
    eta_p_vec = zeros(1, num_points*num_rpms);
    Ct_vec    = zeros(1, num_points*num_rpms);
    Cp_vec    = zeros(1, num_points*num_rpms);

    fprintf('Running Sweep (J = %.2f to J = %.2f) with %d points...\n', ...
            J_min_tested, J_max_tested, num_points*num_rpms);


    % --- 2. Execution Loop ---
    i_end = 0;
    for i = 1:num_points
        v_current = v_vec(i);

        for j = 1:num_rpms
            rpm_current = rpm_vec(j);
            n_rot_current = n_rot_vec(j);

            J_vec(i_end*num_rpms + j) = v_current / (n_rot_current*D);
            try
                % Execute BEMT
                res = runSABEMMT(v_current, rpm_current, Nb, R, c, theta_deg, prop, env);

                % Store results
                eta_p_vec(i_end*num_rpms + j) = res.eta_p;
                Ct_vec(i_end*num_rpms + j)    = res.Ct_std; 
                Cp_vec(i_end*num_rpms + j)    = res.Cp_std;

            catch
                eta_p_vec(i_end*num_rpms + j) = NaN;
                Ct_vec(i_end*num_rpms + j)    = NaN;
                Cp_vec(i_end*num_rpms + j)    = NaN;
            end
        end
        i_end = i;
    end

    % --- 3. Process, Sort, and Scatter Plot (Separate Figures) ---

    % A. Filtering: Identify valid indices (remove NaN and Inf)
    valid_mask = ~isnan(J_vec) & ~isinf(J_vec) & ...
                 ~isnan(eta_p_vec) & ~isinf(eta_p_vec) & ...
                 ~isnan(Ct_vec) & ~isnan(Cp_vec) & ...
                 (J_vec < J_max_tested); 

    % Extract clean vectors
    J_final     = J_vec(valid_mask);
    eta_p_final = eta_p_vec(valid_mask);
    Ct_final    = Ct_vec(valid_mask);
    Cp_final    = Cp_vec(valid_mask);

    % B. Sorting
    [J_sorted, sort_order] = sort(J_final);
    eta_p_sorted = eta_p_final(sort_order);
    Ct_sorted    = Ct_final(sort_order);
    Cp_sorted    = Cp_final(sort_order);

    % Common Design Point Calculation
    D = 2*R;
    n = outputs.rpm/60; % rpm during flight
    J_target = v / (n * D);
    % disp(J_target); % To see if it is within the limits stablished for the plot

    % =======================================================================
    % COMBINED FIGURE OF PERFORMANCE PLOTS
    % =======================================================================
    figure('Name', '7. Performance Clouds', 'Position', [100 100 1600 550], 'Color', 'w');

    % --- SUBPLOT 1: Efficiency ---
    subplot(1, 3, 1);
    % Position Adjustment
    pos = get(gca, 'Position');
    pos(4) = 0.60; pos(2) = 0.15;
    set(gca, 'Position', pos);

    hold on; grid on;
    scatter(J_sorted, eta_p_sorted, 10, [0 0.4470 0.7410], 'filled'); % Blue
    xline(J_target, '--r', 'LineWidth', 2, ...
        'Label', sprintf('Operating condition ($J=%.2f$)', J_target), ...
        'LabelVerticalAlignment', 'bottom', ...
        'LabelOrientation', 'horizontal', ... 
        'Interpreter', 'latex');

    xlabel('Advance Ratio, $J$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Efficiency, $\eta_p$', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{A. Efficiency Cloud}', 'Interpreter', 'latex', 'FontSize', 14);
    xlim([0, floor(J_max_tested * 10) / 10]); ylim([0, 1]);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);


    % --- SUBPLOT 2: Thrust Coefficient ---
    subplot(1, 3, 2);
    % Position Adjustment
    pos = get(gca, 'Position');
    pos(4) = 0.60; pos(2) = 0.15;
    set(gca, 'Position', pos);

    hold on; grid on;
    scatter(J_sorted, Ct_sorted, 10, [0.8500 0.3250 0.0980], 'filled'); % Orange
    xline(J_target, '--r', 'LineWidth', 2, ...
        'Label', sprintf('Operating condition ($J=%.2f$)', J_target), ...
        'LabelVerticalAlignment', 'bottom', ...
        'LabelOrientation', 'horizontal', ...
        'Interpreter', 'latex');

    xlabel('Advance Ratio, $J$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Thrust Coeff., $C_t$', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{B. Thrust Coefficient Cloud}', 'Interpreter', 'latex', 'FontSize', 14);
    xlim([0, floor(J_max_tested * 10) / 10]);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);


    % --- SUBPLOT 3: Power Coefficient ---
    subplot(1, 3, 3);
    % Position Adjustment
    pos = get(gca, 'Position');
    pos(4) = 0.60; pos(2) = 0.15;
    set(gca, 'Position', pos);

    hold on; grid on;
    scatter(J_sorted, Cp_sorted, 10, [0.9290 0.6940 0.1250], 'filled'); % Yellow
    xline(J_target, '--r', 'LineWidth', 2, ...
        'Label', sprintf('Operating condition ($J=%.2f$)', J_target), ...
        'LabelVerticalAlignment', 'bottom', ...
        'LabelOrientation', 'horizontal', ...
        'Interpreter', 'latex');

    xlabel('Advance Ratio, $J$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Power Coeff., $C_p$', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{C. Power Coefficient Cloud}', 'Interpreter', 'latex', 'FontSize', 14);
    xlim([0, floor(J_max_tested * 10) / 10]);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);


    % =======================================================================
    % GLOBAL TITLE (Robust Centering)
    % =======================================================================
    hTitle = annotation('textbox', [0.5, 0.88, 0, 0], ...
        'String', titles, ...
        'Interpreter', 'latex', ...
        'FontSize', 13, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'EdgeColor', 'black', ...
        'LineWidth', 1, ...
        'BackgroundColor', 'white', ...
        'FitBoxToText', 'on');

    drawnow; % Force render to calculate size
    currentPos = get(hTitle, 'Position'); 
    newWidth = currentPos(3);
    newX = 0.5 - (newWidth / 2);
    set(hTitle, 'Position', [newX, currentPos(2), newWidth, currentPos(4)]);

    fprintf('Total points simulated: %d\n', length(J_vec));
    fprintf('Valid points plotted:   %d\n', length(J_sorted));

   



    %% ====================================================================
    %  8. SAVE ALL PLOTS AS JPG TO 'RESULTS' FOLDER
    %  ====================================================================
    % ================= [BOTTOM BLOCK: START] =================
    fprintf('--- Finalizing and Saving Plots ---\n');
    
    % 1. Create results directory if it doesn't exist
    if ~exist('results', 'dir'), mkdir('results'); end
    
    % 2. Find all figures created by the code above
    allFigs = findobj('Type', 'figure');
    
    % 3. Sort them by Figure Number (1, 2, 3...) so they save in order
    if ~isempty(allFigs)
        [~, idx] = sort([allFigs.Number]);
        allFigs = allFigs(idx);
    end

    % 4. SAVE LOOP (Figures remain INVISIBLE here)
    % We use 'exportgraphics' or 'print' because 'saveas' requires visibility.
    fprintf('Saving images (this may take a moment)...\n');
    
    for i = 1:length(allFigs)
        hFig = allFigs(i);
        
        % A. Resize the invisible figure to be large (Full HD)
        % This ensures the saved image has good proportions/layout
        set(hFig, 'Units', 'pixels');
        set(hFig, 'Position', [100 100 1920 1080]); 
        
        % B. Construct Filename
        fName = hFig.Name;
        if isempty(fName), fName = sprintf('Figure_%d', hFig.Number); end
        safeName = regexprep(strrep(fName, ' ', '_'), '[^a-zA-Z0-9_]', '');
        fullPath = fullfile('results', [safeName, '.jpg']);
        
        % C. Save using specific renderers that support invisible figures
        try
            % OPTION 1: Modern MATLAB (R2020a+) - Best Quality
            % exportgraphics handles invisible OpenGL figures perfectly.
            exportgraphics(hFig, fullPath, 'Resolution', 300); 
        catch
            % OPTION 2: Legacy Fallback (R2019b and older)
            % 'print' forces a render even if invisible. 
            % -r150 defines 150 DPI resolution.
            print(hFig, fullPath, '-djpeg', '-r150');
        end
        
        fprintf('Saved: %s\n', fullPath);
    end

    % 5. SHOW LOOP (Restore visibility now that saving is done)
    fprintf('Displaying figures...\n');
    
    % Restore global default
    set(0, 'DefaultFigureVisible', 'on');
    
    for i = length(allFigs):-1:1
        hFig = allFigs(i);
        set(hFig, 'Visible', 'on');
        
        % Maximize the window on screen for the user
        try
            hFig.WindowState = 'maximized';
        catch
            % Fallback for older versions
            set(hFig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        end
        drawnow;
    end
    
    fprintf('Done. Check the ''results'' folder.\n');
    % ================= [BOTTOM BLOCK: END] =================



    
end







% =====================================================================
% Local helpers: filled arrow primitives (same as your finalized style)
% =====================================================================
function hArrow(ax, P0, P1, lw, col, headL, headW)
    P0 = P0(:); P1 = P1(:);
    line(ax, [P0(1) P1(1)], [P0(2) P1(2)], 'Color', col, 'LineWidth', lw);
    hHeadPatch(ax, P0, P1, headL, headW, col);
end

function hHeadPatch(ax, P0, P1, headL, headW, col)
    P0 = P0(:); P1 = P1(:);
    v = P1 - P0;
    nv = norm(v);
    if nv < 1e-12, return; end
    u = v / nv;
    p = [-u(2); u(1)];

    tip  = P1;
    base = P1 - headL*u;
    p1 = base + (headW/2)*p;
    p2 = base - (headW/2)*p;

    patch(ax, [tip(1) p1(1) p2(1)], [tip(2) p1(2) p2(2)], col, ...
        'EdgeColor', col, 'FaceColor', col);
end


% =====================================================================
% Nested helper: render one section in the finalized style
% =====================================================================
function plotAirfoilStressSection(ax, idx, shortTitle, cax, outputs, prop)
    % axes(ax);        % Removed to prevent forcing figure visibility
    hold(ax,'on'); axis(ax,'equal'); box(ax,'on');

    % Extract from SABEMMT the structural results needed
    c = outputs.c;
    theta = outputs.theta;
    theta_deg = rad2deg(theta);
    sigma_total_pts = outputs.sigma_total_pts;
    M_x = outputs.M_x;
    M_y = outputs.M_y;
   
    % Extract the geometric characteristics of the airfoil
    geo = getAirfoilInertias(prop.airfoil_data.x_u, prop.airfoil_data.y_u, ...
                             prop.airfoil_data.x_l, prop.airfoil_data.y_l);
    theta_p = geo.theta_p;
    theta_p_deg = rad2deg(theta_p);


    % Grid style
    grid(ax,'on');
    ax.GridAlpha      = 0.15;
    ax.MinorGridAlpha = 0.08;
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 11;

    % Section rotation
    th = theta(idx);
    Rcw = [cos(-th) -sin(-th); sin(-th) cos(-th)];

    % Airfoil geometry (centroid at origin) scaled by chord
    x_sec = geo.x_pts(:) * c(idx);
    y_sec = geo.y_pts(:) * c(idx);

    XY = Rcw * [x_sec.'; y_sec.'];
    xR = XY(1,:).';
    yR = XY(2,:).';

    % Stress in MPa
    sig_MPa = sigma_total_pts(:, idx) / 1e6;

    % Upper/lower split (by rotated y)
    idx_upper = (yR >= 0);
    idx_lower = ~idx_upper;

    % Limits (consistent per subplot based on this section)
    r_af = max(hypot(xR,yR));
    lim_span = 1.30 * max([r_af, 1e-6]);
    xlim(ax, [-lim_span, lim_span]);
    ylim(ax, [-lim_span, lim_span]);

    % Styles
    axisLW  = 1.5;
    colAxis = [0 0 0];
    colLoc  = [0.15 0.15 0.15];
    colComp = [0.35 0.35 0.35];
    colMom  = [0.05 0.05 0.05];
    dotLW   = 1.2;

    headL = 0.06 * lim_span;
    headW = 0.035 * lim_span;

    % Rotor axes lines
    line(ax, [-lim_span, lim_span], [0 0], 'Color','k', 'LineWidth',axisLW);
    line(ax, [0 0], [-lim_span, lim_span], 'Color','k', 'LineWidth',axisLW);

    % Filled arrowheads on + axes
    hHeadPatch(ax, [0;0], [0.98*lim_span;0], headL, headW, colAxis);
    hHeadPatch(ax, [0;0], [0;0.98*lim_span], headL, headW, colAxis);

    % Axis labels (with y_rotor nudged right)
    text(ax, 0.98*lim_span, 0.01*lim_span, '$x_{\mathrm{rotor}}$', ...
        'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','bottom', ...
        'FontSize',11,'Color','k');
    text(ax, 0.02*lim_span, 0.98*lim_span, '$y_{\mathrm{rotor}}$', ...
        'Interpreter','latex','HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize',11,'Color','k');

    % Local axes (airfoil axes)
    e1 = Rcw * [1;0];   % +x_airfoil
    e2 = Rcw * [0;1];   % +y_airfoil
    axis_len = 0.85*lim_span;

    hArrow(ax, [0;0], [axis_len*e1(1); axis_len*e1(2)], axisLW, colLoc, headL, headW);
    hArrow(ax, [0;0], [axis_len*e2(1); axis_len*e2(2)], axisLW, colLoc, headL, headW);

    text(ax, axis_len*e1(1), axis_len*e1(2), '$x_{\mathrm{airfoil}}$', ...
        'Interpreter','latex','FontSize',11,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',colLoc);
    text(ax, axis_len*e2(1), axis_len*e2(2), '$y_{\mathrm{airfoil}}$', ...
        'Interpreter','latex','FontSize',11,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',colLoc);

    % --- Principal inertia axes (x_p, y_p) ---
    % Compute principal angle from section inertias about centroid.
    % Convention: rotation in the (x_airfoil,y_airfoil) plane.      
    e_xp =  cos(theta_p)*e1 + sin(theta_p)*e2;
    e_yp = -sin(theta_p)*e1 + cos(theta_p)*e2;
    
    colPr = [0 0 1];
    L_xp = 0.55*axis_len;
    L_yp = 0.40*axis_len;
    
    hArrow(ax, [0;0], [L_xp*e_xp(1); L_xp*e_xp(2)], 1.4, colPr, headL, headW);
    hArrow(ax, [0;0], [L_yp*e_yp(1); L_yp*e_yp(2)], 1.4, colPr, headL, headW);
    
    % Label placement: offset slightly normal to each axis to avoid overlap
    off_xp = 0.060*lim_span*[-e_xp(2); e_xp(1)];
    off_yp = 0.060*lim_span*[-e_yp(2); e_yp(1)];
    
    text(ax, L_xp*e_xp(1)+off_xp(1), L_xp*e_xp(2)+off_xp(2), '$x_{\mathrm{princ.inertia}}$', ...
        'Interpreter','latex','FontSize',10,'Color',colPr, ...
        'HorizontalAlignment','left','VerticalAlignment','bottom', ...
        'BackgroundColor','w','Margin',1);
    
    text(ax, L_yp*e_yp(1)+off_yp(1), L_yp*e_yp(2)+off_yp(2), '$y_{\mathrm{princ.inertia}}$', ...
        'Interpreter','latex','FontSize',10,'Color',colPr, ...
        'HorizontalAlignment','left','VerticalAlignment','bottom', ...
        'BackgroundColor','w','Margin',1);
    % Principal Angle Arc
    % The arc must start at the Airfoil X-axis (-theta) and span theta_p
    ang_start = -theta(idx); 
    ang_end   = ang_start + theta_p;
    % Scale radius dynamically to the view limits (e.g., 20% of the view)
    % instead of hardcoding 0.25 meters (which is often way off-screen).
    r_arc_p = 0.2 * lim_span; 
    t_p = linspace(ang_start, ang_end, 30); 
    plot(ax, r_arc_p*cos(t_p), r_arc_p*sin(t_p), 'b', 'LineWidth', 1);
    % Position text slightly beyond the arc
    ang_mid = (ang_start + ang_end) / 2;
    text_r  = r_arc_p * 1.65; % Text at 165% of arc radius
    text(ax, text_r*cos(ang_mid+0.1), text_r*sin(ang_mid+0.1), ...
         sprintf('$\\theta_p=%.1f^\\circ$', theta_p_deg), ...
         'Color', 'b', 'Interpreter', 'latex', 'FontSize', 10, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    % Airfoil outline + stress scatter (same style)
    plot(ax, [xR; xR(1)], [yR; yR(1)], 'k-', 'LineWidth', 1.0);
    scatter(ax, xR(idx_upper), yR(idx_upper), 28, sig_MPa(idx_upper), 'filled', ...
        'MarkerEdgeColor','none','MarkerFaceAlpha',0.95);
    scatter(ax, xR(idx_lower), yR(idx_lower), 28, sig_MPa(idx_lower), 'filled', ...
        'MarkerEdgeColor','none','MarkerFaceAlpha',0.95);

    % Common colormap scale
    caxis(ax, cax);

    % Max stress point marker: colormap-consistent + label sigma_max
    [sig_max_val, iMax] = max(sig_MPa);
    xMax = xR(iMax); yMax = yR(iMax);

    % Map sigma_max to colormap color
    cmap = colormap(ax);
    t = (sig_max_val - cax(1)) / max(eps, (cax(2)-cax(1)));
    t = min(max(t,0),1);
    ci = 1 + round(t*(size(cmap,1)-1));
    colMax = cmap(ci,:);

    plot(ax, xMax, yMax, 'o', 'MarkerSize', 10, 'MarkerFaceColor', colMax, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.0);

    % small label near marker
    sgnX = sign(xMax); if sgnX==0, sgnX=1; end
    sgnY = sign(yMax); if sgnY==0, sgnY=1; end
    off = 0.035*lim_span*[sgnX; sgnY];
    % text(ax, xMax+off(1), yMax+off(2), '$\sigma_{\max}$', ...
    %     'Interpreter','latex','FontSize',10,'Color','k', ...
    %     'HorizontalAlignment','left','VerticalAlignment','bottom', ...
    %     'BackgroundColor','w','Margin',1);

    % Moment resultant + projections (single-headed resultant, projections from origin)
    Msec = [M_x(idx); M_y(idx)];
    vMx = Msec(1) * e1;
    vMy = Msec(2) * e2;
    vM  = vMx + vMy;

    Mmag = norm(vM); if Mmag < 1e-12, Mmag = 1e-12; end
    m_scale = 0.72*lim_span / Mmag;

    VMx = m_scale*vMx;
    VMy = m_scale*vMy;
    VM  = m_scale*vM;

    % Resultant vector (single head)
    hArrow(ax, [0;0], [VM(1); VM(2)], 2.0, colMom, headL, headW);
    text(ax, VM(1), 1.2*VM(2), '$\mathbf{M}_{\mathrm{res}}$', ...
        'Interpreter','latex','FontSize',11,'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom','BackgroundColor','w','Margin',1);

    % Projections from origin
    hArrow(ax, [0;0], [VMx(1); VMx(2)], 1.6, colComp, headL, headW);
    hArrow(ax, [0;0], [VMy(1); VMy(2)], 1.6, colComp, headL, headW);

    % Dotted construction lines to resultant tip
    plot(ax, [VMx(1), VM(1)], [VMx(2), VM(2)], ':', 'Color', colComp, 'LineWidth', dotLW);
    plot(ax, [VMy(1), VM(1)], [VMy(2), VM(2)], ':', 'Color', colComp, 'LineWidth', dotLW);

    % Theta arc (this section’s theta): arc from rotor plane (-x_rotor) to prolongation of -x_airfoil
    v_rot_neg = [-1; 0];
    v_sec_neg = -e1;

    L_ext = 1.35*r_af;
    plot(ax, [0, L_ext*v_sec_neg(1)], [0, L_ext*v_sec_neg(2)], 'k:', 'LineWidth', 1.2);

    arc_R = 1.12*r_af;
    a0 = atan2(v_rot_neg(2), v_rot_neg(1));
    a1 = atan2(v_sec_neg(2), v_sec_neg(1));
    da = atan2(sin(a1-a0), cos(a1-a0));
    t_arc = linspace(a0, a0+da, 120);
    plot(ax, arc_R*cos(t_arc), arc_R*sin(t_arc), 'k-', 'LineWidth', 1.2);

    % theta label (left)
    theta_sec_deg = theta_deg(idx);
    tmid = a0 + 0.55*da;
    posT = 1.08*arc_R*[cos(tmid); sin(tmid)];
    text(ax, posT(1), posT(2), sprintf('$\\theta=%.1f^\\circ$', theta_sec_deg), ...
        'Interpreter','latex','FontSize',11,'HorizontalAlignment','right', ...
        'VerticalAlignment','middle','BackgroundColor','w','Margin',1);

    % Southwest legend (inside axes)
    sigmax = max(sigma_total_pts(:,idx))/1e6;
    % Matlab uses scientific notation only when needed:
    legStr = {
        sprintf('$\\sigma_{max}=%.1f\\,\\mathrm{MPa}$', sigmax)
        sprintf('$M_x=%.4g\\,\\mathrm{Nm}$', Msec(1))
        sprintf('$M_y=%.4g\\,\\mathrm{Nm}$', Msec(2))
    };

    text(ax, -0.90*lim_span, -0.90*lim_span, legStr, ...
        'Interpreter','latex','FontSize',11, ...
        'HorizontalAlignment','left','VerticalAlignment','bottom', ...
        'BackgroundColor','w','EdgeColor',[0 0 0],'Margin',6);

    % Title per subplot (LaTeX; supports multiline via cell array)
    title(ax, shortTitle, 'Interpreter','latex','FontSize',12);

end
