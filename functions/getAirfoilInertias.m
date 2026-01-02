%% ---------- Get Airfoil geometric characteristics and Inertias ----------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: getAirfoilInertias
%   Project : sabemmt
%   Author: Miguel Frade
%   Affiliation at time of publication: Universidad Politecnica de Madrid
%   Date: December 2025
%   License: Creative Commons Attribution-NonCommercial 4.0 (CC BY-NC 4.0)
%
%   Description:
%     Computes normalized area, centroid, and second moments of area of an
%     airfoil cross-section using Green's Theorem.
%  
%   Inputs:
%     x_u, y_u : Upper surface coordinates (LE -> TE)
%     x_l, y_l : Lower surface coordinates (LE -> TE or TE -> LE)
%  
%   Outputs (normalized by chord):
%     geo.A_norm      : Area (A / c^2)
%     geo.I_x_norm    : Ixx about centroid (flapwise, / c^4)
%     geo.I_y_norm    : Iyy about centroid (edgewise, / c^4)
%     geo.I_xy_norm   : Ixy about centroid ( / c^4)
%     geo.xc_norm     : Centroid x-location ( / c )
%     geo.yc_norm     : Centroid y-location ( / c )
%     geo.x_pts       : Centered x coordinates (for stress evaluation)
%     geo.y_pts       : Centered y coordinates
%  
%   Notes:
%     - Automatically normalizes coordinates to chord = 1
%     - Guarantees CCW polygon orientation
%     - Robust against reversed lower-surface ordering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function geo = getAirfoilInertias(x_u, y_u, x_l, y_l)

    % ---------------------------------------------------------------------
    % 1. Normalize coordinates to chord = 1
    % ---------------------------------------------------------------------
    x_all = [x_u(:); x_l(:)];
    min_x = min(x_all);
    max_x = max(x_all);
    c = max_x - min_x;
    
    if abs(min_x) > 1e-6 || abs(c - 1) > 1e-6
        % Shift LE to x = 0
        x_u = x_u - min_x;
        x_l = x_l - min_x;
    
        % Scale by chord
        x_u = x_u / c;   y_u = y_u / c;
        x_l = x_l / c;   y_l = y_l / c;
    end
    
    % ---------------------------------------------------------------------
    % 2. Build closed polygon (upper: LE->TE, lower: TE->LE)
    % ---------------------------------------------------------------------
    x_poly = [x_u(:); flipud(x_l(:))];
    y_poly = [y_u(:); flipud(y_l(:))];
    
    % Close polygon explicitly
    x_poly(end+1) = x_poly(1);
    y_poly(end+1) = y_poly(1);
    
    % ---------------------------------------------------------------------
    % 3. Enforce counter-clockwise orientation
    % ---------------------------------------------------------------------
    signedArea = 0.5 * sum( ...
        x_poly(1:end-1).*y_poly(2:end) - ...
        x_poly(2:end).*y_poly(1:end-1) );
    
    if signedArea < 0
        x_poly = flipud(x_poly);
        y_poly = flipud(y_poly);
    end
    
    % ---------------------------------------------------------------------
    % 4. Green's Theorem Integrals
    % ---------------------------------------------------------------------
    N = length(x_poly) - 1;
    
    A = 0;  Mx = 0;  My = 0;
    Ixx = 0; Iyy = 0; Ixy = 0;
    
    for i = 1:N
        xi  = x_poly(i);     yi  = y_poly(i);
        xj  = x_poly(i+1);   yj  = y_poly(i+1);
    
        a = xi*yj - xj*yi;
    
        A  = A  + a;
        Mx = Mx + (yi + yj) * a;
        My = My + (xi + xj) * a;
    
        Ixx = Ixx + (yi^2 + yi*yj + yj^2) * a;
        Iyy = Iyy + (xi^2 + xi*xj + xj^2) * a;
        Ixy = Ixy + ( ...
              xi*yj + 2*xi*yi + 2*xj*yj + xj*yi ) * a;
    end
    
    % ---------------------------------------------------------------------
    % 5. Finalize integrals
    % ---------------------------------------------------------------------
    A   = A / 2;
    xc  = (My / 6) / A;
    yc  = (Mx / 6) / A;
    
    Ixx = Ixx / 12;
    Iyy = Iyy / 12;
    Ixy = Ixy / 24;
    
    % ---------------------------------------------------------------------
    % 6. Shift to centroid (Parallel Axis Theorem)
    % ---------------------------------------------------------------------
    Ixx_c = Ixx - A * yc^2;
    Iyy_c = Iyy - A * xc^2;
    Ixy_c = Ixy - A * xc * yc;
    
    % ---------------------------------------------------------------------
    % 7. Centered coordinates (useful for stress computations)
    % ---------------------------------------------------------------------
    x_centered = x_poly(1:end-1) - xc;
    y_centered = y_poly(1:end-1) - yc;


    % ---------------------------------------------------------------------
    % 8. Principal inertias (we don't use them in SABEMMT)
    % ---------------------------------------------------------------------
    I_avg = 0.5*(Ixx_c + Iyy_c);
    R = sqrt((0.5*(Ixx_c - Iyy_c))^2 + Ixy_c^2);
    I1 = I_avg + R; % max
    I2 = I_avg - R; % min

    % Calculate Principal Angle (theta_p)
    % tan(2*theta_p) = -2*Ixy / (Ix - Iy)
    numerator = -2 * Ixy_c;
    denominator = Ixx_c - Iyy_c;
    theta_p = 0.5 * atan2(numerator, denominator);
    theta_p_deg = rad2deg(theta_p);
    
    % ---------------------------------------------------------------------
    % 9. Pack output
    % ---------------------------------------------------------------------
    geo.A_norm    = A;
    geo.I_x_norm  = Ixx_c;
    geo.I_y_norm  = Iyy_c;
    geo.I_xy_norm = Ixy_c;
    geo.xc_norm   = xc;
    geo.yc_norm   = yc;

    geo.x_pts = x_centered;
    geo.y_pts = y_centered;
    geo.theta_p = theta_p;
    geo.theta_p_deg = theta_p_deg;

    geo.I1 = I1;
    geo.I2 = I2;

end






% function geo = computeAirfoilInertias(x_u, y_u, x_l, y_l)
% 
%     % COMPUTEAIRFOILGEOMETRY Computes exact geometric properties of the airfoil
%     % using Green's Theorem for a closed polygon.
%     %
%     %   Inputs: Arrays of Upper and Lower surface coordinates.
%     %   Outputs: Struct containing Normalized Area, Inertias, and Centered Points.
%     %
%     %   Note: If the input coordinates are not normalized (e.g. defined in mm 
%     %   or inches), this function will automatically Normalize them to c=1 
%     %   before computing the properties.
% 
%     % 0. Normalization Check & Fix
%     % Find the physical bounds of the input data
%     x_all = [x_u, x_l];
%     min_x = min(x_all);
%     max_x = max(x_all);
%     raw_chord = max_x - min_x;
% 
%     % If data is not roughly 0 to 1, we normalize it.
%     % We use a small tolerance (1e-3) to detect if scaling is needed.
%     if abs(min_x) > 1e-3 || abs(raw_chord - 1) > 1e-3
%         fprintf('  > Airfoil coordinates are not normalized (x from 0 to 1)');
%         fprintf('    Normalizing airfoil coordinates to c=1 (Scale factor: %.4f)...\n', 1/raw_chord);
% 
%         % Shift x to start at 0
%         x_u = x_u - min_x;
%         x_l = x_l - min_x;
% 
%         % Scale both x and y by the raw chord length to preserve aspect ratio
%         x_u = x_u / raw_chord;
%         y_u = y_u / raw_chord;
%         x_l = x_l / raw_chord;
%         y_l = y_l / raw_chord;
%     end
% 
%     % 1. Form a closed polygon (Counter-Clockwise)
%     % x_u goes 0->1. x_l goes 0->1. We need 0->1 (Upper) then 1->0 (Lower).
%     x = [x_u, fliplr(x_l)];
%     y = [y_u, fliplr(y_l)];
% 
%     % Close the loop explicitly if needed
%     if (x(1) ~= x(end)) || (y(1) ~= y(end))
%         x(end+1) = x(1);
%         y(end+1) = y(1);
%     end
% 
%     N = length(x) - 1;
% 
%     % 2. Compute Properties using Green's Theorem (Polygon Formulas)
%     % Initialize
%     A = 0; M_x = 0; M_y = 0;
%     I_xx = 0; I_yy = 0; I_xy = 0;
% 
%     for i = 1:N
%         % Common term a_i = x_i * y_{i+1} - x_{i+1} * y_i
%         a_i = x(i)*y(i+1) - x(i+1)*y(i);
% 
%         % Area
%         A = A + a_i;
% 
%         % First Moments (for Centroid)
%         M_x = M_x + (y(i) + y(i+1)) * a_i;
%         M_y = M_y + (x(i) + x(i+1)) * a_i;
% 
%         % Second Moments (Inertias about origin 0,0)
%         I_xx = I_xx + (y(i)^2 + y(i)*y(i+1) + y(i+1)^2) * a_i;
%         I_yy = I_yy + (x(i)^2 + x(i)*x(i+1) + x(i+1)^2) * a_i;
%         I_xy = I_xy + (x(i)*y(i+1) + 2*x(i)*y(i) + 2*x(i+1)*y(i+1) + x(i+1)*y(i)) * a_i;
%     end
% 
%     % Finalize Integration
%     A = A / 2;
%     x_c = (M_y / 6) / A;
%     y_c = (M_x / 6) / A;
% 
%     I_xx = I_xx / 12;
%     I_yy = I_yy / 12;
%     I_xy = I_xy / 24;
% 
%     % 3. Shift to Centroid (Parallel Axis Theorem)
%     I_x_cent = I_xx - A * y_c^2;
%     I_y_cent = I_yy - A * x_c^2;
%     I_xy_cent = I_xy - A * x_c * y_c;
% 
%     % 4. Center the coordinates (for stress search loop)
%     x_centered = x - x_c;
%     y_centered = y - y_c;
% 
%     % 5. Pack Output
%     geo.A_norm = A;
%     geo.I_x_norm = I_x_cent;   % Flapwise Inertia (about chord)
%     geo.I_y_norm = I_y_cent;   % Edgewise Inertia (about thickness)
%     geo.I_xy_norm = I_xy_cent; % Product of Inertia
%     geo.xc_norm = x_c;         % Centroid location (usually ~0.4)
%     geo.yc_norm = y_c;
% 
%     % Points to iterate over for max stress (remove the duplicate closing point)
%     geo.x_pts = x_centered(1:end-1);
%     geo.y_pts = y_centered(1:end-1);
% end