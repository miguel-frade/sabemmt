%% ----------------------------- SABEMMT ------------------------------- %%
%  ------ (Structures And Blade Element Modified Momentum Theory) ------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: runSABEMMT
%   Project: sabemmt
%   Author: Miguel Frade
%   Affiliation at time of publication: Universidad Politecnica de Madrid
%   Date: December 2025
%   License: MIT
%
%   Description:
%       A robust high-fidelity BEMT solver that implements an Engineering 
%       Modification of the classical Momentum Theory for multi-regime 
%       simulation, considers the effects of Reynolds number at different 
%       sections of the blade and implements Prandtl tip loss function.
%       This code is a general-purpose aerodynamic and structural model
%       for the analysis and optimization of:
%           - Airplane propellers
%           - Helicopter rotors in axial flight
%           - Multirotor propellers in axial flight
%
%       THEORETICAL BASIS:
%       The Momentum Theory used here includes an Engineering Modification
%       to correctly handle the Turbulent Wake State and Vortex Ring State
%       regimes (Cuerva et. al. 2006).
%       The Structural Model is a conservative rectangle approximation of
%       the blade sections that considers the rotation of forces according
%       to the pitch distribution of the blade and computes the total
%       stress at the point of maximum stress of each section.
%
%       Reference:
%       Cuerva, A., Sanz-Andrés, A., Meseguer, J., & Espino, J. L. (2006). 
%       "An Engineering Modification of the Blade Element Momentum Equation 
%       for Vertical Descent". Journal of the American Helicopter Society, 
%       51(4), 349-354.
%
%   Inputs:
%       v, rpm      - Operation point
%       Nb, R       - Number of Blades and Radius
%       c, theta    - Geometry distributions
%       prop_const  - Airfoil polars and material data
%       environment - Atmospheric data
%
%   Outputs:
%       out - Struct containing:
%             * Aero: T, Q, P, eta_p, Ct, Cp, FoM
%             * Distributions: dT/dx, dQ/dx, alpha, lambda along radius
%             * Structural: Centrifugal and Bending Stresses
%
%       The output coefficients Ct_std, Cq_std, Cp_std follow the standard
%       definitions for airplane propellers:
%         Ct = T / (rho*rps^2*D^4);  Where rps = revs per second
%         Cq = Q / (rho*rps^2*D^5);
%         Cp = P / (rho*rps^3*D^5) = 2*pi*Cq;
%       But all the intermediate calculations will use the definitions 
%       of Ct, Cq, Cp as typically used for helicopters (also outputted).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = runSABEMMT(v, rpm, Nb, R, c, theta_deg, prop_const, environment)

    % % Heartbit to see that the code is running:
    % fprintf('Running SABEMMT...\n'); 
    
    % --- 1. UNPACK & PRE-CALCULATE CONSTANTS ---
    rho = environment.rho;
    nu = environment.nu;
    rho_mat = prop_const.rho_mat;
    mesh = prop_const.mesh; 
    dx = prop_const.dx; 
    Fc = prop_const.airfoil.Fc;   % Expecting 2D griddedInterpolant {alpha, Re}
    Fd = prop_const.airfoil.Fd;   % Expecting 2D griddedInterpolant {alpha, Re}
    alpha_min = prop_const.airfoil.alpha_min; 
    alpha_max = prop_const.airfoil.alpha_max; 
    cd_max_val = prop_const.airfoil.cd_max_val;
    cla_data = prop_const.airfoil.cla_data;
    cl0_data = prop_const.airfoil.cl0_data; 

    Area = prop_const.airfoil.A_norm .* (c.^2);
    t = prop_const.airfoil.t_norm .* c;
    I_x = prop_const.airfoil.I_x_norm .* (c.^4); 
    I_y = prop_const.airfoil.I_y_norm .* (c.^4);

    omega = rpm * (2*pi/60);
    theta = theta_deg .* (pi/180); 
    S = pi * (R^2);
    r = R .* mesh;

    omegaR = omega * R;
    sigma = c .* (Nb / (pi*R)); 
    mu = v / omegaR; 
    
    % --- 2. BEMT VECTORIZED SOLVER ---
    
    % Initial Guess
    cte1 = sigma .* cla_data / 8;
    mu_over_mesh = mu ./ mesh;
    cte2 = 1 + (mu_over_mesh.^2);
    D_init = (mu + cte1.*cte2).^2 + 4*cte1.*cte2.*((mesh.*theta - mu) + (cl0_data/cla_data).*mesh);
    D_init = max(D_init, 0);
    lambda = 0.5 .* (-mu - cte1.*cte2 + sqrt(D_init)); 
    
    % Constants (for MMT Model: Álvaro Cuerva et. al. 2006)
    A_coeff = 0.745; 
    B_coeff = 0.447; 
    B_sq_mu_sq = B_coeff^2 * mu^2;
    
    % Pre-calculate inverse for speed
    inv_denom_const = 1 ./ (8 * A_coeff .* mesh); 
    
    TOL = 1e-8;
    MAX_ITER = 100; 
    
    % consts struct
    consts = struct('mu', mu, 'Nb', Nb, 'A', A_coeff, 'B_sq_mu_sq', B_sq_mu_sq, ...
                    'alpha_min', alpha_min, 'alpha_max', alpha_max, ...
                    'cd_max_val', cd_max_val, ...
                    'omegaR', omegaR, 'nu', nu);


    % --- 2.1 ALGEBRAIC VECTORIZED ITERATION ---
    converged_mask = false(size(lambda)); 
    x = mesh;
    
    for k = 1:MAX_ITER
        lam_old = lambda;
        
        % 1. Geometry
        y = mu + lambda;
        
        % Hypotenuse 'h' (Normalized Velocity Magnitude)
        h = sqrt(x.^2 + y.^2); 
        
        % 2. Prandtl Tip Loss
        abs_y = abs(y);
        abs_y(abs_y < 1e-7) = 1e-7; 
        exponent = (0.5 * Nb) .* (x - 1) .* h ./ abs_y;
        arg = exp(exponent);
        arg(arg > 1) = 1; arg(arg < -1) = -1;
        F = (2/pi) * acos(arg);
        F(F < 1e-6) = 1e-6; 
        
        % 3. Aerodynamics (Reynolds Aware)
        phi = atan(y ./ x);
        alpha = theta - phi;
        
        % Reynolds Number Calculation: Re = (V_local * chord) / nu
        % V_local = h * omegaR
        Re_local = (h .* omegaR .* c) ./ nu;
        
        % 2D Interpolation
        cl = Fc(alpha, Re_local);
        cd = Fd(alpha, Re_local);
        
        outside_idx = (alpha < alpha_min) | (alpha > alpha_max);
        if any(outside_idx)
            cl(outside_idx) = 0; 
            cd(outside_idx) = cd_max_val;
        end
        
        % 4. Forces & Update
        force_term = sigma .* h .* (cl .* x - cd .* y);
        D_i = B_sq_mu_sq + y.^2; 
        denom_part = F .* sqrt(D_i);
        denom_part(denom_part < eps) = eps;
        
        lambda_new = force_term ./ denom_part .* inv_denom_const;
        
        % Damping (for convergence)
        lambda = 0.6 * lambda_new + 0.4 * lam_old;
        
        % Check convergence
        diff = abs(lambda - lam_old);
        converged_mask = diff < TOL;
        
        if all(converged_mask)
            break;
        end
    end
    
    % --- 2.2 ROBUST FALLBACK FOR UNCONVERGED SECTIONS USING FZERO IF NECCESSARY ---
    bad_idx = find(~converged_mask | isnan(lambda) | isinf(lambda));
    if ~isempty(bad_idx)
        for i = bad_idx
             a = -10 * abs(mu) - 20;
             b = 10 + 10 * abs(mu);
             try
                 % Solve the flow at one section at a time
                 fa = scalar_fun(a, mesh(i), sigma(i), theta(i), c(i), Fc, Fd, consts);
                 fb = scalar_fun(b, mesh(i), sigma(i), theta(i), c(i), Fc, Fd, consts);
                 if fa * fb < 0
                     lambda(i) = fzero(@(lam) scalar_fun(lam, mesh(i), sigma(i), theta(i), c(i), Fc, Fd, consts), ...
                                       [a, b], optimset('Display','off','TolX',TOL));
                 end
             catch
             end
        end
        % Re-calculate critical variables for post-processing
        y = mu + lambda;
        h = sqrt(x.^2 + y.^2);
        phi = atan(y ./ x);
        alpha = theta - phi;
        
        % Re-calculate Re for final outputs
        Re_local = (h .* omegaR .* c) ./ nu;
        cl = Fc(alpha, Re_local); 
        cd = Fd(alpha, Re_local);
        
        outside = (alpha < alpha_min) | (alpha > alpha_max);
        if any(outside), cl(outside)=0; cd(outside)=cd_max_val; end
    end
    
    % --- 3. FINAL CALCULATIONS ---
    common_h = 0.5 .* sigma .* h;
    dCtdx = common_h .* (cl .* x - cd .* y);
    dCFtdx = common_h .* (cd .* x + cl .* y);
    
    rho_S_omegaR_sq = rho * S * omegaR^2;
    dTdx = rho_S_omegaR_sq .* dCtdx;
    dFtdx = rho_S_omegaR_sq .* dCFtdx;
    
    Ct = sum(dCtdx .* dx);
    dCqdx = dCFtdx .* mesh;
    Cq = sum(dCqdx .* dx);
    
    T = rho_S_omegaR_sq * Ct;
    Q = rho * S * R * omegaR^2 * Cq;
    P = rho * S * omegaR^3 * Cq;
    
    if abs(P) > 1e-7
        eta_p = v * T / P; % Propulsive efficiency
        FoM = Ct*sqrt(0.5*Ct)/Cq; % Figure Of Merit
    else
        warning('P is close to zero (abs(P) <= tolerance_P).');
        eta_p = NaN; FoM = NaN;
    end
    
    % Construct standard coefficients for out struct
    Ct_std = Ct * pi^3 / 4; 
    Cq_std = Cq * pi^3 / 8;
    Cp_std = 2 * pi * Cq;
    
    % Mach number at the tip
    M_tip = sqrt(omegaR^2 + v^2)/environment.sound_speed; % We ingnore v_induced(x)

    % --- 4. STUCTURAL ANALYSIS ---
    dT = dTdx .* dx; dFt = dFtdx .* dx;
    dT_vec = dT(:)'; dFt_vec = dFt(:)'; r_vec = r(:)';
    
    S1_dT = fliplr(cumsum(fliplr(dT_vec)));
    S1_dFt = fliplr(cumsum(fliplr(dFt_vec)));
    Sr_dT = fliplr(cumsum(fliplr(dT_vec .* r_vec)));
    Sr_dFt = fliplr(cumsum(fliplr(dFt_vec .* r_vec)));
    
    V_y_rotor = S1_dT; V_x_rotor = S1_dFt;
    M_x_rotor = Sr_dT - r_vec .* S1_dT;
    M_y_rotor = - (Sr_dFt - r_vec .* S1_dFt);
    
    cos_ntheta = cos(-theta); sin_ntheta = sin(-theta);
    M_x =  M_x_rotor .* cos_ntheta + M_y_rotor .* sin_ntheta;
    M_y = -M_x_rotor .* sin_ntheta + M_y_rotor .* cos_ntheta;
    
    sigma_bend_max_Ix = M_x .* (t/2) ./ I_x;
    sigma_bend_max_Iy = - M_y .* (c/2) ./ I_y;
    
    rho_mat_omega_sq = rho_mat * omega^2;
    integrand = rho_mat_omega_sq .* Area .* r;
    cum_from_root = cumtrapz(r, integrand); 
    F_cf_from_i = cum_from_root(end) - cum_from_root; 
    
    safe_A = Area; safe_A(safe_A < 1e-8) = 1e-8; 
    sigma_cf = F_cf_from_i ./ safe_A; sigma_cf(end) = 0; 
    sigma_total_max = sigma_bend_max_Ix + sigma_bend_max_Iy + sigma_cf;
    
    out = struct('R', R, 'rpm', rpm, 'T', T, 'Q', Q, 'P', P, ...
                 'Ct', Ct, 'Cq', Cq, 'Ct_std', Ct_std, 'Cq_std', Cq_std, ...
                 'Cp_std', Cp_std, 'eta_p', eta_p, 'FoM', FoM, ...
                 'mesh', mesh, 'r', r, ...
                 'c', c, 'theta', theta, 'phi', phi, 'alpha', alpha, ...
                 'lambda', lambda, 'cl', cl, 'cd', cd, 'M_tip', M_tip, ...
                 'dTdx', dTdx,'dFtdx', dFtdx, 'dT', dT, ...
                 'dFt', dFt, 'V_x', V_x_rotor, 'V_y', V_y_rotor, ...
                 'M_x', M_x, 'M_y', M_y, 'sigma_cf', sigma_cf, ...
                 'sigma_bend_max_Ix', sigma_bend_max_Ix, ...
                 'sigma_bend_max_Iy', sigma_bend_max_Iy, ...
                 'sigma_total_max', sigma_total_max);
end

% ------------------ LOCAL HELPER FUNCTION: scalar_fun --------------------
function val = scalar_fun(lambda_val, mesh_i, sigma_i, theta_i, c_i, Fc_handle, Fd_handle, consts)
    y = lambda_val + consts.mu;
    x = mesh_i;
    h = sqrt(x^2 + y^2); 
    phi_val = atan(y / x);
    
    abs_y = abs(y);
    if abs_y < 1e-7, abs_y = 1e-7; end
    exponent = 0.5 * consts.Nb * (x - 1) * h / abs_y;
    
    arg = exp(exponent);
    if arg > 1, arg = 1; elseif arg < -1, arg = -1; end
    Fp = (2/pi) * acos(arg);
    if Fp < 1e-6, Fp = 1e-6; end
    
    D_i = consts.B_sq_mu_sq + y^2; 
    
    alpha_i = theta_i - phi_val;
    
    Re_i = (h * consts.omegaR * c_i) / consts.nu;
    
    cl_i = Fc_handle(alpha_i, Re_i);
    cd_i = Fd_handle(alpha_i, Re_i);
    
    if (alpha_i < consts.alpha_min) || (alpha_i > consts.alpha_max)
        cl_i = 0; cd_i = consts.cd_max_val;
    end
    
    force_term = sigma_i * h * (cl_i * x - cd_i * y);
    val = 8 * consts.A * Fp * lambda_val * sqrt(D_i) * x - force_term;
end





