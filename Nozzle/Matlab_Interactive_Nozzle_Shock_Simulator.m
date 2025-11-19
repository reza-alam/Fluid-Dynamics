function nozzle_shock_sim
   % NOZZLE_SHOCK_SIM
   % A GUI to visualize compressible flow, normal shocks, and velocity vectors
   % in a Converging-Diverging (C-D) nozzle.
   
   % CLEANUP
   close all; 
   clc;

   % =========================================================
   %               USER CONFIGURATION PARAMETERS
   % =========================================================
   N_ARROW_COLS = 21;      % Number of columns of arrows (Keep odd for symmetry)
   ARROW_SCALE  = 0.00025/2; % Length scale of arrows (Decrease to make shorter)
   % =========================================================

   % ---------------------------------------------------------
   % 1. SYSTEM PARAMETERS & GEOMETRY SETUP
   % ---------------------------------------------------------
   gamma = 1.4;            % Ratio of specific heats (Air)
   R = 287;                % Gas constant (J/kg*K)
   T0 = 300;               % Stagnation Temperature (K)
   P0 = 101325;            % Inlet Stagnation Pressure (Pa)
   
   % Define Nozzle Geometry
   % Using 101 points to ensure exact center at index 51
   x_len = 101;            
   x = linspace(0, 1, x_len); 
   
   % Geometry Parameters
   A_inlet = 0.25;
   A_throat = 0.1;
   A_exit  = 0.25;
   
   % Mathematical shape function (Cosine based for single throat at x=0.5)
   A = zeros(1, x_len);
   for i = 1:x_len
       if x(i) < 0.5
           % Converging Section
           xi = x(i) / 0.5;
           shape_factor = 0.5 * (1 + cos(pi * xi));
           A(i) = A_throat + (A_inlet - A_throat) * shape_factor;
       else
           % Diverging Section
           xi = (x(i) - 0.5) / 0.5;
           shape_factor = 0.5 * (1 - cos(pi * xi));
           A(i) = A_throat + (A_exit - A_throat) * shape_factor;
       end
   end
   
   % Calculate Radius
   Y_wall = sqrt(A / pi);

   % ---------------------------------------------------------
   % 2. PRE-CALCULATE VECTOR GRID (SYMMETRIC STATIONS)
   % ---------------------------------------------------------
   % Generate evenly spaced indices based on User Parameter
   viz_indices = round(linspace(1, x_len, N_ARROW_COLS)); 
   
   num_y_points = 9; % Number of vertical arrows per column
   X_field = [];
   Y_field = [];
   Map_Indices = []; % Maps 2D point back to 1D solution index
   
   for k = 1:length(viz_indices)
       idx = viz_indices(k);
       
       % Spread arrows vertically at this station
       % Use 0.85 factor to keep them well inside the wall
       y_local = linspace(-Y_wall(idx)*0.85, Y_wall(idx)*0.85, num_y_points);
       x_local = repmat(x(idx), 1, num_y_points);
       
       X_field = [X_field, x_local];
       Y_field = [Y_field, y_local];
       Map_Indices = [Map_Indices, repmat(idx, 1, num_y_points)];
   end

   % ---------------------------------------------------------
   % 3. PRE-CALCULATION OF ISENTROPIC BRANCHES
   % ---------------------------------------------------------
   M_sub_iso = zeros(1, x_len);
   M_sup_iso = zeros(1, x_len);
   options = optimset('Display','off', 'TolX', 1e-6);
   
   for i = 1:x_len
       area_ratio = A(i) / A_throat;
       func = @(M) (1/M^2) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * M^2))^((gamma+1)/(gamma-1)) - area_ratio^2;
       try, M_sub_iso(i) = fzero(func, [1e-4, 1], options); catch, M_sub_iso(i) = 1; end
       try, M_sup_iso(i) = fzero(func, [1, 10], options); catch, M_sup_iso(i) = 1; end
   end
   
   % Critical Pressures
   [P_ratio_sub_exit, ~] = get_isentropic_ratios(M_sub_iso(end), gamma);
   M_before_shock = M_sup_iso(end);
   M_after_shock = get_normal_shock_M2(M_before_shock, gamma);
   P0_ratio_shock = get_stagnation_pressure_ratio(M_before_shock, gamma);
   [p_static_ratio, ~] = get_isentropic_ratios(M_after_shock, gamma);
   P_ratio_shock_exit = p_static_ratio * P0_ratio_shock; 

   % ---------------------------------------------------------
   % 4. GUI SETUP
   % ---------------------------------------------------------
   hFig = figure('Name', 'Nozzle Flow Simulator', ...
                 'NumberTitle', 'off', ...
                 'Position', [50, 50, 900, 800], ...
                 'Color', 'w');

   % --- AXES 1: Nozzle Geometry (TOP) ---
   ax1 = axes('Parent', hFig, 'Position', [0.1, 0.70, 0.8, 0.25]);
   hold(ax1, 'on');
   axis(ax1, 'normal'); 
   xlim(ax1, [0, 1]);
   ylim(ax1, [-0.3, 0.3]); 
   title(ax1, 'Velocity Field (Color = Magnitude)');
   ylabel(ax1, 'Radius (m)');
   set(ax1, 'XTickLabel', []); 
   
   % 1. Create Patch Object for Color Background
   % We draw a polygon filling the nozzle.
   % Vertices order: Top wall (Left->Right) then Bottom wall (Right->Left)
   patch_x = [x, fliplr(x)];
   patch_y = [Y_wall, fliplr(-Y_wall)];
   
   % Initial Color Data (Will be updated)
   % We need color data for each vertex.
   patch_c = zeros(size(patch_x)); 
   
   hPatch = patch(ax1, patch_x, patch_y, patch_c, 'EdgeColor', 'k', 'LineWidth', 2);
   colormap(ax1, jet); % or parula
   caxis(ax1, [0, 500]); % Fixed velocity range [0, 500 m/s]
   %c = colorbar(ax1);
   %c.Label.String = 'Velocity (m/s)';
   
   % 2. Vector Plot (Black Arrows)
   zeros_vec = zeros(size(X_field));
   hQuiver = quiver(ax1, X_field, Y_field, zeros_vec, zeros_vec, ...
       0, 'k', 'LineWidth', 1.2, 'MaxHeadSize', 0.8, 'AutoScale', 'off'); 
   
   % --- AXES 2: Pressure Distribution (MIDDLE) ---
   ax2 = axes('Parent', hFig, 'Position', [0.1, 0.42, 0.8, 0.20]);
   hold(ax2, 'on'); grid(ax2, 'on');
   xlim(ax2, [0, 1]);
   ylim(ax2, [0, 1.1]);
   ylabel(ax2, 'Pressure Ratio (P/P_0)');
   set(ax2, 'XTickLabel', []); 
   
   hPressLine = plot(ax2, x, ones(1, x_len), 'r-', 'LineWidth', 2);
   yline(ax2, P_ratio_sub_exit, 'g--', 'Subsonic Limit');
   yline(ax2, P_ratio_shock_exit, 'm--', 'Shock Exit Limit');

   % --- AXES 3: Mach Number (BOTTOM) ---
   ax3 = axes('Parent', hFig, 'Position', [0.1, 0.15, 0.8, 0.20]);
   hold(ax3, 'on'); grid(ax3, 'on');
   xlim(ax3, [0, 1]);
   ylim(ax3, [0, 3.0]);
   ylabel(ax3, 'Mach Number');
   xlabel(ax3, 'Position x (m)');
   
   hMachLine = plot(ax3, x, zeros(1, x_len), 'b-', 'LineWidth', 2);
   yline(ax3, 1, 'k--', 'Sonic (M=1)');

   % --- Slider Control ---
   uicontrol('Style', 'text', 'Position', [300, 60, 400, 20], ...
       'String', 'Adjust Back Pressure (Pb/P0):', ...
       'BackgroundColor', 'w', 'FontSize', 10);
   
   slider = uicontrol('Style', 'slider', ...
       'Min', 0.425, 'Max', P_ratio_sub_exit, 'Value', P_ratio_sub_exit, ...
       'Position', [200, 30, 600, 20], ...
       'SliderStep', [0.01, 0.04], ...
       'Callback', @update_flow);
       
   % Initial Update
   update_flow(slider, []);

   % ---------------------------------------------------------
   % 5. CORE SOLVER
   % ---------------------------------------------------------
   function update_flow(src, ~)
       Pb_ratio = get(src, 'Value');
       M_dist = zeros(1, x_len);
       P_dist = zeros(1, x_len);
       mode_str = '';

       % -- SOLVER LOGIC --
       if Pb_ratio >= (P_ratio_sub_exit - 1e-4)
           for k = 1:x_len
               M_dist(k) = M_sub_iso(k); 
               [pr, ~] = get_isentropic_ratios(M_dist(k), gamma);
               P_dist(k) = pr; 
           end
           mode_str = 'Subsonic Isentropic';
           
       elseif Pb_ratio < P_ratio_sub_exit && Pb_ratio > P_ratio_shock_exit
           best_err = 100; shock_idx = x_len;
           start_node = floor(x_len/2);
           for k = start_node : x_len
               M1 = M_sup_iso(k);
               P0_loss = get_stagnation_pressure_ratio(M1, gamma);
               A_star_2 = A_throat * (1/P0_loss);
               AR_exit_2 = A_exit / A_star_2;
               func_sub = @(m) (1/m^2) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * m^2))^((gamma+1)/(gamma-1)) - AR_exit_2^2;
               try, M_ex_calc = fzero(func_sub, [1e-5, 1], options); catch, M_ex_calc = 0.1; end
               [pr_static, ~] = get_isentropic_ratios(M_ex_calc, gamma);
               P_exit_calc = pr_static * P0_loss;
               if abs(P_exit_calc - Pb_ratio) < best_err
                   best_err = abs(P_exit_calc - Pb_ratio);
                   shock_idx = k;
               end
           end
           for k = 1:x_len
               if k < shock_idx
                   if k <= x_len/2, M_dist(k) = M_sub_iso(k); else, M_dist(k) = M_sup_iso(k); end
                   [pr, ~] = get_isentropic_ratios(M_dist(k), gamma);
                   P_dist(k) = pr;
               elseif k == shock_idx
                   M_dist(k) = M_sup_iso(k); 
                   [pr, ~] = get_isentropic_ratios(M_dist(k), gamma);
                   P_dist(k) = pr;
               else
                   M1_shock = M_sup_iso(shock_idx);
                   P0_loss = get_stagnation_pressure_ratio(M1_shock, gamma);
                   A_star_2 = A_throat * (1/P0_loss);
                   AR_local = A(k) / A_star_2;
                   func_sub = @(m) (1/m^2) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * m^2))^((gamma+1)/(gamma-1)) - AR_local^2;
                   try, M_local = fzero(func_sub, [1e-5, 1], options); catch, M_local = 0.1; end
                   M_dist(k) = M_local;
                   [pr_static, ~] = get_isentropic_ratios(M_local, gamma);
                   P_dist(k) = pr_static * P0_loss;
               end
           end
           mode_str = 'Normal Shock';
           
       else
           for k = 1:x_len
               if k <= x_len/2, M_dist(k) = M_sub_iso(k); else, M_dist(k) = M_sup_iso(k); end
               [pr, ~] = get_isentropic_ratios(M_dist(k), gamma);
               P_dist(k) = pr;
           end
           mode_str = 'Fully Supersonic';
       end
       
       % -- UPDATE VISUALS --
       
       T_dist = T0 ./ (1 + (gamma-1)/2 .* M_dist.^2);
       a_dist = sqrt(gamma * R .* T_dist);
       V_dist = M_dist .* a_dist;
       
       % 1. Update Patch Color (Magnitude Background)
       % Vertices are Top [1..101] then Bottom [101..1]
       % V_dist is 1x101. We assume uniform velocity across Y.
       % Map colors to vertices:
       patch_c_new = [V_dist, fliplr(V_dist)];
       set(hPatch, 'CData', patch_c_new);
       
       % 2. Update Quiver (Direction & Magnitude by Length)
       U_field = V_dist(Map_Indices) * ARROW_SCALE;
       set(hQuiver, 'UData', U_field, 'VData', zeros(size(U_field)));
       
       % 3. Plots
       set(hPressLine, 'YData', P_dist);
       title(ax2, sprintf('Pressure Ratio (Regime: %s)', mode_str));
       set(hMachLine, 'YData', M_dist);
   end
end

% ---------------------------------------------------------
% HELPER FUNCTIONS
% ---------------------------------------------------------
function [P_ratio, T_ratio] = get_isentropic_ratios(M, gamma)
   factor = 1 + (gamma-1)/2 * M^2;
   T_ratio = 1 / factor;
   P_ratio = (1 / factor)^(gamma/(gamma-1));
end

function M2 = get_normal_shock_M2(M1, gamma)
   num = 1 + (gamma-1)/2 * M1^2;
   den = gamma * M1^2 - (gamma-1)/2;
   M2 = sqrt(num / den);
end

function P0_ratio = get_stagnation_pressure_ratio(M1, gamma)
   term1 = ((gamma+1)*M1^2 / (2 + (gamma-1)*M1^2)) ^ (gamma/(gamma-1));
   term2 = ((gamma+1) / (2*gamma*M1^2 - (gamma-1))) ^ (1/(gamma-1));
   P0_ratio = term1 * term2;
end
