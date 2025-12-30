%% Interactive Capillary-Gravity Wave Analysis
% Visualizes the dispersion relation and resonance conditions.
% ALL text elements are rendered using the LaTeX interpreter.

clc; clear; close all;

%% 1. USER CONFIGURATION
% ---------------------------------------------------------
Bo = 100;                % Bond Number
kappa_max = 25.0;        % Max range for the plot x-axis

% Slider Configuration for Kappa_0
k0_min = 0;              
k0_max = 25;           
k0_step = 0.1;           

num_points = 4000;       % Resolution
% ---------------------------------------------------------

k_grid = linspace(0, kappa_max, num_points);

%% 2. Setup Figure and UI
fig = figure('Name', 'Interactive Dispersion Relation', 'Color', 'w', ...
    'Position', [100, 100, 900, 600]);

% Create Axes
ax = axes('Parent', fig, 'Position', [0.12, 0.25, 0.8, 0.65]); % Adjusted margins
hold(ax, 'on');
box(ax, 'on');
grid(ax, 'on');

% --- LATEX FORMATTING FOR AXES ---
set(ax, 'FontSize', 18); 
set(ax, 'TickLabelInterpreter', 'latex'); % Makes the numbers LaTeX font
xlabel(ax, '$\kappa$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel(ax, '$\Omega$', 'Interpreter', 'latex', 'FontSize', 18);
%title(ax, sprintf('\\textbf{Dispersion Relation} ($\\mathrm{Bo} = %.1f$)', Bo),'Interpreter', 'latex', 'FontSize', 18);

% --- Plot Objects ---
% 1. Base Curve
h_base = plot(ax, NaN, NaN, 'b-', 'LineWidth', 2, ...
    'DisplayName', '$\Omega(\kappa)$');

% 2. Shifted Curve
h_shift = plot(ax, NaN, NaN, 'r--', 'LineWidth', 2, ...
    'DisplayName', '$\Omega(\kappa - \kappa_0) + \Omega_0$');

% 3. Marker for Kappa0
h_point = plot(ax, NaN, NaN, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, ...
    'HandleVisibility', 'off');

% 4. Intersection Marker
h_intersect = plot(ax, NaN, NaN, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10, ...
    'LineWidth', 2, 'DisplayName', '$\mathrm{Intersection}$');

% Legend with LaTeX
lgd = legend(ax, 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 18);

% --- Slider UI ---
% Label "Adjust Kappa 0"
uicontrol('Parent', fig, 'Style', 'text', ...
    'Position', [150, 65, 600, 20], ...
    'String', 'Adjust Source Wavenumber:', 'BackgroundColor', 'w', ...
    'FontSize', 12, 'FontName', 'Times New Roman'); % Standard font for UI instructions

% Slider Control
slider_range = k0_max - k0_min;
normalized_step = k0_step / slider_range;

h_slider = uicontrol('Parent', fig, 'Style', 'slider', ...
    'Min', k0_min, 'Max', k0_max, 'Value', k0_min, ...
    'SliderStep', [normalized_step, normalized_step*5], ...
    'Position', [150, 35, 600, 20]);

% --- DYNAMIC LATEX LABEL ---
% We use an annotation instead of uicontrol so we can use LaTeX interpreter
dim = [0.82, 0.04, 0.15, 0.1]; % [x y w h]
h_val_label = annotation('textbox', dim, 'String', '', ...
    'Interpreter', 'latex', 'FontSize', 16, ...
    'EdgeColor', 'none', 'FitBoxToText', 'on');

%% 3. Pre-calculate Base Curve
Omega_base = get_omega(k_grid, Bo);
set(h_base, 'XData', k_grid, 'YData', Omega_base);

%% 4. Link Callback and Initialize
h_slider.Callback = @(src, event) update_plot(src, h_val_label, h_shift, h_point, h_intersect, k_grid, Omega_base, Bo);

% Trigger initialization
update_plot(h_slider, h_val_label, h_shift, h_point, h_intersect, k_grid, Omega_base, Bo);

%% 5. Helper Functions

function update_plot(slider, label_annot, h_shift, h_point, h_intersect, k_grid, Omega_base, Bo)
    % Get slider value
    k0 = get(slider, 'Value');
    
    % Update the LaTeX annotation
    % We use \kappa_0 to look math-like
    set(label_annot, 'String', sprintf('$\\kappa_0 = %.2f$', k0));
    
    % Calculate Omega 0
    w0 = get_omega(k0, Bo);
    
    % Update Red Point
    set(h_point, 'XData', k0, 'YData', w0);
    
    % --- Calculate Shifted Curve ---
    valid_indices = k_grid >= k0;
    x_local = k_grid(valid_indices); 
    k_input = x_local - k0; 
    
    Omega_local = get_omega(k_input, Bo);
    Omega_shifted = Omega_local + w0;
    
    set(h_shift, 'XData', x_local, 'YData', Omega_shifted);
    
    % --- Find Intersection ---
    y_base_segment = Omega_base(valid_indices);
    diff_curves = y_base_segment - Omega_shifted;
    
    % Zero crossing detection
    idx_cross = find(diff_curves(1:end-1) .* diff_curves(2:end) <= 0);
    
    % Ignore trivial intersection at start point
    idx_cross = idx_cross(x_local(idx_cross) > k0 + 0.1); 
    
    if ~isempty(idx_cross)
        i_int = idx_cross(1); 
        x1 = x_local(i_int); x2 = x_local(i_int+1);
        y1 = diff_curves(i_int); y2 = diff_curves(i_int+1);
        fraction = abs(y1) / (abs(y1) + abs(y2));
        k_int_exact = x1 + fraction * (x2 - x1);
        
        w_int_exact = interp1(k_grid, Omega_base, k_int_exact);
        
        set(h_intersect, 'XData', k_int_exact, 'YData', w_int_exact, 'Visible', 'on');
    else
        set(h_intersect, 'Visible', 'off');
    end
end

function w = get_omega(k, Bo)
    % Dispersion Relation: Omega^2 = [k + (k^3 / Bo)] * tanh(k)
    term_gravity = k;
    term_capillary = (k.^3) ./ Bo;
    term_depth = tanh(k);
    w_sq = (term_gravity + term_capillary) .* term_depth;
    w = sqrt(w_sq);
end