function couette_flow()
close all
   % Parameters
   h_val = 1.0;  % Channel height
   mu_val = 1.0; % Viscosity
   
   % Create Figure
   fig = figure('Name', 'Couette-Poiseuille Flow (Refined Arrows)', ...
                'NumberTitle', 'off', ...
                'Position', [100, 100, 900, 600], ...
                'Color', 'w');

   % Create Axes
   ax = axes('Parent', fig, 'Position', [0.1, 0.35, 0.8, 0.6]);
   hold(ax, 'on');
   grid(ax, 'on');
   title(ax, 'Velocity Profile u(y)');
   xlabel(ax, 'Velocity u (m/s)');
   ylabel(ax, 'Channel Height y (m)');
   
   % Initialize Plot Objects
   % 1. The Velocity Profile Curve
   h_curve = plot(ax, NaN, NaN, 'b-', 'LineWidth', 2);
   
   % 2. The Arrow Stems (Lines)
   h_stems = plot(ax, NaN, NaN, 'r-', 'LineWidth', 1.5);
   
   % 3. The Arrow Heads (Filled Triangles)
   h_heads = patch(ax, 'XData', [], 'YData', [], 'FaceColor', 'r', 'EdgeColor', 'none');
   
   % Draw plates
   yline(ax, 0, 'k-', 'LineWidth', 3);      % Bottom
   yline(ax, h_val, 'k-', 'LineWidth', 3);  % Top
   
   % Fixed Axes Limits
   xlim(ax, [-15, 15]);
   ylim(ax, [0, h_val]);

   % --- UI Controls (Dual Input) ---

   % 1. Top Plate Velocity (U)
   uicontrol('Style', 'text', 'Parent', fig, ...
       'Position', [50, 130, 150, 20], ...
       'String', 'Top Plate Velocity (U):', ...
       'HorizontalAlignment', 'left', 'BackgroundColor', 'w','fontsize',14);
       
   sld_U = uicontrol('Style', 'slider', 'Parent', fig, ...
       'Position', [200, 130, 400, 20], ...
       'Min', -15, 'Max', 15, 'Value', 0);
       
   txt_U = uicontrol('Style', 'edit', 'Parent', fig, ...
       'Position', [620, 130, 60, 20], ...
       'String', '0', 'BackgroundColor', 'w','fontsize',14);

   % 2. Pressure Gradient (dP/dx)
   uicontrol('Style', 'text', 'Parent', fig, ...
       'Position', [50, 80, 150, 20], ...
       'String', 'Pressure Grad (dP/dx):', ...
       'HorizontalAlignment', 'left', 'BackgroundColor', 'w','fontsize',14);
       
   sld_P = uicontrol('Style', 'slider', 'Parent', fig, ...
       'Position', [200, 80, 400, 20], ...
       'Min', -120, 'Max', 120, 'Value', 0);
       
   txt_P = uicontrol('Style', 'edit', 'Parent', fig, ...
       'Position', [620, 80, 60, 20], ...
       'String', '0', 'BackgroundColor', 'w','fontsize',14);

   % --- Callbacks ---

   addlistener(sld_U, 'Value', 'PostSet', @(src, event) sync_inputs('slider_U'));
   addlistener(sld_P, 'Value', 'PostSet', @(src, event) sync_inputs('slider_P'));
   set(txt_U, 'Callback', @(src, event) sync_inputs('text_U'));
   set(txt_P, 'Callback', @(src, event) sync_inputs('text_P'));

   function sync_inputs(source)
       U_val = get(sld_U, 'Value');
       P_val = get(sld_P, 'Value');

       switch source
           case 'slider_U'
               U_val = get(sld_U, 'Value');
               set(txt_U, 'String', num2str(U_val, '%.2f'));
           case 'slider_P'
               P_val = get(sld_P, 'Value');
               set(txt_P, 'String', num2str(P_val, '%.2f'));
           case 'text_U'
               val = str2double(get(txt_U, 'String'));
               val = max(min(val, 10), -10); 
               U_val = val;
               set(sld_U, 'Value', val);
               set(txt_U, 'String', num2str(val, '%.2f'));
           case 'text_P'
               val = str2double(get(txt_P, 'String'));
               val = max(min(val, 20), -20);
               P_val = val;
               set(sld_P, 'Value', val);
               set(txt_P, 'String', num2str(val, '%.2f'));
       end
       update_plot(U_val, P_val);
   end

   function update_plot(U, dPdx)
       y = linspace(0, h_val, 100);
       term1 = (y ./ h_val) .* U;
       term2 = (h_val^2 / (2 * mu_val)) * (-dPdx) .* (y ./ h_val) .* (1 - (y ./ h_val));
       u = term1 + term2;
       
       set(h_curve, 'XData', u, 'YData', y);
       
       % --- Custom Arrow Drawing ---
       
       % 1. Select points for arrows (Decimate)
       % Using every 6th point creates a nice vertical spacing
       idx = 1:6:length(y); 
       y_q = y(idx);
       u_q = u(idx);
       
       % 2. Draw Stems
       x_stems = zeros(3*length(idx), 1);
       y_stems = zeros(3*length(idx), 1);
       
       % [0, u, NaN, 0, u, NaN...]
       x_stems(1:3:end) = 0;
       x_stems(2:3:end) = u_q;
       x_stems(3:3:end) = NaN;
       
       y_stems(1:3:end) = y_q;
       y_stems(2:3:end) = y_q;
       y_stems(3:3:end) = NaN;
       
       set(h_stems, 'XData', x_stems, 'YData', y_stems);
       
       % 3. Draw Filled Heads
       % Adjusted Dimensions for "Nice Small Size"
       head_len = 0.4;   % Horizontal length of head
       head_wid = 0.015; % Vertical half-width (Total height = 0.03)
                         % Since spacing is approx 0.06, this leaves 50% gap.
       
       dirs = sign(u_q);
       dirs(dirs == 0) = 1;
       
       % X Coordinates (Tip, Top Base, Bottom Base)
       % Note: If velocity is smaller than head_len, we scale the head down
       % to prevent the arrow pointing backwards past zero.
       actual_len = min(head_len, abs(u_q)); 
       
       X_tips = u_q;
       X_base = u_q - (dirs .* actual_len);
       X_patch = [X_tips; X_base; X_base];
       
       % Y Coordinates
       Y_tips = y_q;
       Y_top  = y_q + head_wid;
       Y_bot  = y_q - head_wid;
       Y_patch = [Y_tips; Y_top; Y_bot];
       
       % Hide heads for very small velocities to avoid clutter
       mask = abs(u_q) < 0.05;
       X_patch(:, mask) = NaN;
       Y_patch(:, mask) = NaN;
       
       set(h_heads, 'XData', X_patch, 'YData', Y_patch);
   end

   % Initial Draw
   sync_inputs('slider_U');
end
