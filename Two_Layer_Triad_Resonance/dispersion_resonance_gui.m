function dispersion_resonance_gui
    % Two-Layer Fluid Dispersion: Resonance Analysis & GUI
    % 
    
    clear; clc; close all;
    % =========================================================================
    % 1. CONTROL PANEL (CONFIGURATION)
    % =========================================================================
    
    % --- Initial Values ---
    cfg.k0_init     = 1.5;    % Reference Wavenumber
    cfg.R_init      = 0.85;   % Density Ratio
    cfg.h_init      = 1.0;    % Depth Ratio
    
    cfg.xlim_init   = 8.0;    % X-Axis Limit (+/-)
    cfg.ymax_init   = 4.0;    % Y-Axis Upper Limit
    cfg.ymin_init   = 0.0;   % Y-Axis Lower Limit
    
    cfg.num_points  = 1000;   % Grid resolution
    
    % --- Slider Settings: [Min, Max, Step_Size] ---
    cfg.k0_sets     = [0.1,   5.0,   0.01];   
    cfg.R_sets      = [0.01,  0.99,  0.001]; 
    cfg.h_sets      = [0.1,   5.0,   0.01];  
    
    cfg.xlim_sets   = [2.0,   20.0,  0.1];   
    cfg.ymax_sets   = [0.1,   15.0,  0.1];    
    cfg.ymin_sets   = [-15.0, 0.0,   0.1];    
    
    % --- Visual Styles ---
    cfg.style_orig_surf  = {'Color', 'b', 'LineWidth', 3, 'LineStyle', '-'};
    cfg.style_orig_int   = {'Color', 'r', 'LineWidth', 3, 'LineStyle', '--'};
    cfg.style_shift_surf = {'Color', 'b', 'LineWidth', 0.5, 'LineStyle', '-'};
    cfg.style_shift_int  = {'Color', 'r', 'LineWidth', 1.5, 'LineStyle', ':'}; 
    
    % =========================================================================
    % 2. GUI SETUP
    % =========================================================================
    
    fig = figure('Name', 'Resonance Analysis v7', ...
                 'Color', 'w', 'Units', 'pixels', 'Position', [100, 100, 1200, 800]);
             
    % Control Panel Container
    controlPanel = uipanel('Parent', fig, 'Title', 'Control Panel', ...
                           'FontSize', 12, 'BackgroundColor', 'w', ...
                           'Position', [0.02 0.05 0.25 0.9]);
    
    % --- THE TRICK: Invisible Axes for LaTeX Text ---
    % Standard uicontrols don't support LaTeX. We create an invisible axes
    % spanning the control panel [0,1]x[0,1] to draw text() objects on it.
    lbl_ax = axes('Parent', controlPanel, 'Position', [0 0 1 1], ...
                  'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
                       
    % --- Helper to create Slider + Edit Box + LaTeX Label ---
    function [sld, edt] = create_control_row(parent, label_axes, y_pos, latex_str, settings, init_val, plot_update_func)
        % settings = [min, max, step]
        min_v = settings(1); max_v = settings(2); step_v = settings(3);
        
        % Normalize step
        range_v = max_v - min_v;
        if range_v <= 0, range_v = 1; end
        norm_step = [step_v/range_v, (step_v*10)/range_v];
        
        % 1. LaTeX Label (Using text() on the invisible axes)
        % y_pos is normalized (0 to 1). We offset slightly up for the label.
        % UPDATED: FontSize changed to 12
        text(label_axes, 0.05, y_pos + 0.055, latex_str, ...
             'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold', ...
             'Color', 'k', 'HorizontalAlignment', 'left');
              
        % 2. Edit Box
        edt = uicontrol(parent, 'Style', 'edit', 'String', num2str(init_val), ...
                  'Units', 'normalized', 'Position', [0.70, y_pos+0.032, 0.25, 0.04], ...
                  'BackgroundColor', 'w', 'ForegroundColor', 'k', ...
                  'Callback', {@edit_callback});
              
        % 3. Slider
        sld = uicontrol(parent, 'Style', 'slider', ...
                  'Min', min_v, 'Max', max_v, 'Value', init_val, ...
                  'SliderStep', norm_step, ...
                  'UserData', struct('step', step_v, 'editBox', edt), ... 
                  'Units', 'normalized', 'Position', [0.05, y_pos, 0.9, 0.03], ...
                  'Callback', {@slider_callback});
              
        % --- Local Callbacks ---
        function slider_callback(src, ~)
            data = get(src, 'UserData');
            val = get(src, 'Value');
            val = round(val / data.step) * data.step; 
            set(src, 'Value', val);
            set(data.editBox, 'String', num2str(val));
            plot_update_func();
        end
        function edit_callback(src, ~)
            str = get(src, 'String');
            val = str2double(str);
            if isnan(val), val = get(sld, 'Value'); 
            else
                if val < min_v, val = min_v; end
                if val > max_v, val = max_v; end
            end
            set(sld, 'Value', val);
            set(src, 'String', num2str(val));
            plot_update_func();
        end
    end

    % --- Create Controls ---
    % Now we pass 'lbl_ax' and write full LaTeX strings
    
    [sld_xlim, edt_xlim] = create_control_row(controlPanel, lbl_ax, 0.92, ...
        'Plot Range $X (\pm\kappa)$', cfg.xlim_sets, cfg.xlim_init, @updatePlot);
    
    [sld_ymax, edt_ymax] = create_control_row(controlPanel, lbl_ax, 0.84, ...
        'Y Max ($\Omega$)', cfg.ymax_sets, cfg.ymax_init, @updatePlot);
        
    [sld_ymin, edt_ymin] = create_control_row(controlPanel, lbl_ax, 0.76, ...
        'Y Min ($\Omega$)', cfg.ymin_sets, cfg.ymin_init, @updatePlot);
    
    % Divider
    uicontrol(controlPanel, 'Style', 'frame', 'Units', 'normalized', ...
              'Position', [0.02 0.73 0.96 0.002], 'BackgroundColor', [0.8 0.8 0.8]); 
    
    [sld_k0, edt_k0] = create_control_row(controlPanel, lbl_ax, 0.65, ...
        'Shift Point $(\kappa_0)$', cfg.k0_sets, cfg.k0_init, @updatePlot);
        
    [sld_R,  edt_R]  = create_control_row(controlPanel, lbl_ax, 0.55, ...
        'Density Ratio $(R = \rho_1/\rho_2)$', cfg.R_sets, cfg.R_init, @updatePlot);
        
    [sld_h,  edt_h]  = create_control_row(controlPanel, lbl_ax, 0.45, ...
        'Depth Ratio $(h = h_1/h_2)$', cfg.h_sets, cfg.h_init, @updatePlot);
    
    % Legend Info
    infoStr = sprintf([...
        'Thick Blue: Orig. Surface\n' ...
        'Thick Red:  Orig. Internal\n' ...
        'Thin Blue:  Shifted Surface\n' ...
        'Thin Yell:  Shifted Internal\n' ...
        'Green O:    All Intersections']);
    uicontrol(controlPanel, 'Style', 'text', 'String', infoStr, ...
              'Units', 'normalized', 'Position', [0.05 0.1 0.9 0.25], ...
              'BackgroundColor', 'w', 'HorizontalAlignment', 'left', ...
              'FontAngle', 'italic', 'FontSize', 10);
          
    % =========================================================================
    % 3. PLOTTING AREA
    % =========================================================================
    
    ax = axes('Parent', fig, 'Position', [0.35 0.15 0.6 0.8]);
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    
    xlabel(ax, '$\kappa$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel(ax, '$\Omega$', 'Interpreter', 'latex', 'FontSize', 20);
    set(ax, 'FontSize', 20, 'TickLabelInterpreter', 'latex');
    
    updatePlot();
    
    % =========================================================================
    % 4. UPDATE LOGIC
    % =========================================================================
    function updatePlot(~, ~)
        xl_val   = get(sld_xlim, 'Value');
        ymax_val = get(sld_ymax, 'Value');
        ymin_val = get(sld_ymin, 'Value');
        k0_val   = get(sld_k0, 'Value');
        R_val    = get(sld_R, 'Value');
        h_val    = get(sld_h, 'Value');
        
        if ymin_val >= ymax_val, ymin_val = ymax_val - 0.1; end
        
        % Calculation
        calc_limit = xl_val * 1.1; 
        k_vec = linspace(-calc_limit, calc_limit, cfg.num_points);
        
        [Om_S_orig, Om_I_orig] = solve_dispersion(k_vec, R_val, h_val);
        Orig_Y = {Om_S_orig, -Om_S_orig, Om_I_orig, -Om_I_orig};
        
        [Om_S_0, ~] = solve_dispersion(k0_val, R_val, h_val);
        Omega0 = Om_S_0;
        
        k_shift = k_vec - k0_val;
        [Om_S_shift, Om_I_shift] = solve_dispersion(k_shift, R_val, h_val);
        Shift_Y = { ...
            Omega0 + Om_S_shift, Omega0 - Om_S_shift, ...
            Omega0 + Om_I_shift, Omega0 - Om_I_shift ...
        };
        
        % Plotting
        cla(ax);
        
        % Orig
        plot(ax, k_vec, Orig_Y{1}, cfg.style_orig_surf{:});
        plot(ax, k_vec, Orig_Y{2}, cfg.style_orig_surf{:});
        plot(ax, k_vec, Orig_Y{3}, cfg.style_orig_int{:});
        plot(ax, k_vec, Orig_Y{4}, cfg.style_orig_int{:});
        
        % Shifted
        for i = 1:4
            if i <= 2
                plot(ax, k_vec, Shift_Y{i}, cfg.style_shift_surf{:});
            else
                plot(ax, k_vec, Shift_Y{i}, cfg.style_shift_int{:});
            end
        end
        
        % Point P
        plot(ax, k0_val, Omega0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
        
        % Intersections
        for i = 1:4
            y1 = Orig_Y{i};
            for j = 1:4
                y2 = Shift_Y{j};
                diff_v = y1 - y2;
                sign_chg = (diff_v(1:end-1) .* diff_v(2:end)) <= 0;
                idx_list = find(sign_chg);
                
                for k_idx = 1:length(idx_list)
                    id = idx_list(k_idx);
                    x_a = k_vec(id); x_b = k_vec(id+1);
                    d_a = diff_v(id); d_b = diff_v(id+1);
                    frac = abs(d_a)/(abs(d_a)+abs(d_b));
                    k_int = x_a + frac*(x_b - x_a);
                    
                    val_a = y1(id); val_b = y1(id+1);
                    om_int = val_a + frac*(val_b - val_a);
                    
                    if om_int >= ymin_val && om_int <= ymax_val
                         plot(ax, k_int, om_int, 'go', 'LineWidth', 2, 'MarkerSize', 9);
                    end
                end
            end
        end
        
        xlim(ax, [-xl_val, xl_val]);
        ylim(ax, [ymin_val, ymax_val]);
    end
    
    function [Om_Surf, Om_Int] = solve_dispersion(k_in, R, h)
        k = k_in; 
        k(abs(k) < 1e-9) = 1e-9;
        
        th_k  = tanh(k);
        th_kh = tanh(k * h);
        coth_k  = 1 ./ th_k;
        coth_kh = 1 ./ th_kh;
        
        A = R + (coth_k .* coth_kh);
        B = -k .* (coth_k + coth_kh);
        C = (k.^2) * (1 - R);
        
        discriminant = sqrt(B.^2 - 4 .* A .* C);
        
        Omega2_1 = (-B + discriminant) ./ (2 * A);
        Omega2_2 = (-B - discriminant) ./ (2 * A);
        
        Omega2_1(Omega2_1 < 0) = 0;
        Omega2_2(Omega2_2 < 0) = 0;
        
        Om_Surf = sqrt(Omega2_1);
        Om_Int  = sqrt(Omega2_2);
        
        Om_Surf(abs(k_in) < 1e-8) = 0;
        Om_Int(abs(k_in) < 1e-8)  = 0;
    end
end
