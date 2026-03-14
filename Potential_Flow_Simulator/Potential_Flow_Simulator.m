function PotentialFlowSimulator()
% PotentialFlowSimulator - Interactive 2D potential flow visualization
%
% Superpose fundamental solutions:
%   1) Source/Sink:  u = (m/2pi)(x-x0)/r^2,  v = (m/2pi)(y-y0)/r^2
%      m > 0: source,  m < 0: sink
%   2) Vortex:  u = -(Gamma/2pi)(y-y0)/r^2,  v = (Gamma/2pi)(x-x0)/r^2
%      Gamma > 0: counterclockwise
%   3) Doublet (x-oriented):
%      u = -(kappa/2pi)(y^2-x^2)/r^4,  v = (kappa/2pi)(2xy)/r^4
%   4) Uniform flow:  u = U*cos(alpha),  v = U*sin(alpha)
%
% Controls:
%   - Add/remove elements via buttons
%   - Select element in listbox to edit properties
%   - Strength: text field + slider
%   - Location: (x, y) text fields
%   - Angle: for uniform flow only (slider + text)
%   - Toggle velocity vectors and streamlines independently
%   - Grid density slider

   %% ====== DATA STORAGE ======
   elements = struct('type',{},'strength',{},'x',{},'y',{},'angle',{});
   selectedIdx = 0;

   %% ====== MAIN FIGURE ======
   fig = figure('Name','2D Potential Flow Simulator',...
       'Position',[50 50 1350 750],...
       'MenuBar','none','NumberTitle','off',...
       'Color',[0.94 0.94 0.94],...
       'CloseRequestFcn',@(src,~)delete(src));

   %% ====== LEFT PANEL: CONTROLS ======
   leftPanel = uipanel(fig,'Title','  Flow Elements  ',...
       'Units','pixels','Position',[10 10 370 730],...
       'FontSize',12,'FontWeight','bold');

   % --- Element listbox ---
   uicontrol(leftPanel,'Style','text','String','Current Elements:',...
       'Units','pixels','Position',[10 680 200 18],...
       'HorizontalAlignment','left','FontSize',10,...
       'BackgroundColor',get(leftPanel,'BackgroundColor'));

   hList = uicontrol(leftPanel,'Style','listbox',...
       'Units','pixels','Position',[10 500 345 180],...
       'FontSize',10,'Callback',@onSelectElement);

   % --- Add / Remove buttons ---
   uicontrol(leftPanel,'Style','pushbutton','String','+ Add Element',...
       'Units','pixels','Position',[10 460 168 32],...
       'FontSize',10,'FontWeight','bold',...
       'BackgroundColor',[0.6 0.85 0.6],...
       'Callback',@onAdd);
   uicontrol(leftPanel,'Style','pushbutton','String','- Remove Selected',...
       'Units','pixels','Position',[187 460 168 32],...
       'FontSize',10,...
       'BackgroundColor',[0.95 0.65 0.65],...
       'Callback',@onRemove);

   % --- Properties panel ---
   propPanel = uipanel(leftPanel,'Title','  Next / Selected Element Properties  ',...
       'Units','pixels','Position',[10 10 345 440],...
       'FontSize',10,'FontWeight','bold');

   ypos = 390; dy = 45;

   % Type dropdown
   uicontrol(propPanel,'Style','text','String','Type:',...
       'Units','pixels','Position',[10 ypos 60 20],...
       'HorizontalAlignment','left','FontSize',10,...
       'BackgroundColor',get(propPanel,'BackgroundColor'));
   hType = uicontrol(propPanel,'Style','popupmenu',...
       'String',{'Source/Sink','Vortex','Doublet','Uniform Flow'},...
       'Units','pixels','Position',[80 ypos-2 240 26],...
       'FontSize',10,'Callback',@onTypeChange);
   ypos = ypos - dy;

   % Strength
   uicontrol(propPanel,'Style','text','String','Strength:',...
       'Units','pixels','Position',[10 ypos 70 20],...
       'HorizontalAlignment','left','FontSize',10,...
       'BackgroundColor',get(propPanel,'BackgroundColor'));
   hStrEdit = uicontrol(propPanel,'Style','edit','String','1.00',...
       'Units','pixels','Position',[80 ypos-2 80 26],...
       'FontSize',10,'Callback',@onStrEditChange);
   hStrSlider = uicontrol(propPanel,'Style','slider',...
       'Units','pixels','Position',[170 ypos 150 22],...
       'Min',-20,'Max',20,'Value',1,...
       'Callback',@onStrSliderChange);
   % Try continuous slider update (works R2014b+)
   try
       addlistener(hStrSlider,'ContinuousValueChange',@onStrSliderChange);
   catch
   end
   ypos = ypos - dy;

   % X location
   uicontrol(propPanel,'Style','text','String','X pos:',...
       'Units','pixels','Position',[10 ypos 50 20],...
       'HorizontalAlignment','left','FontSize',10,...
       'BackgroundColor',get(propPanel,'BackgroundColor'));
   hXEdit = uicontrol(propPanel,'Style','edit','String','0.00',...
       'Units','pixels','Position',[60 ypos-2 80 26],...
       'FontSize',10,'Callback',@onLocChange);

   % Y location
   uicontrol(propPanel,'Style','text','String','Y pos:',...
       'Units','pixels','Position',[160 ypos 50 20],...
       'HorizontalAlignment','left','FontSize',10,...
       'BackgroundColor',get(propPanel,'BackgroundColor'));
   hYEdit = uicontrol(propPanel,'Style','edit','String','0.00',...
       'Units','pixels','Position',[210 ypos-2 80 26],...
       'FontSize',10,'Callback',@onLocChange);
   ypos = ypos - dy;

   % Angle (visible only for uniform flow)
   lblAngle = uicontrol(propPanel,'Style','text','String','Angle (deg):',...
       'Units','pixels','Position',[10 ypos 90 20],...
       'HorizontalAlignment','left','FontSize',10,...
       'BackgroundColor',get(propPanel,'BackgroundColor'),'Visible','off');
   hAngleEdit = uicontrol(propPanel,'Style','edit','String','0.0',...
       'Units','pixels','Position',[100 ypos-2 60 26],...
       'FontSize',10,'Callback',@onLocChange,'Visible','off');
   hAngleSlider = uicontrol(propPanel,'Style','slider',...
       'Units','pixels','Position',[170 ypos 150 22],...
       'Min',-180,'Max',180,'Value',0,...
       'Callback',@onAngleSliderChange,'Visible','off');
   try
       addlistener(hAngleSlider,'ContinuousValueChange',@onAngleSliderChange);
   catch
   end
   ypos = ypos - dy - 10;

   % Info text
   infoStr = {
       'Workflow:'
       '  1. Set type/strength/location below'
       '  2. Click "+ Add Element"'
       '  3. Click element in list to edit it'
       ''
       'Conventions:'
       '  Source/Sink: m>0 source, m<0 sink'
       '  Vortex: Gamma>0 counterclockwise'
       '  Doublet: x-oriented'
       '  Uniform: speed + angle'
   };
   uicontrol(propPanel,'Style','text','String',infoStr,...
       'Units','pixels','Position',[10 20 310 ypos-20],...
       'HorizontalAlignment','left','FontSize',9,...
       'BackgroundColor',get(propPanel,'BackgroundColor'));

   %% ====== RIGHT PANEL: PLOT ======
   rightPanel = uipanel(fig,'Title','  Flow Visualization  ',...
       'Units','pixels','Position',[390 10 945 730],...
       'FontSize',12,'FontWeight','bold');

   ax = axes(rightPanel,'Units','normalized',...
       'Position',[0.08 0.12 0.88 0.82]);

   % Toggle checkboxes
   hCBVec = uicontrol(rightPanel,'Style','checkbox',...
       'String','Velocity Vectors','Value',1,...
       'Units','normalized','Position',[0.02 0.065 0.2 0.035],...
       'FontSize',10,'Callback',@(~,~)updatePlot());

   hCBStream = uicontrol(rightPanel,'Style','checkbox',...
       'String','Streamlines','Value',0,...
       'Units','normalized','Position',[0.22 0.065 0.18 0.035],...
       'FontSize',10,'Callback',@(~,~)updatePlot());

   % Grid density
   uicontrol(rightPanel,'Style','text','String','Grid:',...
       'Units','normalized','Position',[0.02 0.025 0.05 0.03],...
       'HorizontalAlignment','right','FontSize',9,...
       'BackgroundColor',get(rightPanel,'BackgroundColor'));
   hGridSlider = uicontrol(rightPanel,'Style','slider',...
       'Units','normalized','Position',[0.08 0.03 0.12 0.025],...
       'Min',10,'Max',60,'Value',30,...
       'Callback',@(~,~)updatePlot());
   try
       addlistener(hGridSlider,'ContinuousValueChange',@(~,~)updatePlot());
   catch
   end

   % Vector scale
   uicontrol(rightPanel,'Style','text','String','Vec scale:',...
       'Units','normalized','Position',[0.38 0.025 0.1 0.03],...
       'HorizontalAlignment','right','FontSize',9,...
       'BackgroundColor',get(rightPanel,'BackgroundColor'));
   hVecScale = uicontrol(rightPanel,'Style','slider',...
       'Units','normalized','Position',[0.49 0.03 0.14 0.025],...
       'Min',0.05,'Max',3.0,'Value',1.0,...
       'Callback',@(~,~)updatePlot());
   try
       addlistener(hVecScale,'ContinuousValueChange',@(~,~)updatePlot());
   catch
   end

   % Streamline density
   uicontrol(rightPanel,'Style','text','String','Stream dens:',...
       'Units','normalized','Position',[0.65 0.025 0.12 0.03],...
       'HorizontalAlignment','right','FontSize',9,...
       'BackgroundColor',get(rightPanel,'BackgroundColor'));
   hStreamDens = uicontrol(rightPanel,'Style','slider',...
       'Units','normalized','Position',[0.78 0.03 0.14 0.025],...
       'Min',0.5,'Max',5,'Value',2,...
       'Callback',@(~,~)updatePlot());
   try
       addlistener(hStreamDens,'ContinuousValueChange',@(~,~)updatePlot());
   catch
   end

   % Initial empty plot
   updatePlot();

   %% ====== CALLBACKS ======
   function onAdd(~,~)
       % Read current control values to create the new element
       n = length(elements) + 1;
       elements(n).type     = get(hType,'Value');
       elements(n).strength = str2double(get(hStrEdit,'String'));
       if isnan(elements(n).strength), elements(n).strength = 1; end
       xv = str2double(get(hXEdit,'String'));
       yv = str2double(get(hYEdit,'String'));
       elements(n).x = xv; if isnan(xv), elements(n).x = 0; end
       elements(n).y = yv; if isnan(yv), elements(n).y = 0; end
       av = str2double(get(hAngleEdit,'String'));
       elements(n).angle = av; if isnan(av), elements(n).angle = 0; end
       % Deselect so controls are free for staging next element
       selectedIdx = 0;
       refreshList();
       resetControls();
       updatePlot();
   end

   function onRemove(~,~)
       if selectedIdx < 1 || selectedIdx > length(elements), return; end
       elements(selectedIdx) = [];
       selectedIdx = 0;
       refreshList();
       resetControls();
       updatePlot();
   end

   function onSelectElement(~,~)
       idx = get(hList,'Value');
       if selectedIdx == 0
           % Hint row is at position 1, so actual elements start at 2
           actualIdx = idx - 1;
           if actualIdx >= 1 && actualIdx <= length(elements)
               selectedIdx = actualIdx;
               loadProperties(selectedIdx);
               refreshList();
           end
       else
           % No hint row, indices match directly
           if idx > 0 && idx <= length(elements)
               if idx == selectedIdx
                   % Click same item again = deselect
                   selectedIdx = 0;
                   refreshList();
                   resetControls();
               else
                   selectedIdx = idx;
                   loadProperties(idx);
                   refreshList();
               end
           end
       end
   end

   function onTypeChange(~,~)
       toggleAngleControls(get(hType,'Value') == 4);
       if selectedIdx < 1 || selectedIdx > length(elements), return; end
       % Editing existing element
       elements(selectedIdx).type = get(hType,'Value');
       refreshList();
       updatePlot();
   end

   function onStrEditChange(~,~)
       val = str2double(get(hStrEdit,'String'));
       if isnan(val), return; end
       set(hStrSlider,'Value', max(-20, min(20, val)));
       if selectedIdx < 1 || selectedIdx > length(elements), return; end
       elements(selectedIdx).strength = val;
       refreshList();
       updatePlot();
   end

   function onStrSliderChange(~,~)
       val = get(hStrSlider,'Value');
       set(hStrEdit,'String',sprintf('%.2f',val));
       if selectedIdx < 1 || selectedIdx > length(elements), return; end
       elements(selectedIdx).strength = val;
       refreshList();
       updatePlot();
   end

   function onAngleSliderChange(~,~)
       val = get(hAngleSlider,'Value');
       set(hAngleEdit,'String',sprintf('%.1f',val));
       if selectedIdx < 1 || selectedIdx > length(elements), return; end
       elements(selectedIdx).angle = val;
       refreshList();
       updatePlot();
   end

   function onLocChange(~,~)
       if selectedIdx < 1 || selectedIdx > length(elements), return; end
       xv = str2double(get(hXEdit,'String'));
       yv = str2double(get(hYEdit,'String'));
       if ~isnan(xv), elements(selectedIdx).x = xv; end
       if ~isnan(yv), elements(selectedIdx).y = yv; end
       if elements(selectedIdx).type == 4
           av = str2double(get(hAngleEdit,'String'));
           if ~isnan(av)
               elements(selectedIdx).angle = av;
               set(hAngleSlider,'Value', max(-180, min(180, av)));
           end
       end
       refreshList();
       updatePlot();
   end

   %% ====== UI HELPERS ======
   function toggleAngleControls(visible)
       vis = 'off';
       if visible, vis = 'on'; end
       set(lblAngle,'Visible',vis);
       set(hAngleEdit,'Visible',vis);
       set(hAngleSlider,'Visible',vis);
   end

   function refreshList()
       typeNames = {'Src/Snk','Vortex','Doublet','Uniform'};
       strs = cell(1, length(elements));
       for i = 1:length(elements)
           e = elements(i);
           t = typeNames{e.type};
           if e.type == 4
               strs{i} = sprintf('#%d  %s   S=%.2f   ang=%.1f deg',...
                   i, t, e.strength, e.angle);
           else
               strs{i} = sprintf('#%d  %s   S=%.2f   (%.1f, %.1f)',...
                   i, t, e.strength, e.x, e.y);
           end
       end
       if isempty(strs)
           strs = {'(set type below, then click Add Element)'};
           set(hList,'String',strs,'Value',1);
       else
           if selectedIdx > 0 && selectedIdx <= length(elements)
               % Highlight selected element with arrow
               strs{selectedIdx} = ['>> ' strs{selectedIdx}];
               set(hList,'String',strs,'Value',selectedIdx);
           else
               % No selection - add hint at top
               strs = [{'-- (click element to edit) --'}, strs];
               set(hList,'String',strs,'Value',1);
           end
       end
   end

   function resetControls()
       % Reset controls to defaults for staging next element
       set(hType,'Value',1);
       set(hStrEdit,'String','1.00');
       set(hStrSlider,'Value',1);
       set(hXEdit,'String','0.00');
       set(hYEdit,'String','0.00');
       set(hAngleEdit,'String','0.0');
       set(hAngleSlider,'Value',0);
       toggleAngleControls(false);
   end

   function loadProperties(idx)
       e = elements(idx);
       set(hType,'Value',e.type);
       set(hStrEdit,'String',sprintf('%.2f',e.strength));
       set(hStrSlider,'Value', max(-20, min(20, e.strength)));
       set(hXEdit,'String',sprintf('%.2f',e.x));
       set(hYEdit,'String',sprintf('%.2f',e.y));
       set(hAngleEdit,'String',sprintf('%.1f',e.angle));
       set(hAngleSlider,'Value', max(-180, min(180, e.angle)));
       toggleAngleControls(e.type == 4);
   end

   %% ====== VELOCITY COMPUTATION ======
   function [u, v] = computeVelocity(X, Y)
       % Superpose contributions from all elements
       u = zeros(size(X));
       v = zeros(size(X));

       for i = 1:length(elements)
           e = elements(i);
           dx = X - e.x;
           dy = Y - e.y;
           r2 = dx.^2 + dy.^2;
           r2(r2 < 1e-8) = 1e-8;  % regularize singularity

           switch e.type
               case 1  % Source/Sink
                   % u = (m / 2*pi) * (x-x0)/r^2
                   % v = (m / 2*pi) * (y-y0)/r^2
                   coeff = e.strength / (2*pi);
                   u = u + coeff * dx ./ r2;
                   v = v + coeff * dy ./ r2;

               case 2  % Vortex (Gamma > 0 = counterclockwise)
                   % u = -(Gamma / 2*pi) * (y-y0)/r^2
                   % v =  (Gamma / 2*pi) * (x-x0)/r^2
                   coeff = e.strength / (2*pi);
                   u = u - coeff * dy ./ r2;
                   v = v + coeff * dx ./ r2;

               case 3  % Doublet (x-oriented)
                   % phi = -(kappa/2*pi) * (x-x0)/r^2
                   % u = -(kappa/2*pi) * ((y-y0)^2 - (x-x0)^2) / r^4
                   % v =  (kappa/2*pi) * 2*(x-x0)*(y-y0) / r^4
                   coeff = e.strength / (2*pi);
                   r4 = r2.^2;
                   r4(r4 < 1e-16) = 1e-16;
                   u = u - coeff * (dy.^2 - dx.^2) ./ r4;
                   v = v + coeff * (2 * dx .* dy) ./ r4;

               case 4  % Uniform flow
                   alpha = e.angle * pi / 180;
                   u = u + e.strength * cos(alpha);
                   v = v + e.strength * sin(alpha);
           end
       end
   end

   %% ====== PLOTTING ======
   function updatePlot()
       cla(ax); hold(ax,'on');

       if isempty(elements)
           title(ax,'Add elements to visualize potential flow',...
               'FontSize',13);
           axis(ax,[-5 5 -5 5]);
           axis(ax,'equal');
           grid(ax,'on');
           xlabel(ax,'x'); ylabel(ax,'y');
           return;
       end

       % --- Auto domain ---
       xc = []; yc = [];
       for i = 1:length(elements)
           if elements(i).type ~= 4
               xc(end+1) = elements(i).x; %#ok<AGROW>
               yc(end+1) = elements(i).y; %#ok<AGROW>
           end
       end
       if isempty(xc)
           xmin = -5; xmax = 5; ymin = -5; ymax = 5;
       else
           xmin = min(xc) - 4; xmax = max(xc) + 4;
           ymin = min(yc) - 4; ymax = max(yc) + 4;
       end
       % Ensure minimum size
       if (xmax - xmin) < 6
           mid = (xmax + xmin)/2;
           xmin = mid - 3; xmax = mid + 3;
       end
       if (ymax - ymin) < 6
           mid = (ymax + ymin)/2;
           ymin = mid - 3; ymax = mid + 3;
       end

       Ngrid = max(10, round(get(hGridSlider,'Value')));

       % --- Velocity field ---
       [X, Y] = meshgrid(linspace(xmin, xmax, Ngrid),...
                         linspace(ymin, ymax, Ngrid));
       [u, v] = computeVelocity(X, Y);

       % --- Velocity vectors ---
       if get(hCBVec,'Value')
           speed = sqrt(u.^2 + v.^2);
           % Clamp at 95th percentile to avoid singularity spikes
           capSpd = prctile(speed(:), 95);
           if capSpd > 0
               clampMask = speed > capSpd;
               scaleFactor = capSpd ./ max(speed, 1e-10);
               uPlot = u;
               vPlot = v;
               uPlot(clampMask) = u(clampMask) .* scaleFactor(clampMask);
               vPlot(clampMask) = v(clampMask) .* scaleFactor(clampMask);
           else
               uPlot = u;
               vPlot = v;
           end
           % Apply user scale factor (slider)
           vecSc = get(hVecScale,'Value');
           quiver(ax, X, Y, uPlot * vecSc, vPlot * vecSc, 0,...
               'Color',[0.15 0.3 0.7],'LineWidth',0.6);
       end

       % --- Streamlines ---
       if get(hCBStream,'Value')
           Nfine = Ngrid * 4;
           [Xf, Yf] = meshgrid(linspace(xmin, xmax, Nfine),...
                               linspace(ymin, ymax, Nfine));
           [uf, vf] = computeVelocity(Xf, Yf);
           dens = get(hStreamDens,'Value');
           hsl = streamslice(ax, Xf, Yf, uf, vf, dens);
           set(hsl,'Color',[0.85 0.15 0.1],'LineWidth',0.9);
       end

       % --- Mark element locations ---
       for i = 1:length(elements)
           e = elements(i);
           if e.type == 4, continue; end  % uniform has no location
           switch e.type
               case 1  % Source/Sink
                   if e.strength >= 0
                       plot(ax, e.x, e.y, 'o','MarkerSize',12,...
                           'MarkerFaceColor',[0.2 0.8 0.2],...
                           'MarkerEdgeColor','k','LineWidth',1.5);
                   else
                       plot(ax, e.x, e.y, 's','MarkerSize',12,...
                           'MarkerFaceColor',[0.9 0.2 0.2],...
                           'MarkerEdgeColor','k','LineWidth',1.5);
                   end
               case 2  % Vortex
                   plot(ax, e.x, e.y, 'd','MarkerSize',12,...
                       'MarkerFaceColor',[0.8 0.2 0.8],...
                       'MarkerEdgeColor','k','LineWidth',1.5);
               case 3  % Doublet
                   plot(ax, e.x, e.y, '^','MarkerSize',12,...
                       'MarkerFaceColor',[0.2 0.8 0.8],...
                       'MarkerEdgeColor','k','LineWidth',1.5);
           end
           % Label
           text(ax, e.x + 0.2, e.y + 0.2,...
               sprintf('#%d',i),'FontSize',9,'FontWeight','bold');
       end

       hold(ax,'off');
       axis(ax,[xmin xmax ymin ymax]);
       axis(ax,'equal');
       grid(ax,'on');
       xlabel(ax,'x','FontSize',11);
       ylabel(ax,'y','FontSize',11);
       title(ax,'2D Potential Flow (superposition)','FontSize',13);
   end

end
