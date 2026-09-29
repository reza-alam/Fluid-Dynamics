clear all; close all; clc;

% NACA 2412: interactive Kutta-condition teaching demo
%
% Three interchangeable singularity representations:
%   1) Constant-strength source panels + one circulation vortex
%   2) Linearly varying vortex panels; circulation constrained directly
%   3) Near-surface doublets + one interior circulation vortex
%
% A common collocation / method-of-fundamental-solutions framework is used
% so the three constructions can be compared interactively.
%
% The freestream is always horizontal.  The airfoil is rotated to alpha.
% Positive Gamma is counterclockwise.

%% User settings
Uinf   = 1.0;
alpha  = 5.0*pi/180;       % aerodynamic angle of attack [rad]
c      = 1.0;
Nside  = 60;               % control panels per surface
Nplot  = 140;              % smoother outline for plotting

% Interior contraction factors for the three representations.
% These are numerical-placement parameters, not aerodynamic parameters.
lambdaSource  = 0.96605;
lambdaVortex  = 0.88462;
lambdaDoublet = 0.98681;

%% Geometry: coarse control boundary and smooth display boundary
[xb0,yb0] = naca2412_closed(Nside,c);
[xbp0,ybp0] = naca2412_closed(Nplot,c);

thetaBody = -alpha;        % rotate body; keep freestream horizontal
xp = 0.25*c;
yp = 0.0;

[xb,yb]   = rotate_xy(xb0,yb0,thetaBody,xp,yp);
[xbp,ybp] = rotate_xy(xbp0,ybp0,thetaBody,xp,yp);

x1 = xb(1:end-1);  y1 = yb(1:end-1);
x2 = xb(2:end);    y2 = yb(2:end);

xc = 0.5*(x1+x2);
yc = 0.5*(y1+y2);

dx = x2-x1;
dy = y2-y1;
L  = hypot(dx,dy);

tx = dx./L;
ty = dy./L;

% Counterclockwise boundary => outward normal is right of tangent.
nx =  ty;
ny = -tx;

Np = numel(L);

% Interior reference point.
xref = mean(xb(1:end-1));
yref = mean(yb(1:end-1));

% Unit normalized circulation corresponds to Gamma = U_inf*c.
GammaUnit = Uinf*c;

% Horizontal freestream.
Un = Uinf*nx(:);
Ut = Uinf*tx(:);

%% Grid used for streamline/vector interpolation
xv = linspace(-0.40*c,1.40*c,360);
yv = linspace(-0.46*c,0.46*c,220);
[X,Y] = meshgrid(xv,yv);

inside = inpolygon(X,Y,xbp,ybp);

% Sample the velocity a small distance OUTSIDE the body.  This is used for
% robust Kutta residuals and stagnation-point detection, especially for
% vortex panels where the panel-limit tangential velocity is discontinuous.
epsSurf = 0.0010*c;
xSurf = xc(:) + epsSurf*nx(:);
ySurf = yc(:) + epsSurf*ny(:);

%% Common interior circulation-vortex influence
[unCirc,utCirc] = point_vortex_surface( ...
    xc,yc,nx,ny,tx,ty,xref,yref,1.0);

[UCircUnit,VCircUnit] = point_vortex_field( ...
    X,Y,xref,yref,GammaUnit);

%% ========================================================================
% METHOD 1: SOURCE PANELS + INTERIOR CIRCULATION VORTEX
% ========================================================================

% Use constant-strength source PANELS on the actual airfoil boundary.
% This fixes the streamline leakage caused by interior point sources.

[AnSP,AtSP,AnVP,AtVP] = panel_influence( ...
    xc,yc,x1,y1,L,tx,ty,nx,ny,8);

% Influence of one interior circulation vortex corresponding to
% Gamma/(U_inf*c) = 1.
[unCirc,utCirc] = point_vortex_surface( ...
    xc,yc,nx,ny,tx,ty,xref,yref,GammaUnit);

% sigma(G*) = sigma0 + G* sigma1
sigma0 = AnSP \ (-Un);
sigma1 = AnSP \ (-unCirc);

Vt0S = Ut + AtSP*sigma0;
Vt1S = AtSP*sigma1 + utCirc;

GkS = kutta_value(Vt0S,Vt1S);

% Flow-field basis from actual source panels.
[Us0,Vs0] = panel_field( ...
    X,Y,x1,y1,L,tx,ty,sigma0,0,8);

[Us1,Vs1] = panel_field( ...
    X,Y,x1,y1,L,tx,ty,sigma1,0,8);

MS = make_method_struct( ...
    'Source panels', ...
    'sourcepanel', ...
    xc(:),yc(:),[],[], ...
    sigma0,sigma1, ...
    Vt0S,Vt1S,GkS, ...
    Uinf+Us0,Vs0, ...
    Us1+UCircUnit,Vs1+VCircUnit, ...
    true,xref,yref);

%% ========================================================================
% METHOD 2: LINEAR VORTEX PANELS
% ========================================================================

% Use linearly varying vortex strength on each boundary panel:
%
%   gamma(s) = gamma_j N_1(s) + gamma_{j+1} N_2(s)
%
% There are Np+1 nodal vortex-sheet strengths.  We impose:
%   (1) Np no-penetration conditions, evaluated just outside the body;
%   (2) one prescribed total-circulation condition.
%
% This avoids the panel-to-panel sawtooth strength pattern produced by the
% old constant-strength vortex-panel solve.

[AnVL,AtVL] = linear_vortex_panel_influence( ...
    xSurf,ySurf,nx,ny,tx,ty,x1,y1,L,tx,ty,18);

% Exact circulation weights:
% Gamma = integral(gamma ds)
%       = sum_j L_j (gamma_j + gamma_{j+1})/2.
circW = zeros(Np+1,1);

for j = 1:Np
    circW(j)   = circW(j)   + 0.5*L(j);
    circW(j+1) = circW(j+1) + 0.5*L(j);
end

AV = [AnVL; circW.'];

% Basis solution for Gamma/(U_inf*c)=0.
gamma0 = AV \ [-Un; 0];

% Basis increment for Gamma/(U_inf*c)=1.
gamma1 = AV \ [zeros(Np,1); GammaUnit];

% The collocation points are already on the exterior side of the sheet, so
% these are the physical outer-surface tangential velocities.
Vt0V = Ut + AtVL*gamma0;
Vt1V = AtVL*gamma1;

GkV = kutta_value(Vt0V,Vt1V);

% Full flow field generated by the linear vortex panels.
[Uv0,Vv0] = linear_vortex_panel_field( ...
    X,Y,x1,y1,L,tx,ty,gamma0,12);

[Uv1,Vv1] = linear_vortex_panel_field( ...
    X,Y,x1,y1,L,tx,ty,gamma1,12);

MV = make_method_struct( ...
    'Vortex panels', ...
    'vortexpanel', ...
    xb(:),yb(:),[],[], ...
    gamma0,gamma1, ...
    Vt0V,Vt1V,GkV, ...
    Uinf+Uv0,Vv0, ...
    Uv1,Vv1, ...
    false,xref,yref);

%% ========================================================================
% METHOD 3: DOUBLET PANELS + INTERIOR CIRCULATION VORTEX
% ========================================================================

% Use true constant-strength DOUBLEt PANELS on the airfoil boundary.
% A constant doublet panel is equivalent to a vortex pair at its endpoints.
% Unlike the old interior point-doublet approximation, this directly
% enforces impermeability on the actual body surface.

[AnD,AtD] = doublet_panel_influence( ...
    xc,yc,nx,ny,tx,ty,x1,y1,x2,y2);

% The doublet system has a constant-strength null mode.  Fix its arbitrary
% additive constant with mean(mu)=0.
mu0 = constrained_ls_linear( ...
    AnD,-Un,ones(Np,1),0);

mu1 = constrained_ls_linear( ...
    AnD,-GammaUnit*unCirc,ones(Np,1),0);

% Full-field doublet contribution.
[Ud0,Vd0] = doublet_panel_field( ...
    X,Y,x1,y1,x2,y2,mu0);

[Ud1,Vd1] = doublet_panel_field( ...
    X,Y,x1,y1,x2,y2,mu1);

% Exterior surface velocity for Kutta/stagnation diagnostics.
[Ud0s,Vd0s] = doublet_panel_field( ...
    xSurf,ySurf,x1,y1,x2,y2,mu0);

[Ud1s,Vd1s] = doublet_panel_field( ...
    xSurf,ySurf,x1,y1,x2,y2,mu1);

[UCs,VCs] = point_vortex_field( ...
    xSurf,ySurf,xref,yref,GammaUnit);

Vt0D = (Uinf + Ud0s(:)).*tx(:) + Vd0s(:).*ty(:);
Vt1D = (Ud1s(:)+UCs(:)).*tx(:) + ...
       (Vd1s(:)+VCs(:)).*ty(:);

GkD = kutta_value(Vt0D,Vt1D);

MD = make_method_struct( ...
    'Doublet panels', ...
    'doublet', ...
    xc(:),yc(:),nx(:),ny(:), ...
    mu0,mu1, ...
    Vt0D,Vt1D,GkD, ...
    Uinf+Ud0,Vd0, ...
    Ud1+UCircUnit,Vd1+VCircUnit, ...
    true,xref,yref);

methods = [MS MV MD];

fprintf('\nKutta circulation values:\n');
for k = 1:numel(methods)
    fprintf('  %-16s Gamma_K/(U_inf*c) = %+8.5f\n', ...
        methods(k).name,methods(k).GammaKstar);
end

%% Figure and controls
fig = figure('Color','w', ...
    'Position',[80 80 1280 720], ...
    'Name','NACA 2412: singularities and the Kutta condition', ...
    'NumberTitle','off');

axPos = [0.07 0.20 0.88 0.69];
ax = axes('Parent',fig,'Position',axPos);

% Method selector
annotation(fig,'textbox',[0.07 0.915 0.13 0.040], ...
    'String','Construction method:', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','middle', ...
    'FontSize',11);

methodPopup = uicontrol(fig,'Style','popupmenu', ...
    'Units','normalized', ...
    'Position',[0.205 0.920 0.21 0.036], ...
    'String',{'Source panels','Vortex panels','Doublet panels'}, ...
    'Value',1, ...
    'FontSize',10, ...
    'Callback',@method_changed);

% Fixed wide normalized circulation range.
Glim = max(3.0,5.0*max(abs([methods.GammaKstar])));
Gmin = -Glim;
Gmax =  Glim;

sliderPos = [0.21 0.118 0.49 0.028];

dGsmall = 0.001;   % left/right arrow click increment
dGlarge = 0.020;   % trough/page click increment
sliderStep = [ ...
    min(1,dGsmall/(Gmax-Gmin)), ...
    min(1,dGlarge/(Gmax-Gmin))];

slider = uicontrol(fig,'Style','slider', ...
    'Units','normalized', ...
    'Position',sliderPos, ...
    'Min',Gmin,'Max',Gmax,'Value',0, ...
    'SliderStep',sliderStep, ...
    'BusyAction','cancel', ...
    'Interruptible','off', ...
    'Callback',@redraw_demo);

annotation(fig,'textbox',[0.055 0.108 0.145 0.040], ...
    'String','$\Gamma/(U_\infty c)$', ...
    'Interpreter','latex', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','middle', ...
    'FontSize',12);

valueText = annotation(fig,'textbox',[0.715 0.108 0.085 0.040], ...
    'String','', ...
    'Interpreter','latex', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','middle', ...
    'FontSize',11);

uicontrol(fig,'Style','pushbutton', ...
    'Units','normalized', ...
    'Position',[0.82 0.110 0.10 0.040], ...
    'String','Set Kutta', ...
    'FontWeight','bold', ...
    'Callback',@set_kutta);

kuttaLabel = annotation(fig,'textbox', ...
    [0.21 0.086 0.07 0.026], ...
    'String','$\Gamma_K$', ...
    'Interpreter','latex', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'Color',[0 0.45 0], ...
    'FontWeight','bold', ...
    'FontSize',11);

% Streamline density
densityMin = 9;
densityMax = 121;
densityDefault = 23;

densitySlider = uicontrol(fig,'Style','slider', ...
    'Units','normalized', ...
    'Position',[0.21 0.067 0.49 0.026], ...
    'Min',densityMin,'Max',densityMax,'Value',densityDefault, ...
    'SliderStep',[1/(densityMax-densityMin) 10/(densityMax-densityMin)], ...
    'BusyAction','cancel', ...
    'Interruptible','off', ...
    'Callback',@redraw_demo);

annotation(fig,'textbox',[0.055 0.058 0.145 0.038], ...
    'String','$N_{\mathrm{streamlines}}$', ...
    'Interpreter','latex', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','middle', ...
    'FontSize',11);

densityText = annotation(fig,'textbox',[0.715 0.058 0.085 0.038], ...
    'String','', ...
    'Interpreter','latex', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','middle', ...
    'FontSize',11);

% Velocity-vector toggle
vectorButton = uicontrol(fig,'Style','togglebutton', ...
    'Units','normalized', ...
    'Position',[0.82 0.018 0.10 0.040], ...
    'String','Velocity vectors', ...
    'Value',0, ...
    'FontWeight','bold', ...
    'Callback',@redraw_demo);

% Velocity-vector density
vectorDensityMin = 0.5;
vectorDensityMax = 3.0;
vectorDensityDefault = 1.0;

vectorDensitySlider = uicontrol(fig,'Style','slider', ...
    'Units','normalized', ...
    'Position',[0.21 0.020 0.49 0.026], ...
    'Min',vectorDensityMin,'Max',vectorDensityMax, ...
    'Value',vectorDensityDefault, ...
    'SliderStep',[0.1/(vectorDensityMax-vectorDensityMin) ...
                  0.5/(vectorDensityMax-vectorDensityMin)], ...
    'BusyAction','cancel', ...
    'Interruptible','off', ...
    'Callback',@redraw_demo);

annotation(fig,'textbox',[0.055 0.011 0.145 0.038], ...
    'String','$N_{V,x}\times N_{V,y}$', ...
    'Interpreter','latex', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','middle', ...
    'FontSize',11);

vectorDensityText = annotation(fig,'textbox',[0.715 0.011 0.095 0.038], ...
    'String','', ...
    'Interpreter','latex', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','middle', ...
    'FontSize',11);

%% Store data
D.ax = ax;
D.axPos = axPos;
D.methods = methods;
D.methodPopup = methodPopup;

D.slider = slider;
D.sliderPos = sliderPos;
D.valueText = valueText;
D.kuttaLabel = kuttaLabel;

D.densitySlider = densitySlider;
D.densityText = densityText;

D.vectorButton = vectorButton;
D.vectorDensitySlider = vectorDensitySlider;
D.vectorDensityText = vectorDensityText;

D.xb = xbp;
D.yb = ybp;
D.xc = xc;
D.yc = yc;
D.tx = tx;
D.ty = ty;

D.X = X;
D.Y = Y;
D.inside = inside;

D.Uinf = Uinf;
D.alpha = alpha;
D.c = c;
D.Gmin = Gmin;
D.Gmax = Gmax;

setappdata(fig,'KuttaDemoData',D);

update_kutta_label(fig);
redraw_demo(slider,[]);

%% ========================================================================
function method_changed(src,~)

fig = ancestor(src,'figure');
update_kutta_label(fig);
redraw_demo(src,[]);

end

%% ========================================================================
function update_kutta_label(fig)

D = getappdata(fig,'KuttaDemoData');
m = get(D.methodPopup,'Value');

Gk = D.methods(m).GammaKstar;

frac = (Gk-D.Gmin)/(D.Gmax-D.Gmin);
frac = max(0,min(1,frac));

p = D.sliderPos;
xCenter = p(1) + p(3)*frac;

set(D.kuttaLabel,'Position', ...
    [xCenter-0.035 0.086 0.07 0.026]);

end

%% ========================================================================
function redraw_demo(src,~)

fig = ancestor(src,'figure');
D = getappdata(fig,'KuttaDemoData');

mIndex = get(D.methodPopup,'Value');
M = D.methods(mIndex);

Gstar = get(D.slider,'Value');

nStreams = max(9,round(get(D.densitySlider,'Value')));
showVectors = logical(get(D.vectorButton,'Value'));

vectorScale = get(D.vectorDensitySlider,'Value');
nVecX = max(6,round(23*vectorScale));
nVecY = max(4,round(13*vectorScale));

strength = M.s0 + Gstar*M.s1;
Vt = M.Vt0 + Gstar*M.Vt1;

kuttaResidual = Vt(1) + Vt(end);

U = M.U0 + Gstar*M.U1;
V = M.V0 + Gstar*M.V1;

U(D.inside) = NaN;
V(D.inside) = NaN;

ax = D.ax;

legend(ax,'off');
delete(allchild(ax));
cla(ax,'reset');
set(ax,'Position',D.axPos);
hold(ax,'on');

%% Streamlines
ySeeds = linspace(-0.43*D.c,0.43*D.c,nStreams);
xSeeds = -0.38*D.c*ones(size(ySeeds));

hs = streamline(ax,D.X,D.Y,U,V,xSeeds,ySeeds);

if ~isempty(hs)
    set(hs, ...
        'Color',[0.20 0.45 0.82], ...
        'LineWidth',0.95, ...
        'HandleVisibility','off');
end

%% Optional velocity-vector overlay
hVectorLegend = [];

if showVectors
    ix = unique(round(linspace(1,size(D.X,2),nVecX)));
    iy = unique(round(linspace(1,size(D.X,1),nVecY)));

    Xq = D.X(iy,ix);
    Yq = D.Y(iy,ix);
    Uq = U(iy,ix);
    Vq = V(iy,ix);

    good = isfinite(Uq) & isfinite(Vq);

    quiver(ax,Xq(good),Yq(good),Uq(good),Vq(good),0.75, ...
        'Color',[0.10 0.10 0.10], ...
        'LineWidth',0.70, ...
        'MaxHeadSize',0.70, ...
        'HandleVisibility','off');

    hVectorLegend = plot(ax,NaN,NaN,'-', ...
        'Color',[0.10 0.10 0.10], ...
        'LineWidth',1.1);
end

%% Airfoil
patch(ax,D.xb,D.yb,[0.93 0.93 0.93], ...
    'EdgeColor',[0 0 0], ...
    'LineWidth',1.7, ...
    'FaceAlpha',0.60, ...
    'HandleVisibility','off');

%% Singularity markers
[hPos,hNeg] = draw_singularities(ax,M,strength,D.c);

% Source and doublet representations use a separate interior circulation
% vortex.  The vortex representation does not: its total circulation is
% the sum of the distributed vortex strengths.
hCirc = [];
if M.hasCircVortex
    circColor = [0.05 0.55 0.15];

    if abs(Gstar) > 1e-12
        plot(ax,M.xCirc,M.yCirc,'o', ...
            'MarkerSize',10, ...
            'MarkerFaceColor','none', ...
            'MarkerEdgeColor',circColor, ...
            'LineWidth',2.0, ...
            'HandleVisibility','off');

        text(ax,M.xCirc,M.yCirc,'$\Gamma$', ...
            'Interpreter','latex', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'Color',circColor, ...
            'FontSize',9, ...
            'FontWeight','bold');
    end

    hCirc = plot(ax,NaN,NaN,'o', ...
        'MarkerSize',8, ...
        'MarkerFaceColor','none', ...
        'MarkerEdgeColor',circColor, ...
        'LineWidth',1.7);
end

%% Surface stagnation points
[xstag,ystag] = surface_stagnation(D.xc,D.yc,Vt);

if ~isempty(xstag)
    scatter(ax,xstag,ystag,70,'d', ...
        'MarkerFaceColor',[0.95 0.25 0.75], ...
        'MarkerEdgeColor',[0.40 0 0.30], ...
        'LineWidth',1.0, ...
        'HandleVisibility','off');
end

%% Trailing-edge surface-velocity vectors
arrowScale = 0.10*D.c/D.Uinf;

vu = Vt(1)*[D.tx(1) D.ty(1)];
vl = Vt(end)*[D.tx(end) D.ty(end)];

quiver(ax,D.xc(1),D.yc(1), ...
    arrowScale*vu(1),arrowScale*vu(2),0, ...
    'Color',[0.10 0.55 0.10], ...
    'LineWidth',1.6, ...
    'MaxHeadSize',0.8, ...
    'HandleVisibility','off');

quiver(ax,D.xc(end),D.yc(end), ...
    arrowScale*vl(1),arrowScale*vl(2),0, ...
    'Color',[0.10 0.55 0.10], ...
    'LineWidth',1.6, ...
    'MaxHeadSize',0.8, ...
    'HandleVisibility','off');

%% Kutta status
tol = 0.03*D.Uinf;

if abs(kuttaResidual) < tol
    teColor = [0.10 0.65 0.10];
    statusLatex = '\mathrm{Kutta\ satisfied}';
else
    teColor = [0.85 0.20 0.10];
    statusLatex = '\mathrm{not\ Kutta}';
end

plot(ax,D.xb(1),D.yb(1),'o', ...
    'MarkerSize',9, ...
    'MarkerFaceColor',teColor, ...
    'MarkerEdgeColor','k', ...
    'HandleVisibility','off');

%% Horizontal freestream arrow
q0x = -0.34*D.c;
q0y =  0.405*D.c;

quiver(ax,q0x,q0y,0.17*D.c,0,0, ...
    'Color',[0.10 0.10 0.10], ...
    'LineWidth',1.5, ...
    'MaxHeadSize',0.8, ...
    'HandleVisibility','off');

text(ax,q0x+0.06*D.c,q0y+0.018*D.c,'$U_\infty$', ...
    'Interpreter','latex', ...
    'FontSize',12, ...
    'HorizontalAlignment','center');

%% Information box
line1 = sprintf([ ...
    '$\\alpha=%.1f^{\\circ},\\quad ' ...
    '\\Gamma/(U_\\infty c)=%+.3f,\\quad ' ...
    '\\Gamma_K/(U_\\infty c)=%+.3f$'], ...
    D.alpha*180/pi,Gstar,M.GammaKstar);

line2 = sprintf([ ...
    '$\\left(V_{t,u}+V_{t,l}\\right)/U_\\infty=%+.3f,\\quad %s$'], ...
    kuttaResidual/D.Uinf,statusLatex);

text(ax,-0.35*D.c,0.355*D.c, ...
    {M.name,line1,line2}, ...
    'Interpreter','latex', ...
    'FontSize',10, ...
    'VerticalAlignment','top', ...
    'BackgroundColor','w', ...
    'EdgeColor',[0.75 0.75 0.75], ...
    'Margin',6);

%% Legend
hStream = plot(ax,NaN,NaN,'-', ...
    'Color',[0.20 0.45 0.82], ...
    'LineWidth',1.2);

hStag = plot(ax,NaN,NaN,'d', ...
    'MarkerFaceColor',[0.95 0.25 0.75], ...
    'MarkerEdgeColor',[0.40 0 0.30]);

[legendHandles,legendLabels] = method_legend( ...
    M,hStream,hPos,hNeg,hStag,hCirc,hVectorLegend,showVectors);

legend(ax,legendHandles,legendLabels, ...
    'Interpreter','latex', ...
    'Location','northeast');

%% Axes
axis(ax,'equal');
xlim(ax,[-0.40 1.40]*D.c);
ylim(ax,[-0.46 0.46]*D.c);

xlabel(ax,'$x/c$','Interpreter','latex');
ylabel(ax,'$y/c$','Interpreter','latex');

grid(ax,'on');
box(ax,'on');

set(ax, ...
    'FontSize',10, ...
    'Layer','top', ...
    'TickLabelInterpreter','latex');

%% Control readouts
set(D.valueText,'String',sprintf('$%+.3f$',Gstar));
set(D.densityText,'String',sprintf('$%d$',nStreams));
set(D.vectorDensityText,'String',sprintf('$%d\\times%d$',nVecX,nVecY));

if showVectors
    set(D.vectorButton,'String','Hide vectors');
else
    set(D.vectorButton,'String','Velocity vectors');
end

drawnow;

end

%% ========================================================================
function set_kutta(src,~)

fig = ancestor(src,'figure');
D = getappdata(fig,'KuttaDemoData');

m = get(D.methodPopup,'Value');
set(D.slider,'Value',D.methods(m).GammaKstar);

redraw_demo(D.slider,[]);

end

%% ========================================================================
function [hPos,hNeg] = draw_singularities(ax,M,s,Dc)

mag = abs(s);
scale = max(mag);

if scale < eps
    ms = 24*ones(size(mag));
else
    ms = 18 + 52*mag/scale;
end

pos = s >= 0;
neg = s < 0;

switch M.type

    case 'sourcepanel'
        scatter(ax,M.xs(pos),M.ys(pos),ms(pos), ...
            'o', ...
            'MarkerFaceColor',[0.90 0.10 0.10], ...
            'MarkerEdgeColor',[0.65 0 0], ...
            'LineWidth',0.7, ...
            'HandleVisibility','off');

        scatter(ax,M.xs(neg),M.ys(neg),ms(neg), ...
            's', ...
            'MarkerFaceColor',[0.05 0.05 0.05], ...
            'MarkerEdgeColor',[0 0 0], ...
            'LineWidth',0.7, ...
            'HandleVisibility','off');

        hPos = plot(ax,NaN,NaN,'o', ...
            'MarkerFaceColor',[0.90 0.10 0.10], ...
            'MarkerEdgeColor',[0.65 0 0]);

        hNeg = plot(ax,NaN,NaN,'s', ...
            'MarkerFaceColor',[0.05 0.05 0.05], ...
            'MarkerEdgeColor',[0 0 0]);

    case 'vortexpanel'
        scatter(ax,M.xs(pos),M.ys(pos),ms(pos), ...
            '^', ...
            'MarkerFaceColor',[0.90 0.10 0.10], ...
            'MarkerEdgeColor',[0.65 0 0], ...
            'LineWidth',0.7, ...
            'HandleVisibility','off');

        scatter(ax,M.xs(neg),M.ys(neg),ms(neg), ...
            'v', ...
            'MarkerFaceColor',[0.05 0.05 0.05], ...
            'MarkerEdgeColor',[0 0 0], ...
            'LineWidth',0.7, ...
            'HandleVisibility','off');

        hPos = plot(ax,NaN,NaN,'^', ...
            'MarkerFaceColor',[0.90 0.10 0.10], ...
            'MarkerEdgeColor',[0.65 0 0]);

        hNeg = plot(ax,NaN,NaN,'v', ...
            'MarkerFaceColor',[0.05 0.05 0.05], ...
            'MarkerEdgeColor',[0 0 0]);

    case 'doublet'
        scatter(ax,M.xs(pos),M.ys(pos),ms(pos), ...
            'd', ...
            'MarkerFaceColor',[0.90 0.10 0.10], ...
            'MarkerEdgeColor',[0.65 0 0], ...
            'LineWidth',0.7, ...
            'HandleVisibility','off');

        scatter(ax,M.xs(neg),M.ys(neg),ms(neg), ...
            'd', ...
            'MarkerFaceColor',[0.05 0.05 0.05], ...
            'MarkerEdgeColor',[0 0 0], ...
            'LineWidth',0.7, ...
            'HandleVisibility','off');

        % Small orientation ticks for the doublet axes.
        tick = 0.010*Dc;
        for j = 1:numel(M.xs)
            plot(ax, ...
                M.xs(j)+tick*[-M.dirx(j) M.dirx(j)], ...
                M.ys(j)+tick*[-M.diry(j) M.diry(j)], ...
                '-', ...
                'Color',[0.30 0.30 0.30], ...
                'LineWidth',0.5, ...
                'HandleVisibility','off');
        end

        hPos = plot(ax,NaN,NaN,'d', ...
            'MarkerFaceColor',[0.90 0.10 0.10], ...
            'MarkerEdgeColor',[0.65 0 0]);

        hNeg = plot(ax,NaN,NaN,'d', ...
            'MarkerFaceColor',[0.05 0.05 0.05], ...
            'MarkerEdgeColor',[0 0 0]);
end

end

%% ========================================================================
function [hh,ll] = method_legend( ...
    M,hStream,hPos,hNeg,hStag,hCirc,hVector,showVectors)

switch M.type
    case 'sourcepanel'
        hh = [hStream hPos hNeg hStag];
        ll = { ...
            '$\mathrm{Streamlines}$', ...
            '$\sigma_j>0\ \mathrm{(source\ panels)}$', ...
            '$\sigma_j<0\ \mathrm{(sink\ panels)}$', ...
            '$\mathrm{Stagnation\ points}$'};

    case 'vortexpanel'
        hh = [hStream hPos hNeg hStag];
        ll = { ...
            '$\mathrm{Streamlines}$', ...
            '$\gamma_j>0\ \mathrm{(CCW)}$', ...
            '$\gamma_j<0\ \mathrm{(CW)}$', ...
            '$\mathrm{Stagnation\ points}$'};

    case 'doublet'
        hh = [hStream hPos hNeg hStag];
        ll = { ...
            '$\mathrm{Streamlines}$', ...
            '$\mu_j>0\ \mathrm{(doublets)}$', ...
            '$\mu_j<0\ \mathrm{(doublets)}$', ...
            '$\mathrm{Stagnation\ points}$'};
end

if M.hasCircVortex
    hh = [hh hCirc];
    ll = [ll {'$\mathrm{Circulation\ vortex}\ \Gamma$'}];
end

if showVectors
    hh = [hh hVector];
    ll = [ll {'$\mathbf{V}(x,y)$'}];
end

end

%% ========================================================================
function Gk = kutta_value(Vt0,Vt1)

R0 = Vt0(1) + Vt0(end);
R1 = Vt1(1) + Vt1(end);

Gk = -R0/R1;

end

%% ========================================================================
function [xstag,ystag] = surface_stagnation(xc,yc,Vt)

% Smooth only along the ordered surface contour.  This suppresses
% panel-scale sign chatter without connecting across the sharp TE.
V = Vt(:);
n = numel(V);

win = nine_odd(min(n,max(7,round(n/18))));
Vsm = movmean(V,win,'Endpoints','shrink');

xc = xc(:);
yc = yc(:);

candX = [];
candY = [];

for i = 1:n-1
    if Vsm(i) == 0
        candX(end+1,1) = xc(i); %#ok<AGROW>
        candY(end+1,1) = yc(i); %#ok<AGROW>

    elseif Vsm(i)*Vsm(i+1) < 0
        f = abs(Vsm(i))/(abs(Vsm(i))+abs(Vsm(i+1)));

        candX(end+1,1) = xc(i) + f*(xc(i+1)-xc(i)); %#ok<AGROW>
        candY(end+1,1) = yc(i) + f*(yc(i+1)-yc(i)); %#ok<AGROW>
    end
end

% Merge roots that are merely adjacent-panel versions of the same
% physical stagnation point.
if isempty(candX)
    xstag = [];
    ystag = [];
    return
end

keep = true(size(candX));
mergeTol = 0.035;

for i = 2:numel(candX)
    if hypot(candX(i)-candX(i-1),candY(i)-candY(i-1)) < mergeTol
        candX(i-1) = 0.5*(candX(i-1)+candX(i));
        candY(i-1) = 0.5*(candY(i-1)+candY(i));
        keep(i) = false;
    end
end

xstag = candX(keep);
ystag = candY(keep);

end

%% ========================================================================
function w = nine_odd(n)

% Return a positive odd integer no larger than n.
w = max(1,round(n));

if mod(w,2) == 0
    w = w-1;
end

w = max(1,w);

end

%% ========================================================================
function M = make_method_struct( ...
    name,type,xs,ys,dirx,diry,s0,s1,Vt0,Vt1,Gk,U0,V0,U1,V1, ...
    hasCircVortex,xCirc,yCirc)

M.name = name;
M.type = type;

M.xs = xs(:);
M.ys = ys(:);

M.dirx = dirx;
M.diry = diry;

M.s0 = s0(:);
M.s1 = s1(:);

M.Vt0 = Vt0(:);
M.Vt1 = Vt1(:);

M.GammaKstar = Gk;

M.U0 = U0;
M.V0 = V0;
M.U1 = U1;
M.V1 = V1;

M.hasCircVortex = hasCircVortex;
M.xCirc = xCirc;
M.yCirc = yCirc;

end

%% ========================================================================
function [An,At] = point_source_influence( ...
    xc,yc,nx,ny,tx,ty,xs,ys)

Np = numel(xc);
Ns = numel(xs);

An = zeros(Np,Ns);
At = zeros(Np,Ns);

for j = 1:Ns
    rx = xc(:)-xs(j);
    ry = yc(:)-ys(j);
    r2 = rx.^2 + ry.^2;

    u = rx./(2*pi*r2);
    v = ry./(2*pi*r2);

    An(:,j) = nx(:).*u + ny(:).*v;
    At(:,j) = tx(:).*u + ty(:).*v;
end

end

%% ========================================================================
function [An,At] = point_vortex_influence( ...
    xc,yc,nx,ny,tx,ty,xs,ys)

Np = numel(xc);
Ns = numel(xs);

An = zeros(Np,Ns);
At = zeros(Np,Ns);

for j = 1:Ns
    rx = xc(:)-xs(j);
    ry = yc(:)-ys(j);
    r2 = rx.^2 + ry.^2;

    u = -ry./(2*pi*r2);
    v =  rx./(2*pi*r2);

    An(:,j) = nx(:).*u + ny(:).*v;
    At(:,j) = tx(:).*u + ty(:).*v;
end

end

%% ========================================================================
function [un,ut] = point_vortex_surface( ...
    xc,yc,nx,ny,tx,ty,x0,y0,Gamma)

rx = xc(:)-x0;
ry = yc(:)-y0;
r2 = rx.^2 + ry.^2;

u = -Gamma*ry./(2*pi*r2);
v =  Gamma*rx./(2*pi*r2);

un = nx(:).*u + ny(:).*v;
ut = tx(:).*u + ty(:).*v;

end

%% ========================================================================
function [An,At] = linear_vortex_panel_influence( ...
    xp,yp,nxp,nyp,txp,typ,x1,y1,L,tx,ty,nQuad)

% Influence of linearly varying vortex panels evaluated at arbitrary
% collocation points.  Unknown strengths live at panel ENDPOINTS, so an
% N-panel closed contour has N+1 unknown nodal strengths; the two TE values
% are intentionally independent unless the Kutta condition selects them.

Nc = numel(xp);
Np = numel(L);

An = zeros(Nc,Np+1);
At = zeros(Nc,Np+1);

[z,w] = gauss_legendre(nQuad);

for j = 1:Np

    s  = 0.5*(z+1)*L(j);
    ww = 0.5*L(j)*w;

    N1 = 1 - s/L(j);
    N2 = s/L(j);

    xq = x1(j) + tx(j)*s;
    yq = y1(j) + ty(j)*s;

    for i = 1:Nc

        rx = xp(i)-xq;
        ry = yp(i)-yq;
        r2 = max(rx.^2 + ry.^2,1e-14);

        u = -ry./(2*pi*r2);
        v =  rx./(2*pi*r2);

        un = nxp(i)*u + nyp(i)*v;
        ut = txp(i)*u + typ(i)*v;

        An(i,j)   = An(i,j)   + sum(ww.*N1.*un);
        An(i,j+1) = An(i,j+1) + sum(ww.*N2.*un);

        At(i,j)   = At(i,j)   + sum(ww.*N1.*ut);
        At(i,j+1) = At(i,j+1) + sum(ww.*N2.*ut);
    end
end

end

%% ========================================================================
function [u,v] = linear_vortex_panel_field( ...
    X,Y,x1,y1,L,tx,ty,gammaNode,nQuad)

u = zeros(size(X));
v = zeros(size(Y));

[z,w] = gauss_legendre(nQuad);

for j = 1:numel(L)

    s  = 0.5*(z+1)*L(j);
    ww = 0.5*L(j)*w;

    N1 = 1 - s/L(j);
    N2 = s/L(j);

    for k = 1:nQuad

        xq = x1(j) + tx(j)*s(k);
        yq = y1(j) + ty(j)*s(k);

        gammaQ = gammaNode(j)*N1(k) + gammaNode(j+1)*N2(k);

        rx = X-xq;
        ry = Y-yq;
        r2 = max(rx.^2 + ry.^2,1e-12);

        u = u + gammaQ*ww(k).*(-ry)./(2*pi*r2);
        v = v + gammaQ*ww(k).*rx./(2*pi*r2);
    end
end

end

%% ========================================================================
function [An,At] = doublet_panel_influence( ...
    xc,yc,nx,ny,tx,ty,x1,y1,x2,y2)

% Velocity induced by a unit constant-strength 2-D doublet panel.
% A doublet panel equals a +vortex at endpoint 2 and -vortex at endpoint 1.

Np = numel(xc);

An = zeros(Np,Np);
At = zeros(Np,Np);

for j = 1:Np

    rx1 = xc(:)-x1(j);
    ry1 = yc(:)-y1(j);
    r21 = max(rx1.^2 + ry1.^2,1e-14);

    rx2 = xc(:)-x2(j);
    ry2 = yc(:)-y2(j);
    r22 = max(rx2.^2 + ry2.^2,1e-14);

    u = (-ry2./r22 + ry1./r21)/(2*pi);
    v = ( rx2./r22 - rx1./r21)/(2*pi);

    An(:,j) = nx(:).*u + ny(:).*v;
    At(:,j) = tx(:).*u + ty(:).*v;
end

end

%% ========================================================================
function [u,v] = doublet_panel_field(X,Y,x1,y1,x2,y2,mu)

u = zeros(size(X));
v = zeros(size(Y));

for j = 1:numel(mu)

    rx1 = X-x1(j);
    ry1 = Y-y1(j);
    r21 = max(rx1.^2 + ry1.^2,1e-12);

    rx2 = X-x2(j);
    ry2 = Y-y2(j);
    r22 = max(rx2.^2 + ry2.^2,1e-12);

    u = u + mu(j)*(-ry2./r22 + ry1./r21)/(2*pi);
    v = v + mu(j)*( rx2./r22 - rx1./r21)/(2*pi);
end

end

%% ========================================================================
function x = constrained_ls_linear(A,b,c,d)

% Solve min ||A*x-b||_2 subject to c'*x=d.
c = c(:);

K = [A.'*A, c; c.', 0];
rhs = [A.'*b; d];

sol = K \ rhs;
x = sol(1:size(A,2));

end

%% ========================================================================
function [u,v] = point_source_field(X,Y,xs,ys,q)

u = zeros(size(X));
v = zeros(size(Y));

for j = 1:numel(xs)
    rx = X-xs(j);
    ry = Y-ys(j);
    r2 = max(rx.^2 + ry.^2,1e-12);

    u = u + q(j)*rx./(2*pi*r2);
    v = v + q(j)*ry./(2*pi*r2);
end

end

%% ========================================================================
function [u,v] = point_vortex_field(X,Y,x0,y0,Gamma)

rx = X-x0;
ry = Y-y0;
r2 = max(rx.^2 + ry.^2,1e-12);

u = -Gamma*ry./(2*pi*r2);
v =  Gamma*rx./(2*pi*r2);

end

%% ========================================================================
function [u,v] = point_vortex_field_many(X,Y,xs,ys,gamma)

u = zeros(size(X));
v = zeros(size(Y));

for j = 1:numel(xs)
    rx = X-xs(j);
    ry = Y-ys(j);
    r2 = max(rx.^2 + ry.^2,1e-12);

    u = u - gamma(j)*ry./(2*pi*r2);
    v = v + gamma(j)*rx./(2*pi*r2);
end

end

%% ========================================================================
function [u,v] = point_doublet_field( ...
    X,Y,xs,ys,dirx,diry,mu)

u = zeros(size(X));
v = zeros(size(Y));

for j = 1:numel(xs)
    rx = X-xs(j);
    ry = Y-ys(j);

    r2 = max(rx.^2 + ry.^2,1e-12);
    rd = dirx(j)*rx + diry(j)*ry;

    u = u - mu(j)*(dirx(j)./r2 - 2*rd.*rx./r2.^2)/(2*pi);
    v = v - mu(j)*(diry(j)./r2 - 2*rd.*ry./r2.^2)/(2*pi);
end

end

%% ========================================================================
function [AnS,AtS,AnV,AtV] = panel_influence( ...
    xc,yc,x1,y1,L,tx,ty,nx,ny,nQuad)

Np = numel(L);

AnS = zeros(Np,Np);
AtS = zeros(Np,Np);
AnV = zeros(Np,Np);
AtV = zeros(Np,Np);

[z,w] = gauss_legendre(nQuad);

for j = 1:Np
    s  = 0.5*(z+1)*L(j);
    ww = 0.5*L(j)*w;

    xq = x1(j) + tx(j)*s;
    yq = y1(j) + ty(j)*s;

    for i = 1:Np
        if i == j
            % Exterior limiting values for a CCW boundary.
            AnS(i,j) = 0.5;
            AtS(i,j) = 0.0;
            AnV(i,j) = 0.0;
            AtV(i,j) = 0.5;
            continue
        end

        rx = xc(i)-xq;
        ry = yc(i)-yq;
        r2 = rx.^2 + ry.^2;

        us = sum(ww.*rx./(2*pi*r2));
        vs = sum(ww.*ry./(2*pi*r2));

        uv = sum(ww.*(-ry)./(2*pi*r2));
        vv = sum(ww.*rx./(2*pi*r2));

        AnS(i,j) = nx(i)*us + ny(i)*vs;
        AtS(i,j) = tx(i)*us + ty(i)*vs;

        AnV(i,j) = nx(i)*uv + ny(i)*vv;
        AtV(i,j) = tx(i)*uv + ty(i)*vv;
    end
end

end

%% ========================================================================
function [u,v] = panel_field( ...
    X,Y,x1,y1,L,tx,ty,sigma,gamma,nQuad)

u = zeros(size(X));
v = zeros(size(Y));

[z,w] = gauss_legendre(nQuad);

for j = 1:numel(L)
    s  = 0.5*(z+1)*L(j);
    ww = 0.5*L(j)*w;

    for k = 1:nQuad
        xq = x1(j) + tx(j)*s(k);
        yq = y1(j) + ty(j)*s(k);

        rx = X-xq;
        ry = Y-yq;
        r2 = max(rx.^2 + ry.^2,1e-12);

        dS = ww(k);

        u = u + sigma(j)*dS.*rx./(2*pi*r2);
        v = v + sigma(j)*dS.*ry./(2*pi*r2);

        if isscalar(gamma)
            gammaj = gamma;
        else
            gammaj = gamma(j);
        end

        u = u + gammaj*dS.*(-ry)./(2*pi*r2);
        v = v + gammaj*dS.*rx./(2*pi*r2);
    end
end

end

%% ========================================================================
function [x,w] = gauss_legendre(n)

k = 1:n-1;
beta = k./sqrt(4*k.^2-1);

T = diag(beta,1) + diag(beta,-1);
[V,D] = eig(T);

x = diag(D);
[x,idx] = sort(x);
V = V(:,idx);

w = 2*(V(1,:)'.^2);

end

%% ========================================================================
function g = vortex_panel_constrained_solve(A,b,L,targetGamma)

% Solve the vortex-panel no-penetration equations while prescribing
% total circulation:
%
%       A*g = b,       sum(g_j L_j) = targetGamma.
%
% For a closed contour A has an (approximately) one-dimensional null
% mode corresponding to circulation.  We use an SVD to isolate that mode,
% which is substantially more stable than interior point vortices.

[U,S,V] = svd(A,'econ');
s = diag(S);

if isempty(s)
    error('Vortex-panel influence matrix is empty.');
end

tol = max(size(A))*eps(max(s))*100;
r = sum(s > tol);

% Minimum-norm particular solution of A*g=b.
if r > 0
    gp = V(:,1:r)*((U(:,1:r).'*b)./s(1:r));
else
    gp = zeros(size(A,2),1);
end

% Numerical nullspace.
Z = V(:,r+1:end);
w = L(:);

if ~isempty(Z)
    a = Z.'*w;
    denom = a.'*a;

    if denom > 1e-20
        delta = targetGamma - w.'*gp;
        g = gp + Z*(a*(delta/denom));
        return
    end
end

% Fallback: exact circulation constraint with least-squares impermeability.
K = [A.'*A, w; w.', 0];
rhs = [A.'*b; targetGamma];
sol = K \ rhs;
g = sol(1:size(A,2));

end

%% ========================================================================
function g = constrained_ls_sum(A,b,targetSum)

% Solve min ||A*g-b||_2 subject to sum(g)=targetSum.
% This enforces total circulation without placing an artificial vortex
% in the middle of the airfoil.
n = size(A,2);
C = ones(1,n);

K = [A.'*A, C.'; C, 0];
rhs = [A.'*b; targetSum];

sol = K \ rhs;
g = sol(1:n);

end

%% ========================================================================
function [xb,yb] = naca2412_closed(Nside,c)

m = 0.02;
p = 0.40;
t = 0.12;

x = linspace(0,1,Nside+1);

% Closed trailing edge.
yt = 5*t*( ...
      0.2969*sqrt(x) ...
    - 0.1260*x ...
    - 0.3516*x.^2 ...
    + 0.2843*x.^3 ...
    - 0.1036*x.^4);

yc = zeros(size(x));
dyc = zeros(size(x));

i1 = x < p;
i2 = ~i1;

yc(i1)  = m/p^2*(2*p*x(i1)-x(i1).^2);
dyc(i1) = 2*m/p^2*(p-x(i1));

yc(i2)  = m/(1-p)^2*((1-2*p)+2*p*x(i2)-x(i2).^2);
dyc(i2) = 2*m/(1-p)^2*(p-x(i2));

theta = atan(dyc);

xu = x - yt.*sin(theta);
yu = yc + yt.*cos(theta);

xl = x + yt.*sin(theta);
yl = yc - yt.*cos(theta);

% Counterclockwise boundary:
% TE -> upper -> LE -> lower -> TE
xb = c*[fliplr(xu) xl(2:end)];
yb = c*[fliplr(yu) yl(2:end)];

end

%% ========================================================================
function [xr,yr] = rotate_xy(x,y,theta,xp,yp)

xr = xp + cos(theta)*(x-xp) - sin(theta)*(y-yp);
yr = yp + sin(theta)*(x-xp) + cos(theta)*(y-yp);

end
