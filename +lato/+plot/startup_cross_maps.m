function [heatmapFig, surfaceFig] = startup_cross_maps(maps, opts)
%STARTUP_CROSS_MAPS Plot finite-time 2-D maps and four-variable surfaces.
%   Surface coordinates are (Omega, scanned parameter, maximum amplitude),
%   while face color is the finite-time peak-envelope growth rate.  That
%   color is a numerical startup diagnostic, not a Floquet exponent.

arguments
    maps (1,1) struct
    opts (1,1) struct = struct()
end
if ~isfield(opts,'visible'), opts.visible = "on"; end
if ~isfield(opts,'make_surface'), opts.make_surface = true; end
if ~isfield(opts,'L_reference_angle_deg'), opts.L_reference_angle_deg = 0; end
if ~isfield(opts,'quantity_mode'), opts.quantity_mode = "maximum"; end
opts.visible = string(opts.visible);
opts.quantity_mode = string(opts.quantity_mode);
mustBeMember(opts.quantity_mode,["maximum","steady"]);

assert(isfield(maps,'kind') && string(maps.kind)=="startup_cross_maps", ...
    'lato:plot:InvalidStartupMaps','Expected startup_cross_maps output.');
maps = ensure_completion_masks(maps);

if opts.quantity_mode == "maximum"
    heatmapFig = figure('Color','w','Visible',opts.visible, ...
        'Units','centimeters','Position',[2 2 18.2 13.8], ...
        'PaperPositionMode','auto');
    tl = tiledlayout(heatmapFig,2,2, ...
        'TileSpacing','compact','Padding','compact');
    amplitudeLimit = amplitude_colour_limit(maps);
    plot_heatmap(nexttile(tl),maps.A_omega,maps, ...
        1e3*maps.A_omega.parameter_values,'Pivot amplitude A (mm)', ...
        "amplitude",amplitudeLimit,opts);
    plot_heatmap(nexttile(tl),maps.A_omega,maps, ...
        1e3*maps.A_omega.parameter_values,'Pivot amplitude A (mm)', ...
        "growth",amplitudeLimit,opts);
    plot_heatmap(nexttile(tl),maps.theta_seed_omega,maps, ...
        maps.theta_seed_omega.parameter_values,'Initial seed angle (deg)', ...
        "amplitude",amplitudeLimit,opts);
    plot_heatmap(nexttile(tl),maps.L_omega,maps, ...
        maps.L_omega.parameter_values,'Length L (m)', ...
        "amplitude",amplitudeLimit,opts);
else
    heatmapFig = figure('Color','w','Visible',opts.visible, ...
        'Units','centimeters','Position',[2 2 27.3 6.9], ...
        'PaperPositionMode','auto');
    tl = tiledlayout(heatmapFig,1,3, ...
        'TileSpacing','compact','Padding','compact');
    amplitudeLimit = steady_amplitude_colour_limit(maps);
    plot_heatmap(nexttile(tl),maps.A_omega,maps, ...
        1e3*maps.A_omega.parameter_values,'Pivot amplitude A (mm)', ...
        "steady",amplitudeLimit,opts);
    plot_heatmap(nexttile(tl),maps.theta_seed_omega,maps, ...
        maps.theta_seed_omega.parameter_values,'Initial seed angle (deg)', ...
        "steady",amplitudeLimit,opts);
    plot_heatmap(nexttile(tl),maps.L_omega,maps, ...
        maps.L_omega.parameter_values,'Length L (m)', ...
        "steady",amplitudeLimit,opts);
end
lato.plot.paper_style(heatmapFig);
set(findall(heatmapFig,'Type','axes'),'FontSize',10.2,'LineWidth',1.0);
set(findall(heatmapFig,'Type','text'),'FontSize',10.2);

surfaceFig = gobjects(0);
if ~opts.make_surface
    return
end
surfaceFig = figure('Color','w','Visible',opts.visible, ...
    'Units','centimeters','Position',[2 2 18.2 7.5], ...
    'PaperPositionMode','auto');
tl2 = tiledlayout(surfaceFig,1,3,'TileSpacing','compact','Padding','loose');
axesList = gobjects(1,3);
axesList(1) = plot_surface(nexttile(tl2),maps.A_omega,maps, ...
    1e3*maps.A_omega.parameter_values,'A (mm)');
axesList(2) = plot_surface(nexttile(tl2),maps.L_omega,maps, ...
    maps.L_omega.parameter_values,'L (m)');
axesList(3) = plot_surface(nexttile(tl2),maps.theta_seed_omega,maps, ...
    maps.theta_seed_omega.parameter_values,'Seed angle (deg)');
muLimit = growth_colour_limit(maps);
for ax = axesList
    clim(ax,[-muLimit,muLimit]);
end
colormap(surfaceFig,blue_white_red(256));
cb = colorbar(axesList(3));
cb.Label.String = 'Finite-time envelope growth rate (s^{-1})';
lato.plot.paper_style(surfaceFig);
end

function plot_heatmap(ax,pair,maps,yValues,yLabel,quantity,amplitudeLimit,opts)
[X,Y] = meshgrid(pair.omega,yValues);
switch quantity
    case "amplitude"
        z = min(pair.max_amplitude_deg,90);
        z(pair.out_of_scope) = 90;
        draw_nonuniform_field(ax,X,Y,z,pair.valid | pair.out_of_scope);
        cb = colorbar(ax);
        cb.Label.String = 'Maximum half-opening (deg)';
        clim(ax,[0,amplitudeLimit]);
        colormap(ax,turbo(256));
    case "growth"
        z = pair.envelope_growth_rate_per_s;
        muLimit = max(abs(z(pair.valid)),[],'omitnan');
        if isempty(muLimit) || ~isfinite(muLimit) || muLimit==0
            muLimit=1;
        end
        % Out-of-scope points are known samples, not missing data. Folding
        % them into the interpolated field removes NaN cracks and the false
        % blocky white boundary produced by a separate categorical surface.
        z(pair.out_of_scope) = muLimit;
        draw_nonuniform_field(ax,X,Y,z,pair.valid | pair.out_of_scope);
        cb = colorbar(ax);
        cb.Label.String = 'Finite-time envelope rate (s^{-1})';
        clim(ax,[-muLimit,muLimit]);
        colormap(ax,blue_white_red(256));
    case "steady"
        z = min(pair.terminal_steady_amplitude_deg,90);
        valid = pair.valid & pair.terminal_steady & isfinite(z);
        z(pair.out_of_scope) = 90;
        valid = valid | pair.out_of_scope;
        draw_nonuniform_field(ax,X,Y,z,valid);
        cb = colorbar(ax);
        cb.Label.String = 'Steady half-opening (deg)';
        clim(ax,[0,amplitudeLimit]);
        colormap(ax,turbo(256));
end
set(ax,'YDir','normal','Color','white');
hold(ax,'on');
draw_masks(ax,X,Y,pair);
draw_two_to_one(ax,pair,maps,yValues,opts);
draw_floquet_onset(ax,pair,maps);
xlabel(ax,'Drive frequency \Omega (rad s^{-1})');
ylabel(ax,yLabel);
xlim(ax,[min(pair.omega),max(pair.omega)]);
ylim(ax,[min(yValues),max(yValues)]);
box(ax,'on');
end

function handle = draw_nonuniform_field(ax,X,Y,z,valid)
% Use physical X/Y vertex coordinates; imagesc would silently treat the
% locally refined axes as uniform and create a false kink at step changes.
z(~valid) = NaN;
handle = surf(ax,X,Y,zeros(size(z)),z,'EdgeColor','none', ...
    'FaceColor','interp');
view(ax,2);
axis(ax,'tight');
end

function ax = plot_surface(ax,pair,maps,yValues,yLabel)
[X,Y] = meshgrid(pair.omega,yValues);
Z = min(pair.max_amplitude_deg,90);
C = pair.envelope_growth_rate_per_s;
Z(~pair.valid) = NaN;
C(~pair.valid) = NaN;
surf(ax,X,Y,Z,C,'EdgeColor','none','FaceColor','interp');
hold(ax,'on');
draw_surface_scope_points(ax,X,Y,pair);
draw_two_to_one_3d(ax,pair,maps,yValues);
xlabel(ax,'\Omega (rad s^{-1})');
ylabel(ax,yLabel);
zlabel(ax,'Maximum half-opening (deg)');
zUpper = pair_amplitude_limit(pair);
if any(pair.out_of_scope,'all'), zUpper = 90; end
zlim(ax,[0,zUpper]);
if zUpper >= 90
    zticks(ax,[0,45,90]);
else
    zticks(ax,unique([0,0.5*zUpper,zUpper]));
end
view(ax,[-38,28]);
grid(ax,'on'); box(ax,'on');
end

function draw_masks(ax,X,Y,pair)
% Startup maps intentionally retain steady/transient flags in the data but
% do not draw their noisy finite-window boundary over the response field.
% The >=90-degree scope mask is rendered at the 90-degree color ceiling;
% point markers are reserved for incomplete integrations and failures.
if any(pair.incomplete,'all')
    scatter(ax,X(pair.incomplete),Y(pair.incomplete),22, ...
        [0.25 0.25 0.25],'+','LineWidth',1.1);
end
if any(pair.failed,'all')
    scatter(ax,X(pair.failed),Y(pair.failed),24,[0 0 0],'x','LineWidth',1.2);
end
end

function draw_surface_scope_points(ax,X,Y,pair)
if any(pair.out_of_scope,'all')
    scatter3(ax,X(pair.out_of_scope),Y(pair.out_of_scope), ...
        90*ones(sum(pair.out_of_scope,'all'),1),24,[.75 .05 .05], ...
        'x','LineWidth',1.1);
end
if any(pair.incomplete,'all')
    scatter3(ax,X(pair.incomplete),Y(pair.incomplete), ...
        zeros(sum(pair.incomplete,'all'),1),24,[.25 .25 .25], ...
        '+','LineWidth',1.1);
end
if any(pair.failed,'all')
    scatter3(ax,X(pair.failed),Y(pair.failed), ...
        zeros(sum(pair.failed,'all'),1),24,[0 0 0],'x','LineWidth',1.1);
end
end

function draw_two_to_one(ax,pair,maps,yValues,opts)
switch string(pair.parameter)
    case "L"
        omega = nonlinear_damped_drive_frequency(maps.g, ...
            pair.parameter_values,maps.beta,opts.L_reference_angle_deg);
        plot(ax,omega,yValues,'k--','LineWidth',1.35, ...
            'HandleVisibility','off');
    case "theta_seed_deg"
        omega = nonlinear_damped_drive_frequency(maps.g, ...
            maps.options.L_fixed,maps.beta,pair.parameter_values);
        plot(ax,omega,yValues,'k--','LineWidth',1.35, ...
            'HandleVisibility','off');
    otherwise
        omega = nonlinear_damped_drive_frequency(maps.g, ...
            maps.options.L_fixed,maps.beta,0);
        xline(ax,omega,'k--','LineWidth',1.35, ...
            'HandleVisibility','off');
end
end

function draw_floquet_onset(ax,pair,maps)
% Transform the averaged F2 boundary into the physical coordinates of F7.
beta = 0.01;
if isfield(maps,'beta'), beta = maps.beta; end
g = maps.g;
omega = pair.omega;
switch string(pair.parameter)
    case "A"
        L = maps.options.L_fixed;
        omega0 = sqrt(g/L);
        x = omega/omega0;
        hBoundary = 2*sqrt((x-2).^2+(beta/omega0).^2);
        amplitudeBoundaryMm = 1e3*g*hBoundary./omega.^2;
        plot(ax,omega,amplitudeBoundaryMm,'--','Color',[0.80 0.08 0.08], ...
            'LineWidth',1.45, ...
            'HandleVisibility','off');
end
end

function omegaDrive = nonlinear_damped_drive_frequency(g,L,beta,thetaDeg)
% Exact undamped finite-amplitude period (elliptic K), followed by the
% weak-linear-damping frequency correction for theta''+beta*theta'+...=0.
% A closed-form exact period for the damped nonlinear pendulum does not
% exist, so the damping part is explicitly an approximation.
parameter = sin(deg2rad(thetaDeg)/2).^2;
[K,~] = ellipke(parameter);
omegaUndamped = (pi./(2*K)).*sqrt(g./L);
omegaDamped = sqrt(max(omegaUndamped.^2-(beta/2).^2,0));
omegaDrive = 2*omegaDamped;
end

function draw_two_to_one_3d(ax,pair,maps,yValues)
switch string(pair.parameter)
    case "L"
        plot3(ax,2*sqrt(maps.g./pair.parameter_values),yValues, ...
            zeros(size(yValues)),'k--','LineWidth',1.15);
    otherwise
        plot3(ax,2*maps.omega0_reference*ones(size(yValues)),yValues, ...
            zeros(size(yValues)),'k--','LineWidth',1.15);
end
end

function limit = growth_colour_limit(maps)
values = [maps.A_omega.envelope_growth_rate_per_s(:); ...
    maps.L_omega.envelope_growth_rate_per_s(:); ...
    maps.theta_seed_omega.envelope_growth_rate_per_s(:)];
values = abs(values(isfinite(values)));
if isempty(values)
    limit = 1;
else
    limit = max(values);
    if limit == 0, limit = 1; end
end
end

function limit = amplitude_colour_limit(maps)
if has_out_of_scope(maps)
    limit = 90;
    return
end
values = [maps.A_omega.max_amplitude_deg(:); ...
    maps.L_omega.max_amplitude_deg(:); ...
    maps.theta_seed_omega.max_amplitude_deg(:)];
values = values(isfinite(values) & values <= 90);
if isempty(values)
    limit = 90;
else
    limit = min(90,max(5,5*ceil(max(values)/5)));
end
end

function limit = steady_amplitude_colour_limit(maps)
if has_out_of_scope(maps)
    limit = 90;
    return
end
values = [maps.A_omega.terminal_steady_amplitude_deg(:); ...
    maps.L_omega.terminal_steady_amplitude_deg(:); ...
    maps.theta_seed_omega.terminal_steady_amplitude_deg(:)];
values = values(isfinite(values) & values <= 90);
if isempty(values)
    limit = 90;
else
    limit = min(90,max(5,5*ceil(max(values)/5)));
end
end

function tf = has_out_of_scope(maps)
tf = any(maps.A_omega.out_of_scope,'all') || ...
    any(maps.L_omega.out_of_scope,'all') || ...
    any(maps.theta_seed_omega.out_of_scope,'all');
end

function limit = pair_amplitude_limit(pair)
values = pair.max_amplitude_deg(:);
values = values(isfinite(values) & values <= 90);
if isempty(values)
    limit = 90;
else
    limit = min(90,max(5,5*ceil(max(values)/5)));
end
end

function cmap = blue_white_red(n)
x = linspace(0,1,n).';
anchors = [0.05 0.20 0.65; 0.92 0.97 1.00; ...
    1.00 0.96 0.92; 0.70 0.05 0.05];
position = [0,0.47,0.53,1];
cmap = interp1(position,anchors,x,'linear');
end

function maps = ensure_completion_masks(maps)
% Cached quick maps made before completion masks were added remain plottable.
pairs = {'A_omega','L_omega','theta_seed_omega'};
for k = 1:numel(pairs)
    pair = maps.(pairs{k});
    if ~isfield(pair,'incomplete')
        pair.incomplete = false(size(pair.failed));
    end
    if ~isfield(pair,'completed_requested_interval')
        pair.completed_requested_interval = ...
            ~pair.failed & ~pair.out_of_scope & ~pair.incomplete;
    end
    if ~isfield(pair,'terminal_steady')
        pair.terminal_steady = false(size(pair.steady));
    end
    if ~isfield(pair,'terminal_steady_amplitude_deg')
        pair.terminal_steady_amplitude_deg = nan(size(pair.max_amplitude_deg));
    end
    if ~isfield(pair,'terminal_steady_relative_mean_deviation')
        pair.terminal_steady_relative_mean_deviation = ...
            nan(size(pair.max_amplitude_deg));
    end
    pair.valid = pair.completed_requested_interval & ...
        ~pair.out_of_scope & ~pair.failed & isfinite(pair.max_amplitude_deg);
    pair.transient = ~pair.steady & pair.valid;
    maps.(pairs{k}) = pair;
end
end
