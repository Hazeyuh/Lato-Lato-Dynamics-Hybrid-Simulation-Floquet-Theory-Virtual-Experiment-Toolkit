function [fig, graphics] = startup_phase_diagram(result, options)
%STARTUP_PHASE_DIAGRAM Plot the small-angle Floquet onset parameter map.
%   The stable region is gray.  Color in the unstable region is the
%   normalized numerical growth rate mu_max/omega0.  The manuscript view
%   keeps the averaged-theory boundary dashed and omits the redundant
%   numerical zero line by default because the color mask already shows it.

arguments
    result (1, 1) struct
    options.Visible (1, 1) string = "on"
    options.ShowNumericBoundary (1, 1) logical = false
    options.ShowTheoryBoundary (1, 1) logical = true
end

required = {'x', 'h', 'mu_max_normalized', 'numeric_boundary', 'theory'};
for k = 1:numel(required)
    if ~isfield(result, required{k})
        error('lato:plot:InvalidStartupResult', ...
            'Missing startup-map field: %s.', required{k});
    end
end

fig = figure('Name', 'F2: small-angle onset map', ...
    'NumberTitle', 'off', 'Color', 'white', ...
    'Visible', char(options.Visible), 'Units','centimeters', ...
    'Position', [2, 2, 18.2, 11.8], 'PaperPositionMode','auto');
ax = axes(fig);
hold(ax, 'on');
ax.Color = [0.91, 0.91, 0.91];

mu = result.mu_max_normalized;
unstable = mu > result.solver.StabilityToleranceNormalized;
image_handle = imagesc(ax, result.x, result.h, mu);
image_handle.AlphaData = double(unstable);
image_handle.HandleVisibility = 'off';
set(ax, 'YDir', 'normal');

positive_mu = mu(unstable & isfinite(mu));
if isempty(positive_mu)
    color_limits = [0, 1];
else
    color_limits = [0, max(positive_mu)];
    if color_limits(2) <= eps
        color_limits(2) = 1;
    end
end
clim(ax, color_limits);
colormap(ax, parula(256));
cb = colorbar(ax);
cb.Label.String = '\mu_{max}/\omega_0';
cb.Label.Interpreter = 'tex';

numeric_valid = isfinite(result.numeric_boundary.h);
numeric_line = gobjects(0);
if options.ShowNumericBoundary
    numeric_line = plot(ax, result.numeric_boundary.x(numeric_valid), ...
        result.numeric_boundary.h(numeric_valid), '-', ...
        'Color', [0.05, 0.05, 0.05], 'LineWidth', 1.5, ...
        'DisplayName', 'Numerical boundary');
end
theory_line = gobjects(0);
if options.ShowTheoryBoundary
    theory_line = plot(ax, result.theory.x, result.theory.h_boundary, '--', ...
        'Color', [0.05, 0.05, 0.05], 'LineWidth', 2.1, ...
        'DisplayName', 'Averaged onset boundary');
end

xlim(ax, [min(result.x), max(result.x)]);
ylim(ax, [min(result.h), max(result.h)]);
xlabel(ax, 'Drive-frequency ratio, \Omega/\omega_0', 'Interpreter', 'tex');
ylabel(ax, 'Acceleration modulation, h = a_0/g', 'Interpreter', 'tex');
grid(ax, 'on');
ax.GridAlpha = 0.12;
ax.Layer = 'top';
ax.FontName = 'Times New Roman';
ax.FontSize = 10.5;
ax.LineWidth = 1.0;
box(ax, 'on');

legend_handles = [numeric_line, theory_line];
if ~isempty(legend_handles)
    legend(ax, legend_handles, ...
        'Location', 'southwest', 'Box', 'off', ...
        'FontName', 'Times New Roman','FontSize',9.5);
end
hold(ax, 'off');

graphics = struct('axes', ax, 'image', image_handle, ...
    'colorbar', cb, 'numeric_boundary', numeric_line, ...
    'theory_boundary', theory_line);

end
