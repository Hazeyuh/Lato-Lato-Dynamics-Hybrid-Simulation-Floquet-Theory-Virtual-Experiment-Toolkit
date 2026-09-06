function fig = mechanism_timeline(result, output_base)
%MECHANISM_TIMELINE Plot the parameter-modulation mechanism and response.
t = result.time(:);
d = result.drive;
swing_deg = rad2deg(0.5 * abs(result.theta2(:) - result.theta1(:)));
fig = figure('Color', 'white', 'Position', [80 80 1080 780]);
layout = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(layout);
plot(ax1, t, d.y, 'k-', 'LineWidth', 1.1);
ylabel(ax1, 'y_p (m)');
title(ax1, 'Pivot motion');

ax2 = nexttile(layout);
plot(ax2, t, d.a, 'Color', [0.10 0.35 0.75], 'LineWidth', 1.0, ...
    'DisplayName', 'Pivot acceleration');
hold(ax2, 'on');
plot(ax2, t, result.parameters.g + d.a, 'Color', [0.75 0.20 0.12], ...
    'LineWidth', 1.0, 'DisplayName', 'g + y_p''''');
hold(ax2, 'off');
ylabel(ax2, 'm s^{-2}');
legend(ax2, 'Location', 'best');
title(ax2, 'Parametric modulation of effective gravity');

ax3 = nexttile(layout);
plot(ax3, t, swing_deg, 'Color', [0.55 0.15 0.65], 'LineWidth', 1.1);
yline(ax3, 45, '--', 'Large-amplitude threshold', 'Color', [0.30 0.30 0.30]);
yline(ax3, 90, ':', 'Out-of-scope threshold', 'Color', [0.60 0.10 0.10]);
ylabel(ax3, 'Half-opening (deg)');
title(ax3, 'Relative angular response');

ax4 = nexttile(layout);
if isfield(result, 'energy') && isfield(result.energy, 'total_mechanical')
    plot(ax4, t, result.energy.total_mechanical, 'Color', [0.10 0.55 0.35], 'LineWidth', 1.0);
    ylabel(ax4, 'Energy (J)');
else
    stairs(ax4, t, double(result.slack1 | result.slack2), 'k-', 'LineWidth', 1.0);
    ylabel(ax4, 'Slack state');
end
xlabel(ax4, 'Time (s)');
title(ax4, sprintf('Hybrid response (%d ball-ball collisions)', ...
    result.metrics.collision_count));
linkaxes([ax1 ax2 ax3 ax4], 'x');
title(layout, sprintf('Lato-Lato mechanism: A = %.3f m, \\Omega = %.3f rad s^{-1}', ...
    result.parameters.drive.amplitude, result.parameters.drive.omega));
lato.plot.paper_style(fig);
if nargin >= 2 && ~isempty(output_base)
    exportgraphics(fig, output_base + ".png", 'Resolution', 300);
    exportgraphics(fig, output_base + ".pdf", 'ContentType', 'vector');
end
end
