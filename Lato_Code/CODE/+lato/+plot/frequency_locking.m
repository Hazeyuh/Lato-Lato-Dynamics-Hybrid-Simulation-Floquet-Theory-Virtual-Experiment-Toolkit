function figs = frequency_locking(result, output_dir, prefix)
%FREQUENCY_LOCKING Plot spectrum, phase mismatch and Poincare samples.
if nargin < 3, prefix = "locking"; end
if nargin < 2 || isempty(output_dir), output_dir = pwd; end
if ~isfolder(output_dir), mkdir(output_dir); end
freq = lato.estimate_response_frequency(result);
phase = lato.estimate_phase_locking(result);
samples = lato.poincare_samples(result);

figs.spectrum = figure('Color', 'white', 'Position', [100 100 760 480]);
plot(freq.omega_axis, freq.spectrum / max(freq.spectrum + eps), 'k-', 'LineWidth', 1.2);
xline(result.parameters.drive.omega / 2, '--', '\\Omega/2', 'Color', [0.75 0.20 0.12]);
xline(result.parameters.drive.omega, ':', '\\Omega', 'Color', [0.35 0.35 0.35]);
xlabel('Angular frequency (rad s^{-1})'); ylabel('Normalized amplitude');
title(sprintf('Response spectrum: \\omega_{peak} = %.3f rad s^{-1}', freq.omega_response));
xlim([0, 1.25*result.parameters.drive.omega]);
ylim([0, 1.05]);
lato.plot.paper_style(figs.spectrum);
exportgraphics(figs.spectrum, fullfile(output_dir, prefix + "_spectrum.png"), 'Resolution', 300);
exportgraphics(figs.spectrum, fullfile(output_dir, prefix + "_spectrum.pdf"), 'ContentType', 'vector');

figs.phase = figure('Color', 'white', 'Position', [120 120 900 420]);
layout = tiledlayout(figs.phase, 1, 2, 'TileSpacing', 'compact');
ax1 = nexttile(layout);
plot(ax1, phase.time, phase.psi, '-', 'Color', [0.10 0.35 0.75], 'LineWidth', 0.9);
xlabel(ax1, 'Time (s)'); ylabel(ax1, '\\psi (rad)');
title(ax1, sprintf('2:1 phase mismatch, R = %.3f', phase.psi_resultant));
ax2 = nexttile(layout);
plot(ax2, rad2deg(samples.theta), samples.omega, 'o', 'MarkerSize', 4, ...
    'MarkerFaceColor', [0.75 0.20 0.12], 'MarkerEdgeColor', 'none');
xlabel(ax2, 'Relative angle (deg)'); ylabel(ax2, 'Relative angular velocity (rad s^{-1})');
title(ax2, 'Drive-period Poincare samples');
lato.plot.paper_style(figs.phase);
exportgraphics(figs.phase, fullfile(output_dir, prefix + "_phase_poincare.png"), 'Resolution', 300);
exportgraphics(figs.phase, fullfile(output_dir, prefix + "_phase_poincare.pdf"), 'ContentType', 'vector');
end
