function make_F2_startup_phase_diagram()
%MAKE_F2_STARTUP_PHASE_DIAGRAM Small-angle Floquet onset phase diagram.
% Edit only the short USER SETTINGS block, then press Run in MATLAB.

script_path = mfilename('fullpath');
project_root = fileparts(fileparts(script_path));
addpath(project_root);
addpath(fullfile(project_root, 'configs'));

%% USER SETTINGS ---------------------------------------------------------
parameter_file = @paper_parameters;   % function in configs/
nx = 121;                             % preview: 61; final: 121
nh = 101;                             % preview: 61; final: 101
x_range = [1.4, 2.6];                 % Omega/omega0
h_range = [0.0, 0.8];                 % h = A*Omega^2/g
show_figure = true;
% ------------------------------------------------------------------------

if nx < 50 || nh < 50
    error('F2 requires at least 50 grid points in each dimension.');
end

p = parameter_file();
x_values = linspace(x_range(1), x_range(2), nx);
h_values = linspace(h_range(1), h_range(2), nh);

fprintf('Computing F2 on a %d-by-%d grid with ode45...\n', nx, nh);
result = lato.theory.startup_phase_diagram(p, ...
    'XValues', x_values, 'HValues', h_values, ...
    'RelTol', 2e-9, 'AbsTol', 2e-11, ...
    'MaxStepFraction', 1/64, ...
    'StabilityToleranceNormalized', 2e-8, ...
    'RefineBoundary', true, 'BoundaryTolH', 2e-6, ...
    'ShowProgress', true);

results_dir = fullfile(project_root, 'results', 'F2_startup_phase_diagram');
figures_dir = fullfile(project_root, 'figures');
if ~isfolder(results_dir), mkdir(results_dir); end
if ~isfolder(figures_dir), mkdir(figures_dir); end

save(fullfile(results_dir, 'F2_startup_phase_diagram_data.mat'), ...
    'result', 'p', '-v7.3');
writetable(result.grid_table, ...
    fullfile(results_dir, 'F2_startup_phase_diagram_grid.csv'));
boundary_table = table(result.numeric_boundary.x(:), ...
    result.numeric_boundary.h(:), ...
    result.theory.h_boundary(:), ...
    result.numeric_boundary.residual_mu_normalized(:), ...
    result.numeric_boundary.status(:), ...
    'VariableNames', {'omega_over_omega0', 'numeric_boundary_a0_over_g', ...
    'averaged_boundary_a0_over_g', 'numeric_residual_mu_over_omega0', ...
    'refinement_status'});
writetable(boundary_table, ...
    fullfile(results_dir, 'F2_startup_phase_diagram_boundaries.csv'));

visibility = "off";
if show_figure, visibility = "on"; end
[fig, ~] = lato.plot.startup_phase_diagram(result, 'Visible', visibility);
exportgraphics(fig, fullfile(figures_dir, 'F2_startup_phase_diagram.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, fullfile(figures_dir, 'F2_startup_phase_diagram.pdf'), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');

fprintf('F2 complete in %.1f s (grid only).\n', result.runtime.grid_elapsed_s);
fprintf('Refined numerical boundary at %d/%d frequency columns.\n', ...
    result.numeric_boundary.num_refined, numel(result.x));
fprintf('Data:   %s\n', results_dir);
fprintf('Figure: %s\n', fullfile(figures_dir, 'F2_startup_phase_diagram.png'));

end
