function result = startup_phase_diagram(p, options)
%STARTUP_PHASE_DIAGRAM Compute the unique small-angle onset parameter map.
%   RESULT = lato.theory.startup_phase_diagram(P) scans
%   x = Omega/omega0 and h = a0/g.  Every pixel is a one-period ODE45
%   Floquet calculation; the principal numerical zero boundary is refined
%   independently of the plot grid.

arguments
    p (1, 1) struct
    options.XValues (1, :) double {mustBeFinite, mustBePositive} = linspace(1.4, 2.6, 61)
    options.HValues (1, :) double {mustBeFinite, mustBeNonnegative} = linspace(0, 0.8, 61)
    options.RelTol (1, 1) double {mustBeFinite, mustBePositive} = 2e-9
    options.AbsTol (1, 1) double {mustBeFinite, mustBePositive} = 2e-11
    options.MaxStepFraction (1, 1) double {mustBeFinite, mustBePositive} = 1/64
    options.StabilityToleranceNormalized (1, 1) double ...
        {mustBeFinite, mustBeNonnegative} = 2e-8
    options.RefineBoundary (1, 1) logical = true
    options.BoundaryTolH (1, 1) double {mustBeFinite, mustBePositive} = 2e-6
    options.ShowProgress (1, 1) logical = false
end

if numel(options.XValues) < 50 || numel(options.HValues) < 50
    error('lato:theory:InsufficientGrid', ...
        'The onset map requires at least 50 points in both dimensions.');
end
if any(diff(options.XValues) <= 0) || any(diff(options.HValues) <= 0)
    error('lato:theory:NonmonotonicGrid', ...
        'XValues and HValues must be strictly increasing.');
end

s = lato.theory.startup_parameters(p);
x_values = double(options.XValues(:).');
h_values = double(options.HValues(:).');
nx = numel(x_values);
nh = numel(h_values);

mu_s_inv = nan(nh, nx);
mu_normalized = nan(nh, nx);
spectral_radius = nan(nh, nx);
rho_abs_min = nan(nh, nx);
rho_abs_max = nan(nh, nx);
determinant_error = nan(nh, nx);
stability_code = zeros(nh, nx, 'int8');
accepted_time_points = zeros(nh, nx, 'uint16');

start_time = tic;
for ix = 1:nx
    omega = x_values(ix) * s.omega0_rad_s;
    for ih = 1:nh
        point = lato.theory.floquet_growth_rate(p, omega, h_values(ih), ...
            'RelTol', options.RelTol, ...
            'AbsTol', options.AbsTol, ...
            'MaxStepFraction', options.MaxStepFraction, ...
            'StabilityToleranceNormalized', options.StabilityToleranceNormalized);
        multiplier_abs = sort(abs(point.multipliers));
        mu_s_inv(ih, ix) = point.mu_max_s_inv;
        mu_normalized(ih, ix) = point.mu_max_normalized;
        spectral_radius(ih, ix) = point.spectral_radius;
        rho_abs_min(ih, ix) = multiplier_abs(1);
        rho_abs_max(ih, ix) = multiplier_abs(end);
        determinant_error(ih, ix) = point.determinant_relative_error;
        stability_code(ih, ix) = point.stability_code;
        accepted_time_points(ih, ix) = point.solver.accepted_time_points;
    end
    if options.ShowProgress && (ix == 1 || mod(ix, 5) == 0 || ix == nx)
        fprintf('Floquet onset map: %d/%d frequency columns (%.1f%%)\n', ...
            ix, nx, 100 * ix / nx);
    end
end
grid_elapsed_s = toc(start_time);

theory = lato.theory.averaged_startup_boundary(p, x_values, h_values);
if options.RefineBoundary
    boundary = lato.theory.refine_instability_boundary(p, x_values, ...
        h_values, mu_normalized, ...
        'RelTol', options.RelTol, 'AbsTol', options.AbsTol, ...
        'MaxStepFraction', options.MaxStepFraction, ...
        'TolH', options.BoundaryTolH);
else
    boundary = struct('x', x_values, 'h', nan(size(x_values)), ...
        'num_refined', 0);
end

[X, H] = meshgrid(x_values, h_values);
omega_grid = X * s.omega0_rad_s;
a0_grid = H * s.g_m_s2;
amplitude_grid = a0_grid ./ omega_grid.^2;

result.model = "linear small-angle taut-string Floquet model";
result.x = x_values;
result.h = h_values;
result.X = X;
result.H = H;
result.omega_rad_s = omega_grid;
result.a0_m_s2 = a0_grid;
result.equivalent_amplitude_m = amplitude_grid;
result.mu_max_s_inv = mu_s_inv;
result.mu_max_normalized = mu_normalized;
result.spectral_radius = spectral_radius;
result.rho_abs_min = rho_abs_min;
result.rho_abs_max = rho_abs_max;
result.stability_code = stability_code;
result.determinant_relative_error = determinant_error;
result.accepted_time_points = accepted_time_points;
result.theory = theory;
result.numeric_boundary = boundary;
result.parameters = s;
result.current_point = struct('x', s.current_x, 'h', s.current_h);
result.solver = struct('name', 'ode45', ...
    'RelTol', options.RelTol, 'AbsTol', options.AbsTol, ...
    'MaxStepFraction', options.MaxStepFraction, ...
    'StabilityToleranceNormalized', options.StabilityToleranceNormalized, ...
    'BoundaryTolH', options.BoundaryTolH);
result.runtime = struct('grid_elapsed_s', grid_elapsed_s, ...
    'grid_size', [nh, nx], 'computed_at', datetime('now'));

result.grid_table = table(X(:), H(:), omega_grid(:), a0_grid(:), ...
    amplitude_grid(:), mu_s_inv(:), mu_normalized(:), ...
    spectral_radius(:), rho_abs_min(:), rho_abs_max(:), ...
    stability_code(:), determinant_error(:), double(accepted_time_points(:)), ...
    'VariableNames', {'omega_over_omega0', 'a0_over_g', ...
    'omega_rad_s', 'a0_m_s2', 'equivalent_amplitude_m', ...
    'mu_max_s_inv', 'mu_max_over_omega0', 'spectral_radius', ...
    'rho_abs_min', 'rho_abs_max', 'stability_code', ...
    'determinant_relative_error', 'accepted_time_points'});

end
