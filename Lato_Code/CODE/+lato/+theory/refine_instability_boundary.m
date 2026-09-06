function boundary = refine_instability_boundary(p, x_values, h_values, ...
    mu_normalized, options)
%REFINE_INSTABILITY_BOUNDARY Refine the first principal-tongue zero in h.
%   Each frequency column is bracketed with the supplied grid, then the
%   zero of the numerical Floquet exponent is refined with FZERO.  NaN
%   means that no onset boundary was bracketed inside the requested h span.

arguments
    p (1, 1) struct
    x_values (1, :) double {mustBeFinite, mustBePositive}
    h_values (1, :) double {mustBeFinite, mustBeNonnegative}
    mu_normalized (:, :) double {mustBeFinite}
    options.RelTol (1, 1) double {mustBeFinite, mustBePositive} = 2e-9
    options.AbsTol (1, 1) double {mustBeFinite, mustBePositive} = 2e-11
    options.MaxStepFraction (1, 1) double {mustBeFinite, mustBePositive} = 1/64
    options.TolH (1, 1) double {mustBeFinite, mustBePositive} = 2e-6
end

if ~isequal(size(mu_normalized), [numel(h_values), numel(x_values)])
    error('lato:theory:GridSizeMismatch', ...
        'mu_normalized must have size numel(h_values)-by-numel(x_values).');
end

s = lato.theory.startup_parameters(p);
h_values = h_values(:);
boundary_h = nan(size(x_values));
bracket_low = nan(size(x_values));
bracket_high = nan(size(x_values));
residual = nan(size(x_values));
status = strings(size(x_values));

ode_options = struct('RelTol', options.RelTol, 'AbsTol', options.AbsTol, ...
    'MaxStepFraction', options.MaxStepFraction, ...
    'StabilityTolerance', 0);
root_options = optimset('Display', 'off', 'TolX', options.TolH);

for ix = 1:numel(x_values)
    column = mu_normalized(:, ix);
    crossing = find(column(1:end-1) <= 0 & column(2:end) >= 0, 1, 'first');
    if isempty(crossing)
        status(ix) = "not bracketed";
        continue;
    end

    h_low = h_values(crossing);
    h_high = h_values(crossing + 1);
    bracket_low(ix) = h_low;
    bracket_high(ix) = h_high;
    omega = x_values(ix) * s.omega0_rad_s;
    objective = @(h) mu_at_h(p, omega, h, ode_options);

    if column(crossing) == 0
        h_root = h_low;
        status(ix) = "grid zero";
    elseif column(crossing + 1) == 0
        h_root = h_high;
        status(ix) = "grid zero";
    else
        try
            h_root = fzero(objective, [h_low, h_high], root_options);
            status(ix) = "refined";
        catch
            h_root = h_low - column(crossing) * (h_high - h_low) / ...
                (column(crossing + 1) - column(crossing));
            status(ix) = "linear fallback";
        end
    end

    boundary_h(ix) = h_root;
    residual(ix) = objective(h_root);
end

boundary.x = x_values;
boundary.h = boundary_h;
boundary.bracket_low_h = bracket_low;
boundary.bracket_high_h = bracket_high;
boundary.residual_mu_normalized = residual;
boundary.status = status;
boundary.tolerance_h = options.TolH;
boundary.num_refined = nnz(isfinite(boundary_h));

end

function value = mu_at_h(p, omega, h, ode_options)
point = lato.theory.floquet_growth_rate(p, omega, h, ...
    'RelTol', ode_options.RelTol, ...
    'AbsTol', ode_options.AbsTol, ...
    'MaxStepFraction', ode_options.MaxStepFraction, ...
    'StabilityTolerance', 0);
value = point.mu_max_normalized;
end
