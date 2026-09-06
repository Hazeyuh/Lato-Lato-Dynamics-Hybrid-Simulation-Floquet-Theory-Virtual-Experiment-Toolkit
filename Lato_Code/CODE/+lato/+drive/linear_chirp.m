function d = linear_chirp(t, cfg)
%LINEAR_CHIRP Constant-amplitude, phase-continuous linear frequency sweep.
% Required fields: omega_start, sweep_rate (rad/s^2), amplitude, phase0,
% time_origin. Optional omega_min/omega_max clamp the instantaneous rate.

tau = t - cfg.time_origin;
omega0 = cfg.omega_start;
k = cfg.sweep_rate;

if isfield(cfg, 'omega_min') || isfield(cfg, 'omega_max')
    lower = -inf;
    upper = inf;
    if isfield(cfg, 'omega_min'), lower = cfg.omega_min; end
    if isfield(cfg, 'omega_max'), upper = cfg.omega_max; end
    omega = min(max(omega0 + k .* tau, lower), upper);
    phase = cfg.phase0 + arrayfun(@(x) integrated_clamped_phase( ...
        x, omega0, k, lower, upper), tau);
else
    omega = omega0 + k .* tau;
    phase = cfg.phase0 + omega0 .* tau + 0.5 .* k .* tau.^2;
end

% A constant displacement amplitude gives exact, self-consistent y/v/a.
d.y = cfg.amplitude .* cos(phase);
d.v = -cfg.amplitude .* omega .* sin(phase);
d.a = -cfg.amplitude .* (omega.^2 .* cos(phase) + ...
    k .* active_ramp(omega, cfg) .* sin(phase));
d.phase = phase;
d.omega = omega;
d.A = cfg.amplitude + zeros(size(t));

end

function active = active_ramp(omega, cfg)
active = ones(size(omega));
if isfield(cfg, 'omega_min')
    active(omega <= cfg.omega_min) = 0;
end
if isfield(cfg, 'omega_max')
    active(omega >= cfg.omega_max) = 0;
end
end

function integral = integrated_clamped_phase(tau, omega0, k, lower, upper)
if abs(k) <= eps
    integral = min(max(omega0, lower), upper) * tau;
    return;
end
% Piecewise integration is evaluated by splitting at all clamp crossings.
points = [0, tau];
crossings = [(lower - omega0) / k, (upper - omega0) / k];
lo = min(0, tau);
hi = max(0, tau);
points = [points, crossings(crossings > lo & crossings < hi)]; %#ok<AGROW>
points = sort(points);
integral = 0;
for j = 1:(numel(points) - 1)
    a = points(j);
    b = points(j + 1);
    mid = 0.5 * (a + b);
    if omega0 + k * mid < lower
        contribution = lower * (b - a);
    elseif omega0 + k * mid > upper
        contribution = upper * (b - a);
    else
        contribution = omega0 * (b - a) + 0.5 * k * (b^2 - a^2);
    end
    integral = integral + contribution;
end
if tau < 0
    integral = -integral;
end
end
