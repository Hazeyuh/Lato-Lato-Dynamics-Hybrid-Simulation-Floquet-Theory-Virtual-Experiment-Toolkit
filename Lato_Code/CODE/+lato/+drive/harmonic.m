function d = harmonic(t, cfg)
%HARMONIC Analytic sinusoidal pivot motion with upward-positive y.
%   y = A*cos(phase), phase = phase0 + omega*(t-time_origin).

tau = t - cfg.time_origin;
phase = cfg.phase0 + cfg.omega .* tau;
d.y = cfg.amplitude .* cos(phase);
d.v = -cfg.amplitude .* cfg.omega .* sin(phase);
d.a = -cfg.amplitude .* cfg.omega.^2 .* cos(phase);
d.phase = phase;
d.omega = cfg.omega + zeros(size(t));
d.A = cfg.amplitude + zeros(size(t));

end
