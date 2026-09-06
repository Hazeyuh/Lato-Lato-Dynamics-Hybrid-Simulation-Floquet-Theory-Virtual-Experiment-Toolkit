function d = arbitrary(t, cfg)
%ARBITRARY Adapter for a differentiable user-defined pivot trajectory.
% Required function handles: displacement, velocity, acceleration.
% Optional phase and instantaneous_omega handles support cycle analysis.
% Supplying analytic derivatives is deliberate: numerical differentiation at
% joins can create non-physical acceleration spikes.

required = {'displacement', 'velocity', 'acceleration'};
for k = 1:numel(required)
    assert(isfield(cfg, required{k}) && isa(cfg.(required{k}), 'function_handle'), ...
        'lato:InvalidDrive', 'arbitrary drive requires cfg.%s.', required{k});
end
d.y = cfg.displacement(t);
d.v = cfg.velocity(t);
d.a = cfg.acceleration(t);
if isfield(cfg, 'phase') && isa(cfg.phase, 'function_handle')
    d.phase = cfg.phase(t);
else
    d.phase = NaN(size(t));
end
if isfield(cfg, 'instantaneous_omega') && ...
        isa(cfg.instantaneous_omega, 'function_handle')
    d.omega = cfg.instantaneous_omega(t);
else
    d.omega = NaN(size(t));
end
if isfield(cfg, 'amplitude')
    d.A = cfg.amplitude + zeros(size(t));
elseif isfield(cfg, 'A')
    d.A = cfg.A + zeros(size(t));
else
    d.A = NaN(size(t));
end

end
