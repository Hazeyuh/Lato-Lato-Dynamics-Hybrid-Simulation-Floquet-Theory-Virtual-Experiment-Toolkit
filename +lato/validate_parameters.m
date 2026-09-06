function p = validate_parameters(p)
%VALIDATE_PARAMETERS Validate and normalize the public simulator contract.
%   P = lato.validate_parameters(P) checks physical and numerical inputs,
%   fills optional fields, and supplies a few read-only compatibility aliases.

arguments
    p (1,1) struct
end

required_positive = {'g', 'mass', 'radius', 'ell'};
for k = 1:numel(required_positive)
    name = required_positive{k};
    assert(isfield(p, name) && isscalar(p.(name)) && ...
        isfinite(p.(name)) && p.(name) > 0, ...
        'lato:InvalidParameter', '%s must be a positive finite scalar.', name);
end
assert(p.radius < p.ell, 'lato:InvalidGeometry', ...
    'radius must be smaller than ell.');

p = default_field(p, 'beta', 0.010);
p = default_field(p, 'e_ball', 0.822);
p = default_field(p, 'e_restring', 0.90);
assert(isfinite(p.beta) && p.beta >= 0, 'lato:InvalidParameter', ...
    'beta must be finite and non-negative.');
validate_restitution(p.e_ball, 'e_ball');
validate_restitution(p.e_restring, 'e_restring');

p = default_field(p, 'contact_angle', asin(p.radius / p.ell));
p = default_field(p, 'initial_opening_offset', deg2rad(2.0));

assert(isfield(p, 'drive') && isstruct(p.drive), ...
    'lato:MissingDrive', 'p.drive must be a struct.');
p.drive = default_field(p.drive, 'function', @lato.drive.harmonic);
if isfield(p.drive, 'A') && ~isempty(p.drive.A)
    p.drive.amplitude = p.drive.A;
else
    p.drive = default_field(p.drive, 'amplitude', 0.020);
end
p.drive.A = p.drive.amplitude;
p.drive = default_field(p.drive, 'omega', 15.0);
p.drive = default_field(p.drive, 'phase0', 0.0);
p.drive = default_field(p.drive, 'time_origin', 0.0);
p.drive = default_field(p.drive, 'breakpoints', []);
assert(isa(p.drive.function, 'function_handle'), 'lato:InvalidDrive', ...
    'p.drive.function must be a function handle.');
assert(isempty(p.drive.breakpoints) || ...
    (isvector(p.drive.breakpoints) && all(isfinite(p.drive.breakpoints))), ...
    'lato:InvalidDrive', 'p.drive.breakpoints must be a finite vector.');
p.drive.breakpoints = unique(p.drive.breakpoints(:).');

p = default_field(p, 'air', struct());
p.air = default_field(p.air, 'enabled', true);
p.air = default_field(p.air, 'density', 1.225);
p.air = default_field(p.air, 'dynamic_viscosity', 1.81e-5);
p.air = default_field(p.air, 'model', "schiller-naumann");
p.air = default_field(p.air, 'constant_cd', 0.47);
assert(p.air.density > 0 && p.air.dynamic_viscosity > 0, ...
    'lato:InvalidAirProperties', 'Air properties must be positive.');
assert(any(strcmpi(string(p.air.model), ["schiller-naumann", "constant"])), ...
    'lato:InvalidDragModel', 'Unknown drag model "%s".', p.air.model);

p = default_field(p, 'events', struct());
p.events = default_field(p.events, 'enable_collision', true);
p.events = default_field(p.events, 'enable_slack', true);
p.events = default_field(p.events, 'slack_tolerance', 5e-4);
p.events = default_field(p.events, 'restring_capture_speed', 1e-4);
p.events = default_field(p.events, 'collision_tolerance', 1e-9);
p.events = default_field(p.events, 'collision_capture_speed', 1e-3);
p.events = default_field(p.events, 'maximum_count', 20000);
assert(p.events.slack_tolerance >= 0 && ...
    p.events.collision_tolerance >= 0 && ...
    p.events.collision_capture_speed >= 0, 'lato:InvalidEventTolerance', ...
    'Event tolerances must be non-negative.');

p = default_field(p, 'solver', struct());
p.solver = default_field(p.solver, 't_start', 0.0);
p.solver = default_field(p.solver, 't_end', 20.0);
if isfield(p.solver, 't_final') && ~isempty(p.solver.t_final)
    p.solver.t_end = p.solver.t_start + p.solver.t_final;
end
p.solver.t_final = p.solver.t_end - p.solver.t_start;
p.solver = default_field(p.solver, 'output_dt', 2e-3);
p.solver = default_field(p.solver, 'relative_tolerance', 1e-7);
p.solver = default_field(p.solver, 'absolute_tolerance', 1e-9);
p.solver = default_field(p.solver, 'maximum_step', 1e-2);
p.solver = default_field(p.solver, 'stop_at_angle_limit', true);
p.solver = default_field(p.solver, 'angle_limit', deg2rad(90));
assert(p.solver.t_end > p.solver.t_start, 'lato:InvalidTimeSpan', ...
    'solver.t_end must be greater than solver.t_start.');
assert(all([p.solver.output_dt, p.solver.relative_tolerance, ...
    p.solver.absolute_tolerance, p.solver.maximum_step] > 0), ...
    'lato:InvalidSolverSetting', 'ODE45 tolerances and steps must be positive.');

p = default_field(p, 'stability', struct());
p.stability = default_field(p.stability, 'window_cycles', 5);
p.stability = default_field(p.stability, 'relative_tolerance', 0.02);
p.stability = default_field(p.stability, 'required_consecutive_passes', 3);
p.stability = default_field(p.stability, 'minimum_cycles', 20);
p.stability = default_field(p.stability, 'maximum_cycles', 150);

p = default_field(p, 'analysis', struct());
p.analysis = default_field(p.analysis, 'large_angle', deg2rad(45));
p.analysis = default_field(p.analysis, 'maximum_interpretable_angle', ...
    deg2rad(90));
assert(p.analysis.large_angle < p.analysis.maximum_interpretable_angle, ...
    'lato:InvalidAnalysisLimit', ...
    'The large-angle threshold must be below the interpretation limit.');

if ~isfield(p, 'initial_state') || isempty(p.initial_state)
    p.initial_state.theta = [-1, 1] .* ...
        (p.contact_angle + p.initial_opening_offset);
    p.initial_state.omega = [0, 0];
    p.initial_state.slack = [false, false];
end
p.initial_state = normalise_initial_state(p.initial_state);

% Compatibility aliases used by the audited predecessor library. New code
% should use beta/e_ball and p.drive.* directly.
p.mechanical_damping_beta = p.beta;
p.collision_restitution = p.e_ball;
p.restring_restitution = p.e_restring;
p.drive_amplitude = p.drive.amplitude;
p.drive_omega = p.drive.omega;
p.drive_phase = p.drive.phase0;
p.t_end = p.solver.t_end;

end

function s = normalise_initial_state(s)
if isfield(s, 'ball1') && isfield(s, 'ball2')
    names = {'theta', 'omega', 'slack', 'x', 'y', 'vx', 'vy'};
    for k = 1:numel(names)
        name = names{k};
        s.(name) = [value_or(s.ball1, name, NaN), ...
            value_or(s.ball2, name, NaN)];
    end
end
s = default_field(s, 'theta', [NaN, NaN]);
s = default_field(s, 'omega', [0, 0]);
s = default_field(s, 'slack', [false, false]);
s = default_field(s, 'contact', false);
s = default_field(s, 'x', [NaN, NaN]);
s = default_field(s, 'y', [NaN, NaN]);
s = default_field(s, 'vx', [NaN, NaN]);
s = default_field(s, 'vy', [NaN, NaN]);
names = {'theta', 'omega', 'slack', 'x', 'y', 'vx', 'vy'};
for k = 1:numel(names)
    name = names{k};
    assert(numel(s.(name)) == 2, 'lato:InvalidInitialState', ...
        'initial_state.%s must contain one value per ball.', name);
    s.(name) = reshape(s.(name), 1, 2);
end
s.slack = logical(s.slack);
s.contact = logical(s.contact);
assert(isscalar(s.contact), 'lato:InvalidInitialState', ...
    'initial_state.contact must be scalar.');
assert(~s.contact || ~any(s.slack), 'lato:InvalidInitialState', ...
    'A captured ball-ball contact requires both strings to be taut.');
assert(all(isfinite(s.theta)) && all(isfinite(s.omega)), ...
    'lato:InvalidInitialState', 'Initial theta and omega must be finite.');
for ball = 1:2
    if s.slack(ball)
        assert(all(isfinite([s.x(ball), s.y(ball), s.vx(ball), s.vy(ball)])), ...
            'lato:InvalidInitialState', ...
            'A slack restart requires finite x, y, vx and vy.');
    end
end
end

function value = value_or(s, name, fallback)
if isfield(s, name)
    value = s.(name);
else
    value = fallback;
end
end

function validate_restitution(value, name)
assert(isscalar(value) && isfinite(value) && value >= 0 && value <= 1, ...
    'lato:InvalidRestitution', '%s must lie in [0, 1].', name);
end

function s = default_field(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end
