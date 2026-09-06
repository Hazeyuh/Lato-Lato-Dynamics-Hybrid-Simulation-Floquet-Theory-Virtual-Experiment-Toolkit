function result = simulate_hybrid(p, drive_function)
%SIMULATE_HYBRID ODE45 hybrid simulation of two driven Lato pendulums.
%   RESULT = lato.simulate_hybrid(P) integrates taut and slack motion with
%   event-localized string release, re-tension, ball collision and the
%   90-degree interpretation limit. The pivot drive is P.drive.function.
%
%   RESULT = lato.simulate_hybrid(P, F) temporarily uses F, where
%   D = F(T, P.drive) returns y, v, a, phase and omega fields.

arguments
    p (1,1) struct
    drive_function = []
end

p = lato.validate_parameters(p);
if ~isempty(drive_function)
    assert(isa(drive_function, 'function_handle'), 'lato:InvalidDrive', ...
        'The second input must be a drive function handle.');
    p.drive.function = drive_function;
end
drive_function = p.drive.function;

t_start = p.solver.t_start;
t_end = p.solver.t_end;
d0 = lato.drive.evaluate(drive_function, t_start, p.drive);
[z, modes] = initialise_state(p.initial_state, p, d0);
contact_active = logical(p.initial_state.contact);

time = t_start;
state_history = z(:).';
mode_history = modes(:).';
contact_history = contact_active;
event_time = zeros(0,1);
event_type = strings(0,1);
event_ball = zeros(0,1);
event_value = zeros(0,1);
termination_reason = "completed";

options = odeset('RelTol', p.solver.relative_tolerance, ...
    'AbsTol', p.solver.absolute_tolerance, ...
    'MaxStep', p.solver.maximum_step);
event_counter = 0;
t_current = t_start;
breakpoints = p.drive.breakpoints;
if isfield(p.drive, 'times')
    breakpoints = unique([breakpoints(:); p.drive.times(:)]).';
end

while t_current < t_end - time_epsilon(t_end)
    if event_counter >= p.events.maximum_count
        termination_reason = "maximum_event_count";
        break;
    end

    future_breaks = breakpoints(breakpoints > t_current + time_epsilon(t_current) & ...
        breakpoints < t_end - time_epsilon(t_end));
    if isempty(future_breaks)
        segment_end = t_end;
    else
        segment_end = future_breaks(1);
    end
    tspan = output_times(t_current, segment_end, ...
        p.solver.output_dt, t_start);
    segment_options = odeset(options, 'Events', ...
        @(t, x) hybrid_events(t, x, modes, contact_active, p, drive_function));
    [ts, zs, te, ze, ie] = ode45( ...
        @(t, x) hybrid_rhs(t, x, modes, contact_active, p, drive_function), ...
        tspan, z, segment_options);

    if numel(ts) > 1
        time = [time; ts(2:end)]; %#ok<AGROW>
        state_history = [state_history; zs(2:end,:)]; %#ok<AGROW>
        mode_history = [mode_history; ...
            repmat(modes(:).', numel(ts)-1, 1)]; %#ok<AGROW>
        contact_history = [contact_history; ...
            repmat(contact_active, numel(ts)-1, 1)]; %#ok<AGROW>
    end
    z = zs(end,:).';

    if isempty(ie)
        t_current = segment_end;
        z = project_redundant_state(z, modes, p, ...
            lato.drive.evaluate(drive_function, t_current, p.drive));
        if segment_end < t_end - time_epsilon(t_end)
            [event_time, event_type, event_ball, event_value] = append_event( ...
                event_time, event_type, event_ball, event_value, ...
                t_current, "drive_breakpoint", 0, NaN);
        end
        continue;
    end

    % ODE45 stops at the earliest terminal event. Process every event index
    % reported at that same time, with the angle limit taking priority.
    t_event = te(end);
    z = ze(end,:).';
    indices = unique(ie(abs(te - t_event) <= 20*time_epsilon(t_event)));
    if any(indices == 7)
        [event_time, event_type, event_ball, event_value] = append_event( ...
            event_time, event_type, event_ball, event_value, ...
            t_event, "angle_limit", 0, p.solver.angle_limit);
        termination_reason = "angle_limit";
        break;
    end

    d_event = lato.drive.evaluate(drive_function, t_event, p.drive);
    for index = reshape(indices, 1, [])
        switch index
            case 1
                [z, closing_speed, contact_active] = resolve_collision( ...
                    z, modes, p, d_event);
                [event_time, event_type, event_ball, event_value] = append_event( ...
                    event_time, event_type, event_ball, event_value, ...
                    t_event, "ball_collision", 0, closing_speed);
            case {2, 3}
                ball = index - 1;
                z = taut_to_slack(z, ball, p, d_event);
                modes(ball) = true;
                contact_active = false;
                [event_time, event_type, event_ball, event_value] = append_event( ...
                    event_time, event_type, event_ball, event_value, ...
                    t_event, "string_release", ball, 0);
            case {4, 5}
                ball = index - 3;
                [z, captured, radial_speed] = resolve_restring( ...
                    z, ball, p, d_event);
                modes(ball) = ~captured;
                label = "string_rebound";
                if captured, label = "string_captured"; end
                [event_time, event_type, event_ball, event_value] = append_event( ...
                    event_time, event_type, event_ball, event_value, ...
                    t_event, label, ball, radial_speed);
            case 6
                contact_active = false;
                [event_time, event_type, event_ball, event_value] = append_event( ...
                    event_time, event_type, event_ball, event_value, ...
                    t_event, "contact_release", 0, 0);
        end
        event_counter = event_counter + 1;
    end
    z = project_redundant_state(z, modes, p, d_event);
    t_current = t_event;

    % Store the post-event state at the same physical time. This explicitly
    % represents velocity jumps without smearing them over an output step.
    time = [time; t_event]; %#ok<AGROW>
    state_history = [state_history; z(:).']; %#ok<AGROW>
    mode_history = [mode_history; modes(:).']; %#ok<AGROW>
    contact_history = [contact_history; contact_active]; %#ok<AGROW>
end

[states, diagnostic] = reconstruct_states( ...
    time, state_history, mode_history, p, drive_function);
states.contact = logical(contact_history);
drive = lato.drive.evaluate(drive_function, time, p.drive);
events = struct('time', event_time, 'type', event_type, ...
    'ball', event_ball, 'value', event_value, ...
    'count', numel(event_time));

result.parameters = p;
result.time = time;
result.drive = drive;
result.states = states;
result.events = events;
result.termination = struct('reason', termination_reason, ...
    'time', time(end), 'completed_requested_interval', ...
    termination_reason == "completed");
result.diagnostics = diagnostic;
result.energy = lato.compute_energy(result, p);
result = lato.compute_metrics(result, p);
result.final_state = make_final_state(states);

% Flat aliases keep predecessor plotting scripts and scan helpers usable.
result.pivot_y = drive.y;
result.pivot_velocity = drive.v;
result.pivot_acceleration = drive.a;
result.theta1 = states.theta(:,1);
result.theta2 = states.theta(:,2);
result.omega1 = states.omega(:,1);
result.omega2 = states.omega(:,2);
result.x1 = states.x(:,1);
result.y1 = states.y(:,1);
result.x2 = states.x(:,2);
result.y2 = states.y(:,2);
result.slack1 = states.slack(:,1);
result.slack2 = states.slack(:,2);
result.reynolds1 = diagnostic.reynolds(:,1);
result.reynolds2 = diagnostic.reynolds(:,2);
result.drag_coefficient1 = diagnostic.drag_coefficient(:,1);
result.drag_coefficient2 = diagnostic.drag_coefficient(:,2);
result.swing_angle = result.response.half_opening_angle;
result.midpoint_angle = result.response.midline_angle;
result.cycle_time = result.cycle.mid_time;
result.cycle_peak = result.cycle.peak_half_opening;

end

function dz = hybrid_rhs(t, z, modes, contact_active, p, drive_function)
d = lato.drive.evaluate(drive_function, t, p.drive);
dz = zeros(12,1);
if contact_active && ~any(modes)
    alpha = zeros(1,2);
    omega_common = 0.5*(z(2)+z(8));
    for ball = 1:2
        idx = (ball-1)*6+(1:6);
        alpha(ball) = taut_acceleration(z(idx(1)), omega_common, d, p);
    end
    alpha_common = mean(alpha);
    dz(1) = omega_common;
    dz(2) = alpha_common;
    dz(7) = omega_common;
    dz(8) = alpha_common;
    return;
end
for ball = 1:2
    offset = (ball - 1) * 6;
    idx = offset + (1:6);
    q = z(idx);
    if ~modes(ball)
        theta = q(1);
        omega = q(2);
        alpha = taut_acceleration(theta, omega, d, p);
        dz(idx(1)) = omega;
        dz(idx(2)) = alpha;
    else
        velocity = q(5:6).';
        force = air_force(velocity, p);
        dz(idx(3)) = q(5);
        dz(idx(4)) = q(6);
        dz(idx(5)) = force(1)/p.mass;
        dz(idx(6)) = force(2)/p.mass - p.g;
    end
end
end

function [value, terminal, direction] = hybrid_events( ...
        t, z, modes, contact_active, p, drive_function)
d = lato.drive.evaluate(drive_function, t, p.drive);
value = ones(7,1);
terminal = ones(7,1);
direction = [-1; -1; -1; 1; 1; 1; 1];

[position, ~, theta, omega] = physical_state(z, modes, p, d);
if p.events.enable_collision && ~contact_active
    value(1) = norm(position(2,:) - position(1,:)) - 2*p.radius;
end
for ball = 1:2
    if ~modes(ball) && p.events.enable_slack
        radial = [sin(theta(ball)), -cos(theta(ball))];
        force = air_force(physical_velocity(theta(ball), omega(ball), d.v, p), p);
        tension = p.mass*(p.ell*omega(ball)^2 + ...
            (p.g+d.a)*cos(theta(ball))) + dot(force, radial);
        value(1 + ball) = tension;
    elseif modes(ball)
        distance = norm(position(ball,:) - [0, d.y]);
        value(3 + ball) = distance - (p.ell + p.events.slack_tolerance);
    end
end
if contact_active && ~any(modes)
    alpha1 = taut_acceleration(theta(1), omega(1), d, p);
    alpha2 = taut_acceleration(theta(2), omega(2), d, p);
    value(6) = alpha2-alpha1;
end
value(7) = max(abs(theta)) - p.solver.angle_limit;
if ~p.solver.stop_at_angle_limit
    value(7) = 1;
end
end

function [z, modes] = initialise_state(s, p, d)
modes = logical(s.slack(:).');
z = zeros(12,1);
for ball = 1:2
    idx = (ball - 1)*6 + (1:6);
    theta = s.theta(ball);
    omega = s.omega(ball);
    if modes(ball)
        z(idx) = [theta; omega; s.x(ball); s.y(ball); ...
            s.vx(ball); s.vy(ball)];
    else
        x = p.ell*sin(theta);
        y = d.y-p.ell*cos(theta);
        vx = p.ell*omega*cos(theta);
        vy = d.v+p.ell*omega*sin(theta);
        z(idx) = [theta; omega; x; y; vx; vy];
    end
end
end

function z = project_redundant_state(z, modes, p, d)
[position, velocity, theta, omega] = physical_state(z, modes, p, d);
for ball = 1:2
    idx = (ball - 1)*6 + (1:6);
    z(idx) = [theta(ball); omega(ball); position(ball,1); ...
        position(ball,2); velocity(ball,1); velocity(ball,2)];
end
end

function z = taut_to_slack(z, ball, p, d)
idx = (ball - 1)*6 + (1:6);
theta = z(idx(1));
omega = z(idx(2));
z(idx(3:6)) = [p.ell*sin(theta); d.y-p.ell*cos(theta); ...
    p.ell*omega*cos(theta); d.v+p.ell*omega*sin(theta)];
end

function [z, captured, radial_speed] = resolve_restring(z, ball, p, d)
idx = (ball - 1)*6 + (1:6);
position = z(idx(3:4)).';
velocity = z(idx(5:6)).';
relative_position = position - [0, d.y];
radial = relative_position / max(norm(relative_position), eps);
tangent = [-radial(2), radial(1)];
relative_velocity = velocity - [0, d.v];
radial_speed = dot(relative_velocity, radial);
tangential_speed = dot(relative_velocity, tangent);
radial_after = -p.e_restring * max(radial_speed, 0);
captured = abs(radial_after) <= p.events.restring_capture_speed;
position = [0, d.y] + p.ell*radial;
theta = atan2(radial(1), -radial(2));
omega = tangential_speed/p.ell;
if captured
    velocity = [0, d.v] + tangential_speed*tangent;
else
    velocity = [0, d.v] + radial_after*radial + tangential_speed*tangent;
end
z(idx) = [theta; omega; position(:); velocity(:)];
end

function [z, closing_speed, captured] = resolve_collision(z, modes, p, d)
[position, velocity, ~, ~] = physical_state(z, modes, p, d);
delta = position(2,:) - position(1,:);
normal = delta/max(norm(delta), eps);
relative_normal = dot(velocity(2,:) - velocity(1,:), normal);
closing_speed = max(0, -relative_normal);
captured = false;
if relative_normal >= 0
    return;
end

inverse_effective = zeros(1,2);
tangent = zeros(2,2);
for ball = 1:2
    if modes(ball)
        inverse_effective(ball) = 1/p.mass;
    else
        theta = z((ball-1)*6+1);
        tangent(ball,:) = [cos(theta), sin(theta)];
        inverse_effective(ball) = dot(normal, tangent(ball,:))^2/p.mass;
    end
end
impulse = -(1+p.e_ball)*relative_normal / max(sum(inverse_effective), eps);
for ball = 1:2
    sign_ball = 2*ball-3; % ball 1 receives -J*n; ball 2 receives +J*n
    idx = (ball-1)*6+(1:6);
    if modes(ball)
        z(idx(5:6)) = z(idx(5:6)) + sign_ball*impulse/p.mass*normal(:);
    else
        delta_tangent_speed = sign_ball*impulse/p.mass * ...
            dot(normal, tangent(ball,:));
        z(idx(2)) = z(idx(2)) + delta_tangent_speed/p.ell;
    end
end
z = project_redundant_state(z, modes, p, d);
if ~any(modes) && p.e_ball*closing_speed <= p.events.collision_capture_speed
    contact_angle = asin(p.radius/p.ell);
    midpoint = 0.5*(z(1)+z(7));
    omega_common = 0.5*(z(2)+z(8));
    z(1) = midpoint-contact_angle;
    z(7) = midpoint+contact_angle;
    z(2) = omega_common;
    z(8) = omega_common;
    captured = true;
    z = project_redundant_state(z, modes, p, d);
end
end

function alpha = taut_acceleration(theta, omega, d, p)
velocity = physical_velocity(theta, omega, d.v, p);
force = air_force(velocity, p);
tangent = [cos(theta), sin(theta)];
alpha = -(p.g+d.a)/p.ell*sin(theta) - 2*p.beta*omega + ...
    dot(force,tangent)/(p.mass*p.ell);
end

function [states, diagnostic] = reconstruct_states( ...
        time, z_history, mode_history, p, drive_function)
n = numel(time);
states.theta = zeros(n,2);
states.omega = zeros(n,2);
states.x = zeros(n,2);
states.y = zeros(n,2);
states.vx = zeros(n,2);
states.vy = zeros(n,2);
states.slack = logical(mode_history);
diagnostic.reynolds = zeros(n,2);
diagnostic.drag_coefficient = zeros(n,2);
diagnostic.tension = NaN(n,2);
for k = 1:n
    d = lato.drive.evaluate(drive_function, time(k), p.drive);
    [position, velocity, theta, omega] = physical_state( ...
        z_history(k,:).', mode_history(k,:), p, d);
    states.theta(k,:) = theta;
    states.omega(k,:) = omega;
    states.x(k,:) = position(:,1).';
    states.y(k,:) = position(:,2).';
    states.vx(k,:) = velocity(:,1).';
    states.vy(k,:) = velocity(:,2).';
    for ball = 1:2
        [force, re, cd] = air_force(velocity(ball,:), p);
        diagnostic.reynolds(k,ball) = re;
        diagnostic.drag_coefficient(k,ball) = cd;
        if ~mode_history(k,ball)
            radial = [sin(theta(ball)), -cos(theta(ball))];
            diagnostic.tension(k,ball) = p.mass*(p.ell*omega(ball)^2 + ...
                (p.g+d.a)*cos(theta(ball))) + dot(force, radial);
        end
    end
end
end

function [position, velocity, theta, omega] = physical_state(z, modes, p, d)
position = zeros(2,2);
velocity = zeros(2,2);
theta = zeros(1,2);
omega = zeros(1,2);
for ball = 1:2
    idx = (ball-1)*6+(1:6);
    q = z(idx);
    if ~modes(ball)
        theta(ball) = q(1);
        omega(ball) = q(2);
        position(ball,:) = [p.ell*sin(q(1)), d.y-p.ell*cos(q(1))];
        velocity(ball,:) = physical_velocity(q(1), q(2), d.v, p);
    else
        position(ball,:) = q(3:4).';
        velocity(ball,:) = q(5:6).';
        relative = position(ball,:) - [0, d.y];
        distance = max(norm(relative), eps);
        radial = relative/distance;
        tangent = [-radial(2), radial(1)];
        theta(ball) = atan2(radial(1), -radial(2));
        omega(ball) = dot(velocity(ball,:) - [0,d.v], tangent)/distance;
    end
end
end

function velocity = physical_velocity(theta, omega, pivot_velocity, p)
velocity = [p.ell*omega*cos(theta), ...
    pivot_velocity+p.ell*omega*sin(theta)];
end

function [force, reynolds, cd] = air_force(velocity, p)
speed = norm(velocity);
if speed <= eps
    force = [0,0]; reynolds = 0; cd = 0; return;
end
reynolds = p.air.density*speed*(2*p.radius)/p.air.dynamic_viscosity;
switch lower(string(p.air.model))
    case "schiller-naumann"
        if reynolds <= 1000
            cd = 24/reynolds*(1+0.15*reynolds^0.687);
        else
            cd = 0.44;
        end
    case "constant"
        cd = p.air.constant_cd;
end
if p.air.enabled
    area = pi*p.radius^2;
    force = -0.5*p.air.density*cd*area*speed*velocity;
else
    force = [0,0];
end
end

function times = output_times(t0, t1, dt, origin)
first_index = ceil((t0-origin)/dt - 10*eps);
last_index = floor((t1-origin)/dt + 10*eps);
grid = origin + (first_index:last_index)*dt;
grid = grid(grid > t0 + time_epsilon(t0) & grid < t1-time_epsilon(t1));
times = unique([t0, grid, t1]);
if numel(times) < 2
    times = [t0, t1];
end
end

function e = time_epsilon(t)
e = max(1e-12, 50*eps(max(1,abs(t))));
end

function [times, types, balls, values] = append_event( ...
        times, types, balls, values, time, type, ball, value)
times(end+1,1) = time;
types(end+1,1) = type;
balls(end+1,1) = ball;
values(end+1,1) = value;
end

function state = make_final_state(states)
names = {'theta','omega','slack','x','y','vx','vy'};
for k = 1:numel(names)
    name = names{k};
    state.(name) = reshape(states.(name)(end,:), 1, 2);
end
state.contact = logical(states.contact(end));
state.ball1 = ball_state(state, 1);
state.ball2 = ball_state(state, 2);
end

function b = ball_state(state, index)
names = {'theta','omega','slack','x','y','vx','vy'};
for k = 1:numel(names)
    b.(names{k}) = state.(names{k})(index);
end
end
