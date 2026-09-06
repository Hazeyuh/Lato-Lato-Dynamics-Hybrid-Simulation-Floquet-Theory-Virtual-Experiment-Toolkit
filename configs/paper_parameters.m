function p = paper_parameters()
%PAPER_PARAMETERS Editable manuscript parameter set for the Lato simulator.
% Keep this file as an ordinary MATLAB function so that parameters can be
% selected near the top of a run script with:
%   p = paper_parameters();
%   p.drive.omega = 15.5;
%   result = lato.simulate_hybrid(p);

p.name = "paper-default";

%% Physical parameters (SI units)
p.g = 9.81;                     % gravitational acceleration, m/s^2
p.mass = 0.100;                 % mass of each ball, kg
p.radius = 0.012;               % ball radius, m
p.ell = 0.186;                  % pivot-to-ball-center distance, m
p.beta = 0.160;                 % measured amplitude decay rate, 1/s
p.e_ball = 0.822;               % ball-ball normal restitution coefficient
p.e_restring = 0.90;            % radial restitution at string re-tension

%% Initial state: contact angle plus two degrees of opening on each side
p.contact_angle = asin(p.radius / p.ell);
p.initial_opening_offset = deg2rad(2.0);
p.initial_state.theta = [-1, 1] .* ...
    (p.contact_angle + p.initial_opening_offset);
p.initial_state.omega = [0, 0];
p.initial_state.slack = [false, false];
p.initial_state.contact = false; % true only when both taut balls are captured
% Cartesian fields may be set manually for a slack-string restart. NaN asks
% validate_parameters/simulate_hybrid to construct them from theta/omega.
p.initial_state.x = [NaN, NaN];
p.initial_state.y = [NaN, NaN];
p.initial_state.vx = [NaN, NaN];
p.initial_state.vy = [NaN, NaN];

%% Pivot drive
% Every drive function must return a struct containing y, v, a, phase and
% omega. Replace the function handle here to use a chirp, a step schedule,
% or a user-defined differentiable pivot motion.
p.drive.function = @lato.drive.harmonic;
p.drive.amplitude = 0.020;      % one-sided displacement amplitude A, m
p.drive.omega = 15.0;           % angular frequency Omega, rad/s
p.drive.phase0 = 0.0;           % phase at t = time_origin, rad
p.drive.time_origin = 0.0;      % s
p.drive.breakpoints = [];       % known derivative-discontinuity times, s

%% Aerodynamic drag
p.air.enabled = true;
p.air.density = 1.225;          % kg/m^3
p.air.dynamic_viscosity = 1.81e-5; % Pa s
p.air.model = "schiller-naumann";  % or "constant"
p.air.constant_cd = 0.47;

%% Hybrid events and ODE45 controls
p.events.enable_collision = true;
p.events.enable_slack = true;
p.events.slack_tolerance = 5.0e-4; % outward distance used to arm re-tension
p.events.restring_capture_speed = 1.0e-4; % m/s; capture if radial rebound is tiny
p.events.collision_tolerance = 1.0e-9; % m
p.events.collision_capture_speed = 1.0e-3; % m/s; suppress inelastic chatter
p.events.maximum_count = 20000;

p.solver.t_start = 0.0;
p.solver.t_end = 20.0;
p.solver.output_dt = 2.0e-3;    % saved output spacing; not an Euler time step
p.solver.relative_tolerance = 1.0e-7;
p.solver.absolute_tolerance = 1.0e-9;
p.solver.maximum_step = 1.0e-2;
p.solver.stop_at_angle_limit = true;
p.solver.angle_limit = deg2rad(90.0);

%% Cycle-stability and interpretation rules
p.stability.window_cycles = 5;
p.stability.relative_tolerance = 0.02;
p.stability.required_consecutive_passes = 3;
p.stability.minimum_cycles = 20;
p.stability.maximum_cycles = 150;

p.analysis.large_angle = deg2rad(45.0);
p.analysis.maximum_interpretable_angle = deg2rad(90.0);

end
