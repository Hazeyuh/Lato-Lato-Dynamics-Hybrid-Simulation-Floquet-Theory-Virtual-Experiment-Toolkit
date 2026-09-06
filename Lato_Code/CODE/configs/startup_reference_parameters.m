function p = startup_reference_parameters()
%STARTUP_REFERENCE_PARAMETERS Collision-suppressed 2:1 reference.
% This parameter file uses the same ODE45 core and drive protocol as the
% full model, but suppresses collision and slack events so that the reduced
% parametric-growth mechanism can be isolated.  Its small seed is below the
% physical contact half-angle, so it is not a geometrically valid two-ball
% pre-collision trajectory and must not be interpreted as one.

p = paper_parameters();
p.name = "startup-reference";

% Small antisymmetric seed.  This is a mechanism reference, not the full
% two-ball release condition used by paper_parameters.m.
p.initial_state.theta = deg2rad([-0.5, 0.5]);
p.initial_state.omega = [0, 0];
p.initial_state.contact = false;
p.initial_state.slack = [false, false];

p.events.enable_collision = false;
p.events.enable_slack = false;

% Main Mathieu tongue: Omega approximately equals 2*omega0.  The modest
% displacement gives visible growth in 20 s without intentionally driving
% the response beyond the 90-degree interpretation limit.
p.drive.amplitude = 0.004;
p.drive.omega = 2*sqrt(p.g/p.ell);
p.drive.phase0 = 0;
p.solver.t_end = 20;

end
