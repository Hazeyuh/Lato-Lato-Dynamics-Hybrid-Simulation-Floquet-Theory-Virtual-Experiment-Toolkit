function d = sweep_drive(t, drive)
%SWEEP_DRIVE Analytic phase-continuous frequency/amplitude sweep.
%   The half-cosine frequency ramp has zero slope at both endpoints. For
%   fixed_h/fixed_a0 operation, A(t)=a0/omega(t)^2 and its derivatives are
%   included analytically in pivot velocity and acceleration.

s = drive.sweep;
tau = t - s.t_start;
T = s.duration;
delta = s.omega_end - s.omega_start;

omega = zeros(size(t));
omegaDot = zeros(size(t));
omegaDDot = zeros(size(t));
phase = zeros(size(t));

before = tau <= 0;
after = tau >= T;
inside = ~(before | after);

omega(before) = s.omega_start;
phase(before) = s.phase0 + s.omega_start .* tau(before);

u = tau(inside);
switch s.profile
    case "linear"
        omega(inside) = s.omega_start + delta .* u ./ T;
        omegaDot(inside) = delta ./ T;
        phase(inside) = s.phase0 + s.omega_start .* u + 0.5 .* delta .* u.^2 ./ T;
        phaseEnd = s.phase0 + 0.5 * (s.omega_start + s.omega_end) * T;
    otherwise % half_cosine
        omega(inside) = s.omega_start + 0.5 .* delta .* (1-cos(pi.*u./T));
        omegaDot(inside) = 0.5 .* delta .* pi ./ T .* sin(pi.*u./T);
        omegaDDot(inside) = 0.5 .* delta .* pi^2 ./ T^2 .* cos(pi.*u./T);
        phase(inside) = s.phase0 + s.omega_start .* u + ...
            0.5 .* delta .* (u - T./pi .* sin(pi.*u./T));
        phaseEnd = s.phase0 + 0.5 * (s.omega_start + s.omega_end) * T;
end

omega(after) = s.omega_end;
phase(after) = phaseEnd + s.omega_end .* (tau(after)-T);

switch s.control
    case "fixed_A"
        A = s.A .* ones(size(t));
        ADot = zeros(size(t));
        ADDot = zeros(size(t));
    otherwise
        A = s.a0 ./ omega.^2;
        ADot = -2 .* s.a0 .* omega.^-3 .* omegaDot;
        ADDot = 6 .* s.a0 .* omega.^-4 .* omegaDot.^2 - ...
            2 .* s.a0 .* omega.^-3 .* omegaDDot;
end

% The cosine convention matches lato.drive.harmonic.
y = A .* cos(phase);
v = ADot .* cos(phase) - A .* omega .* sin(phase);
a = ADDot .* cos(phase) - 2 .* ADot .* omega .* sin(phase) - ...
    A .* omegaDot .* sin(phase) - A .* omega.^2 .* cos(phase);

d = struct('y',y,'v',v,'a',a,'phase',phase,'omega',omega, ...
    'A',A,'omega_dot',omegaDot);
end
