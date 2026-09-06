function p = configure_point(p, omega, A, L, tFinal)
%CONFIGURE_POINT Set one scan point using the public parameter layout.
%   The aliases retained here make the scan layer tolerant of early
%   versions of the core while p.drive and p.geometry remain canonical.

arguments
    p (1,1) struct
    omega (1,1) double {mustBeFinite,mustBePositive}
    A (1,1) double {mustBeFinite,mustBeNonnegative}
    L (1,1) double {mustBeFinite,mustBePositive}
    tFinal (1,1) double {mustBeFinite,mustBePositive}
end

p.drive.omega = omega;
p.drive.A = A;
p.drive.amplitude = A;
p.ell = L;
p.solver.t_final = tFinal;

% These aliases are harmless and support old exploratory parameter files.
if isfield(p, 'omega'), p.omega = omega; end
if isfield(p, 'Omega'), p.Omega = omega; end
if isfield(p, 'A'), p.A = A; end
if isfield(p, 'L'), p.L = L; end
if isfield(p, 'geometry'), p.geometry.L = L; end
if isfield(p, 'simulation')
    p.simulation.t_final = tFinal;
    if isfield(p.simulation, 'duration'), p.simulation.duration = tFinal; end
    if isfield(p.simulation, 't_end'), p.simulation.t_end = tFinal; end
end
if isfield(p.solver, 'duration'), p.solver.duration = tFinal; end
if isfield(p.solver, 't_end'), p.solver.t_end = tFinal; end
end
