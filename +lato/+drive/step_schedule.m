function d = step_schedule(t, cfg)
%STEP_SCHEDULE Phase-continuous piecewise-constant frequency schedule.
% cfg.times and cfg.omega_values have equal length; times(1) is the origin.
% amplitude_values may be scalar or one value per interval. Amplitude jumps
% imply displacement jumps and should therefore be avoided in physical runs.

times = cfg.times(:).';
omegas = cfg.omega_values(:).';
assert(numel(times) == numel(omegas) && issorted(times), ...
    'lato:InvalidDrive', 'times and omega_values must be sorted and equal length.');
amplitudes = cfg.amplitude_values;
if isscalar(amplitudes)
    amplitudes = amplitudes + zeros(size(omegas));
else
    amplitudes = amplitudes(:).';
end
assert(numel(amplitudes) == numel(omegas), 'lato:InvalidDrive', ...
    'amplitude_values must be scalar or match omega_values.');

phase0 = cfg.phase0;
phase_at_step = zeros(size(times));
phase_at_step(1) = phase0;
for k = 2:numel(times)
    phase_at_step(k) = phase_at_step(k - 1) + ...
        omegas(k - 1) * (times(k) - times(k - 1));
end

shape = size(t);
tcol = t(:);
y = zeros(size(tcol)); v = y; a = y; phase = y; omega = y;
A_history = y;
for j = 1:numel(tcol)
    idx = find(times <= tcol(j), 1, 'last');
    if isempty(idx), idx = 1; end
    phase(j) = phase_at_step(idx) + omegas(idx) * (tcol(j) - times(idx));
    omega(j) = omegas(idx);
    A_history(j) = amplitudes(idx);
    y(j) = amplitudes(idx) * cos(phase(j));
    v(j) = -amplitudes(idx) * omega(j) * sin(phase(j));
    a(j) = -amplitudes(idx) * omega(j)^2 * cos(phase(j));
end
d.y = reshape(y, shape);
d.v = reshape(v, shape);
d.a = reshape(a, shape);
d.phase = reshape(phase, shape);
d.omega = reshape(omega, shape);
d.A = reshape(A_history, shape);

end
