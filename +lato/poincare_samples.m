function samples = poincare_samples(result)
%POINCARE_SAMPLES Sample the signed relative state once per drive cycle.
t = result.time(:);
theta = 0.5 * (result.theta2(:) - result.theta1(:));
omega = 0.5 * (result.omega2(:) - result.omega1(:));
if isfield(result, 'drive') && isfield(result.drive, 'phase')
    phase = result.drive.phase(:);
else
    phase = result.parameters.drive.omega * t;
end
cycle = floor((phase - phase(1)) / (2*pi));
crossings = find(diff(cycle) > 0) + 1;
samples.time = t(crossings);
samples.theta = theta(crossings);
samples.omega = omega(crossings);
samples.cycle = cycle(crossings);
end
