function analysis = estimate_phase_locking(result, fraction)
%ESTIMATE_PHASE_LOCKING Estimate 2:1 phase mismatch without toolboxes.
if nargin < 2 || isempty(fraction)
    fraction = 0.5;
end
t = result.time(:);
signal = 0.5 * (result.theta2(:) - result.theta1(:));
start_index = max(1, floor((1 - fraction) * numel(t)));
x = detrend(signal(start_index:end), 0);
analytic_signal = local_analytic_signal(x);
phase_ball = unwrap(angle(analytic_signal));
if isfield(result, 'drive') && isfield(result.drive, 'phase')
    drive_phase = result.drive.phase(start_index:end);
else
    omega = result.parameters.drive.omega;
    drive_phase = omega * t(start_index:end);
end
psi_unwrapped = 2 * phase_ball - drive_phase(:);
psi = mod(psi_unwrapped + pi, 2*pi) - pi;
analysis.time = t(start_index:end);
analysis.phase_ball = phase_ball;
analysis.phase_drive = drive_phase(:);
analysis.psi = psi;
analysis.psi_resultant = abs(mean(exp(1i * psi)));
analysis.psi_circular_std = sqrt(max(0, -2 * log(max(analysis.psi_resultant, eps))));
end

function z = local_analytic_signal(x)
n = numel(x);
X = fft(x);
h = zeros(n, 1);
if rem(n, 2) == 0
    h([1, n/2+1]) = 1;
    h(2:n/2) = 2;
else
    h(1) = 1;
    h(2:(n+1)/2) = 2;
end
z = ifft(X .* h);
end
