function analysis = estimate_response_frequency(result, fraction)
%ESTIMATE_RESPONSE_FREQUENCY Estimate the relative-mode angular frequency.
if nargin < 2 || isempty(fraction)
    fraction = 0.5;
end
t = result.time(:);
if isfield(result, 'swing_signed')
    signal = result.swing_signed(:);
else
    signal = 0.5 * (result.theta2(:) - result.theta1(:));
end
start_index = max(1, floor((1 - fraction) * numel(t)));
t_used = t(start_index:end);
x = signal(start_index:end);
x = detrend(x, 0);
dt = median(diff(t_used));
n = numel(x);
window = 0.5 - 0.5 * cos(2 * pi * (0:n-1).' / max(n-1, 1));
spectrum = abs(fft(x .* window));
frequency_hz = (0:n-1).' / (n * dt);
keep = 1:max(2, floor(n/2));
spectrum = spectrum(keep);
frequency_hz = frequency_hz(keep);
if numel(spectrum) > 1
    spectrum(1) = 0;
end
[~, peak_index] = max(spectrum);
analysis.omega_response = 2 * pi * frequency_hz(peak_index);
analysis.frequency_hz = frequency_hz;
analysis.omega_axis = 2 * pi * frequency_hz;
analysis.spectrum = spectrum;
analysis.time_used = t_used;
analysis.signal_used = x;
end
