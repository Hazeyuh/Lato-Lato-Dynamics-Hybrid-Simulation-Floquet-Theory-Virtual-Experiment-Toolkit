function summary = summarize_result(result, varargin)
%SUMMARIZE_RESULT Convert one hybrid simulation to scan-safe observables.
%   Unconverged simulations retain diagnostic extrema, but their steady
%   amplitude is NaN. This prevents later plotting or fitting from silently
%   treating a last-window value as a steady state.

ip = inputParser;
ip.addParameter('Omega', NaN, @(x)isnumeric(x) && isscalar(x));
ip.addParameter('A', NaN, @(x)isnumeric(x) && isscalar(x));
ip.addParameter('L', NaN, @(x)isnumeric(x) && isscalar(x));
ip.parse(varargin{:});

theta = extract_theta(result);
thetaDeg = abs(rad2deg(theta));
maxAmplitudeDeg = extract_metric(result, ...
    {'overall_amplitude_deg','maximum_response_deg'},NaN);
if isnan(maxAmplitudeDeg)
    maxAmplitudeDeg = max(thetaDeg, [], 'all', 'omitnan');
end
if isempty(maxAmplitudeDeg), maxAmplitudeDeg = NaN; end

stable = extract_stable(result);
steadyAmplitudeDeg = extract_metric(result, ...
    {'half_amplitude_deg','half_opening_amplitude_deg', ...
     'steady_amplitude_deg','amplitude_deg'}, NaN);

if isnan(steadyAmplitudeDeg) && stable
    steadyAmplitudeDeg = terminal_cycle_amplitude(result, theta);
end
if ~stable
    steadyAmplitudeDeg = NaN;
end

summary = struct();
summary.omega = ip.Results.Omega;
summary.A = ip.Results.A;
summary.L = ip.Results.L;
summary.h = ip.Results.A .* ip.Results.Omega.^2 ./ 9.80665;
summary.steady = logical(stable);
summary.steady_amplitude_deg = steadyAmplitudeDeg;
summary.max_amplitude_deg = maxAmplitudeDeg;
summary.amplitude_class = lato.scan.classify_amplitude(maxAmplitudeDeg);
summary.out_of_scope = maxAmplitudeDeg > 90 || analysis_limit_reached(result);
summary.large_amplitude = maxAmplitudeDeg > 45 && maxAmplitudeDeg <= 90;
summary.collision_rate_hz = extract_metric(result, ...
    {'collision_rate_hz','collision_rate','impact_rate_hz'}, NaN);
summary.slack_fraction = extract_metric(result, ...
    {'slack_fraction','slack_time_fraction'}, NaN);
summary.midline_offset_deg = extract_metric(result, ...
    {'midline_offset_deg','mean_angle_deg'}, NaN);
summary.mean_energy_J = extract_metric(result, ...
    {'mean_energy_J','mean_energy','energy_mean'}, NaN);
summary.status = "nonsteady";
if summary.steady, summary.status = "steady"; end
if summary.out_of_scope, summary.status = "out_of_scope"; end
end

function tf=analysis_limit_reached(result)
tf=false;
if isfield(result,'metrics') && isfield(result.metrics,'exceeds_analysis_limit')
    tf=logical(result.metrics.exceeds_analysis_limit);
elseif isfield(result,'termination') && isfield(result.termination,'reason')
    tf=string(result.termination.reason)=="angle_limit";
end
end

function theta = extract_theta(result)
if isfield(result, 'theta1') && isfield(result, 'theta2')
    theta = [result.theta1(:), result.theta2(:)];
elseif isfield(result, 'states') && ~isempty(result.states)
    x = result.states;
    if size(x, 2) >= 3
        theta = x(:, [1 3]);
    elseif size(x, 2) >= 2
        theta = x(:, 1:2);
    else
        theta = x(:, 1);
    end
else
    theta = NaN;
end
end

function stable = extract_stable(result)
stable = false;
if isfield(result, 'stability') && isstruct(result.stability)
    candidates = {'steady','is_steady','converged','isStable'};
    for k = 1:numel(candidates)
        name = candidates{k};
        if isfield(result.stability, name)
            stable = logical(result.stability.(name));
            return
        end
    end
elseif isfield(result, 'stable')
    stable = logical(result.stable);
end
end

function value = extract_metric(result, names, fallback)
value = fallback;
containers = {};
if isfield(result, 'metrics') && isstruct(result.metrics)
    containers{end+1} = result.metrics; %#ok<AGROW>
end
if isfield(result, 'cycle') && isstruct(result.cycle)
    containers{end+1} = result.cycle; %#ok<AGROW>
end
containers{end+1} = result;
for c = 1:numel(containers)
    s = containers{c};
    for k = 1:numel(names)
        if isfield(s, names{k}) && isnumeric(s.(names{k})) && ~isempty(s.(names{k}))
            v = s.(names{k});
            value = v(end);
            return
        end
    end
end
end

function amplitudeDeg = terminal_cycle_amplitude(result, theta)
amplitudeDeg = NaN;
if ~isfield(result, 'time') || isempty(result.time) || all(isnan(theta), 'all')
    return
end
t = result.time(:);
duration = t(end) - t(1);
if duration <= 0, return; end
mask = t >= t(end) - 0.2 * duration;
tail = theta(mask, :);
if isempty(tail), return; end
perBall = 0.5 * (max(tail, [], 1) - min(tail, [], 1));
amplitudeDeg = max(rad2deg(perBall), [], 'omitnan');
end
