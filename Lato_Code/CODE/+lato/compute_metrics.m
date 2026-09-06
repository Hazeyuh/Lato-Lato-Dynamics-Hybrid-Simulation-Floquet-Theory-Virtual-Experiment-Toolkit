function result = compute_metrics(result, p)
%COMPUTE_METRICS Add response, per-cycle, stability and scalar summaries.

arguments
    result (1,1) struct
    p (1,1) struct
end

t = result.time(:);
s = result.states;
half_opening = abs(s.theta(:,1)-s.theta(:,2))/2;
midline = 0.5*(s.theta(:,1)+s.theta(:,2));
result.response.half_opening_angle = half_opening;
result.response.midline_angle = midline;

cycle = cycle_metrics(t, half_opening, result.drive);
result.cycle = cycle;
result.stability = lato.assess_stability(cycle.peak_half_opening, p);

duration = max(t(end)-t(1), eps);
collision_count = sum(result.events.type == "ball_collision");
maximum_angle = max(half_opening, [], 'omitnan');
maximum_ball_angle = max(abs(s.theta), [], 'all', 'omitnan');
valid = half_opening < p.analysis.maximum_interpretable_angle;

metrics.overall_amplitude_rad = maximum_angle;
metrics.overall_amplitude_deg = rad2deg(maximum_angle);
metrics.maximum_ball_angle_deg = rad2deg(maximum_ball_angle);
metrics.large_angle_reached = maximum_angle >= p.analysis.large_angle;
metrics.exceeds_analysis_limit = maximum_angle >= ...
    p.analysis.maximum_interpretable_angle || ...
    result.termination.reason == "angle_limit";
metrics.analysis_valid_fraction = mean(valid);
metrics.amplitude_class = "small_or_moderate";
if metrics.large_angle_reached
    metrics.amplitude_class = "large";
end
if metrics.exceeds_analysis_limit
    metrics.amplitude_class = "out_of_scope";
end
metrics.is_stable = result.stability.is_stable;
metrics.stable_cycle = result.stability.stable_cycle;
metrics.steady_amplitude_rad = result.stability.steady_peak;
metrics.steady_amplitude_deg = rad2deg(result.stability.steady_peak);
metrics.half_amplitude_deg = metrics.steady_amplitude_deg;
metrics.completed_cycles = numel(cycle.index);
metrics.collision_count = collision_count;
metrics.collision_rate_hz = collision_count/duration;
metrics.slack_fraction_ball1 = time_fraction(t, s.slack(:,1));
metrics.slack_fraction_ball2 = time_fraction(t, s.slack(:,2));
metrics.slack_fraction = mean([metrics.slack_fraction_ball1, ...
    metrics.slack_fraction_ball2]);
metrics.midline_offset_rms_deg = rad2deg(sqrt(time_average( ...
    t, midline.^2)));
metrics.midline_offset_deg = metrics.midline_offset_rms_deg;
metrics.maximum_reynolds_number = max(result.diagnostics.reynolds, ...
    [], 'all', 'omitnan');
positive_cd = result.diagnostics.drag_coefficient( ...
    result.diagnostics.drag_coefficient > 0);
if isempty(positive_cd)
    metrics.median_drag_coefficient = 0;
else
    metrics.median_drag_coefficient = median(positive_cd, 'omitnan');
end
metrics.mean_energy_J = time_average(t, result.energy.total_mechanical);
metrics.final_energy_change_J = result.energy.total_mechanical(end) - ...
    result.energy.total_mechanical(1);
result.metrics = metrics;

end

function cycle = cycle_metrics(t, response, drive)
phase = drive.phase(:);
omega = drive.omega(:);
A = drive.A(:);
if any(~isfinite(phase)) || numel(phase) ~= numel(t)
    cycle = empty_cycle();
    return;
end
phase = unwrap(phase);
base = phase(1);
indices = floor((phase-base)/(2*pi));
candidate = unique(indices(indices >= 0), 'stable');

cycle.index = zeros(0,1);
cycle.start_time = zeros(0,1);
cycle.end_time = zeros(0,1);
cycle.mid_time = zeros(0,1);
cycle.mean_omega = zeros(0,1);
cycle.mean_A = zeros(0,1);
cycle.peak_half_opening = zeros(0,1);
for value = reshape(candidate,1,[])
    mask = indices == value;
    tk = t(mask);
    if numel(tk) < 2
        continue;
    end
    phase_span = phase(mask);
    % Exclude the terminal partial cycle unless its sampled phase span is
    % within one output interval of 2*pi.
    expected_end = base + (value+1)*2*pi;
    if phase(end) < expected_end - max(abs(omega(mask)))* ...
            max(diff(unique(t(mask))))
        continue;
    end
    cycle.index(end+1,1) = value; %#ok<AGROW>
    cycle.start_time(end+1,1) = tk(1); %#ok<AGROW>
    cycle.end_time(end+1,1) = tk(end); %#ok<AGROW>
    cycle.mid_time(end+1,1) = 0.5*(tk(1)+tk(end)); %#ok<AGROW>
    cycle.mean_omega(end+1,1) = time_average(tk, omega(mask)); %#ok<AGROW>
    cycle.mean_A(end+1,1) = time_average(tk, A(mask)); %#ok<AGROW>
    cycle.peak_half_opening(end+1,1) = max(response(mask)); %#ok<AGROW>
end
end

function cycle = empty_cycle()
names = {'index','start_time','end_time','mid_time','mean_omega', ...
    'mean_A','peak_half_opening'};
for k = 1:numel(names)
    cycle.(names{k}) = zeros(0,1);
end
end

function value = time_fraction(t, indicator)
value = time_average(t, double(indicator));
end

function value = time_average(t, x)
if numel(t) < 2 || t(end) <= t(1)
    value = mean(x, 'omitnan');
else
    [tu, keep] = unique(t, 'stable');
    xu = x(keep);
    value = trapz(tu, xu)/(tu(end)-tu(1));
end
end
