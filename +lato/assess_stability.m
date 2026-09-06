function stability = assess_stability(cycle_peak, p)
%ASSESS_STABILITY Confirm convergence using two rolling cycle windows.
%   The default compares the medians of two consecutive five-cycle windows
%   at 2% relative tolerance for three consecutive tests, with no testing
%   before cycle 20.  Within-window range is retained as a diagnostic but is
%   not a rejection criterion: a bounded non-phase-locked response can have
%   repeatable beating across drive cycles.  A non-converged run never
%   receives a steady-amplitude value.

arguments
    cycle_peak (:,1) double
    p (1,1) struct
end

if isfield(p, 'stability')
    cfg = p.stability;
else
    cfg = p;
end
window = cfg.window_cycles;
tolerance = cfg.relative_tolerance;
required = cfg.required_consecutive_passes;
minimum_cycles = cfg.minimum_cycles;
maximum_cycles = cfg.maximum_cycles;

n = numel(cycle_peak);
pass_flags = false(n,1);
relative_change = NaN(n,1);
relative_range = NaN(n,1);
first_test = max(minimum_cycles, 2*window);
for cycle = first_test:n
    previous = cycle_peak(cycle-2*window+1:cycle-window);
    current = cycle_peak(cycle-window+1:cycle);
    old_level = median(previous, 'omitnan');
    new_level = median(current, 'omitnan');
    scale = max([abs(old_level), abs(new_level), sqrt(eps)]);
    relative_change(cycle) = abs(new_level-old_level)/scale;
    relative_range(cycle) = (max(current)-min(current))/ ...
        max(max(abs(current)), sqrt(eps));
    pass_flags(cycle) = all(isfinite([old_level,new_level])) && ...
        relative_change(cycle) <= tolerance;
end

stable_cycle = NaN;
run_length = 0;
for cycle = first_test:n
    if pass_flags(cycle)
        run_length = run_length+1;
        if run_length >= required
            stable_cycle = cycle;
            break;
        end
    else
        run_length = 0;
    end
end

is_stable = isfinite(stable_cycle);
recent_count = min(window,n);
if recent_count > 0
    recent_peak = median(cycle_peak(end-recent_count+1:end), 'omitnan');
else
    recent_peak = NaN;
end
steady_peak = NaN;
if is_stable
    steady_peak = recent_peak;
end

stability.is_stable = is_stable;
stability.steady = is_stable;
stability.is_steady = is_stable;
stability.converged = is_stable;
stability.isStable = is_stable;
stability.stable_cycle = stable_cycle;
stability.pass_flags = pass_flags;
stability.relative_change = relative_change;
stability.relative_range = relative_range;
stability.steady_peak = steady_peak;
stability.representative_peak = steady_peak;
stability.recent_peak_diagnostic = recent_peak;
stability.cycles_available = n;
stability.reached_maximum_cycles = n >= maximum_cycles && ~is_stable;

end
