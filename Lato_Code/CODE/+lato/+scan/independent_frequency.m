function scan = independent_frequency(p, opts)
%INDEPENDENT_FREQUENCY Independent fixed-parameter frequency scan.
%   SCAN = lato.scan.independent_frequency(P, OPTS) restarts from the same
%   initial condition at every frequency. OPTS.control is fixed_A (default),
%   fixed_h, or fixed_a0. In formal mode at least 50 frequency points are
%   required. Unconverged points retain diagnostics but have NaN steady
%   amplitude.

defaults=struct('omega',linspace(10,20,61),'control',"fixed_A", ...
    'A',NaN,'h',NaN,'a0',NaN,'L',NaN,'min_cycles',20, ...
    'max_cycles',150,'mode',"quick",'keep_runs',false, ...
    'stop_out_of_scope',false);
opts=lato.scan.options(opts,defaults);
opts.control=string(opts.control); opts.mode=string(opts.mode);
mustBeMember(opts.control,["fixed_A","fixed_h","fixed_a0"]);
mustBeMember(opts.mode,["quick","formal"]);

if opts.max_cycles < opts.min_cycles
    error('lato:scan:InvalidCycles', 'max_cycles must be >= min_cycles.');
end
if opts.mode == "formal" && numel(opts.omega) < 50
    error('lato:scan:FormalResolution', ...
        'Formal frequency scans require at least 50 points.');
end

L = opts.L;
if isnan(L), L = get_parameter(p, {'ell','L'}, 0.35); end
baseA = opts.A;
if isnan(baseA), baseA = get_drive_parameter(p, 'amplitude', NaN); end
if isnan(baseA), baseA = get_drive_parameter(p, 'A', 0.02); end
baseA0 = opts.a0;
if isnan(baseA0), baseA0 = baseA * opts.omega(1)^2; end
baseH = opts.h;
if isnan(baseH), baseH = baseA0 / get_parameter(p, {'g'}, 9.80665); end
g = get_parameter(p, {'g'}, 9.80665);

n = numel(opts.omega);
items = repmat(empty_summary(), n, 1);
runs = cell(n, 1);

for k = 1:n
    omega = opts.omega(k);
    switch opts.control
        case "fixed_A"
            A = baseA;
        case "fixed_h"
            A = baseH * g / omega^2;
        case "fixed_a0"
            A = baseA0 / omega^2;
    end
    duration = opts.max_cycles * 2*pi / omega;
    pk = lato.scan.configure_point(p, omega, A, L, duration);
    try
        result = lato.simulate_hybrid(pk);
        items(k) = lato.scan.summarize_result(result, ...
            'Omega', omega, 'A', A, 'L', L);
        if opts.keep_runs, runs{k} = result; end
    catch ME
        items(k) = empty_summary();
        items(k).omega = omega;
        items(k).A = A;
        items(k).L = L;
        items(k).h = A * omega^2 / g;
        items(k).status = "solver_error";
        if opts.keep_runs, runs{k} = ME; end
        warning('lato:scan:PointFailed', ...
            'omega=%g rad/s failed: %s', omega, ME.message);
    end
    if opts.stop_out_of_scope && items(k).out_of_scope
        items = items(1:k);
        runs = runs(1:k);
        break
    end
end

scan = struct();
scan.kind = "independent_frequency";
scan.control = opts.control;
scan.mode = opts.mode;
scan.options = opts;
scan.items = items;
scan.table = lato.scan.summary_table(items);
if opts.keep_runs, scan.runs = runs; end
end

function value = get_parameter(p, names, fallback)
value = fallback;
for k = 1:numel(names)
    if isfield(p, names{k}) && isnumeric(p.(names{k}))
        value = p.(names{k});
        return
    end
end
end

function value = get_drive_parameter(p, name, fallback)
value = fallback;
if isfield(p, 'drive') && isfield(p.drive, name)
    value = p.drive.(name);
end
end

function s = empty_summary()
s = struct('omega',NaN,'A',NaN,'L',NaN,'h',NaN,'steady',false, ...
    'steady_amplitude_deg',NaN,'max_amplitude_deg',NaN, ...
    'amplitude_class',"unavailable",'out_of_scope',false, ...
    'large_amplitude',false,'collision_rate_hz',NaN, ...
    'slack_fraction',NaN,'midline_offset_deg',NaN, ...
    'mean_energy_J',NaN,'status',"not_run");
end
