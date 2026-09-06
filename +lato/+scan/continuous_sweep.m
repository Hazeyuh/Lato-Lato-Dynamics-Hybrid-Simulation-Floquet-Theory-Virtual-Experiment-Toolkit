function sweep = continuous_sweep(p, opts)
%CONTINUOUS_SWEEP Run one phase-continuous up- or down-frequency sweep.
%   Cycle results are process quantities and are never labeled steady.

defaults=struct('omega_start',8,'omega_end',24,'duration',60, ...
    'control',"fixed_A",'A',NaN,'h',NaN,'a0',NaN,'phase0',0, ...
    'profile',"half_cosine");
opts=lato.scan.options(opts,defaults);
opts.control=string(opts.control); opts.profile=string(opts.profile);
mustBeMember(opts.control,["fixed_A","fixed_h","fixed_a0"]);
mustBeMember(opts.profile,["half_cosine","linear"]);
if any([opts.omega_start,opts.omega_end,opts.duration] <= 0)
    error('lato:scan:InvalidSweep','Frequencies and duration must be positive.');
end
if opts.control~="fixed_A" && opts.profile=="linear"
    error('lato:scan:DiscontinuousAmplitudeDerivative', ...
        'fixed_h/fixed_a0 sweeps require the half_cosine profile.');
end

A = opts.A;
if isnan(A), A = get_drive_value(p,'amplitude',NaN); end
if isnan(A), A = get_drive_value(p,'A',0.02); end
g = get_value(p,'g',9.80665);
a0 = opts.a0;
if isnan(a0)
    if ~isnan(opts.h), a0 = opts.h*g;
    else, a0 = A*opts.omega_start^2;
    end
end

pk = p;
pk.drive.function = @lato.scan.sweep_drive;
pk.drive.sweep = struct('omega_start',opts.omega_start, ...
    'omega_end',opts.omega_end,'duration',opts.duration, ...
    't_start',0,'phase0',opts.phase0,'profile',opts.profile, ...
    'control',opts.control,'A',A,'a0',a0);
pk.drive.A = A;
pk.drive.amplitude = A;
pk.drive.omega = opts.omega_start;
pk.solver.t_final = opts.duration;
pk.solver.t_end = pk.solver.t_start + opts.duration;

result = lato.simulate_hybrid(pk);
cycle = lato.scan.cycle_trace(result);
cycle.steady = false(height(cycle),1);
cycle.status = repmat("process_quantity",height(cycle),1);
cycle.status(cycle.out_of_scope) = "out_of_scope";
if isfield(result,'termination') && result.termination.reason=="angle_limit" && ~isempty(cycle)
    cycle.out_of_scope(end)=true;
    cycle.status(end)="out_of_scope";
end

sweep = struct('kind',"continuous_sweep",'control',opts.control, ...
    'direction',direction_name(opts.omega_start,opts.omega_end), ...
    'options',opts,'result',result,'cycle',cycle);
end

function value = get_value(s,name,fallback)
value = fallback;
if isfield(s,name), value = s.(name); end
end

function value = get_drive_value(p,name,fallback)
value = fallback;
if isfield(p,'drive') && isfield(p.drive,name), value = p.drive.(name); end
end

function name = direction_name(a,b)
if b>a, name="up"; elseif b<a, name="down"; else, name="hold"; end
end
