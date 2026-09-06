function sweep = conditional_step_sweep(p, opts)
%CONDITIONAL_STEP_SWEEP Advance frequency only after cycle convergence.
%   Each plateau uses two five-cycle windows, 2% relative tolerance, three
%   consecutive passes, and a 20--150 cycle residence limit by default.
%   Frequency changes use a short analytic half-cosine transition so that
%   fixed-h operation does not introduce displacement/velocity impulses.

defaults=struct('omega',linspace(8,24,51),'control',"fixed_A", ...
    'A',NaN,'h',NaN,'a0',NaN,'phase0',0,'min_cycles',20, ...
    'max_cycles',150,'window_cycles',5,'tolerance',0.02, ...
    'consecutive_passes',3,'transition_cycles',2,'mode',"quick", ...
    'keep_segments',true);
opts=lato.scan.options(opts,defaults);
opts.control=string(opts.control); opts.mode=string(opts.mode);
mustBeMember(opts.control,["fixed_A","fixed_h","fixed_a0"]);
mustBeMember(opts.mode,["quick","formal"]);
if any(opts.omega<=0), error('lato:scan:InvalidSweep','Frequencies must be positive.'); end

if opts.mode == "formal" && numel(opts.omega) < 50
    error('lato:scan:FormalResolution', ...
        'Formal conditional sweeps require at least 50 frequency levels.');
end
if opts.max_cycles < opts.min_cycles
    error('lato:scan:InvalidCycles','max_cycles must be >= min_cycles.');
end

g = get_value(p,'g',9.80665);
A0 = opts.A;
if isnan(A0), A0 = get_drive_value(p,'amplitude',NaN); end
if isnan(A0), A0 = get_drive_value(p,'A',0.02); end
a0 = opts.a0;
if isnan(a0)
    if ~isnan(opts.h), a0=opts.h*g;
    else, a0=A0*opts.omega(1)^2;
    end
end

pk = p;
phase0 = opts.phase0;
globalTime = 0;
globalCycle = 0;
segments = {};
allCycles = table();
levels = repmat(level_template(),numel(opts.omega),1);
nCompleted = 0;
termination = "completed";

for level = 1:numel(opts.omega)
    omega = opts.omega(level);
    Alevel = amplitude_for(omega,opts.control,A0,a0);
    levelCycles = table();
    residence = 0;
    stable = false;

    while residence < opts.max_cycles
        if residence < opts.min_cycles
            requestedCycles = opts.min_cycles-residence;
        else
            requestedCycles = 1;
        end
        requestedCycles = min(requestedCycles,opts.max_cycles-residence);
        duration = requestedCycles*2*pi/omega;
        protocol = make_protocol(omega,omega,duration,phase0, ...
            opts.control,Alevel,a0);
        [result, pk] = run_segment(pk,protocol,duration);
        c = lato.scan.cycle_trace(result);
        if isempty(c)
            termination = "empty_cycle_trace";
            break
        end
        c.level = repmat(level,height(c),1);
        c.level_cycle = (residence+(1:height(c))).';
        c.global_cycle = globalCycle+(1:height(c)).';
        c.start_time = c.start_time+globalTime;
        c.end_time = c.end_time+globalTime;
        c.segment_type = repmat("plateau",height(c),1);
        c.steady = false(height(c),1);
        levelCycles = [levelCycles;c]; %#ok<AGROW>
        allCycles = [allCycles;c]; %#ok<AGROW>
        residence = height(levelCycles);
        globalCycle = globalCycle+height(c);
        if opts.keep_segments
            segments{end+1} = segment_record(result,globalTime,"plateau",level); %#ok<AGROW>
        end
        globalTime = globalTime+duration;
        phase0 = phase0+omega*duration;

        if isfield(result,'termination') && ...
                ~result.termination.completed_requested_interval
            if result.termination.reason=="angle_limit"
                termination="out_of_scope";
            else
                termination="solver_terminated_early";
            end
            break
        end

        stable = converged(levelCycles.amplitude_deg,opts);
        if stable || any(levelCycles.out_of_scope), break; end
    end

    levels(level).omega = omega;
    levels(level).A = Alevel;
    levels(level).h = Alevel*omega^2/g;
    levels(level).cycles = residence;
    levels(level).steady = stable;
    levels(level).steady_amplitude_deg = NaN;
    if stable
        levels(level).steady_amplitude_deg = median( ...
            levelCycles.amplitude_deg(max(1,end-opts.window_cycles+1):end),'omitnan');
        allCycles.steady(end) = true;
        levels(level).status = "steady";
    else
        levels(level).status = "nonsteady_timeout";
    end
    levels(level).max_amplitude_deg = max(levelCycles.max_angle_deg,[],'omitnan');
    levels(level).amplitude_class = lato.scan.classify_amplitude( ...
        levels(level).max_amplitude_deg);
    if any(levelCycles.out_of_scope) || termination=="out_of_scope"
        levels(level).status = "out_of_scope";
        termination = "out_of_scope";
    elseif termination=="solver_terminated_early"
        levels(level).status=termination;
    elseif termination == "empty_cycle_trace"
        levels(level).status = termination;
    elseif ~stable
        termination = "nonsteady_timeout";
    end
    nCompleted = level;

    if termination ~= "completed", break; end
    if level == numel(opts.omega), continue; end

    nextOmega = opts.omega(level+1);
    if opts.transition_cycles > 0
        meanOmega = 0.5*(omega+nextOmega);
        duration = opts.transition_cycles*2*pi/meanOmega;
        protocol = make_protocol(omega,nextOmega,duration,phase0, ...
            opts.control,Alevel,a0);
        [result,pk] = run_segment(pk,protocol,duration);
        c = lato.scan.cycle_trace(result);
        if ~isempty(c)
            c.level = repmat(level,height(c),1);
            c.level_cycle = nan(height(c),1);
            c.global_cycle = globalCycle+(1:height(c)).';
            c.start_time = c.start_time+globalTime;
            c.end_time = c.end_time+globalTime;
            c.segment_type = repmat("transition",height(c),1);
            c.steady = false(height(c),1);
            allCycles = [allCycles;c]; %#ok<AGROW>
            globalCycle = globalCycle+height(c);
        end
        if opts.keep_segments
            segments{end+1} = segment_record(result,globalTime,"transition",level); %#ok<AGROW>
        end
        globalTime = globalTime+duration;
        phase0 = phase0+meanOmega*duration;
    end
end

sweep = struct();
sweep.kind = "conditional_step_sweep";
sweep.control = opts.control;
sweep.direction = direction_name(opts.omega(1),opts.omega(end));
sweep.options = opts;
sweep.levels = struct2table(levels(1:nCompleted));
sweep.cycle = allCycles;
sweep.termination = termination;
if opts.keep_segments, sweep.segments = segments; end
end

function [result,p] = run_segment(p,protocol,duration)
p.drive.function = @lato.scan.sweep_drive;
p.drive.sweep = protocol;
p.drive.A = protocol.A;
p.drive.amplitude = protocol.A;
p.drive.omega = protocol.omega_start;
p.solver.t_final = duration;
p.solver.t_end = p.solver.t_start + duration;
result = lato.simulate_hybrid(p);
if isfield(result,'final_state')
    p.initial_state = result.final_state;
else
    error('lato:scan:MissingFinalState', ...
        'simulate_hybrid must return result.final_state for path continuation.');
end
end

function protocol = make_protocol(w1,w2,duration,phase0,control,A,a0)
protocol = struct('omega_start',w1,'omega_end',w2, ...
    'duration',duration,'t_start',0,'phase0',phase0, ...
    'profile',"half_cosine",'control',control,'A',A,'a0',a0);
end

function tf = converged(amplitude,opts)
cfg.stability.window_cycles=opts.window_cycles;
cfg.stability.relative_tolerance=opts.tolerance;
cfg.stability.required_consecutive_passes=opts.consecutive_passes;
cfg.stability.minimum_cycles=opts.min_cycles;
cfg.stability.maximum_cycles=opts.max_cycles;
assessment=lato.assess_stability(amplitude(:),cfg);
tf=assessment.is_stable;
end

function A = amplitude_for(omega,control,A0,a0)
if control == "fixed_A", A=A0; else, A=a0/omega^2; end
end

function s = segment_record(result,t0,type,level)
s = struct('result',result,'time_offset',t0,'type',type,'level',level);
end

function value = get_value(s,name,fallback)
value=fallback; if isfield(s,name), value=s.(name); end
end

function value = get_drive_value(p,name,fallback)
value=fallback;
if isfield(p,'drive') && isfield(p.drive,name), value=p.drive.(name); end
end

function name = direction_name(a,b)
if b>a, name="up"; elseif b<a, name="down"; else, name="hold"; end
end

function s = level_template()
s=struct('omega',NaN,'A',NaN,'h',NaN,'cycles',0,'steady',false, ...
    'steady_amplitude_deg',NaN,'max_amplitude_deg',NaN, ...
    'amplitude_class',"unavailable",'status',"not_run");
end
