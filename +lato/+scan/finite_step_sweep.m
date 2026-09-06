function sweep = finite_step_sweep(p, opts)
%FINITE_STEP_SWEEP Phase-continuous fixed-residence staircase sweep.
%   Each frequency level is held for a prescribed number of drive cycles.
%   The returned amplitudes are finite-time process quantities, not steady
%   amplitudes.  State and drive phase are continued between levels.

defaults=struct('omega',linspace(8,24,51),'control',"fixed_A", ...
    'A',NaN,'h',NaN,'a0',NaN,'phase0',0,'residence_cycles',12, ...
    'transition_cycles',1,'summary_cycles',5,'mode',"quick", ...
    'keep_segments',false);
opts=lato.scan.options(opts,defaults);
opts.control=string(opts.control); opts.mode=string(opts.mode);
mustBeMember(opts.control,["fixed_A","fixed_h","fixed_a0"]);
mustBeMember(opts.mode,["quick","formal"]);
if any(opts.omega<=0)
    error('lato:scan:InvalidSweep','Frequencies must be positive.');
end
if opts.mode=="formal" && numel(opts.omega)<51
    error('lato:scan:FormalResolution', ...
        'Formal finite-step sweeps require at least 51 frequency levels.');
end
if opts.residence_cycles<1 || opts.summary_cycles<1 || ...
        opts.transition_cycles<0
    error('lato:scan:InvalidCycles','Cycle counts must be positive.');
end

g=get_value(p,'g',9.80665);
A0=opts.A;
if isnan(A0), A0=get_drive_value(p,'amplitude',NaN); end
if isnan(A0), A0=get_drive_value(p,'A',0.02); end
a0=opts.a0;
if isnan(a0)
    if ~isnan(opts.h), a0=opts.h*g;
    else, a0=A0*opts.omega(1)^2;
    end
end

pk=p;
phase0=opts.phase0;
globalTime=0;
globalCycle=0;
allCycles=table();
levels=repmat(level_template(),numel(opts.omega),1);
segments={};
nCompleted=0;
termination="completed";

for level=1:numel(opts.omega)
    omega=opts.omega(level);
    Alevel=amplitude_for(omega,opts.control,A0,a0);
    duration=opts.residence_cycles*2*pi/omega;
    protocol=make_protocol(omega,omega,duration,phase0, ...
        opts.control,Alevel,a0);
    [result,pk]=run_segment(pk,protocol,duration);
    c=lato.scan.cycle_trace(result);
    if isempty(c)
        termination="empty_cycle_trace";
        break
    end
    c.level=repmat(level,height(c),1);
    c.level_cycle=(1:height(c)).';
    c.global_cycle=globalCycle+(1:height(c)).';
    c.start_time=c.start_time+globalTime;
    c.end_time=c.end_time+globalTime;
    c.segment_type=repmat("plateau",height(c),1);
    c.steady=false(height(c),1);
    c.status=repmat("finite_time_process",height(c),1);
    c.status(c.out_of_scope)="out_of_scope";
    allCycles=[allCycles;c]; %#ok<AGROW>
    globalCycle=globalCycle+height(c);
    if opts.keep_segments
        segments{end+1}=segment_record(result,globalTime,"plateau",level); %#ok<AGROW>
    end

    tail=max(1,height(c)-opts.summary_cycles+1):height(c);
    levels(level).omega=omega;
    levels(level).A=Alevel;
    levels(level).h=Alevel*omega^2/g;
    levels(level).cycles=height(c);
    levels(level).process_amplitude_deg=median(c.amplitude_deg(tail),'omitnan');
    levels(level).max_amplitude_deg=max(c.max_angle_deg,[],'omitnan');
    levels(level).amplitude_class=lato.scan.classify_amplitude( ...
        levels(level).max_amplitude_deg);
    levels(level).steady=false;
    levels(level).steady_amplitude_deg=NaN;
    levels(level).status="finite_time_process";

    globalTime=globalTime+duration;
    phase0=phase0+omega*duration;
    nCompleted=level;
    if any(c.out_of_scope) || ...
            (isfield(result,'termination') && result.termination.reason=="angle_limit")
        levels(level).status="out_of_scope";
        termination="out_of_scope";
        break
    elseif isfield(result,'termination') && ...
            ~result.termination.completed_requested_interval
        levels(level).status="solver_terminated_early";
        termination="solver_terminated_early";
        break
    end
    if level==numel(opts.omega), continue; end

    nextOmega=opts.omega(level+1);
    if opts.transition_cycles>0
        meanOmega=0.5*(omega+nextOmega);
        transitionDuration=opts.transition_cycles*2*pi/meanOmega;
        protocol=make_protocol(omega,nextOmega,transitionDuration,phase0, ...
            opts.control,Alevel,a0);
        [result,pk]=run_segment(pk,protocol,transitionDuration);
        ct=lato.scan.cycle_trace(result);
        if ~isempty(ct)
            ct.level=repmat(level,height(ct),1);
            ct.level_cycle=nan(height(ct),1);
            ct.global_cycle=globalCycle+(1:height(ct)).';
            ct.start_time=ct.start_time+globalTime;
            ct.end_time=ct.end_time+globalTime;
            ct.segment_type=repmat("transition",height(ct),1);
            ct.steady=false(height(ct),1);
            ct.status=repmat("transition_process",height(ct),1);
            ct.status(ct.out_of_scope)="out_of_scope";
            allCycles=[allCycles;ct]; %#ok<AGROW>
            globalCycle=globalCycle+height(ct);
        end
        if opts.keep_segments
            segments{end+1}=segment_record(result,globalTime,"transition",level); %#ok<AGROW>
        end
        globalTime=globalTime+transitionDuration;
        phase0=phase0+meanOmega*transitionDuration;
        if isfield(result,'termination') && ...
                ~result.termination.completed_requested_interval
            if result.termination.reason=="angle_limit"
                termination="out_of_scope";
            else
                termination="solver_terminated_early";
            end
            break
        end
    end
end

sweep=struct('kind',"finite_step_sweep",'control',opts.control, ...
    'direction',direction_name(opts.omega(1),opts.omega(end)), ...
    'options',opts,'levels',struct2table(levels(1:nCompleted)), ...
    'cycle',allCycles,'termination',termination);
if opts.keep_segments, sweep.segments=segments; end
end

function [result,p]=run_segment(p,protocol,duration)
p.drive.function=@lato.scan.sweep_drive;
p.drive.sweep=protocol;
p.drive.A=protocol.A;
p.drive.amplitude=protocol.A;
p.drive.omega=protocol.omega_start;
p.solver.t_final=duration;
p.solver.t_end=p.solver.t_start+duration;
result=lato.simulate_hybrid(p);
if ~isfield(result,'final_state')
    error('lato:scan:MissingFinalState', ...
        'simulate_hybrid must return final_state for path continuation.');
end
p.initial_state=result.final_state;
end

function protocol=make_protocol(w1,w2,duration,phase0,control,A,a0)
protocol=struct('omega_start',w1,'omega_end',w2, ...
    'duration',duration,'t_start',0,'phase0',phase0, ...
    'profile',"half_cosine",'control',control,'A',A,'a0',a0);
end

function A=amplitude_for(omega,control,A0,a0)
if control=="fixed_A", A=A0; else, A=a0/omega^2; end
end

function value=get_value(s,name,fallback)
value=fallback; if isfield(s,name), value=s.(name); end
end

function value=get_drive_value(p,name,fallback)
value=fallback;
if isfield(p,'drive') && isfield(p.drive,name), value=p.drive.(name); end
end

function name=direction_name(a,b)
if b>a, name="up"; elseif b<a, name="down"; else, name="hold"; end
end

function s=segment_record(result,t0,type,level)
s=struct('result',result,'time_offset',t0,'type',type,'level',level);
end

function s=level_template()
s=struct('omega',NaN,'A',NaN,'h',NaN,'cycles',0, ...
    'process_amplitude_deg',NaN,'max_amplitude_deg',NaN, ...
    'steady',false,'steady_amplitude_deg',NaN, ...
    'amplitude_class',"unavailable",'status',"not_run");
end
