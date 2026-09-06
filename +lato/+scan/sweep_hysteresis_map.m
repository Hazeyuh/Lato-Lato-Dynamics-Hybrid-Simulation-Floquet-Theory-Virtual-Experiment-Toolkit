function map = sweep_hysteresis_map(p, opts)
%SWEEP_HYSTERESIS_MAP Compare up/down paths over a control-parameter family.
%   The map quantifies finite-time dynamic lag. It must not be interpreted
%   as quasistatic hysteresis unless independent steady-state tests support
%   that stronger claim.

omega0=sqrt(p.g/p.ell);
defaults=struct('control',"fixed_A", ...
    'control_values',linspace(0.002,0.006,51), ...
    'omega_limits',omega0*[1.6 2.4], ...
    'omega_grid',linspace(1.6,2.4,81)*omega0, ...
    'omega_levels',linspace(1.6,2.4,51)*omega0, ...
    'continuous_duration',80,'residence_cycles',12, ...
    'transition_cycles',1,'summary_cycles',5,'profile',"half_cosine", ...
    'mode',"quick",'keep_sweeps',false,'run_continuous',true, ...
    'run_step',true,'checkpoint_file',"",'resume',true);
opts=lato.scan.options(opts,defaults);
opts.control=string(opts.control); opts.mode=string(opts.mode);
mustBeMember(opts.control,["fixed_A","fixed_h"]);
mustBeMember(opts.mode,["quick","formal"]);
if numel(opts.omega_limits)~=2 || any(opts.omega_limits<=0) || ...
        opts.omega_limits(1)>=opts.omega_limits(2)
    error('lato:scan:InvalidSweep','omega_limits must be increasing.');
end
if opts.mode=="formal" && (numel(opts.control_values)<51 || ...
        numel(opts.omega_grid)<51 || ...
        (opts.run_step && numel(opts.omega_levels)<51))
    error('lato:scan:FormalResolution', ...
        'Formal sweep maps require at least 51 values on every varied axis.');
end
validate_axis(opts.control_values,'control_values');
validate_axis(opts.omega_grid,'omega_grid');
if opts.run_step, validate_axis(opts.omega_levels,'omega_levels'); end

values=opts.control_values(:);
n=numel(values);
ng=numel(opts.omega_grid);
nl=numel(opts.omega_levels);
continuousUp=nan(n,ng); continuousDown=nan(n,ng);
stepUp=nan(n,nl); stepDown=nan(n,nl);
continuousMetric=nan(n,1); continuousSigned=nan(n,1);
stepMetric=nan(n,1); stepSigned=nan(n,1);
continuousCoverage=nan(n,2); stepCoverage=nan(n,2);
continuousStop=nan(n,2); stepStop=nan(n,2);
records=table();
sweeps=struct();
if opts.keep_sweeps
    sweeps.continuous_up=cell(n,1); sweeps.continuous_down=cell(n,1);
    sweeps.step_up=cell(n,1); sweeps.step_down=cell(n,1);
end

signature=checkpoint_signature(opts,values,p);
startIndex=1;
checkpointFile=string(opts.checkpoint_file);
if checkpointFile~="" && logical(opts.resume) && isfile(checkpointFile)
    loaded=load(checkpointFile,'checkpoint');
    if ~isfield(loaded,'checkpoint') || ...
            ~isequaln(loaded.checkpoint.signature,signature)
        error('lato:scan:CheckpointMismatch', ...
            'Checkpoint axes or sweep options do not match this request.');
    end
    checkpoint=loaded.checkpoint;
    continuousUp=checkpoint.continuousUp;
    continuousDown=checkpoint.continuousDown;
    stepUp=checkpoint.stepUp;
    stepDown=checkpoint.stepDown;
    continuousMetric=checkpoint.continuousMetric;
    continuousSigned=checkpoint.continuousSigned;
    stepMetric=checkpoint.stepMetric;
    stepSigned=checkpoint.stepSigned;
    continuousCoverage=checkpoint.continuousCoverage;
    stepCoverage=checkpoint.stepCoverage;
    continuousStop=checkpoint.continuousStop;
    stepStop=checkpoint.stepStop;
    records=checkpoint.records;
    if opts.keep_sweeps && isfield(checkpoint,'sweeps')
        sweeps=checkpoint.sweeps;
    end
    startIndex=checkpoint.last_completed+1;
end

for k=startIndex:n
    value=values(k);
    if opts.run_continuous
        upOpts=continuous_options(opts,value,opts.omega_limits(1), ...
            opts.omega_limits(2));
        downOpts=continuous_options(opts,value,opts.omega_limits(2), ...
            opts.omega_limits(1));
        up=lato.scan.continuous_sweep(p,upOpts);
        down=lato.scan.continuous_sweep(p,downOpts);
        continuousUp(k,:)=continuous_row(up,opts.omega_grid);
        continuousDown(k,:)=continuous_row(down,opts.omega_grid);
        [continuousMetric(k),continuousSigned(k),continuousCoverage(k,:)]= ...
            lag_metrics(continuousUp(k,:),continuousDown(k,:),opts.omega_grid);
        continuousStop(k,:)=[last_omega(up),last_omega(down)];
        records=[records;cycle_records(up.cycle,value,"continuous_up")]; %#ok<AGROW>
        records=[records;cycle_records(down.cycle,value,"continuous_down")]; %#ok<AGROW>
        if opts.keep_sweeps
            sweeps.continuous_up{k}=up; sweeps.continuous_down{k}=down;
        end
    end

    if opts.run_step
        upOpts=step_options(opts,value,opts.omega_levels);
        downOpts=step_options(opts,value,fliplr(opts.omega_levels));
        up=lato.scan.finite_step_sweep(p,upOpts);
        down=lato.scan.finite_step_sweep(p,downOpts);
        stepUp(k,:)=step_row(up,opts.omega_levels);
        stepDown(k,:)=step_row(down,opts.omega_levels);
        [stepMetric(k),stepSigned(k),stepCoverage(k,:)]= ...
            lag_metrics(stepUp(k,:),stepDown(k,:),opts.omega_levels);
        stepStop(k,:)=[last_level_omega(up),last_level_omega(down)];
        records=[records;level_records(up.levels,value,"step_up")]; %#ok<AGROW>
        records=[records;level_records(down.levels,value,"step_down")]; %#ok<AGROW>
        if opts.keep_sweeps
            sweeps.step_up{k}=up; sweeps.step_down{k}=down;
        end
    end
    if checkpointFile~=""
        checkpoint=make_checkpoint(signature,k,continuousUp,continuousDown, ...
            stepUp,stepDown,continuousMetric,continuousSigned,stepMetric, ...
            stepSigned,continuousCoverage,stepCoverage,continuousStop, ...
            stepStop,records,sweeps);
        save_checkpoint(checkpointFile,checkpoint);
    end
end

map=struct();
map.kind="finite_time_sweep_path_map";
map.claim="dynamic_lag_not_quasistatic_hysteresis";
map.control=opts.control;
map.control_values=values;
map.omega0=omega0;
map.omega_grid=opts.omega_grid(:).';
map.omega_levels=opts.omega_levels(:).';
map.continuous=struct('up_deg',continuousUp,'down_deg',continuousDown, ...
    'difference_deg',continuousUp-continuousDown, ...
    'mean_abs_lag_deg',continuousMetric, ...
    'mean_signed_gap_deg',continuousSigned, ...
    'coverage_fraction',continuousCoverage,'stop_omega',continuousStop);
map.step=struct('up_deg',stepUp,'down_deg',stepDown, ...
    'difference_deg',stepUp-stepDown,'mean_abs_lag_deg',stepMetric, ...
    'mean_signed_gap_deg',stepSigned,'coverage_fraction',stepCoverage, ...
    'stop_omega',stepStop);
map.records=records;
map.options=opts;
if opts.keep_sweeps, map.sweeps=sweeps; end
end

function signature=checkpoint_signature(opts,values,p)
validated=lato.validate_parameters(p);
signature=struct('control',opts.control,'control_values',values, ...
    'omega_limits',opts.omega_limits(:),'omega_grid',opts.omega_grid(:), ...
    'omega_levels',opts.omega_levels(:),'continuous_duration', ...
    opts.continuous_duration,'residence_cycles',opts.residence_cycles, ...
    'transition_cycles',opts.transition_cycles,'summary_cycles', ...
    opts.summary_cycles,'profile',string(opts.profile),'mode',opts.mode, ...
    'run_continuous',logical(opts.run_continuous), ...
    'run_step',logical(opts.run_step), ...
    'keep_sweeps',logical(opts.keep_sweeps),'parameters',validated);
end

function checkpoint=make_checkpoint(signature,lastCompleted,continuousUp, ...
        continuousDown,stepUp,stepDown,continuousMetric,continuousSigned, ...
        stepMetric,stepSigned,continuousCoverage,stepCoverage, ...
        continuousStop,stepStop,records,sweeps)
checkpoint=struct('signature',signature,'last_completed',lastCompleted, ...
    'continuousUp',continuousUp,'continuousDown',continuousDown, ...
    'stepUp',stepUp,'stepDown',stepDown, ...
    'continuousMetric',continuousMetric, ...
    'continuousSigned',continuousSigned,'stepMetric',stepMetric, ...
    'stepSigned',stepSigned,'continuousCoverage',continuousCoverage, ...
    'stepCoverage',stepCoverage,'continuousStop',continuousStop, ...
    'stepStop',stepStop,'records',records,'sweeps',sweeps);
end

function save_checkpoint(filename,checkpoint)
folder=fileparts(filename);
if folder~="" && ~isfolder(folder), mkdir(folder); end
temporary=filename+".tmp.mat";
save(temporary,'checkpoint','-v7.3');
[ok,message]=movefile(temporary,filename,'f');
if ~ok
    error('lato:scan:CheckpointWrite','Could not write checkpoint: %s',message);
end
end

function out=continuous_options(opts,value,w1,w2)
out=struct('omega_start',w1,'omega_end',w2, ...
    'duration',opts.continuous_duration,'control',opts.control, ...
    'A',NaN,'h',NaN,'a0',NaN,'phase0',0,'profile',opts.profile);
if opts.control=="fixed_A", out.A=value; else, out.h=value; end
end

function out=step_options(opts,value,omega)
out=struct('omega',omega,'control',opts.control,'A',NaN,'h',NaN, ...
    'a0',NaN,'phase0',0,'residence_cycles',opts.residence_cycles, ...
    'transition_cycles',opts.transition_cycles, ...
    'summary_cycles',opts.summary_cycles,'mode',opts.mode, ...
    'keep_segments',false);
if opts.control=="fixed_A", out.A=value; else, out.h=value; end
end

function validate_axis(values,name)
values=values(:);
if isempty(values) || any(~isfinite(values)) || ...
        numel(unique(values))~=numel(values) || any(diff(values)<=0)
    error('lato:scan:InvalidAxis', ...
        '%s must be finite, unique and strictly increasing.',name);
end
end

function y=continuous_row(sweep,grid)
y=nan(size(grid));
if isfield(sweep,'result') && isfield(sweep.result,'termination')
    term=sweep.result.termination;
    if ~term.completed_requested_interval && term.reason~="angle_limit"
        return
    end
end
y=interpolate_cycles(sweep.cycle,grid);
end

function y=step_row(sweep,grid)
y=nan(size(grid));
if any(sweep.termination==["solver_terminated_early","empty_cycle_trace"])
    return
end
y=align_levels(sweep.levels,grid);
end

function y=interpolate_cycles(c,grid)
y=nan(size(grid));
if isempty(c), return; end
valid=isfinite(c.omega_mean) & isfinite(c.amplitude_deg) & ...
    c.amplitude_deg<=90 & ~c.out_of_scope;
if nnz(valid)<2, return; end
[x,ia]=unique(c.omega_mean(valid),'sorted');
a=c.amplitude_deg(valid); a=a(ia);
if numel(x)<2, return; end
inside=grid>=x(1) & grid<=x(end);
y(inside)=interp1(x,a,grid(inside),'linear');
end

function y=align_levels(levels,grid)
y=nan(size(grid));
if isempty(levels), return; end
valid=isfinite(levels.omega) & isfinite(levels.process_amplitude_deg) & ...
    levels.process_amplitude_deg<=90 & levels.status~="out_of_scope";
if ~any(valid), return; end
[tf,loc]=ismembertol(levels.omega(valid),grid,1e-10,'DataScale',1);
values=levels.process_amplitude_deg(valid);
y(loc(tf))=values(tf);
end

function [meanAbs,signed,coverage]=lag_metrics(up,down,omega)
valid=isfinite(up) & isfinite(down);
coverage=[mean(isfinite(up)),mean(isfinite(down))];
if nnz(valid)<2
    meanAbs=NaN; signed=NaN; return
end
x=omega(valid); delta=up(valid)-down(valid);
span=max(x)-min(x);
if span<=0
    meanAbs=mean(abs(delta),'omitnan'); signed=mean(delta,'omitnan');
else
    meanAbs=trapz(x,abs(delta))/span;
    signed=trapz(x,delta)/span;
end
end

function value=last_omega(sweep)
value=NaN;
if ~isempty(sweep.cycle), value=sweep.cycle.omega_mean(end); end
end

function value=last_level_omega(sweep)
value=NaN;
if ~isempty(sweep.levels), value=sweep.levels.omega(end); end
end

function t=cycle_records(c,value,protocol)
if isempty(c), t=table(); return; end
t=table(repmat(value,height(c),1),repmat(string(protocol),height(c),1), ...
    c.omega_mean,c.amplitude_deg,c.out_of_scope,c.status, ...
    'VariableNames',{'control_value','protocol','omega','amplitude_deg', ...
    'out_of_scope','status'});
end

function t=level_records(levels,value,protocol)
if isempty(levels), t=table(); return; end
t=table(repmat(value,height(levels),1), ...
    repmat(string(protocol),height(levels),1),levels.omega, ...
    levels.process_amplitude_deg,levels.status=="out_of_scope",levels.status, ...
    'VariableNames',{'control_value','protocol','omega','amplitude_deg', ...
    'out_of_scope','status'});
end
