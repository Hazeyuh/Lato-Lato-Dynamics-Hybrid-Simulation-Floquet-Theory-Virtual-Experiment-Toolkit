function scan = independent_parameter(p,parameter,values,opts)
%INDEPENDENT_PARAMETER One-dimensional A, L, or omega restart scan.

if nargin<4, opts=struct(); end
defaults=struct('omega',p.drive.omega,'A',p.drive.amplitude,'L',p.ell, ...
    'min_cycles',20,'max_cycles',150,'mode',"quick",'keep_runs',false);
opts=lato.scan.options(opts,defaults);
parameter=string(parameter); opts.mode=string(opts.mode);
mustBeMember(parameter,["A","L","omega"]);
if opts.mode=="formal" && numel(values)<50
    error('lato:scan:FormalResolution', ...
        'Formal one-dimensional scans require at least 50 points.');
end

n=numel(values); items=repmat(template(),n,1); runs=cell(n,1);
for k=1:n
    A=opts.A; L=opts.L; omega=opts.omega;
    switch parameter
        case "A", A=values(k);
        case "L", L=values(k);
        case "omega", omega=values(k);
    end
    pk=lato.scan.configure_point(p,omega,A,L,opts.max_cycles*2*pi/omega);
    try
        r=lato.simulate_hybrid(pk);
        items(k)=lato.scan.summarize_result(r,'Omega',omega,'A',A,'L',L);
        if opts.keep_runs, runs{k}=r; end
    catch ME
        items(k)=template(); items(k).omega=omega; items(k).A=A; items(k).L=L;
        items(k).h=A*omega^2/p.g; items(k).status="solver_error";
        if opts.keep_runs, runs{k}=ME; end
        warning('lato:scan:PointFailed','%s=%g failed: %s',parameter,values(k),ME.message);
    end
end
scan=struct('kind',"independent_parameter",'parameter',parameter, ...
    'values',values,'mode',opts.mode,'items',items, ...
    'table',lato.scan.summary_table(items),'options',opts);
if opts.keep_runs, scan.runs=runs; end
end

function s=template()
s=struct('omega',NaN,'A',NaN,'L',NaN,'h',NaN,'steady',false, ...
    'steady_amplitude_deg',NaN,'max_amplitude_deg',NaN, ...
    'amplitude_class',"unavailable",'out_of_scope',false, ...
    'large_amplitude',false,'collision_rate_hz',NaN, ...
    'slack_fraction',NaN,'midline_offset_deg',NaN, ...
    'mean_energy_J',NaN,'status',"not_run");
end
