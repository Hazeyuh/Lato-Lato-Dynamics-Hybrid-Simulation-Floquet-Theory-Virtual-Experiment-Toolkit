function map = operating_map(p, opts)
%OPERATING_MAP A-L-omega nonlinear operating map or one pairwise slice.
%   Full factorial formal grids require at least 50 samples per dimension
%   and are an explicitly advanced, high-cost option. For manuscript work,
%   use lato.scan.operating_slices instead.

defaults=struct('A',linspace(0.005,0.035,51), ...
    'L',linspace(0.14,0.24,51),'omega',linspace(8,24,61), ...
    'mode',"quick",'min_cycles',20,'max_cycles',150, ...
    'checkpoint_file',"",'continue_on_error',true);
opts=lato.scan.options(opts,defaults);
opts.mode=string(opts.mode); opts.checkpoint_file=string(opts.checkpoint_file);
mustBeMember(opts.mode,["quick","formal"]);

counts=[numel(opts.A),numel(opts.L),numel(opts.omega)];
isSlice=sum(counts==1)==1 && all(counts(counts>1)>=50);
if opts.mode == "formal" && ~(all(counts>=50) || isSlice)
    error('lato:scan:FormalResolution', ...
        'Formal maps require >=50 points on every varied dimension.');
end
if opts.max_cycles < opts.min_cycles
    error('lato:scan:InvalidCycles', 'max_cycles must be >= min_cycles.');
end

nA = numel(opts.A); nL = numel(opts.L); nW = numel(opts.omega);
shape = [nA, nL, nW];
if prod(shape)>1e4
    warning('lato:scan:HighCostMap', ...
        ['This request contains %d hybrid simulations. Use ', ...
         'lato.scan.operating_slices for the manuscript default.'],prod(shape));
end
steadyAmplitude = nan(shape);
maxAmplitude = nan(shape);
collisionRate = nan(shape);
slackFraction = nan(shape);
midlineOffset = nan(shape);
meanEnergy = nan(shape);
steady = false(shape);
large = false(shape);
outOfScope = false(shape);
status = strings(shape);
g = get_gravity(p);

rows = repmat(row_template(), prod(shape), 1);
q = 0;
for j = 1:nL
    for i = 1:nA
        for k = 1:nW
            q = q + 1;
            omega = opts.omega(k); A = opts.A(i); L = opts.L(j);
            duration = opts.max_cycles * 2*pi / omega;
            pk = lato.scan.configure_point(p, omega, A, L, duration);
            try
                result = lato.simulate_hybrid(pk);
                s = lato.scan.summarize_result(result, ...
                    'Omega', omega, 'A', A, 'L', L);
            catch ME
                s = row_template();
                s.omega = omega; s.A = A; s.L = L;
                s.h = A * omega^2 / g;
                s.status = "solver_error";
                if ~opts.continue_on_error, rethrow(ME); end
                warning('lato:scan:PointFailed', ...
                    'A=%g, L=%g, omega=%g failed: %s', A, L, omega, ME.message);
            end
            rows(q) = s;
            steadyAmplitude(i,j,k) = s.steady_amplitude_deg;
            maxAmplitude(i,j,k) = s.max_amplitude_deg;
            collisionRate(i,j,k) = s.collision_rate_hz;
            slackFraction(i,j,k) = s.slack_fraction;
            midlineOffset(i,j,k) = s.midline_offset_deg;
            meanEnergy(i,j,k) = s.mean_energy_J;
            steady(i,j,k) = s.steady;
            large(i,j,k) = s.large_amplitude;
            outOfScope(i,j,k) = s.out_of_scope;
            status(i,j,k) = s.status;
        end
    end
    if strlength(opts.checkpoint_file) > 0
        checkpoint = assemble_map(); %#ok<NASGU>
        save(opts.checkpoint_file, 'checkpoint', '-v7.3');
    end
end

map = assemble_map();

    function value = assemble_map()
        value = struct();
        value.kind = "A_L_omega_operating_map";
        value.mode = opts.mode;
        value.A = opts.A;
        value.L = opts.L;
        value.omega = opts.omega;
        value.h = reshape(opts.A,[],1,1) .* reshape(opts.omega,1,1,[]).^2 / g;
        value.steady_amplitude_deg = steadyAmplitude;
        value.max_amplitude_deg = maxAmplitude;
        value.steady = steady;
        value.large_amplitude = large;
        value.out_of_scope = outOfScope;
        value.collision_rate_hz = collisionRate;
        value.slack_fraction = slackFraction;
        value.midline_offset_deg = midlineOffset;
        value.mean_energy_J = meanEnergy;
        value.status = status;
        value.table = lato.scan.summary_table(rows(1:q));
        value.options = opts;
    end
end

function g = get_gravity(p)
g = 9.80665;
if isfield(p, 'g'), g = p.g; end
end

function s = row_template()
s = struct('omega',NaN,'A',NaN,'L',NaN,'h',NaN,'steady',false, ...
    'steady_amplitude_deg',NaN,'max_amplitude_deg',NaN, ...
    'amplitude_class',"unavailable",'out_of_scope',false, ...
    'large_amplitude',false,'collision_rate_hz',NaN, ...
    'slack_fraction',NaN,'midline_offset_deg',NaN, ...
    'mean_energy_J',NaN,'status',"not_run");
end
