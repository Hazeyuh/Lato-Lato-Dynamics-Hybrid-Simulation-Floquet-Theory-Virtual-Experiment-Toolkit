function slices = operating_slices(p,opts)
%OPERATING_SLICES Three pairwise parameter maps for the F3 main result.
%   Produces A-Omega at fixed L, L-Omega at fixed A, and A-L at a
%   representative Omega. This avoids an unnecessary 51x51x61 full cube.

defaults=struct('A',linspace(.005,.040,51), ...
    'L',linspace(.75*p.ell,1.30*p.ell,51), ...
    'omega',linspace(8,24,61),'A_fixed',p.drive.amplitude, ...
    'L_fixed',p.ell,'omega_fixed',p.drive.omega,'mode',"quick", ...
    'min_cycles',20,'max_cycles',150,'checkpoint_dir',"", ...
    'continue_on_error',true);
opts=lato.scan.options(opts,defaults);
opts.mode=string(opts.mode);
if opts.mode=="formal" && any([numel(opts.A),numel(opts.L),numel(opts.omega)]<50)
    error('lato:scan:FormalResolution', ...
        'Formal pairwise slices require >=50 points on every varied axis.');
end

base=struct('mode',opts.mode,'min_cycles',opts.min_cycles, ...
    'max_cycles',opts.max_cycles,'checkpoint_file',"", ...
    'continue_on_error',opts.continue_on_error);
if strlength(opts.checkpoint_dir)>0 && ~isfolder(opts.checkpoint_dir)
    mkdir(opts.checkpoint_dir);
end

ao=base; ao.A=opts.A; ao.L=opts.L_fixed; ao.omega=opts.omega;
if strlength(opts.checkpoint_dir)>0
    ao.checkpoint_file=fullfile(opts.checkpoint_dir,'F3_A_omega_checkpoint.mat');
end
slices.A_omega=lato.scan.operating_map(p,ao);

lo=base; lo.A=opts.A_fixed; lo.L=opts.L; lo.omega=opts.omega;
if strlength(opts.checkpoint_dir)>0
    lo.checkpoint_file=fullfile(opts.checkpoint_dir,'F3_L_omega_checkpoint.mat');
end
slices.L_omega=lato.scan.operating_map(p,lo);

al=base; al.A=opts.A; al.L=opts.L; al.omega=opts.omega_fixed;
if strlength(opts.checkpoint_dir)>0
    al.checkpoint_file=fullfile(opts.checkpoint_dir,'F3_A_L_checkpoint.mat');
end
slices.A_L=lato.scan.operating_map(p,al);
slices.kind="pairwise_operating_slices";
slices.mode=opts.mode;
slices.options=opts;
end
