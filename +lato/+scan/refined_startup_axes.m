function axes = refined_startup_axes(baseMaps,p)
%REFINED_STARTUP_AXES Targeted F7 axes that preserve every formal point.

arguments
    baseMaps (1,1) struct
    p (1,1) struct
end

axes = struct();
axes.A = refine_existing_and_extend( ...
    baseMaps.A_omega.parameter_values,0.001,0.015);
axes.L = refine_existing_and_extend( ...
    baseMaps.L_omega.parameter_values,0.10,0.24);
axes.theta_seed_deg = refine_existing( ...
    baseMaps.theta_seed_omega.parameter_values);
axes.omega_A = refine_window(baseMaps.A_omega.omega,[13,16]);
axes.omega_seed = refine_window( ...
    baseMaps.theta_seed_omega.omega,[13,16]);
omegaLBase = extend_upper(baseMaps.L_omega.omega, ...
    2*sqrt(p.g/min(axes.L))+1);
centreRange = [2*sqrt(p.g/max(axes.L))-1, ...
    2*sqrt(p.g/min(axes.L))+1];
axes.omega_L = refine_window(omegaLBase,centreRange);
end

function values = refine_existing(values)
values = reshape(values,1,[]);
midpoints = 0.5*(values(1:end-1)+values(2:end));
values = unique([values,midpoints],'sorted');
end

function values = refine_existing_and_extend(values,newMinimum,newMaximum)
old = reshape(values,1,[]);
oldStep = median(diff(old));
values = refine_existing(old);
lower = old(1)-oldStep:-oldStep:newMinimum;
upper = old(end)+oldStep:oldStep:newMaximum;
values = unique([values,lower,upper,newMinimum,newMaximum],'sorted');
values = values(values>=newMinimum-10*eps & values<=newMaximum+10*eps);
end

function values = refine_window(values,window)
values = reshape(values,1,[]);
midpoints = 0.5*(values(1:end-1)+values(2:end));
midpoints = midpoints(midpoints>=window(1) & midpoints<=window(2));
values = unique([values,midpoints],'sorted');
end

function values = extend_upper(values,newMaximum)
values = reshape(values,1,[]);
step = median(diff(values));
extra = values(end)+step:step:newMaximum;
if isempty(extra) || extra(end)<newMaximum-10*eps
    extra(end+1) = newMaximum;
end
values = unique([values,extra],'sorted');
end
