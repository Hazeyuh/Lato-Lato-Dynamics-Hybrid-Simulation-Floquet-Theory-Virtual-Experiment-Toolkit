function figures = sweep_hysteresis_maps(map, opts)
%SWEEP_HYSTERESIS_MAPS Plot finite-time up/down sweep maps and diagnostics.

if nargin<2, opts=struct(); end
if ~isfield(opts,'visible'), opts.visible="on"; end
opts.visible=string(opts.visible);

xContinuous=map.omega_grid/map.omega0;
xStep=map.omega_levels/map.omega0;
[controlDisplay,controlLabel]=display_control(map);
continuousAmplitudeLimit=positive_limit([map.continuous.up_deg(:); ...
    map.continuous.down_deg(:)]);
stepAmplitudeLimit=positive_limit([map.step.up_deg(:);map.step.down_deg(:)]);
continuousDifferenceLimit=symmetric_limit(map.continuous.difference_deg(:));
stepDifferenceLimit=symmetric_limit(map.step.difference_deg(:));

figures.maps=figure('Color','w','Visible',opts.visible, ...
    'Position',[60 60 1400 760]);
tl=tiledlayout(figures.maps,2,3,'TileSpacing','compact','Padding','compact');
draw_amplitude(nexttile(tl),xContinuous,controlDisplay, ...
    map.continuous.up_deg,'Continuous up',continuousAmplitudeLimit);
draw_amplitude(nexttile(tl),xContinuous,controlDisplay, ...
    map.continuous.down_deg,'Continuous down',continuousAmplitudeLimit);
draw_difference(nexttile(tl),xContinuous,controlDisplay, ...
    map.continuous.difference_deg,'Continuous up - down', ...
    continuousDifferenceLimit);
draw_amplitude(nexttile(tl),xStep,controlDisplay,map.step.up_deg, ...
    'Finite-residence steps up',stepAmplitudeLimit);
draw_amplitude(nexttile(tl),xStep,controlDisplay,map.step.down_deg, ...
    'Finite-residence steps down',stepAmplitudeLimit);
draw_difference(nexttile(tl),xStep,controlDisplay,map.step.difference_deg, ...
    'Step up - down',stepDifferenceLimit);
for ax=findall(figures.maps,'Type','axes').'
    xlabel(ax,'Cycle-mean \Omega / \omega_0'); ylabel(ax,controlLabel);
end
modeLabel=upper(char(map.options.mode));
controlLabelShort=char(map.control);
title(tl,sprintf(['F8a  Collision-suppressed startup paths — %s, %s ' ...
    '(row-wise scales; white = missing / out of scope)'], ...
    modeLabel,controlLabelShort));

figures.metrics=figure('Color','w','Visible',opts.visible, ...
    'Position',[90 90 1250 800]);
tl2=tiledlayout(figures.metrics,2,2,'TileSpacing','compact','Padding','compact');
representative_paths(nexttile(tl2),xContinuous,controlDisplay, ...
    map.continuous.up_deg,map.continuous.down_deg,'Continuous paths');
representative_paths(nexttile(tl2),xStep,controlDisplay, ...
    map.step.up_deg,map.step.down_deg,'Finite-residence step paths');

ax=nexttile(tl2); hold(ax,'on');
plot(ax,controlDisplay,map.continuous.mean_abs_lag_deg,'o-', ...
    'LineWidth',1.4,'DisplayName','continuous');
plot(ax,controlDisplay,map.step.mean_abs_lag_deg,'s-', ...
    'LineWidth',1.4,'DisplayName','steps');
truncated=any(map.step.coverage_fraction<1-1e-10,2);
if any(truncated)
    plot(ax,controlDisplay(truncated),map.step.mean_abs_lag_deg(truncated), ...
        'rx','LineWidth',1.5,'MarkerSize',7, ...
        'DisplayName','steps on truncated overlap');
end
xlabel(ax,controlLabel); ylabel(ax,'Mean absolute up/down gap (deg)');
title(ax,'Dynamic-lag metric'); grid(ax,'on'); box(ax,'on'); legend(ax);

ax=nexttile(tl2); hold(ax,'on');
plot(ax,controlDisplay,map.continuous.coverage_fraction(:,1),'o-', ...
    'LineWidth',1.25,'DisplayName','continuous up');
plot(ax,controlDisplay,map.continuous.coverage_fraction(:,2),'o--', ...
    'LineWidth',1.25,'DisplayName','continuous down');
plot(ax,controlDisplay,map.step.coverage_fraction(:,1),'s-', ...
    'LineWidth',1.25,'DisplayName','step up');
plot(ax,controlDisplay,map.step.coverage_fraction(:,2),'s--', ...
    'LineWidth',1.25,'DisplayName','step down');
xlabel(ax,controlLabel); ylabel(ax,'Frequency-axis coverage'); ylim(ax,[0 1.05]);
title(ax,'Coverage before stop / 90^\circ limit');
grid(ax,'on'); box(ax,'on'); legend(ax,'Location','best');
title(tl2,sprintf(['F8b  Collision-suppressed dynamic path dependence — %s, %s; ' ...
    'not a quasistatic hysteresis claim'],modeLabel,controlLabelShort));
end

function limit=positive_limit(values)
limit=max(values,[],'omitnan');
if isempty(limit) || ~isfinite(limit), limit=1; end
limit=max(limit,1);
end

function limit=symmetric_limit(values)
limit=max(abs(values),[],'omitnan');
if isempty(limit) || ~isfinite(limit) || limit==0, limit=1; end
end

function draw_amplitude(ax,x,y,z,label,limit)
h=imagesc(ax,x,y,z); set(h,'AlphaData',isfinite(z));
set(ax,'YDir','normal','Color',[0.94 0.94 0.94]);
colormap(ax,turbo(256)); clim(ax,[0 limit]); colorbar(ax);
hold(ax,'on');
finiteValues=z(isfinite(z));
if any(finiteValues>=45) && any(finiteValues<45)
    contour(ax,x,y,z,[45 45],'k--','LineWidth',1.1);
end
title(ax,label); box(ax,'on');
end

function draw_difference(ax,x,y,z,label,limit)
h=imagesc(ax,x,y,z); set(h,'AlphaData',isfinite(z));
set(ax,'YDir','normal','Color',[0.94 0.94 0.94]);
colormap(ax,diverging_map(257)); clim(ax,[-limit limit]); colorbar(ax);
title(ax,label); box(ax,'on');
end

function representative_paths(ax,x,control,up,down,label)
hold(ax,'on');
indices=unique(round(linspace(1,numel(control),min(3,numel(control)))));
colors=lines(numel(indices));
for j=1:numel(indices)
    k=indices(j);
    plot(ax,x,up(k,:),'LineWidth',1.5,'Color',colors(j,:), ...
        'DisplayName',sprintf('%g up',control(k)));
    plot(ax,x,down(k,:),'--','LineWidth',1.5,'Color',colors(j,:), ...
        'DisplayName',sprintf('%g down',control(k)));
end
maximum=max([up(:);down(:)],[],'omitnan');
if isempty(maximum) || ~isfinite(maximum), maximum=1; end
if maximum>=45
    yline(ax,45,'k:','45^\circ','HandleVisibility','off');
    ylim(ax,[0 max(50,1.08*maximum)]);
else
    ylim(ax,[0 max(1,1.12*maximum)]);
    text(ax,0.98,0.94,'All shown paths remain below 45^\circ', ...
        'Units','normalized','HorizontalAlignment','right', ...
        'Color',[0.35 0.35 0.35]);
end
xlabel(ax,'Cycle-mean \Omega / \omega_0'); ylabel(ax,'Process amplitude (deg)');
title(ax,label); grid(ax,'on'); box(ax,'on'); legend(ax,'Location','best');
end

function [values,label]=display_control(map)
if map.control=="fixed_A"
    values=1e3*map.control_values;
    label='Fixed displacement A (mm)';
else
    values=map.control_values;
    label='Fixed acceleration ratio h = a_0/g';
end
end

function cmap=diverging_map(n)
if nargin<1, n=257; end
half=ceil(n/2);
blue=[0.18 0.37 0.68]; white=[1 1 1]; red=[0.70 0.16 0.20];
c1=[linspace(blue(1),white(1),half).', ...
    linspace(blue(2),white(2),half).',linspace(blue(3),white(3),half).'];
c2=[linspace(white(1),red(1),n-half+1).', ...
    linspace(white(2),red(2),n-half+1).', ...
    linspace(white(3),red(3),n-half+1).'];
cmap=[c1;c2(2:end,:)];
end
