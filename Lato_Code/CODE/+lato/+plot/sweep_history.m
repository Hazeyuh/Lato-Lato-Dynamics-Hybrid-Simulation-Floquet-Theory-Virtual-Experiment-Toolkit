function fig = sweep_history(collection, opts)
%SWEEP_HISTORY Plot F6 process histories from continuous/step sweeps.

if nargin<2, opts=struct(); end
if ~isfield(opts,'visible'), opts.visible="on"; end
opts.visible=string(opts.visible);

names = fieldnames(collection);
fig = figure('Color','w','Visible',opts.visible,'Position',[100 100 1180 820]);
tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
ax1=nexttile(tl); hold(ax1,'on');
ax2=nexttile(tl); hold(ax2,'on');
ax3=nexttile(tl); hold(ax3,'on');
colors = lines(max(1,numel(names)));
maximumAmplitude = 0;

for k=1:numel(names)
    item=collection.(names{k});
    if ~isfield(item,'cycle') || isempty(item.cycle), continue; end
    c=item.cycle;
    maximumAmplitude=max(maximumAmplitude,max(c.amplitude_deg,[],'omitnan'));
    tmid=0.5*(c.start_time+c.end_time);
    label=strrep(names{k},'_',' ');
    plot(ax1,tmid,c.omega_mean,'LineWidth',1.35,'Color',colors(k,:), ...
        'DisplayName',label);
    plot(ax2,tmid,c.amplitude_deg,'LineWidth',1.35,'Color',colors(k,:), ...
        'DisplayName',label);
    plot(ax3,c.omega_mean,c.amplitude_deg,'LineWidth',1.15, ...
        'Color',colors(k,:),'DisplayName',label);
    if any(c.out_of_scope)
        scatter(ax2,tmid(c.out_of_scope),c.amplitude_deg(c.out_of_scope), ...
            28,colors(k,:),'x','HandleVisibility','off');
        scatter(ax3,c.omega_mean(c.out_of_scope),c.amplitude_deg(c.out_of_scope), ...
            28,colors(k,:),'x','HandleVisibility','off');
    end
end

ylabel(ax1,'Cycle-mean \Omega (rad s^{-1})');
ylabel(ax2,'Cycle amplitude (deg)');
ylabel(ax3,'Cycle amplitude (deg)'); xlabel(ax3,'Cycle-mean \Omega (rad s^{-1})');
xlabel(ax1,'Time (s)'); xlabel(ax2,'Time (s)');
if maximumAmplitude >= 45
    yline(ax2,45,'k--','45^\circ threshold','HandleVisibility','off');
    yline(ax3,45,'k--','HandleVisibility','off');
end
if maximumAmplitude >= 90
    yline(ax2,90,'r--','90^\circ scope limit','HandleVisibility','off');
    yline(ax3,90,'r--','HandleVisibility','off');
end
upper=max(2,1.15*maximumAmplitude);
if maximumAmplitude>=45, upper=max(50,upper); end
if maximumAmplitude>=90, upper=max(95,upper); end
ylim(ax2,[0 upper]); ylim(ax3,[0 upper]);
if maximumAmplitude < 45
    text(ax2,0.99,0.92,'All responses remain below 45^\circ', ...
        'Units','normalized','HorizontalAlignment','right','Color',[.3 .3 .3]);
end
for ax=[ax1 ax2 ax3], grid(ax,'on'); box(ax,'on'); end
legend(ax1,'Location','eastoutside');
title(tl,'F6  Phase-continuous and condition-triggered sweep histories');
end
