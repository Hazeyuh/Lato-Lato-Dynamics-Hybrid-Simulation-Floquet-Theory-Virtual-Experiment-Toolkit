function fig = operating_maps(map, opts)
%OPERATING_MAPS Plot F3 from an A-L-omega scan without interpolation.

if nargin<2, opts=struct(); end
if ~isfield(opts,'visible'), opts.visible="on"; end
opts.visible=string(opts.visible);

if isfield(map,'kind') && string(map.kind)=="pairwise_operating_slices"
    fig=plot_pairwise_slices(map,opts.visible);
    return
end

if ~isfield(opts,'L_index'), opts.L_index=ceil(numel(map.L)/2); end
if ~isfield(opts,'A_index'), opts.A_index=ceil(numel(map.A)/2); end

j = min(opts.L_index,numel(map.L));
i = min(opts.A_index,numel(map.A));
fig = figure('Color','w','Visible',opts.visible,'Position',[100 100 1180 780]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl);
z1 = squeeze(map.steady_amplitude_deg(:,j,:));
imagesc(ax1,map.omega,map.A*1e3,z1);
set(ax1,'YDir','normal');
xlabel(ax1,'Drive frequency \Omega (rad s^{-1})');
ylabel(ax1,'Pivot amplitude A (mm)');
title(ax1,sprintf('Steady response at L = %.3f m',map.L(j)));
decorate_map(ax1,z1,squeeze(map.steady(:,j,:)), ...
    squeeze(map.out_of_scope(:,j,:)),map.omega,map.A*1e3);

ax2 = nexttile(tl);
z2 = squeeze(map.steady_amplitude_deg(i,:,:));
imagesc(ax2,map.omega,map.L,z2);
set(ax2,'YDir','normal');
xlabel(ax2,'Drive frequency \Omega (rad s^{-1})');
ylabel(ax2,'Length L (m)');
title(ax2,sprintf('Steady response at A = %.1f mm',1e3*map.A(i)));
decorate_map(ax2,z2,squeeze(map.steady(i,:,:)), ...
    squeeze(map.out_of_scope(i,:,:)),map.omega,map.L);

ax3 = nexttile(tl);
[W,AA] = meshgrid(map.omega,map.A);
h = AA.*W.^2/9.80665;
scatter(ax3,W(:),h(:),28,z1(:),'filled');
hold(ax3,'on');
nonsteady = ~squeeze(map.steady(:,j,:)) & ~squeeze(map.out_of_scope(:,j,:));
scatter(ax3,W(nonsteady),h(nonsteady),18,[.35 .35 .35],'x');
oos = squeeze(map.out_of_scope(:,j,:));
scatter(ax3,W(oos),h(oos),24,[.75 .1 .1],'x','LineWidth',1.2);
xlabel(ax3,'Drive frequency \Omega (rad s^{-1})');
ylabel(ax3,'Excitation strength h = A\Omega^2/g');
title(ax3,'Same slice in physical coordinates');
grid(ax3,'on'); box(ax3,'on');

ax4 = nexttile(tl);
validFraction = squeeze(mean(map.steady,1));
largeFraction = squeeze(mean(map.large_amplitude,1));
plot(ax4,map.omega,100*validFraction(j,:),'-o','LineWidth',1.4, ...
    'DisplayName','steady');
hold(ax4,'on');
plot(ax4,map.omega,100*largeFraction(j,:),'-s','LineWidth',1.4, ...
    'DisplayName','45^\circ < amplitude \leq 90^\circ');
xlabel(ax4,'Drive frequency \Omega (rad s^{-1})');
ylabel(ax4,'Fraction across A (%)');
title(ax4,'Data-quality and amplitude classes');
ylim(ax4,[0 100]); grid(ax4,'on'); box(ax4,'on'); legend(ax4,'Location','best');

cb = colorbar(ax2);
cb.Label.String = 'Steady half-amplitude (deg)';
clim(ax1,[0 90]); clim(ax2,[0 90]); clim(ax3,[0 90]);
colormap(fig,turbo(256));
title(tl,'F3  Nonlinear A-L-\Omega operating map');
end

function fig=plot_pairwise_slices(slices,visible)
ao=slices.A_omega; lo=slices.L_omega; al=slices.A_L;
zAO=squeeze(ao.steady_amplitude_deg(:,1,:));
zLO=squeeze(lo.steady_amplitude_deg(1,:,:));
zAL=squeeze(al.steady_amplitude_deg(:,:,1));
fig=figure('Color','w','Visible',visible,'Position',[100 100 1180 780]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax1=nexttile(tl); imagesc(ax1,ao.omega,1e3*ao.A,zAO); set(ax1,'YDir','normal');
xlabel(ax1,'Drive frequency \Omega (rad s^{-1})'); ylabel(ax1,'A (mm)');
title(ax1,sprintf('A-\\Omega at L = %.3f m',ao.L));
decorate_map(ax1,zAO,squeeze(ao.steady(:,1,:)), ...
    squeeze(ao.out_of_scope(:,1,:)),ao.omega,1e3*ao.A);

ax2=nexttile(tl); imagesc(ax2,lo.omega,lo.L,zLO); set(ax2,'YDir','normal');
xlabel(ax2,'Drive frequency \Omega (rad s^{-1})'); ylabel(ax2,'L (m)');
title(ax2,sprintf('L-\\Omega at A = %.1f mm',1e3*lo.A));
decorate_map(ax2,zLO,squeeze(lo.steady(1,:,:)), ...
    squeeze(lo.out_of_scope(1,:,:)),lo.omega,lo.L);

ax3=nexttile(tl); imagesc(ax3,al.L,1e3*al.A,zAL); set(ax3,'YDir','normal');
xlabel(ax3,'L (m)'); ylabel(ax3,'A (mm)');
title(ax3,sprintf('A-L at \\Omega = %.2f rad s^{-1}',al.omega));
decorate_map(ax3,zAL,squeeze(al.steady(:,:,1)), ...
    squeeze(al.out_of_scope(:,:,1)),al.L,1e3*al.A);

ax4=nexttile(tl); [W,AA]=meshgrid(ao.omega,ao.A);
h=AA.*W.^2/9.80665;
scatter(ax4,W(:),h(:),25,zAO(:),'filled'); hold(ax4,'on');
nonsteady=~squeeze(ao.steady(:,1,:)) & ~squeeze(ao.out_of_scope(:,1,:));
scatter(ax4,W(nonsteady),h(nonsteady),16,[.3 .3 .3],'x');
oos=squeeze(ao.out_of_scope(:,1,:));
scatter(ax4,W(oos),h(oos),22,[.8 .1 .1],'x','LineWidth',1.1);
xlabel(ax4,'Drive frequency \Omega (rad s^{-1})');
ylabel(ax4,'h = A\Omega^2/g'); title(ax4,'A-\Omega slice in physical coordinates');
grid(ax4,'on'); box(ax4,'on');

for ax=[ax1 ax2 ax3 ax4], clim(ax,[0 90]); end
colormap(fig,turbo(256)); cb=colorbar(ax2); cb.Label.String='Steady half-amplitude (deg)';
title(tl,'F3  Pairwise nonlinear operating maps');
end

function decorate_map(ax,z,steady,oos,x,y)
hold(ax,'on');
[X,Y] = meshgrid(x,y);
nonsteady = ~steady & ~oos;
scatter(ax,X(nonsteady),Y(nonsteady),15,[.3 .3 .3],'x');
scatter(ax,X(oos),Y(oos),20,[.8 .1 .1],'x','LineWidth',1.2);
if any(isfinite(z),'all') && min(z,[],'all','omitnan') <= 45 && max(z,[],'all','omitnan') >= 45
    contour(ax,X,Y,z,[45 45],'k--','LineWidth',1.1);
end
box(ax,'on');
end
