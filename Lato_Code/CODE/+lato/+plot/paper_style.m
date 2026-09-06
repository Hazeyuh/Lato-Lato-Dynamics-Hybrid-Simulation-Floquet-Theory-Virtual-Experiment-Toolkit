function paper_style(fig)
%PAPER_STYLE Apply the common paper figure appearance.
if nargin < 1 || isempty(fig)
    fig = gcf;
end
set(fig, 'Color', 'white');
axes_handles = findall(fig, 'Type', 'axes');
for ax = reshape(axes_handles, 1, [])
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 10, ...
        'LineWidth', 0.8, 'Box', 'on', 'Layer', 'top');
    grid(ax, 'on');
    ax.GridAlpha = 0.12;
end
text_handles = findall(fig, '-property', 'FontName');
set(text_handles, 'FontName', 'Times New Roman');
end
