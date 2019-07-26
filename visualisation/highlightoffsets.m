function [] = highlightoffsets( offsets, Tp, hfig)

if ~exist('hfig', 'var')
    hfig = gcf;
end

% highlight pattern
D_MAX = 20;
patt_locs = offsets;
p_length = Tp;
hold on

for i = 1:numel(patt_locs)
    start_loc = patt_locs(i);
    xpos = start_loc + p_length + D_MAX;
    ylim = get( gca, 'Ylim' );
    patch([start_loc xpos xpos start_loc ], [ 0 0 ylim(2) ylim(2) ], [0.6 0.4 0.9], 'FaceAlpha',0.5, 'EdgeColor','none')
end

drawnow;
hold off


end