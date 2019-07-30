function [] = highlightoffsets( offsets, Tp, Np, hfig)

if ~exist('hfig', 'var')
    hfig = gcf;
end

ylim = get( gca, 'Ylim' );
if ~exist('Np', 'var')
    Np = ylim(2);
end

% highlight pattern
D_MAX = 20;
patt_locs = offsets;
p_length = Tp;
hold on

for i = 1:numel(patt_locs)
    start_loc = patt_locs(i);
    xpos = start_loc + p_length + D_MAX;
    
    patch([start_loc xpos xpos start_loc ], [ ylim(1) ylim(1) Np Np ], [0.6 0.4 0.9], 'FaceAlpha',0.5, 'EdgeColor','none')
end

drawnow;
hold off


end