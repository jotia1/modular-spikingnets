function [] = visualiseaccuracy(net, out, hfig)

if ~exist('hfig', 'var')
    hfig = gcf;
end

ylim = get( gca, 'Ylim' );
if ~exist('Np', 'var')
    Np = ylim(2);
end

patt_locs = out.offsets;
D_MAX = net.delay_max;
GREEN = [0.1 0.9 0.1];
RED = [0.8 0.1 0.1];
PURPLE = [0.6 0.4 0.9];
YELLOW = [0.95 0.5 0];
p_length = net.Tp;
stt = out.spike_time_trace(out.spike_time_trace(:, 2) == 2001, :);
hold on

last_slot = 0;

for i = 1:numel(patt_locs)
    start_loc = patt_locs(i);
    xpos = start_loc + p_length + D_MAX;
    
    output_spikes = stt(stt(:, 1) > start_loc & stt(:, 1) < xpos, :);
    
    if numel(output_spikes) > 0
        patch([start_loc xpos xpos start_loc ], [ ylim(1) ylim(1) Np Np ], GREEN, 'FaceAlpha',0.5, 'EdgeColor','none')
    else
        patch([start_loc xpos xpos start_loc ], [ ylim(1) ylim(1) Np Np ], RED, 'FaceAlpha',0.5, 'EdgeColor','none')        
    end
    
    % Also colour gaps between patterns
    output_spikes = stt(stt(:, 1) > last_slot & stt(:, 1) < start_loc, :);
    gap_xs = [last_slot, start_loc, start_loc, last_slot];
    if numel(output_spikes) > 0  % Mistake
        patch(gap_xs, [ ylim(1) ylim(1) Np Np ], YELLOW, 'FaceAlpha',0.5, 'EdgeColor','none')
    %else
        %patch(gap_xs, [ ylim(1) ylim(1) Np Np ], GREEN, 'FaceAlpha',0.1, 'EdgeColor','none')        
    end
    last_slot = xpos;
end

drawnow;
hold off

end