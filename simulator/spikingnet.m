function [ out ] = spikingnet( net )
%SPIKINGNET Computational engine of a spiking network
%   Given a network description run the simulation and compute the final
%   state variables of the network
%   In connectivity matrices, rows are the pre-synaptic (from) neurons and
%   columns are the post-synaptic (to) neurons.

out = struct();
out.timing_info.init_time = tic;
[res, missing] = validatenetwork(net);
assert(res, missing);

if net.rand_seed ~= -1
    rng(net.rand_seed); 
end

ms_per_sec = 1000;

delays = sparse(net.delays);
delayst = zeros(numel(net.delays_to_save), net.N, net.sim_time_sec * ms_per_sec);
variance = sparse(net.variance);
vart = zeros(numel(net.variance_to_save), net.N, net.sim_time_sec * ms_per_sec);
%v_threst = zeros(numel(net.v_thres_to_save), net.sim_time_sec * ms_per_sec);
w = sparse(net.w);

N = net.N;
N_inp = net.group_sizes(1);
v = ones(N, 1) * net.v_rest;
vt = zeros(numel(net.voltages_to_save), net.sim_time_sec * ms_per_sec);
last_spike_time = zeros(net.N, 1) * -Inf;
% Dynamic threshold stuff
v_thres = ones(size(v)) * net.v_thres;
net.thres_rise = net.thres_rise * net.dynamic_threshold; % zero if false

conns = w ~= 0;
dApre = sparse(zeros(size(w)));
dApost = sparse(zeros(size(w)));
STDPdecaypre = exp(-1/net.taupre);
STDPdecaypost = exp(-1/net.taupost);
active_spikes = cell(net.delay_max, 1);  % To track when spikes arrive
active_idx = 1;
%I0 = net.fgi;
%If = I0 * 0.596;
%Tf = 31.4;

% output variables
out.timing_info.init_toc = toc(out.timing_info.init_time);
out.timing_info.sim_sec_tics = uint64(zeros(net.sim_time_sec, 1));
out.timing_info.sim_sec_tocs = zeros(net.sim_time_sec, 1);
out.timing_info.full_sec_tocs = zeros(net.sim_time_sec, 1);
out.timing_info.plotting_tics = uint64([]);
out.timing_info.plotting_tocs = [];

out.spike_time_trace = [];

debug = [];
disp('starting simulation');
for sec = 1 : net.sim_time_sec
    out.timing_info.sim_sec_times(sec) = tic;
    spike_time_trace = [];
    
    % Trim data into seconds to speed searching later
    idxs = net.ts > (sec - 1) * 1000 & net.ts <= (sec * 1000);
    inp_trimmed = net.inp(idxs);
    ts_trimmed = net.ts(idxs);
    
    for ms = 1 : ms_per_sec
        
        time = (sec - 1) * ms_per_sec + ms;
        
        % Simulated annealing
%         if time < Tf * ms_per_sec
%             m = (I0 - If) / (Tf * ms_per_sec);
%             net.fgi = net.fgi - m;
%         end
        
        %% Calculate input
        Iapp = zeros(size(v));
        t0 = time - last_spike_time;
        t0_negu = t0 - round(delays);
        p = net.fgi ./ sqrt(2 * pi * variance);
        g = p .* exp(- (t0_negu .^ 2) ./ (2 * variance));
        gaussian_values = w .* g;
        gaussian_values(isnan(gaussian_values)) = 0;
        Iapp = sum(gaussian_values, 1)';

        %% Update membrane voltages  
        v = v + (net.v_rest + Iapp - v) / net.neuron_tau;
        vt(:, time) = v(net.voltages_to_save);
        
        %% Deal with neurons that just spiked
        fired_naturally = find(v >= v_thres);
        fired_pixels = inp_trimmed(ts_trimmed == time);

        if net.supervising
            %  if      supervised                and       it fired recently         or     is firing now.
            if numel(find(fired_pixels == 4)) > 0 && (time - last_spike_time(4) < 30 || numel(find(fired_naturally == 4)) > 0)   % Only supervise if we havent seen a pixel 4 fire recently.
                fired_pixels(find(fired_pixels == 4)) = [];
                fprintf('Remove supervision %d\n', sec*1000 + ms);
            end
                    %  if      supervised                and       it fired recently         or     is firing now.
            if numel(find(fired_pixels == 5)) > 0 && (time - last_spike_time(5) < 30 || numel(find(fired_naturally == 5)) > 0)   % Only supervise if we havent seen a pixel 5 fire recently.
                fired_pixels(find(fired_pixels == 5)) = [];
                fprintf('Remove supervision %d\n', sec*1000 + ms);
            end
        end
        fired = [fired_naturally; fired_pixels';]; % unique in case of SDVL supervision later
        spike_time_trace = [spike_time_trace; time*ones(length(fired),1), fired];
        last_spike_time(fired) = time; 
        
        v(fired) = net.v_reset;
        v_thres(fired) = v_thres(fired) + net.thres_rise;
        %debug = [debug; v_thres'];
        
        
        % lateral inhibition
        %if sum(fired > N_inp) > 0
        %    v(4:5) = net.v_reset;
        %end
        
        %% STDP
        % Any pre-synaptics weights should be increased
        w(:, fired) = w(:, fired) + (dApost(:, fired) .* conns(:, fired));
        % Any post synaptic weights decrease
        w(fired, :) = w(fired, :) + (dApre(fired, :) .* conns(fired, :));
        dApost(fired, :) = dApost(fired, :) + net.Apost;
        dApre(:, fired) = dApre(:, fired) + net.Apre;
        
        %% SDVL
        % TODO - consider if this should be done at the start of the ms or
        % the end (here).
        %[pre_idxs, post_idxs] = ind2sub(size(conns), find(conns(:, fired))); 
        t0 = repmat(time - last_spike_time, 1, numel(fired));
        t0_negu = t0 - delays(:, fired);  % TODO - check this works with multiple fired
        abst0_negu = abs(t0_negu);
        k = (variance(:, fired) + 0.9) .^2;
        shifts = sign(t0_negu) .* k .* net.nu;
        
        % Update SDVL mean
        du = zeros(size(t0_negu));
        du(t0 >= net.a2) = -k(t0 >= net.a2) .* net.nu;
        du(abst0_negu >= net.a1) = shifts(abst0_negu >= net.a1);% TODO: verify this line is correct, made an edit without checkign the maths.

        delays(:, fired) = delays(:, fired) + (du .* conns(:, fired));
        delays(conns) = max(1, min(net.delay_max, delays(conns)));

        % Update SDVL variance
        dv = zeros(size(t0_negu));
        dv(abst0_negu <= net.b2) = -k(abst0_negu <= net.b2) .* net.nv;
        dv(abst0_negu >= net.b1) = k(abst0_negu >= net.b1) .* net.nv;

        variance(:, fired) = variance(:, fired) + (dv .* conns(:, fired));
        variance(conns) = max(net.variance_min, min(net.variance_max, variance(conns)));

        % STDP decay
        dApre = dApre * STDPdecaypre;
        dApost = dApost * STDPdecaypost;
        active_idx = mod(active_idx, net.delay_max) + 1;
    
        % Intrinsic plasticity threshold decay
        v_thres = v_thres - (net.thres_rise * net.thres_freq / ms_per_sec);
        %v_threst(:, time) = v_thres(:, net.v_thres_to_save)';
        
        % Synaptic bounding - limit w to [0, w_max]
        w = max(0, min(net.w_max, w)); 
        
        % Save variables
        % NOTE: variances saved have rotated to what is usual
        % that is, vart(1, :, :) is the variances of the first variance to
        % save. (as opposed to vart(:, 1, :) which might seem sensible).
        delayst(:, :, time) = delays(:, net.delays_to_save)';
        vart(:, :, time) = variance(:, net.variance_to_save)';
        
    end  % of ms for loop

    
out.spike_time_trace = [out.spike_time_trace; spike_time_trace]; % TODO - optimise for speed if necessary
out.timing_info.full_sec_tocs(sec) = toc(out.timing_info.sim_sec_times(sec));
if mod(sec, 1) == 0 && net.print_progress
    fprintf('Sim sec: %d, Real sec: %.3f \n', sec, toc(out.timing_info.init_time));
    clf;
    % TODO : Fix plotting
    %rasterspiketimes(spike_time_trace, 2000, 1);
    drawnow;
end

end    % of seconds for loop
if net.print_progress
    fprintf('Real seconds per simulated: %.3f \n', toc(out.timing_info.init_time) / sec);
end
out.vt = vt;
out.vart = vart;
out.delayst = delayst;
%out.v_threst = v_threst;
out.w = w;
out.variance = variance;
out.delays = delays;
out.debug = debug;

end















