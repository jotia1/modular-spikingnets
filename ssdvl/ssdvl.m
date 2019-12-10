function [ out ] = ssdvl( net )
%SSDVL A single neuron with ssdvl on the input synapses
%   Takes a network struct and returns an output object


out = struct();
out.timing_info.init_time = tic;
[isvalid, missing] = validatenetwork(net);
assert(isvalid, 'Network missing field...');

assert( net.nv == net.nu, 'ssdvl needs nu == nv');

if net.rand_seed ~= -1
    rng(net.rand_seed); 
end

ms_per_sec = 1000;
sim_time_ms = net.sim_time_sec * ms_per_sec;

% SDVL variables
delays = net.delays;
delayst = zeros(numel(net.delays_to_save), net.N, sim_time_ms);
variance = net.variance;
vart = zeros(numel(net.variance_to_save), net.N, sim_time_ms);
v_threst = zeros(numel(net.v_thres_to_save), sim_time_ms);
iappt = zeros(numel(net.iapp_to_save), sim_time_ms);

N = net.N;
N_inp = net.group_sizes(1);
v = ones(N, 1) * net.v_rest;
vt = zeros(numel(net.voltages_to_save), sim_time_ms);
last_spike_time = zeros(net.N, 1) * -Inf;

% TODO : current_steps should be based off step size and max delay etc. or 
% we can just make it a hyperparam...???
current_steps = 40;  
upcoming_current = zeros(N, current_steps);
upcur_idx = 1;

% Dynamic threshold parameters
v_thres = ones(size(v)) * net.v_thres;
net.thres_rise = net.thres_rise * net.dynamic_threshold; % zero if false

% Izhikevich params
dt = 1;
d = ones(N, 1) * 8;
a = ones(N, 1) * 0.02;
u = ones(N, 1) * -14;
do_rounding = true;

% STDP variables
% HACK FOR SINGLE NEURON, should be NxN but is quite slow an unnecessary
% right now
assert(net.group_sizes(2) == 1, 'Simulator only supports 1 output neuron');
dApre = zeros(N, 1);
dApost = zeros(1, N);
%dApre = zeros(size(w));
%dApost = zeros(size(w));
STDPdecaypre = exp(-1/net.taupre);
STDPdecaypost = exp(-1/net.taupost);
w = sparse(net.w);
conns = w ~= 0;

%Simulated annealing parameters
I0 = net.fgi;
If = net.If;
Tf = net.Tf;
anneal_gradient = (I0 - If) / (Tf * ms_per_sec);
anneal_freq_ms = 500;

if net.fixed_integrals
    variance_precision = 0.01;
    % TODO : Rename ptable to something more representative
    [ ptable ] = buildlookuptable(net.variance_min : variance_precision : net.variance_max, ...
                                1 : net.delay_max, ...
                                1 : current_steps, ...
                                net.fgi);
end

% output variables
out.timing_info.init_toc = toc(out.timing_info.init_time);
out.timing_info.sim_sec_tics = uint64(zeros(net.sim_time_sec, 1));
out.timing_info.sim_sec_tocs = zeros(net.sim_time_sec, 1);
out.timing_info.full_sec_tocs = zeros(net.sim_time_sec, 1);
out.timing_info.plotting_tics = uint64([]);
out.timing_info.plotting_tocs = [];
out.timing_info.profiling_tocs = zeros(10, 10 * 1000);  % Total of tocs

out.spike_time_trace = [];
out.offsets = [];
debug = zeros(sim_time_ms, 0);

if net.print_progress
    disp('Starting simulation');
end

for sec = 1 : net.sim_time_sec
    out.timing_info.sim_sec_times(sec) = tic;
    spike_time_trace = [];
    
    if isa(net.data_generator,'function_handle')
        [inp_trimmed, ts_trimmed, ~, ~, offsets] = net.data_generator();
        ts_trimmed = ts_trimmed + (sec -1) * 1000;
        out.offsets = [out.offsets, offsets + (sec -1) * 1000];
    else
        % Trim data into seconds to speed searching later
        idxs = net.ts > (sec - 1) * 1000 & net.ts <= (sec * 1000);
        inp_trimmed = net.inp(idxs);
        ts_trimmed = net.ts(idxs);
    end
    
    for ms = 1 : ms_per_sec
        
        ms_tic = tic;
        
        time = (sec - 1) * ms_per_sec + ms;
        
        % Simulated annealing
        if net.use_simulated_annealing && mod(time, anneal_freq_ms) == 0 && time < (Tf * ms_per_sec)
            net.fgi = net.fgi - (anneal_gradient * anneal_freq_ms);
            [ ptable ] = buildlookuptable(net.variance_min : variance_precision : net.variance_max, ...
                                1 : net.delay_max, ...
                                1 : current_steps, ...
                                net.fgi);       
        end
        
        %% Calculate input
        Iapp = upcoming_current(:, upcur_idx);
        %debug(time, 4) = Iapp(sum(net.group_sizes));
        
        %% TIMER
        out.timing_info.profiling_tocs(1, time) = toc(ms_tic);
        
        %% Update membrane voltages
        if net.use_izhikevich
            dv = (0.04 * v + 5) .* v + 140 - u;
            v = v + (dv + Iapp) * dt;
            du = a .* (0.2 * v - u);
            u = u + dt * du;
        else
            v = v + ((net.v_rest - v) / net.neuron_tau) + Iapp;
        end
        vt(:, time) = v(net.voltages_to_save);
        
        %% TIMER 
        out.timing_info.profiling_tocs(2, time) = toc(ms_tic);
        
        %% Deal with neurons that just spiked
        fired_naturally = find(v >= v_thres);
        fired_pixels = inp_trimmed(ts_trimmed == time);
        
        if net.use_izhikevich
            % TODO : consider, should be be applied too pixels (inputs) too?
            vt(fired_naturally, time) = 35;
            u(fired_naturally) = u(fired_naturally) + d(fired_naturally);
        end

        if net.supervising
            % TODO : Supervision needs to be made general.
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
        
        %% TIMER
        out.timing_info.profiling_tocs(3, time) = toc(ms_tic);
        
        %% Update upcoming current based on who spiked
        if do_rounding
            fired_delays = round(delays(fired, :));
        else
            fired_delays = delays(fired, :);
        end
        
        %% TIMER
        out.timing_info.profiling_tocs(4, time) = toc(ms_tic);
        fired_values = [];
        
            fired_var = round(variance(fired, :), 2);
            fired_conns = find(fired_delays);
            % TODO : hardcoded offset (-9) for variable table look up
            conn_variances = round((round(fired_var(fired_conns), 2) / variance_precision) - 9);
            conn_delays = fired_delays(fired_conns);
            fired_w = w(fired, :);
            conn_w = fired_w(fired_conns);
            
            
        if numel(conn_variances) > 0
            fired_idxs = sub2ind(size(ptable), ...
                        repmat(conn_variances, current_steps, 1), ...
                        repmat(conn_delays, current_steps, 1), ...
                        reshape(repmat(1:current_steps, numel(conn_delays), 1), [], 1));

            fired_values = ptable(fired_idxs);
        end
        
        stepped_current = conn_w .* reshape(fired_values, numel(fired_conns), []);
        output_current = sum(stepped_current);
        
        weighted_gauss_samples = zeros(N, current_steps);
        weighted_gauss_samples(N, :) = output_current;
        weighted_gauss_samples(isnan(weighted_gauss_samples)) = 0;
        
        upcoming_current(:, upcur_idx) = 0;
        upcur_idx = mod(upcur_idx, current_steps) + 1;
        idx_diff = current_steps - upcur_idx + 1;
        upcoming_current(:, upcur_idx:end) = upcoming_current(:, upcur_idx:end) + weighted_gauss_samples(:, 1 : idx_diff);
        upcoming_current(:, 1 : upcur_idx - 1) = upcoming_current(:, 1 : upcur_idx - 1) + weighted_gauss_samples(:, idx_diff + 1 : end);
   
        v(fired) = net.v_reset;
        v_thres(fired) = v_thres(fired) + net.thres_rise;
        %debug = [debug; v_thres'];
        
        %% TIMER
        out.timing_info.profiling_tocs(5, time) = toc(ms_tic);
        
        %% lateral inhibition
        if numel(intersect(fired, net.lateral_inhibition)) > 0 
            v(net.lateral_inhibition) = net.v_reset;         
        end
        
        % if testing
        if sec < (net.sim_time_sec - net.test_seconds)
        
            %% STDP
            % Any pre-synaptics weights should be increased
            w(:, fired) = w(:, fired) + (repmat(dApost(fired), net.N, 1) .* conns(:, fired));
            % Any post synaptic weights decrease
            w(fired, :) = w(fired, :) + (repmat(dApre(fired), 1, net.N) .* conns(fired, :));
            dApost(fired) = dApost(fired) + net.Apost;
            dApre(fired) = dApre(fired) + net.Apre;

            % STDP decay
            dApre = dApre * STDPdecaypre;
            dApost = dApost * STDPdecaypost;
        
        end
        
        %% TIMER
        out.timing_info.profiling_tocs(6, time) = toc(ms_tic);
        
        %% SDVL
        % if testing
        if sec < (net.sim_time_sec - net.test_seconds)
            
            t0 = repmat(time - last_spike_time, 1, numel(fired));
            t0_negu = t0 - delays(:, fired);
            abst0_negu = abs(t0_negu);
            k = (variance(:, fired) + 0.9) .^2;
            shifts = sign(t0_negu) .* k .* net.nu;

            % Update SDVL mean
            du = zeros(size(t0_negu));
            %du(t0 >= net.a2) = -k(t0 >= net.a2) .* net.nu;
            du(t0 >= 1) = -k(t0 >= 1) .* net.nu;
            du(abst0_negu >= net.a1) = shifts(abst0_negu >= net.a1);

            delays(:, fired) = delays(:, fired) + (du .* conns(:, fired));
            delays(conns) = max(1, min(net.delay_max, delays(conns)));

            % Update SDVL variance
            dvar = (ones(size(t0_negu)) .* -k) .* net.nv;  % TODO : check k is the right size and this operation is happening nicely. 
            %dvar(abst0_negu <= net.b2) = -k(abst0_negu <= net.b2) .* net.nv;
            dvar(abst0_negu >= net.b1) = k(abst0_negu >= net.b1) .* net.nv;

            variance(:, fired) = variance(:, fired) + (dvar .* conns(:, fired));
            variance(conns) = max(net.variance_min, min(net.variance_max, variance(conns)));
        
        end
            
        %% TIMER
        out.timing_info.profiling_tocs(7, time) = toc(ms_tic); 
        
        %% TIMER
        out.timing_info.profiling_tocs(8, time) = toc(ms_tic);
    
        %% Intrinsic plasticity threshold decay
        v_thres(N_inp + 1 :end) = v_thres(N_inp + 1 :end) - (net.thres_rise * net.thres_freq / ms_per_sec);
        v_threst(:, time) = v_thres(net.v_thres_to_save)';
        
        %% Synaptic bounding - limit w to [0, w_max]
        w = max(0, min(net.w_max, w)); 
        
        %% Save variables
        % NOTE: variances saved have rotated to what is usual
        % that is, vart(1, :, :) is the variances of the first variance to
        % save. (as opposed to vart(:, 1, :) which might seem sensible).
        delayst(:, :, time) = delays(:, net.delays_to_save)';
        vart(:, :, time) = variance(:, net.variance_to_save)';
        iappt(:, time) = Iapp(net.iapp_to_save);
        
        %% TIMER
        out.timing_info.profiling_tocs(9, time) = toc(ms_tic);
        
        %% TIMER
        out.timing_info.profiling_tocs(10, time) = toc(ms_tic);
        

        
    end  % of ms for loop
    
    out.spike_time_trace = [out.spike_time_trace; spike_time_trace];
    out.timing_info.full_sec_tocs(sec) = toc(out.timing_info.sim_sec_times(sec));
    if mod(sec, 2) == 0 && net.print_progress
        fprintf('Sim sec: %d, Total real: %.3f, %.3f sec/ss \n', sec, toc(out.timing_info.init_time), toc(out.timing_info.sim_sec_times(1))/sec);
        % Can do optional plotting here later
    end
end    % of seconds for loop

if net.print_progress
    fprintf('Sim sec: %d, Total real: %.3f, %.3f sec/ss + %.3f seconds of intialisation\n', sec, toc(out.timing_info.init_time), toc(out.timing_info.sim_sec_times(1))/sec, out.timing_info.init_toc);
end

out.vt = vt;
out.vart = vart;
out.delayst = delayst;
out.v_threst = v_threst;
out.iappt = iappt;
out.w = w;
out.variance = variance;
out.delays = delays;
out.debug = debug;

end



function [ p ] = fixedintegrals(net, vars, sample_starts)
    % This method was updated and then deprecated for the much faster
    % look up table for getting integrals.
    assert(false, 'deprecated');
    %% Fixed Integrals: Update peak values for any that fired
    p = net.fgi ./ sqrt(2 * pi * vars);  
    epsilon = 0.1;
    do = true;
    small_peaks = 0;
    big_peaks = 0;
    while do
        p = p + (0.005 .* small_peaks);
        full_integral = zeros(size(sample_starts));      
        p = p - (0.005 .* big_peaks);     % Also do bigs
        
        for j = 1 : 40
            % Note: Converting exp( ... ) to a large constant matrix may be
            % faster than looping. (esp. twice...)
            step_integral = p .* exp(- ((sample_starts + j) .^ 2) ./ (2 * vars));
            full_integral = full_integral + step_integral;
        end
        small_peaks = full_integral < (net.fgi - epsilon);
        big_peaks = full_integral > (net.fgi + epsilon);
        do = any(small_peaks(:)) | any(big_peaks(:));
    end
    p(isinf(p)) = 0;

end



function [ ptable ] = buildlookuptable(var_range, delays_range, steps_range, fgi)
%% BUILDLOOKUPTABLE - Table for postsynaptic currents for given delay and variance
%   Build a 3D table, where rows are a range of variances, columns are a
%   range of delays and depths are an amount of current to deliver at that
%   time step for the given delay and variance. 
%   
%   Params:
%       The ranges to use for each parameter
%
%   Example usage:
%       [ ptable ] = buildlookuptable(0.1 : 0.01 : 10, 1 : 1 : 20, 1 : 40)
%

    accuracy = 1e-4;
    ptable_filename = [erase(sprintf('ptable_%.4f_%.0e.mat', fgi, accuracy), '.'), '.mat'];
    
    if exist(ptable_filename, 'file') == 2
        load(ptable_filename, 'ptable');
        fprintf('Loading ptable: %s\n', ptable_filename);
        return
    else
        fprintf('Creating ptable: %s\n', ptable_filename);
    end

    steps = repmat(reshape(steps_range, 1, 1, []), numel(var_range), numel(delays_range), 1);
    delays = repmat(reshape(delays_range, 1, [], 1), numel(var_range), 1, numel(steps_range));
    variances = repmat(reshape(var_range, [], 1, 1), 1, numel(delays_range), numel(steps_range));

    exptable = exp(- ((steps - delays) .^ 2) ./ (2 * variances));
    sample_starts = squeeze(-round(delays(:, :, 1)));
    vars = squeeze(variances(:, :, 1));

    %% Fixed Integrals: Update peak values for any that fired
    peaks = [];
    p = fgi ./ sqrt(2 * pi * vars);  
    epsilon = fgi * accuracy;
    do = true;
    small_peaks = 0;
    big_peaks = 0;
    adjustment_term = fgi * accuracy * 0.05;
    count = 0;
    show_plot = false;
    while do
        count = count + 1;
        p = p + (adjustment_term .* small_peaks);
        full_integral = zeros(size(sample_starts));      
        p = p - (adjustment_term .* big_peaks);     % Do bigs in parallel

        if show_plot && mod(count, 100) == 0
           %[mean(p(:)), sum(small_peaks(:)), sum(big_peaks(:))]
           peaks = [peaks; sum(small_peaks(:)), sum(big_peaks(:))];
           plot(peaks);
           legend({'small', 'big'});
           drawnow();
        end

        for j = 1 : 40
            % Note: Converting exp( ... ) to a large constant matrix may be
            % faster than looping. (esp. twice...)
            step_integral = p .* exp(- ((sample_starts + j) .^ 2) ./ (2 * vars));
            full_integral = full_integral + step_integral;
        end
        small_peaks = full_integral < (fgi - epsilon);
        big_peaks = full_integral > (fgi + epsilon);
        do = any(small_peaks(:)) | any(big_peaks(:));
    end
    %p(isinf(p)) = 0;

    ptable = repmat(p, 1, 1, 40) .* exptable;
    try  
        save(ptable_filename, 'ptable');
        fprintf('Saved ptable: %s\n', ptable_filename);
    catch % In case two threads both try saving
        fprintf('Error saving ptable: %s \n', ptable_filename);
    end
end






