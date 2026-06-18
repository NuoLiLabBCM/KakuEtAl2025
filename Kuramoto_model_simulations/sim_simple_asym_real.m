% simulate a kuramotor model with 2 separate coupling terms
% K12: coupling from oscillator 1 to oscillator 2
% K21: coupling from oscillator 2 to oscillator 1
%
% simulate nosiy smoothed frequency - mimicking empirical data
%
% Simulation 1 - bidirectional coupling
% Simulation 2 - unidirectional breathing->licking coupling
% Simulation 3 - unidirectional breathing<-licking coupling
%
% Analyze: 1) PRC curve; 2) volcano plots and associated quantifications
%


clear all
close all

for i_sim = 1:3


    if i_sim == 1       % bidirectional
        K12 = 10;                       % influence of oscillator 1 on 2
        K21 = 10;                       % influence of oscillator 2 on 1
        w1 = 1.5*(2*pi);                % frequency of oscillator 1 - convert freq. in Hz to angular freq.
        w2 = 7*(2*pi);                  % frequency of oscillator 2
        sigma1 = 2*(2*pi);              % std dev of cycle-to-cycle frequency noise
        sigma2 = 3*(2*pi);              % std dev of cycle-to-cycle frequency noise
        tau_omega = 2;                  % smoothing time constant for frequency changes
        theta0 = [0; 2*rand(1)*pi];     %initial phases [theta1_0; theta2_0]


    elseif i_sim == 2   % Breathin -> Lick
        K12 = 10;                       % influence of oscillator 1 on 2
        K21 = 0;                        % influence of oscillator 2 on 1
        w1 = 1.5*(2*pi);                % frequency of oscillator 1 - convert freq. in Hz to angular freq.
        w2 = 7*(2*pi);                  % frequency of oscillator 2
        sigma1 = 2*(2*pi);              % std dev of cycle-to-cycle frequency noise
        sigma2 = 3*(2*pi);              % std dev of cycle-to-cycle frequency noise
        tau_omega = 2;                  % smoothing time constant for frequency changes
        theta0 = [0; 2*rand(1)*pi];     %initial phases [theta1_0; theta2_0]

    elseif i_sim == 3   % Lick -> Breath
        K12 = 0;                        % influence of oscillator 1 on 2
        K21 = 10;                       % influence of oscillator 2 on 1
        w1 = 1.5*(2*pi);                % frequency of oscillator 1 - convert freq. in Hz to angular freq.
        w2 = 7*(2*pi);                  % frequency of oscillator 2
        sigma1 = 2*(2*pi);              % std dev of cycle-to-cycle frequency noise
        sigma2 = 3*(2*pi);              % std dev of cycle-to-cycle frequency noise
        tau_omega = 2;                  % smoothing time constant for frequency changes
        theta0 = [0; 2*rand(1)*pi];     %initial phases [theta1_0; theta2_0]

    end

    tspan = linspace(0,500,50000);

    [t, theta, omega_hist, omega_target_hist] = ...
        kuramoto_asymmetric_noisy_smooth( ...
        K12, ...        % K12: 1 -> 2
        K21, ...        % K21: 2 -> 1
        w1, w2, ...     % mean frequencies
        sigma1, sigma2, ... % sigma1, sigma2
        5, ...      % tau_omega
        tspan, ...
        theta0);


    % Wrap phases to [0, 2*pi)
    theta_wrapped = mod(theta, 2*pi);

    % Plot
    figure;
    subplot(2,1,1);
    plot(t, theta_wrapped(:,1), 'LineWidth', 1.5); hold on;
    plot(t, theta_wrapped(:,2), 'LineWidth', 1.5);

    xlabel('Time');
    ylabel('Phase');
    if i_sim==1
        title('Bidirectional Kuramoto Oscillators (Wrapped Phase)');
    elseif i_sim==2
        title('Breath->Lick Kuramoto Oscillators (Wrapped Phase)');
    elseif i_sim==3;
        title('Breath<-Lick Kuramoto Oscillators (Wrapped Phase)');
    end
    legend('\theta_1','\theta_2');

    xlim([0 5])
    yticks([0 pi 2*pi]);
    yticklabels({'0','\pi','2\pi'});




    %% Compute frequency
    % Detect peak times from wrap events
    % A peak occurs just before wrapped phase jumps from near 2*pi back to 0.
    % Detect large negative jumps in wrapped phase.

    dtheta1 = diff(theta_wrapped(:,1));
    dtheta2 = diff(theta_wrapped(:,2));

    % Threshold for wrap detection
    thr = -pi;   % any big negative jump counts as a wrap

    idx1 = find(dtheta1 < thr);
    idx2 = find(dtheta2 < thr);

    % Peak times: use sample just before the wrap
    tpeak1 = t(idx1);
    tpeak2 = t(idx2);

    % Peak phase values, optional
    peak1 = theta_wrapped(idx1,1);
    peak2 = theta_wrapped(idx2,2);

    subplot(2,3,4);
    [y x] = hist(1./diff(tpeak1),0:.5:10); plot(x,y,'LineWidth',1.5,'Color',[0 .6 0]);  hold on
    [y x] = hist(1./diff(tpeak2),0:.5:10); plot(x,y,'LineWidth',1.5,'Color',[.6 0 0]);  hold on
    xlim([0 10]);


    %% PRC
    % oscillator 1
    E = event_triggered_prc_from_theta_v2(t, theta, ...
        'driver_idx', 2, ...
        'target_idx', 1, ...
        'transient_time', 20, ...
        'baseline_mode', 'coupled_mean_period', ...
        'event_selection', 'all');

    subplot(2,3,5); hold on
    plot_prc_with_running_sem(E.target_phase, E.dphase_shift);

    xlabel('Phase of oscillator 1 at oscillator 2 event');
    ylabel('Phase shift of oscillator 1');
    title('Event-triggered PRC-like curve');
    xticks([0 pi 2*pi]);
    xticklabels({'0','\pi','2\pi'});
    ylim([-1.8 1.8])
    grid on;


    % oscillator 2
    E = event_triggered_prc_from_theta_v2(t, theta, ...
        'driver_idx', 1, ...
        'target_idx', 2, ...
        'transient_time', 20, ...
        'baseline_mode', 'coupled_mean_period', ...
        'event_selection', 'all');

    subplot(2,3,6); hold on
    plot_prc_with_running_sem(E.target_phase, E.dphase_shift);

    xlabel('Phase of oscillator 2 at oscillator 1 event');
    ylabel('Phase shift of oscillator 2');
    title('Event-triggered PRC-like curve');
    xticks([0 pi 2*pi]);
    xticklabels({'0','\pi','2\pi'});
    ylim([-1.8 1.8])
    grid on;


    %% volcano plots for aligning phase = pi time of timeseries 2 (lick) relative to peak time of timeseries 1 (breath), sorted by period duration of timeseries1

    % Find oscillator 2 phase = pi crossing times
    tpi2 = crossing_times_phase(t, theta(:,2), pi);

    iti_peak1 = diff(tpeak1);
    [~, i_sorted_peak1] = sort(iti_peak1,'descend');
    i_sorted_peak1 = i_sorted_peak1 + 1;

    dt_peak1_aligned = [];      % [period1# peak1 time]
    dt_pi2_aligned   = [];      % [period1# osc2 pi time]

    i_cycle = 0;

    for k = i_sorted_peak1'

        i_cycle = i_cycle + 1;

        % Plot oscillator 1 cycle boundaries relative to oscillator 1 peak
        peak1_time_tmp = [tpeak1(k-1); tpeak1(k)] - tpeak1(k-1);

        dt_peak1_aligned = cat(1, dt_peak1_aligned, ...
            [repmat(i_cycle, size(peak1_time_tmp,1), 1), peak1_time_tmp]);

        % Find oscillator 2 phase = pi events around this oscillator 1 peak
        i_pi2_tmp = find(tpi2 > (tpeak1(k-1)-1) & tpi2 < (tpeak1(k-1)+2));

        pi2_time_tmp = tpi2(i_pi2_tmp) - tpeak1(k-1);

        dt_pi2_aligned = cat(1, dt_pi2_aligned, ...
            [repmat(i_cycle, size(pi2_time_tmp,1), 1), pi2_time_tmp]);

    end

    %% Plot aligned times
    figure;
    subplot(2,2,1);
    plot(dt_peak1_aligned(:,2), dt_peak1_aligned(:,1), 'xr'); hold on
    plot(dt_pi2_aligned(:,2),   dt_pi2_aligned(:,1),   '.k')

    xlabel('Time relative to oscillator 1 peak');
    ylabel('Oscillator 1 cycles sorted by duration');
    title('Oscillator 2 Phase = \pi Times Aligned to Oscillator 1 Peaks');
    ylim([0 max([dt_pi2_aligned(:,1); dt_peak1_aligned(:,1)])+1])
    xlim([-.4 1])



    %% volcano plots for aligning phase = pi time of timeseries 1 (breath) relative to peak time of timeseries 2 (lick), sorted by period duration of timeseries2

    % Find oscillator 2 phase = pi crossing times
    tpi2 = crossing_times_phase(t, theta(:,2), pi);

    iti_pi2 = diff(tpi2);
    [~, i_sorted_pi2] = sort(iti_pi2, 'descend');
    i_sorted_pi2 = i_sorted_pi2 + 1;

    dt_peak1_aligned = [];   % [cycle# osc1 peak time relative to osc2 pi]
    dt_pi2_aligned   = [];   % [cycle# osc2 pi times relative to osc2 pi]

    i_cycle = 0;

    for k = i_sorted_pi2'

        i_cycle = i_cycle + 1;

        % Align to oscillator 2 phase = pi
        t_ref = tpi2(k-1);

        % Plot current and next oscillator-2 phase = pi times
        pi2_time_tmp = [tpi2(k-1); tpi2(k)] - t_ref;
        dt_pi2_aligned = cat(1, dt_pi2_aligned, ...
            [repmat(i_cycle, size(pi2_time_tmp,1), 1), pi2_time_tmp]);

        % Find oscillator 1 peaks near oscillator 2 phase = pi
        i_peak1_tmp = find(tpeak1 > (t_ref - 0.1) & tpeak1 < (t_ref + 0.4));

        peak1_time_tmp = tpeak1(i_peak1_tmp) - t_ref;

        dt_peak1_aligned = cat(1, dt_peak1_aligned, ...
            [repmat(i_cycle, size(peak1_time_tmp,1), 1), peak1_time_tmp]);
    end

    %% Plot aligned peak times
    subplot(2,2,2);
    plot(dt_pi2_aligned(:,2), dt_pi2_aligned(:,1), 'xr'); hold on
    plot(dt_peak1_aligned(:,2), dt_peak1_aligned(:,1), '.k')

    xlabel('Time relative to oscillator 2 phase = \pi');
    ylabel('Oscillator 2 cycles sorted by pi-to-pi duration');
    title('Oscillator 1 Peak Times Aligned to Oscillator 2 Phase = \pi');
    ylim([0 max([dt_pi2_aligned(:,1); dt_peak1_aligned(:,1)])+1])
    xlim([-0.1 0.4])


    %% Additional volcano-derived analyses
    % Event definitions:
    %   oscillator 1 event = peak of oscillator 1
    %   oscillator 2 event = phase = pi of oscillator 2

    %% 1) Histogram of oscillator 1 peak timing relative to oscillator 2 phase=pi
    %    Separate slow vs fast oscillator 2 cycles based on volcano sorting.

    % Make sure oscillator 2 phase = pi times exist
    tpi2 = crossing_times_phase(t, theta(:,2), pi);

    % Sort oscillator 2 cycles by pi-to-pi duration
    iti_pi2 = diff(tpi2);
    [~, i_sorted_pi2] = sort(iti_pi2, 'descend');
    i_sorted_pi2 = i_sorted_pi2 + 1;

    dt_pi2_aligned_hist   = [];     % [sorted cycle# osc2 pi times relative to osc2 pi]
    dt_peak1_aligned_hist = [];     % [sorted cycle# osc1 peak times relative to osc2 pi]

    i_cycle = 0;

    for k = i_sorted_pi2'

        i_cycle = i_cycle + 1;

        % Reference event: oscillator 2 phase = pi
        t_ref = tpi2(k-1);

        % Oscillator 2 pi-to-pi cycle boundaries
        pi2_time_tmp = [tpi2(k-1); tpi2(k)] - t_ref;
        dt_pi2_aligned_hist = cat(1, dt_pi2_aligned_hist, ...
            [repmat(i_cycle, size(pi2_time_tmp,1), 1), pi2_time_tmp]);

        % Oscillator 1 peaks around this oscillator 2 phase=pi event
        i_peak1_tmp = find(tpeak1 > (t_ref - 0.3) & tpeak1 < (t_ref + 0.6));

        peak1_time_tmp = tpeak1(i_peak1_tmp) - t_ref;
        dt_peak1_aligned_hist = cat(1, dt_peak1_aligned_hist, ...
            [repmat(i_cycle, size(peak1_time_tmp,1), 1), peak1_time_tmp]);
    end

    % Define slow and fast oscillator 2 cycles using sorted volcano rows
    n_cycle = max(dt_pi2_aligned_hist(:,1));

    slow_rows = 1 : round(n_cycle*0.5);       % slower osc2 cycles
    fast_rows = round(n_cycle*0.5)+1 : round(n_cycle);     % faster osc2 cycles

    peak1_slow = dt_peak1_aligned_hist(ismember(dt_peak1_aligned_hist(:,1), slow_rows), 2);
    peak1_fast = dt_peak1_aligned_hist(ismember(dt_peak1_aligned_hist(:,1), fast_rows), 2);

    subplot(2,2,4); hold on

    bin_edges = -0.3:0.005:0.6;
    [y_slow, x_slow] = hist(peak1_slow, bin_edges);
    [y_fast, x_fast] = hist(peak1_fast, bin_edges);

    plot(x_slow, y_slow, 'Color', [.7 .7 .7], 'LineWidth', 1.5);
    plot(x_fast, y_fast, 'b', 'LineWidth', 1.5);

    xlabel('Oscillator 1 peak time relative to oscillator 2 phase=\pi');
    ylabel('Count');
    title('Oscillator 1 peak timing for slow vs fast oscillator 2 cycles');
    legend('Slow oscillator 2 cycles', 'Fast oscillator 2 cycles');
    xlim([-.15 .3]);


    %% 2) Oscillator 2 event frequency vs oscillator 1 frequency
    % Event definitions:
    %   oscillator 1 event = peak of oscillator 1
    %   oscillator 2 event = phase = pi of oscillator 2
    %
    % This is analogous to the old subplot(1,3,2):
    %   x = oscillator 1 frequency
    %   y = oscillator 2 event frequency
    %   then sort by x and bin-average.

    % Make sure oscillator 2 phase = pi times exist
    tpi2 = crossing_times_phase(t, theta(:,2), pi);

    data_all = [];
    % columns:
    % 1: oscillator 1 cycle index
    % 2: number of oscillator 2 phase=pi events within oscillator 1 cycle
    % 3: oscillator 1 interval
    % 4: oscillator 2 phase=pi interval
    % 5: oscillator 1 frequency
    % 6: oscillator 2 event frequency

    for k = 2:length(tpeak1)

        t_start = tpeak1(k-1);
        t_end   = tpeak1(k);

        % Oscillator 1 cycle duration
        iti_peak1_tmp = t_end - t_start;

        % Oscillator 2 phase=pi events inside this oscillator 1 cycle
        i_pi2_in = find(tpi2 > t_start & tpi2 < t_end);
        n_pi2_tmp = numel(i_pi2_in);

        % Need at least one oscillator 2 event to estimate an interval.
        if n_pi2_tmp == 0
            continue

        elseif n_pi2_tmp == 1
            % If only one osc2 event occurs inside the osc1 cycle,
            % use the interval from the previous osc2 event.
            idx = i_pi2_in(1);

            if idx < 2
                continue
            end

            iti_pi2_tmp = tpi2(idx) - tpi2(idx-1);

        else
            % If multiple osc2 events occur inside the osc1 cycle,
            % use their mean interval.
            iti_pi2_tmp = mean(diff(tpi2(i_pi2_in)));
        end

        freq_peak1_tmp = 1 / iti_peak1_tmp;
        freq_pi2_tmp   = 1 / iti_pi2_tmp;

        data_all(end+1,:) = [k, n_pi2_tmp, iti_peak1_tmp, iti_pi2_tmp, ...
            freq_peak1_tmp, freq_pi2_tmp];
    end


    % Smoothed frequency-frequency plot with event-count purity filter
    subplot(2,2,3); hold on

    x = data_all(:,5);        % oscillator 1 frequency
    y = data_all(:,6);        % oscillator 2 phase=pi event frequency
    n_events = data_all(:,2); % number of osc2 events per osc1 cycle

    % Remove invalid values
    valid = isfinite(x) & isfinite(y) & isfinite(n_events);
    x = x(valid);
    y = y(valid);
    n_events = n_events(valid);

    % Sort by oscillator 1 frequency
    [~, i_sort] = sort(x);
    x = x(i_sort);
    y = y(i_sort);
    n_events = n_events(i_sort);

    % Bin-average with event-count purity criterion
    bin_size = 10;
    purity_threshold = 0.8;   % keep bins where dominant event count is >= 8/10

    x_new = [];
    y_new = [];
    n_new = [];
    purity_new = [];

    for i = 1:bin_size:length(x)-bin_size+1

        idx_bin = i:i+bin_size-1;

        x_bin = x(idx_bin);
        y_bin = y(idx_bin);
        n_bin = n_events(idx_bin);

        unique_n = unique(n_bin);

        % Find dominant event-count category in this bin
        counts = zeros(size(unique_n));
        for j = 1:length(unique_n)
            counts(j) = sum(n_bin == unique_n(j));
        end

        [max_count, i_max] = max(counts);
        dominant_n = unique_n(i_max);
        purity = max_count / bin_size;

        % Skip mixed bins without a dominant event-count category
        if purity < purity_threshold
            continue
        end

        x_new(end+1,1) = mean(x_bin);
        y_new(end+1,1) = mean(y_bin);
        n_new(end+1,1) = dominant_n;
        purity_new(end+1,1) = purity;
    end

    unique_n_plot = unique(n_new)';

    for n = unique_n_plot

        idx = n_new == n;

        plot(x_new(idx), y_new(idx), '.k', 'MarkerSize', 18);

    end

    xlabel('Oscillator 1 frequency');
    ylabel('Oscillator 2 phase=\pi event frequency');
    title('Oscillator 2 event frequency vs oscillator 1 frequency');


    xlim([1 4])
    ylim([5.5 8])


end
