% simulate an kuramotor model with 2 separate coupling terms
% simulate nosiy smoothed frequency - mimicking empirical data, for paper

% bidirectional coupling, volcano plots come out naturally

clear all
close all


% bidirectional
K12 = 10;                       % influence of oscillator 1 on 2
K21 = 10;                       % influence of oscillator 2 on 1
w1 = 1.5*(2*pi);                % frequency of oscillator 1 - convert freq. in Hz to angular freq.
w2 = 7*(2*pi);                  % frequency of oscillator 2
sigma1 = 2*(2*pi);              % std dev of cycle-to-cycle frequency noise
sigma2 = 3*(2*pi);              % std dev of cycle-to-cycle frequency noise
tau_omega = 2;                  % smoothing time constant for frequency changes
theta0 = [0; 2*rand(1)*pi];     %initial phases [theta1_0; theta2_0]


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
title('Bidirectional Kuramoto Oscillators (Wrapped Phase)');
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
subplot(1,2,1);
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
subplot(1,2,2);
plot(dt_pi2_aligned(:,2), dt_pi2_aligned(:,1), 'xr'); hold on
plot(dt_peak1_aligned(:,2), dt_peak1_aligned(:,1), '.k')

xlabel('Time relative to oscillator 2 phase = \pi');
ylabel('Oscillator 2 cycles sorted by pi-to-pi duration');
title('Oscillator 1 Peak Times Aligned to Oscillator 2 Phase = \pi');
ylim([0 max([dt_pi2_aligned(:,1); dt_peak1_aligned(:,1)])+1])
xlim([-0.1 0.4])




%% Conditional phase/hazard analysis for Kuramoto model
% Model-comparison analysis for a continuously coupled oscillator model.
%
% In the Kuramoto model, the two oscillators have time-varying intrinsic
% frequencies, but their phase evolution is continuously coupled. Thus, the
% probability or timing of lick events can depend on the ongoing breathing
% phase, not only on discrete breath or lick events.
%
% To compare this model fairly with the stall/lockout model, we apply the
% same post-breath exclusion window used for the stall model analysis.

t_breath = tpeak1(:);   % breathing events, oscillator 1 phase reset
t_lick   = tpi2(:);     % lick events, oscillator 2 phase = pi

phi_breath_t = mod(theta(:,1), 2*pi);
phi_lick_t   = mod(theta(:,2), 2*pi);

opts = struct;
opts.n_phase_bins = 18;          % no smoothing

% Exclude the first 75 ms after each breath, matching the exclusion used
% for the stall/lockout model. This removes the immediate post-breath
% window from both models so that any remaining phase dependence reflects
% coupling outside the event-triggered lockout/release period.
opts.exclude_post_breath = 0.075;

opts.lick_mode = 'all';
opts.bout_gap = 0.5;

R_kuramoto = plot_lick_hazard_decomposition( ...
    t(:), ...
    t_lick, ...
    phi_lick_t(:), ...
    phi_breath_t(:), ...
    t_breath, ...
    opts, ...
    'Kuramoto model');
