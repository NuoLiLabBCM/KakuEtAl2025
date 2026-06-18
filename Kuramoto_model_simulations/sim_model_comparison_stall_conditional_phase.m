%% Simple lick-breath reciprocal stall simulation

clear all; 
close all;

% rng(42);  % reproducible random seed

LICK_FREQ   = 8.0;
BREATH_FREQ = 2.75;
CV_LICK     = 0.12;
CV_BREATH   = 0.24;

STALL    = 0.050;   % seconds
DURATION = 300.0;   % seconds

% Initial event times
tl = period_rand(1 / LICK_FREQ, CV_LICK, 0.02);
tb = period_rand(1 / BREATH_FREQ, CV_BREATH, 0.05);

licks = [];
breaths = [];

while min(tl, tb) < DURATION

    if tl <= tb
        licks(end+1) = tl; %#ok<SAGROW>

        % If breath would occur too soon after lick, delay it
        if tb < tl + STALL
            tb = jitter_time(tl + STALL);
        end

        tl = tl + period_rand(1 / LICK_FREQ, CV_LICK, 0.02);

    else
        breaths(end+1) = tb; %#ok<SAGROW>

        % If lick would occur too soon after breath, delay it
        if tl < tb + STALL
            tl = jitter_time(tb + STALL);
        end

        tb = tb + period_rand(1 / BREATH_FREQ, CV_BREATH, 0.05);
    end
end

licks = licks(:);
breaths = breaths(:);

%% Plot rasters

figure('Color', 'w', 'Position', [100 100 1200 700]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Breaths aligned to lick onset
nexttile;
plot_raster(licks, breaths, ...
    'Breaths aligned to lick onset', ...
    'Time rel. lick onset (ms)', ...
    'ILI');

% Licks aligned to breath onset
nexttile;
plot_raster(breaths, licks, ...
    'Licks aligned to breath onset', ...
    'Time rel. breath onset (ms)', ...
    'IBI');

%% Save figure

% outDir = 'D:\Dropbox\code\analysis\lickingBreathingSim';
% if ~exist(outDir, 'dir')
%     mkdir(outDir);
% end
% 
% outFile = fullfile(outDir, 'lick_breath_raster_simple.png');
% exportgraphics(gcf, outFile, 'Resolution', 150);


%% Convert discrete event times into continuous phase traces

dt = 0.001;  % 1 ms resolution
tvec = (0:dt:DURATION)';

phi_lick   = event_times_to_phase(tvec, licks);
phi_breath = event_times_to_phase(tvec, breaths);

% Optional: convert phase to oscillator-like traces
lick_osc   = sin(phi_lick);
breath_osc = sin(phi_breath);

%% Plot phase traces

figure('Color', 'w', 'Position', [100 100 1000 500]);

plot(tvec, phi_lick, 'LineWidth', 1);
hold on;
plot(tvec, phi_breath, 'LineWidth', 1);

xlim([0 5]);  % show first 5 seconds
ylim([0 2*pi]);

yticks([0 pi 2*pi]);
yticklabels({'0', '\pi', '2\pi'});

xlabel('Time (s)');
ylabel('Phase');
legend({'Lick phase', 'Breath phase'});
box off;

% %% Plot oscillator-like traces
% 
% figure('Color', 'w', 'Position', [100 100 1000 400]);
% 
% plot(tvec, lick_osc, 'LineWidth', 1);
% hold on;
% plot(tvec, breath_osc, 'LineWidth', 1);
% 
% xlim([0 5]);
% xlabel('Time (s)');
% ylabel('sin(phase)');
% legend({'Lick oscillator', 'Breath oscillator'});
% box off;
% 


%% Conditional phase/hazard analysis for stall model
% Model-comparison analysis for an event-gated stall/lockout model.
%
% In the stall model, the two oscillators have independent intrinsic timing,
% but expressed events are reciprocally constrained by a short lockout window:
% a breath transiently prevents licking, and a lick transiently prevents breathing.
% Thus, any dependence of lick timing on breathing in this model should arise
% from the event-triggered lockout/release rule, rather than from continuous
% dependence on breathing phase.

dt = 0.005;
tvec = (0:dt:DURATION)';

phi_lick_t   = event_times_to_phase(tvec, licks);
phi_breath_t = event_times_to_phase(tvec, breaths);

opts = struct;
opts.n_phase_bins = 18;

% Exclude the post-breath lockout/release window.
% STALL is 50 ms. We add a 25 ms buffer because licks that would have occurred
% during the lockout are pushed to approximately breath_time + STALL + jitter.
% If only 0–50 ms is excluded, the first time bin after 50 ms can contain an
% artificial pile-up of released licks, producing a spurious peak in lick rate
% at an early breathing phase.
opts.exclude_post_breath = STALL + 0.025;

opts.lick_mode = 'all';
opts.bout_gap = 0.5;

R_stall = plot_lick_hazard_decomposition( ...
    tvec, ...
    licks(:), ...
    phi_lick_t(:), ...
    phi_breath_t(:), ...
    breaths(:), ...
    opts, ...
    'Stall model');



%% Local functions

function x = period_rand(meanPeriod, cv, minPeriod)
    x = max(meanPeriod + randn * cv * meanPeriod, minPeriod);
end

function t = jitter_time(t)
    t = t + randn * 0.005;
end

function [xs, ys, ivl, rank] = make_raster(ref, evts)

    % Inter-event interval in ms
    ivl = diff(ref) * 1000;

    % Rank trials by interval, with shortest interval at the top
    [~, order] = sort(ivl, 'descend');
    rank = zeros(size(ivl));
    rank(order) = 0:numel(ivl)-1;

    xs = [];
    ys = [];

    for i = 1:numel(ivl)
        rel = (evts - ref(i)) * 1000;
        m = rel >= -500 & rel <= 500;

        xs = [xs; rel(m)]; %#ok<AGROW>
        ys = [ys; repmat(rank(i), sum(m), 1)]; %#ok<AGROW>
    end
end

function plot_raster(ref, evts, plotTitle, xLabel, intervalName)

    [x, y, ivl, rank] = make_raster(ref, evts);

    scatter(x, y, 2, [0.58 0.44 0.86], 'filled', ...
        'MarkerFaceAlpha', 0.8, ...
        'MarkerEdgeAlpha', 0);
    hold on;

    % Plot the next reference interval for each trial
    in_win = ivl >= -500 & ivl <= 500;
    scatter(ivl(in_win), rank(in_win), 10, 'k', 'filled');

    xline(0, 'k--', 'LineWidth', 1.2, 'Alpha', 0.8);

    xlim([-500 500]);

    xlabel(xLabel);
    ylabel(sprintf('Trial sorted by %s, shortest at top', intervalName));
    title(plotTitle);

    box off;
end


function phi = event_times_to_phase(t, events)

    t = t(:);
    events = events(:);

    phi = nan(size(t));

    for i = 1:numel(events)-1
        idx = t >= events(i) & t < events(i+1);

        phi(idx) = 2*pi * ...
            (t(idx) - events(i)) / ...
            (events(i+1) - events(i));
    end

    phi = mod(phi, 2*pi);
end
