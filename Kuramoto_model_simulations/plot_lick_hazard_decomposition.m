function R = plot_lick_hazard_decomposition(t, lick_events, phi_lick_t, ...
    phi_breath_t, breath_events, opts, plot_title)

    t = t(:);
    lick_events = lick_events(:);
    breath_events = breath_events(:);
    phi_lick_t = phi_lick_t(:);
    phi_breath_t = phi_breath_t(:);

    if ~isfield(opts, 'n_phase_bins')
        opts.n_phase_bins = 18;
    end

    if ~isfield(opts, 'exclude_post_breath')
        opts.exclude_post_breath = 0.050;
    end

    if ~isfield(opts, 'lick_mode')
        opts.lick_mode = 'all';
    end

    if ~isfield(opts, 'bout_gap')
        opts.bout_gap = 0.5;
    end

    dt = median(diff(t));

    %% Select lick events

    switch lower(opts.lick_mode)

        case 'all'
            target_licks = lick_events;

        case 'initiation'
            is_init = [true; diff(lick_events) > opts.bout_gap];
            target_licks = lick_events(is_init);

        otherwise
            error('opts.lick_mode must be all or initiation');
    end

    %% Binary lick-event vector

    y = false(size(t));

    idx_event = round((target_licks - t(1)) / dt) + 1;
    idx_event = idx_event(idx_event >= 1 & idx_event <= numel(t));

    y(idx_event) = true;

    %% Valid time bins

    valid = isfinite(phi_lick_t) & isfinite(phi_breath_t);

    valid = valid & t >= min(breath_events) & t <= max(breath_events);
    valid = valid & t >= min(lick_events)   & t <= max(lick_events);

    % Exclude post-breath window to remove immediate lockout/release artifact.
    % This is a time-based exclusion, not a phase-bin exclusion.
    if isfield(opts, 'exclude_post_breath') && opts.exclude_post_breath > 0

        ts_breath = time_since_previous_event(t, breath_events);

        in_post_breath_exclusion = ...
            ts_breath >= 0 & ts_breath < opts.exclude_post_breath;

        valid = valid & ~in_post_breath_exclusion;
    end

    % For bout-initiation analysis, only include time bins in which a new
    % lick bout could occur.
    if strcmpi(opts.lick_mode, 'initiation')

        ts_lick = time_since_previous_event(t, lick_events);

        at_risk = isnan(ts_lick) | ts_lick >= opts.bout_gap;

        % Keep the event bins themselves.
        valid = valid & (at_risk | y);
    end

    %% Phase bins

    nbin = opts.n_phase_bins;

    edges = linspace(0, 2*pi, nbin + 1);
    centers = edges(1:end-1) + diff(edges)/2;

    lick_phase_bin = discretize(mod(phi_lick_t, 2*pi), edges);
    breath_phase_bin = discretize(mod(phi_breath_t, 2*pi), edges);

    %% 1. Raw lick rate as a function of breathing phase

    count_breath = zeros(nbin, 1);
    occ_breath = zeros(nbin, 1);

    for j = 1:nbin

        idx = valid & breath_phase_bin == j;

        count_breath(j) = sum(y(idx));
        occ_breath(j) = sum(idx) * dt;
    end

    raw_rate_by_breath = count_breath ./ occ_breath;

    %% 2. Own-phase effect: lick rate as a function of lick phase

    count_lick = zeros(nbin, 1);
    occ_lick = zeros(nbin, 1);

    for i = 1:nbin

        idx = valid & lick_phase_bin == i;

        count_lick(i) = sum(y(idx));
        occ_lick(i) = sum(idx) * dt;
    end

    own_rate_by_lick = count_lick ./ occ_lick;

    %% 3. Expected breathing-phase rate from own-phase effect alone

    % Build 2D occupancy matrix: own lick phase x breathing phase.
    occ_2d = zeros(nbin, nbin);
    count_2d = zeros(nbin, nbin);

    for i = 1:nbin
        for j = 1:nbin

            idx = valid & lick_phase_bin == i & breath_phase_bin == j;

            occ_2d(i,j) = sum(idx) * dt;
            count_2d(i,j) = sum(y(idx));
        end
    end


    % For each breathing phase bin, ask what lick rate would be expected
    % if only the distribution of own lick phase mattered.
    expected_count_from_own = zeros(nbin, 1);
    expected_rate_from_own = zeros(nbin, 1);

    % convert the own-lick-phase effect into the lick rate you would expect
    % at each breathing phase, assuming breathing phase has no additional effect.
    for j = 1:nbin

        % For each lick phase bin i:
        %   expected count in joint bin (i,j)
        %       = time spent in joint bin (i,j)
        %         × lick rate expected from lick phase bin i.
        % sum across lick phase bins to get the total expected number of licks in breathing phase bin j.
        expected_count_from_own(j) = sum(occ_2d(:,j) .* own_rate_by_lick, ...
            'omitnan');

        expected_rate_from_own(j) = expected_count_from_own(j) / ...
            sum(occ_2d(:,j), 'omitnan');
    end

    residual_rate_by_breath = raw_rate_by_breath - expected_rate_from_own;

    ratio_rate_by_breath = raw_rate_by_breath ./ expected_rate_from_own;

    %% Store outputs

    R = struct;

    R.phase_centers = centers(:);

    R.raw_rate_by_breath = raw_rate_by_breath;
    R.count_breath = count_breath;
    R.occ_breath = occ_breath;

    R.own_rate_by_lick = own_rate_by_lick;
    R.count_lick = count_lick;
    R.occ_lick = occ_lick;

    R.expected_rate_from_own = expected_rate_from_own;
    R.residual_rate_by_breath = residual_rate_by_breath;
    R.ratio_rate_by_breath = ratio_rate_by_breath;

    R.occ_2d = occ_2d;
    R.count_2d = count_2d;

    R.target_licks = target_licks;
    R.valid = valid;
    R.opts = opts;

    %% Plot

    figure;
    plot(centers, ratio_rate_by_breath, 'ko-', ...
        'MarkerFaceColor', 'k', ...
        'LineWidth', 1.2);
    yline(0, 'k--');
    xlim([0 2*pi]);
    ylim([0.5 1.5]);
    line([0 2*pi],[1 1],'color','k','linestyle','--');
    xticks([0 pi 2*pi]);
    xticklabels({'0','\pi','2\pi'});
    xlabel('Breathing phase');
    ylabel('Normalized lick rate');
    title('Lick rate relative to the mean lick rate');
    box off;

    sgtitle(plot_title);
end