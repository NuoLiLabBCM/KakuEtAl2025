function E = event_triggered_prc_from_theta_v2(t, theta, varargin)
% EVENT_TRIGGERED_PRC_FROM_THETA
% Treat each cycle of oscillator 2 as a perturbation event acting on oscillator 1.
%
% Inputs
%   t         : Nx1 time vector
%   theta     : Nx2 matrix of UNWRAPPED phases
%
% Optional name-value pairs
%   'driver_idx'       : oscillator generating perturbation events (default = 2)
%   'target_idx'       : oscillator whose response is measured (default = 1)
%   'transient_time'   : ignore events before this time (default = 0)
%   'baseline_mode'    : 'coupled_mean_period' or 'local_previous_period'
%                        default = 'coupled_mean_period'
%   'event_selection'  : 'all', 'first', or 'last'
%                        default = 'all'
%
% Output struct E fields
%   E.event_times          : times of driver cycle events
%   E.target_phase         : phase of target oscillator at each event, in [0, 2*pi)
%   E.target_next_peak     : observed next peak time of target oscillator
%   E.expected_next_peak   : expected next peak time of target oscillator
%   E.dt_shift             : observed - expected next peak time
%   E.dphase_shift         : PRC-like phase shift = -2*pi*dt_shift/T_target
%   E.target_period_mean   : mean period of target oscillator
%   E.driver_peak_times    : cycle times of driver oscillator
%   E.target_peak_times    : cycle times of target oscillator
%
% Notes
%   - theta must be unwrapped phase, not mod(theta,2*pi)
%   - a "peak" means crossing an integer multiple of 2*pi upward
%   - this is PRC-like / event-triggered, not a classical externally imposed PRC

    p = inputParser;
    addParameter(p, 'driver_idx', 2);
    addParameter(p, 'target_idx', 1);
    addParameter(p, 'transient_time', 0);
    addParameter(p, 'baseline_mode', 'coupled_mean_period');
    addParameter(p, 'event_selection', 'all');   % NEW
    parse(p, varargin{:});

    driver_idx = p.Results.driver_idx;
    target_idx = p.Results.target_idx;
    transient_time = p.Results.transient_time;
    baseline_mode = p.Results.baseline_mode;
    event_selection = p.Results.event_selection;

    t = t(:);
    if size(theta,1) ~= numel(t)
        error('theta must have same number of rows as length(t).');
    end
    if size(theta,2) ~= 2
        error('theta must be Nx2.');
    end
    if driver_idx == target_idx
        error('driver_idx and target_idx must be different.');
    end

    th_driver = theta(:, driver_idx);
    th_target = theta(:, target_idx);

    % Detect upward crossings of multiples of 2*pi
    driver_peak_times = crossing_times_2pi(t, th_driver);
    target_peak_times = crossing_times_2pi(t, th_target);

    % Apply transient cutoff
    driver_peak_times = driver_peak_times(driver_peak_times >= transient_time);
    target_peak_times = target_peak_times(target_peak_times >= transient_time);

    if numel(target_peak_times) < 3
        error('Not enough target cycles after transient to estimate response.');
    end

    % Mean target period
    T_target = mean(diff(target_peak_times));

    % ------------------------------------------------------------
    % NEW: optionally subsample driver events within each target cycle
    % ------------------------------------------------------------
    switch lower(event_selection)
        case 'all'
            % keep current behavior, do nothing

        case {'first','last'}
            driver_keep = [];

            % For each target cycle interval [target_peak_times(j), target_peak_times(j+1))
            % find driver events in that interval and keep first or last.
            for j = 1:(numel(target_peak_times)-1)
                t_start = target_peak_times(j);
                t_stop = target_peak_times(j+1);

                idx_in = find(driver_peak_times >= t_start & driver_peak_times < t_stop);

                if isempty(idx_in)
                    continue
                end

                if strcmpi(event_selection, 'first')
                    driver_keep(end+1,1) = idx_in(1); %#ok<AGROW>
                else
                    driver_keep(end+1,1) = idx_in(end); %#ok<AGROW>
                end
            end

            driver_peak_times = driver_peak_times(driver_keep);

        otherwise
            error('Unknown event_selection. Use ''all'', ''first'', or ''last''.');
    end

    nE = numel(driver_peak_times);

    event_times = nan(nE,1);
    target_phase = nan(nE,1);
    target_next_peak = nan(nE,1);
    expected_next_peak = nan(nE,1);
    dt_shift = nan(nE,1);
    dphase_shift = nan(nE,1);

    for k = 1:nE
        te = driver_peak_times(k);

        % Find phase of target oscillator at event time
        th_te = interp1(t, th_target, te, 'linear', 'extrap');
        target_phase(k) = mod(th_te, 2*pi);
        event_times(k) = te;

        % Observed next target peak after event
        idx_next = find(target_peak_times > te, 1, 'first');
        if isempty(idx_next)
            continue
        end
        t_next_obs = target_peak_times(idx_next);
        target_next_peak(k) = t_next_obs;

        % Expected next peak
        switch baseline_mode
            case 'coupled_mean_period'
                idx_prev = find(target_peak_times <= te, 1, 'last');
                if isempty(idx_prev)
                    continue
                end
                t_prev = target_peak_times(idx_prev);
                t_next_exp = t_prev + T_target;

            case 'local_previous_period'
                idx_prev = find(target_peak_times <= te, 1, 'last');
                if isempty(idx_prev) || idx_prev < 2
                    continue
                end
                t_prev = target_peak_times(idx_prev);
                T_local = target_peak_times(idx_prev) - target_peak_times(idx_prev-1);
                t_next_exp = t_prev + T_local;

            otherwise
                error('Unknown baseline_mode.');
        end

        expected_next_peak(k) = t_next_exp;

        % Shift: positive dt means delayed next target peak
        dt_shift(k) = t_next_obs - t_next_exp;

        % Convert to phase shift: positive dphase means advance
        dphase_shift(k) = -2*pi * dt_shift(k) / T_target;
    end

    % Remove NaN entries where next peak could not be assigned
    valid = ~isnan(dphase_shift);

    E = struct();
    E.event_times = event_times(valid);
    E.target_phase = target_phase(valid);
    E.target_next_peak = target_next_peak(valid);
    E.expected_next_peak = expected_next_peak(valid);
    E.dt_shift = dt_shift(valid);
    E.dphase_shift = dphase_shift(valid);
    E.target_period_mean = T_target;
    E.driver_peak_times = driver_peak_times;
    E.target_peak_times = target_peak_times;
    E.driver_idx = driver_idx;
    E.target_idx = target_idx;
    E.transient_time = transient_time;
    E.baseline_mode = baseline_mode;
    E.event_selection = event_selection;   % NEW
end


function tcross = crossing_times_2pi(t, theta_unwrapped)
% Find times where unwrapped phase crosses an integer multiple of 2*pi upward

    n0 = floor(theta_unwrapped(1) / (2*pi));
    n1 = floor(theta_unwrapped(end) / (2*pi));

    tcross = [];

    for n = (n0+1):n1
        level = 2*pi*n;

        idx = find(theta_unwrapped(1:end-1) < level & theta_unwrapped(2:end) >= level, 1, 'first');

        if ~isempty(idx)
            % Linear interpolation
            t1 = t(idx);
            t2 = t(idx+1);
            y1 = theta_unwrapped(idx);
            y2 = theta_unwrapped(idx+1);

            alpha = (level - y1) / (y2 - y1);
            tc = t1 + alpha*(t2 - t1);
            tcross(end+1,1) = tc; %#ok<AGROW>
        end
    end
end