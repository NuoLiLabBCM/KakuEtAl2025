function tcross = crossing_times_phase(t, theta_unwrapped, target_phase)
% Find times where unwrapped phase crosses target_phase mod 2*pi upward.
%
% target_phase = pi gives mid-cycle phase crossings.

    t = t(:);
    theta_unwrapped = theta_unwrapped(:);

    n0 = floor((theta_unwrapped(1) - target_phase) / (2*pi));
    n1 = floor((theta_unwrapped(end) - target_phase) / (2*pi));

    tcross = [];

    for n = n0:n1
        level = target_phase + 2*pi*n;

        idx = find(theta_unwrapped(1:end-1) < level & ...
                   theta_unwrapped(2:end) >= level, 1, 'first');

        if ~isempty(idx)
            t1 = t(idx);
            t2 = t(idx+1);
            y1 = theta_unwrapped(idx);
            y2 = theta_unwrapped(idx+1);

            alpha = (level - y1) / (y2 - y1);
            tc = t1 + alpha * (t2 - t1);

            tcross(end+1,1) = tc; %#ok<AGROW>
        end
    end
end