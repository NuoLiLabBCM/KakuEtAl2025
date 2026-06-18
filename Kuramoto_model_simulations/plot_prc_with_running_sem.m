function plot_prc_with_running_sem(x, y)
% Plot raw PRC dots plus running/binned mean ± SEM.
%
% x: phase in [0, 2*pi)
% y: phase shift

    x = x(:);
    y = y(:);

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    % Raw dots: small gray
    scatter(x, y, 8, ...
        'MarkerFaceColor', [0.75 0.75 0.75], ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.35);

    % Sort by phase
    [x_sort, idx] = sort(x);
    y_sort = y(idx);

    % Running window parameters
    win_size = 30; %min(15, round(numel(x_sort) * 0.08));   % 8% of points, minimum 15
    step_size = max(1, round(win_size / 4));

    x_mean = [];
    y_mean = [];
    y_sem  = [];

    for i = 1:step_size:(numel(x_sort) - win_size + 1)

        idx_win = i:(i + win_size - 1);

        x_tmp = x_sort(idx_win);
        y_tmp = y_sort(idx_win);

        x_mean(end+1,1) = mean(x_tmp);
        y_mean(end+1,1) = mean(y_tmp);
        y_sem(end+1,1)  = std(y_tmp) / sqrt(numel(y_tmp));
    end

    % SEM shading
    x_patch = [x_mean; flipud(x_mean)];
    y_patch = [y_mean - y_sem; flipud(y_mean + y_sem)];

    patch(x_patch, y_patch, [0 0 0], ...
        'FaceAlpha', 0.18, ...
        'EdgeColor', 'none');

    % Running mean line
    plot(x_mean, y_mean, 'k-', 'LineWidth', 2.5);

    xlim([0 2*pi])
end