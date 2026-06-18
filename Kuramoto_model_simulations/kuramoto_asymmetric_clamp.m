function [t, theta] = kuramoto_asymmetric_clamp(K12, K21, w1, w2, tspan, theta0, t1, t2, clamp_idx)
% Asymmetric 2-oscillator Kuramoto model with hard clamping perturbation
%
% K12: coupling from oscillator 1 to oscillator 2
% K21: coupling from oscillator 2 to oscillator 1
% t1, t2: perturbation interval [t1, t2]
% clamp_idx: which oscillator to clamp (1 or 2)
%
% During [t1, t2], the chosen oscillator is held at theta = 0.

    if nargin < 5 || isempty(tspan)
        tspan = linspace(0,20,10000);
    end

    if nargin < 6 || isempty(theta0)
        theta0 = [0; pi/2];
    end

    if nargin < 7 || isempty(t1)
        t1 = 5;
    end

    if nargin < 8 || isempty(t2)
        t2 = 8;
    end

    if nargin < 9 || isempty(clamp_idx)
        clamp_idx = 1;
    end

    if t2 <= t1
        error('t2 must be greater than t1');
    end

    tspan = tspan(:);

    % Split requested time vector into three parts
    idx_pre  = tspan < t1;
    idx_mid  = tspan >= t1 & tspan <= t2;
    idx_post = tspan > t2;

    t_pre  = tspan(idx_pre);
    t_mid  = tspan(idx_mid);
    t_post = tspan(idx_post);

    theta_all = [];
    t_all = [];

    opts = odeset('RelTol',1e-9,'AbsTol',1e-11);

    %% Segment 1: before perturbation
    theta_init = theta0(:);

    if ~isempty(t_pre)
        ode_pre = @(t,theta) [ ...
            w1 + K21*sin(theta(2) - theta(1)); ...
            w2 + K12*sin(theta(1) - theta(2))];

        [tA, thetaA] = ode45(ode_pre, t_pre, theta_init, opts);

        t_all = [t_all; tA];
        theta_all = [theta_all; thetaA];

        theta_init = thetaA(end,:)';
    else
        % integrate to t1 if needed
        if tspan(1) < t1
            ode_pre = @(t,theta) [ ...
                w1 + K21*sin(theta(2) - theta(1)); ...
                w2 + K12*sin(theta(1) - theta(2))];
            [~, theta_tmp] = ode45(ode_pre, [tspan(1) t1], theta_init, opts);
            theta_init = theta_tmp(end,:)';
        end
    end

    %% Segment 2: during perturbation
    % One oscillator is clamped to 0, the other evolves normally.
    if ~isempty(t_mid)
        switch clamp_idx
            case 1
                % theta1 = 0, only evolve theta2
                theta2_init = theta_init(2);

                ode_mid = @(t,theta2) ...
                    w2 + K12*sin(0 - theta2);

                [tB, theta2B] = ode45(ode_mid, t_mid, theta2_init, opts);
                thetaB = [zeros(length(tB),1), theta2B];

            case 2
                % theta2 = 0, only evolve theta1
                theta1_init = theta_init(1);

                ode_mid = @(t,theta1) ...
                    w1 + K21*sin(0 - theta1);

                [tB, theta1B] = ode45(ode_mid, t_mid, theta1_init, opts);
                thetaB = [theta1B, zeros(length(tB),1)];

            otherwise
                error('clamp_idx must be 1 or 2');
        end

        % Avoid duplicate boundary time if needed
        if ~isempty(t_all) && ~isempty(tB) && tB(1) == t_all(end)
            tB = tB(2:end);
            thetaB = thetaB(2:end,:);
        end

        t_all = [t_all; tB];
        theta_all = [theta_all; thetaB];

        theta_init = thetaB(end,:)';
    else
        % if needed, propagate state to t2 through perturbation
        if tspan(end) > t1 && tspan(1) < t2
            switch clamp_idx
                case 1
                    ode_mid = @(t,theta2) w2 + K12*sin(0 - theta2);
                    [~, theta2_tmp] = ode45(ode_mid, [max(tspan(1),t1) t2], theta_init(2), opts);
                    theta_init = [0; theta2_tmp(end)];
                case 2
                    ode_mid = @(t,theta1) w1 + K21*sin(0 - theta1);
                    [~, theta1_tmp] = ode45(ode_mid, [max(tspan(1),t1) t2], theta_init(1), opts);
                    theta_init = [theta1_tmp(end); 0];
            end
        end
    end

    %% Segment 3: after perturbation
    if ~isempty(t_post)
        ode_post = @(t,theta) [ ...
            w1 + K21*sin(theta(2) - theta(1)); ...
            w2 + K12*sin(theta(1) - theta(2))];

        [tC, thetaC] = ode45(ode_post, t_post, theta_init, opts);

        if ~isempty(t_all) && ~isempty(tC) && tC(1) == t_all(end)
            tC = tC(2:end);
            thetaC = thetaC(2:end,:);
        end

        t_all = [t_all; tC];
        theta_all = [theta_all; thetaC];
    end

    t = t_all;
    theta = theta_all;
end