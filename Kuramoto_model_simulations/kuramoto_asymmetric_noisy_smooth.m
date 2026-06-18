function [t, theta, omega_hist, omega_target_hist] = kuramoto_asymmetric_noisy_smooth( ...
    K12, K21, w1, w2, sigma1, sigma2, tau_omega, tspan, theta0, ...
    omega1_min, omega1_max, omega2_min, omega2_max)
% 2-oscillator Kuramoto with cycle-to-cycle noisy target frequencies
% and smooth frequency relaxation.
%
% Inputs
%   K12, K21           coupling strength
%   w1, w2             mean intrinsic frequencies - angular frequency  [actual frequency in Hz is w/(2*pi)]
%   sigma1, sigma2     std dev of cycle-to-cycle frequency noise
%   tau_omega          smoothing time constant for frequency changes
%   tspan              time vector, e.g. linspace(0,200,50000)
%   theta0             initial phases [theta1_0; theta2_0]
% Bounded frequency draws fixes an issue where theta takes on negative value
% stuck on slow frequency
%
% Outputs
%   t                  time vector
%   theta              unwrapped phases, size length(t)-by-2
%   omega_hist         actual smooth frequencies used at each time point
%   omega_target_hist  target frequencies after cycle redraws


    if nargin < 8 || isempty(tspan)
        tspan = linspace(0,20,10000);
    end

    if nargin < 9 || isempty(theta0)
        theta0 = [0; pi/2];
    end

    if nargin < 10 || isempty(omega1_min)
        omega1_min = 0.1;
    end
    if nargin < 11 || isempty(omega1_max)
        omega1_max = Inf;
    end
    if nargin < 12 || isempty(omega2_min)
        omega2_min = 0.1;
    end
    if nargin < 13 || isempty(omega2_max)
        omega2_max = Inf;
    end

    t = tspan(:);
    nT = length(t);

    theta = zeros(nT,2);
    omega_hist = zeros(nT,2);
    omega_target_hist = zeros(nT,2);

    theta(1,:) = theta0(:)';

    % Initial target intrinsic frequencies
    omega1_target = bounded_randn(w1, sigma1, omega1_min, omega1_max);
    omega2_target = bounded_randn(w2, sigma2, omega2_min, omega2_max);

    % Initial actual intrinsic frequencies
    omega1 = omega1_target;
    omega2 = omega2_target;

    omega_hist(1,:) = [omega1, omega2];
    omega_target_hist(1,:) = [omega1_target, omega2_target];

    cycle1_prev = floor(theta(1,1)/(2*pi));
    cycle2_prev = floor(theta(1,2)/(2*pi));

    for k = 1:nT-1
        dt = t(k+1) - t(k);

        th1 = theta(k,1);
        th2 = theta(k,2);

        % Smooth relaxation toward targets
        omega1 = omega1 + dt*(omega1_target - omega1)/tau_omega;
        omega2 = omega2 + dt*(omega2_target - omega2)/tau_omega;

        f = @(x1,x2,o1,o2) [ ...
            o1 + K21*sin(x2 - x1); ...
            o2 + K12*sin(x1 - x2)];

        k1 = f(th1, th2, omega1, omega2);
        k2 = f(th1 + 0.5*dt*k1(1), th2 + 0.5*dt*k1(2), omega1, omega2);
        k3 = f(th1 + 0.5*dt*k2(1), th2 + 0.5*dt*k2(2), omega1, omega2);
        k4 = f(th1 + dt*k3(1),     th2 + dt*k3(2),     omega1, omega2);

        theta(k+1,:) = theta(k,:) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)';

        cycle1_new = floor(theta(k+1,1)/(2*pi));
        cycle2_new = floor(theta(k+1,2)/(2*pi));

        if cycle1_new > cycle1_prev
            omega1_target = bounded_randn(w1, sigma1, omega1_min, omega1_max);
            cycle1_prev = cycle1_new;
        end

        if cycle2_new > cycle2_prev
            omega2_target = bounded_randn(w2, sigma2, omega2_min, omega2_max);
            cycle2_prev = cycle2_new;
        end

        omega_hist(k+1,:) = [omega1, omega2];
        omega_target_hist(k+1,:) = [omega1_target, omega2_target];
    end
end


function x = bounded_randn(mu, sigma, xmin, xmax)
    x = mu + sigma*randn;
    while x < xmin || x > xmax
        x = mu + sigma*randn;
    end
end