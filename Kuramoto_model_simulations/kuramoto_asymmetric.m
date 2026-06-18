function [t, theta] = kuramoto_asymmetric(K12, K21, w1, w2, tspan, theta0)
% Asymmetric 2-oscillator Kuramoto model
%
% Inputs
% K12: coupling from oscillator 1 to oscillator 2
% K21: coupling from oscillator 2 to oscillator 1
%   w1, w2          mean intrinsic frequencies - angular frequency  [actual frequency in Hz is w/(2*pi)]
%   tspan           vector of time points, e.g. linspace(0,200,20000)
%   theta0          initial phases [theta1_0; theta2_0]
%
% Outputs
%   t               time vector
%   theta           unwrapped phases, size length(t)-by-2

    if nargin < 5 || isempty(tspan)
        tspan = linspace(0,20,10000);
    end

    if nargin < 6 || isempty(theta0)
        theta0 = [0; pi/2];
    end

    ode = @(t,theta) [
        w1 + K21*sin(theta(2) - theta(1));
        w2 + K12*sin(theta(1) - theta(2))
    ];

    opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
    [t, theta] = ode45(ode, tspan, theta0, opts);
end