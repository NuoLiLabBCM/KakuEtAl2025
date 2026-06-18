% simulate a toy kuramotor model with 2 separate coupling terms:
% K12: coupling from oscillator 1 to oscillator 2
% K21: coupling from oscillator 2 to oscillator 1
%
% Simulation a range of coupling strength
%
% Analyze: 1) PRC curve; 2) effect of clamping one oscillator
%

clear all
close all

for K_strength = [0 6 10 25]

    K12 = 10;               % influence of oscillator 1 on 2
    K21 = K_strength;       % influence of oscillator 2 on 1
    w1 = 1.5*(2*pi);        % frequency of oscillator 1 - convert freq. in Hz to angular freq.; breathing
    w2 = 7*(2*pi);          % frequency of oscillator 2; licking


    tspan = linspace(0,200,20000);
    [t, theta] = kuramoto_asymmetric(K12, K21, w1, w2, tspan, [0; rand(1)*pi]);


    % Wrap phases to [0, 2*pi)
    theta_wrapped = mod(theta, 2*pi);

    % Plot
    figure;
    subplot(1,5,1);
    plot(t, theta_wrapped(:,1), 'LineWidth', 1.5); hold on;
    plot(t, theta_wrapped(:,2), 'LineWidth', 1.5);
    % plot(t, theta(:,1), 'LineWidth', 1.5); hold on;
    % plot(t, theta(:,2), 'LineWidth', 1.5);

    xlabel('Time');
    ylabel('Wrapped Phase');
    title(['Coupled Kuramoto Oscillators, K12:',num2str(K12),' K21:',num2str(K21)]);
    legend('1-breathing','2-licking');

    % ylim([0 2*pi]);
    xlim([0 3.5])
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

    subplot(1,5,2);
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
    subplot(1,5,4);
    scatter(E.target_phase, E.dphase_shift, 40, 'filled');
    xlabel('Phase of oscillator 1 at oscillator 2 event');
    ylabel('Phase shift of oscillator 1');
    title('Event-triggered PRC-like curve');
    xticks([0 pi 2*pi]);
    xticklabels({'0','\pi','2\pi'});
    ylim([-1 1])
    xlim([0 2*pi])
    grid on;


    % oscillator 2
    E = event_triggered_prc_from_theta_v2(t, theta, ...
        'driver_idx', 1, ...
        'target_idx', 2, ...
        'transient_time', 20, ...
        'baseline_mode', 'coupled_mean_period', ...
        'event_selection', 'all');
    subplot(1,5,5);
    scatter(E.target_phase, E.dphase_shift, 40, 'filled');
    xlabel('Phase of oscillator 2 at oscillator 1 event');
    ylabel('Phase shift of oscillator 2');
    title('Event-triggered PRC-like curve');
    xticks([0 pi 2*pi]);
    xticklabels({'0','\pi','2\pi'});
    xlim([0 2*pi])
    ylim([-1 1])
    grid on;





    %% Effect of clamping one oscillator
    [t, theta] = kuramoto_asymmetric_clamp( ...
        K12, ...      % K12: 1 -> 2
        K21, ...      % K21: 2 -> 1
        w1, w2, ...   % w1, w2
        tspan, ...
        [rand(1)*2*pi; rand(1)*2*pi], ...
        1, 3, ...   % clamp interval [t1, t2]
        2);         % clamp licking oscillator

    % Wrap phases to [0, 2*pi)
    theta_wrapped = mod(theta, 2*pi);

    subplot(1,5,3);
    plot(t, theta_wrapped(:,1), 'LineWidth', 1.5); hold on;
    plot(t, theta_wrapped(:,2), 'LineWidth', 1.5);

    xlabel('Time');
    ylabel('Phase');
    title('Clamp oscillator 1');

    xlim([0 5])
    yticks([0 pi 2*pi]);
    yticklabels({'0','\pi','2\pi'});


end
