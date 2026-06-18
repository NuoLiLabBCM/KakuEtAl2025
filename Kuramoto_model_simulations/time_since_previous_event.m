function time_since = time_since_previous_event(t, events)

    t = t(:);
    events = events(:);

    prev_event = interp1(events, events, t, 'previous', NaN);

    time_since = t - prev_event;
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