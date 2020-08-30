function [r_x, r_y] = createFMRIInput(trials)
%CREATERANDOMINPUT Create N random visual stimulus points (dictated by
%trials)
    % Ensure points are no closer than 3 units to the boundary
    r_x = randi(100, trials, 1);
    while any(47 < r_x & r_x < 53)
        r_x(47 < r_x & r_x < 53) = randi(100); % TODO: Check if this assigns different values for specified indices.
    end
    r_y = randi(100, trials, 1);
end