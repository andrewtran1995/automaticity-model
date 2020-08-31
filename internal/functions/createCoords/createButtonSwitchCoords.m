function [xs, ys] = createButtonSwitchCoords(trials)
%Create N random visual stimulus points (dictated by trials). Ensure points
%are no closer than three units to the boundary.
    xs = zeros(trials,1);
    for i=1:trials
        xs(i) = randExcludingMiddle();
    end
    ys = randi(100, trials, 1);
end

function [val] = randExcludingMiddle()
    set = setdiff(1:100, 47:53);
    val = set(randi(numel(set)));
end