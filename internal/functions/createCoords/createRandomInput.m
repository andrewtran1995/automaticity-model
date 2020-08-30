function [visualInputMatrix] = createRandomInput(trials)
%CREATERANDOMINPUT Create N random visual stimulus points (dictated by
%trials)
    visualInputMatrix = randi(100, trials, 2);
end

