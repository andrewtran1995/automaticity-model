function [ visualInputMatrix ] = createWallisInput(distance)
%CREATEWALLISINPUT Summary of this function goes here
%   Detailed explanation goes here
    x_diff = [-distance, distance];
    
    visualInputMatrix = zeros(6000,3);
    
    for row=1:size(visualInputMatrix,1)
        x=randi(2);
        y=randi(100);
        visualInputMatrix(row,:) = [50-x_diff(x), y, 0];
    end
    
    return

end

