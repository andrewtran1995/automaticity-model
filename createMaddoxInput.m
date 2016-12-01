function [visualInputMatrix] = createMaddoxInput()
%   Create a 6000X3 matrix where the rows are representative of each
%   trial, and the columns are:
%       (1) X coordinate
%       (2) Y coordinate
%       (3) "Group" that coordinate belongs to (i.e., short, medium, or
%       long)
    x_values = [5, 23, 41, 59, 77, 95];
    y_values = x_values;
    grouping = ['L', 'M', 'S', 'S', 'M', 'L'];
    
%     [X, Y] = meshgrid(x_values, y_values);
%     all_points = reshape(cat(2, X', Y'), [], 2);
%     scatter(all_points(:, 1), all_points(:, 2));
%     axis([0 100 0 100]);

    visualInputMatrix = zeros(6000,3);
    
    for row=1:size(visualInputMatrix, 1)
        x = randi(6);
        y = randi(3)*2 - mod(x,2);
        disp(y)
        visualInputMatrix(row,:) = [x_values(x), y_values(y), grouping(x)];
    end
    
    scatter(visualInputMatrix(:,1), visualInputMatrix(:,2));
    axis([0 100 0 100]);
    
    return;
end

