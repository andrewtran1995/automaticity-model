classdef Coord
    %COORD Represents a coordinate of a visual stimulus.
    
    properties
        x (1,1) {mustBeNonnegative}
        y (1,1) {mustBeNonnegative}
        group
    end
    
    methods
        function obj = Coord(x, y, group)
            if nargin > 0
                obj.x = x;
                obj.y = y;
            end
            if nargin > 2
                obj.group = group;
            end
        end
    end
end

