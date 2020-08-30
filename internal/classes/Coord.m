classdef Coord
    %COORD Represents a coordinate of a visual stimulus.
    
    properties
        x (1,1) {mustBeNonnegative}
        y (1,1) {mustBeNonnegative}
    end
    
    methods
        function obj = Coord(x, y)
            if nargin > 0
                obj.x = x;
                obj.y = y;
            end
        end
    end
end

