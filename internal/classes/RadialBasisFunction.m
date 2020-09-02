classdef RadialBasisFunction
    % Radial Basis Function
    %   Utility class that applies a Radial Basis Function to a visual
    %   stimulus input
    properties
        rbv % grid to keep track of radial basis vector output
    end
    
    properties (Access = private)
        X
        Y
        stimulus_weight
    end
    
    properties (Constant)
        RADIUS = 0.8
    end
    
    methods
        function obj = RadialBasisFunction(GRID_SIZE, stimulus_weight)
            if nargin > 0
                [X, Y] = meshgrid(1:GRID_SIZE, 1:GRID_SIZE);
                obj.X = X;
                obj.Y = Y;
                obj.stimulus_weight = stimulus_weight;
                obj.rbv = zeros(GRID_SIZE);
            end
        end
        
        function obj = resolvestimulus(obj, coord)
            x = coord.x;
            y = coord.y;
            obj.rbv = exp( -(sqrt((y-obj.Y).^2 + (x-obj.X).^2))/obj.RADIUS ) * obj.stimulus_weight;
        end
        
        function dispStimulus(obj, coord)
            colormap('hot');
            imagesc(obj.rbv(ModelConfig.BORDER_SIZE:end-ModelConfig.BORDER_SIZE-1, ModelConfig.BORDER_SIZE:end-ModelConfig.BORDER_SIZE-1));
            title(sprintf('Stimulus: %d, %d', coord.x, coord.y));
        end
    end
    
end