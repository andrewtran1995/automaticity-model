classdef ModelConfigDualTask < ModelConfig
    % Corresponds to the Zeithamova, D. paper.
    % Observes dual task interference in category learning.
    
    properties (Constant)
        name = "ModelConfigDualTask"
        isFROSTEnabled = true
        isCOVISEnabled = true
        isMCLearningEnabled = false
        hasCriterialNoise = true
        hasStroopInterference = true
    end

    methods (Static)
        function vals = hebbianValues()
            vals = Hebbian(1e-9, 400);
        end
    end
    
    methods
        function obj = ModelConfigDualTask()
            obj@ModelConfig();
        end
        
        function [x_coords, y_coords, coord_groups] = loadCoords(~)
            coords = createDualTaskCoords();
            x_coords = Coord.xs(coords);
            y_coords = Coord.ys(coords);
            coord_groups = Coord.groups(coords);
        end
        
        function config = doPreprocessing(~)
            return
        end
        
        function dispResults(config, neurons)
            dispResults@ModelConfig(config, neurons);
        end
    end
end

