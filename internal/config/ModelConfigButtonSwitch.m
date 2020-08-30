classdef ModelConfigButtonSwitch < ModelConfig
    % Corresponds to the Helie, S. paper.
    % Observes a button switch effect with some late-stage dual-task
    % constraints.
    properties (Constant)
        isFROSTEnabled = false
        isCOVISEnabled = true
        isMCLearningEnabled = true
        hasCriterialNoise = false
    end
    
    methods
        function obj = ModelConfigButtonSwitch()
            obj@ModelConfig();
            obj.meta = struct( ...
            'trialsAfterSwitch', 600, ...
            'PMC_A_weights', ones(ModelConfig.GRID_SIZE,ModelConfig.GRID_SIZE,1,4), ...
            'PMC_B_weights', ones(ModelConfig.GRID_SIZE,ModelConfig.GRID_SIZE, 1,4), ...
            'optimization', struct( ...
                'NUM_TRIALS', 11520, 'GROUP_RUN', 0, ...
                'SES_1',      1:480, 'SES_4',    1681:2160, ...
                'SES_10', 5161:5640, 'SES_20', 11041:11520 ...
            ));
        end
        
        function [x_coords, y_coords, coord_groups] = loadCoords(~)
            loaded_input = load('data/buttonSwitch/coords.mat');
            x_coords = loaded_input.x_coordinates;
            y_coords = loaded_input.y_coordinates;
            coord_groups = zeros(length(x_coords), 1);
        end
        
        function config = doPreprocessing(~)
            return
        end
    end
end

