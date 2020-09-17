classdef ModelConfigButtonSwitch < ModelConfig
    % Corresponds to the Helie, S. paper.
    % Observes a button switch effect with some late-stage dual-task
    % constraints.
    properties (Constant)
        name = "ModelConfigButtonSwitch"
        isFROSTEnabled = false
        isCOVISEnabled = true
        isMCLearningEnabled = true
        hasCriterialNoise = false
        hasStroopInterference = false
    end

    methods (Static)
        function vals = hebbianValues()
            vals = HebbianConst(1e-9, 450);
        end
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
            coord_groups = repmat(Category.NONE, length(x_coords), 1);
        end
        
        function config = setTrials(config, preLearningTrials, learningTrials, postLearningTrials)
            config.preLearningTrials = preLearningTrials;
            config.learningTrials = learningTrials;
            config.postLearningTrials = postLearningTrials;
            config.trials = preLearningTrials + learningTrials + postLearningTrials + config.meta.trialsAfterSwitch;
            config.accuracy = zeros(config.trials,1);
        end
        
        function config = doPreprocessing(~)
            return
        end
        
        function dispResults(config, neurons)
            dispResults@ModelConfig(config, neurons);
        end
        
        function dispButtonSwitchLine(config)
            hold on;
            plot([config.trials - config.meta.trialsAfterSwitch; config.trials - config.meta.trialsAfterSwitch], get(gca,'ylim'), 'r');
        end
    end
end

