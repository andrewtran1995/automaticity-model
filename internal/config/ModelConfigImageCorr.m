classdef ModelConfigImageCorr < ModelConfig
    % Corresponds to the Wallis, Jonathan D. paper.
    % Observes image correlation effects.
    
    properties (Constant)
        isFROSTEnabled = false
        isCOVISEnabled = false
        isMCLearningEnabled = false
        hasCriterialNoise = false
    end
    
    methods
        function obj = ModelConfigImageCorr()
            obj@ModelConfig();
        end
        
        function [xs, ys, groups] = loadCoords(~)
            loaded_input = load('data/imageCoor/coords.mat');
            xs = loaded_input.visualInput.x;
            ys = loaded_input.visualInput.y;
            groups = loaded_input.visualInput.groups;
        end
        
        function config = doPreprocessing(~)
            return
        end
    end
end

