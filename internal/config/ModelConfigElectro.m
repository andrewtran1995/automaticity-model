classdef ModelConfigElectro < ModelConfig
    % Corresponds to the Wallis, Jonathan D. paper.
    % Observes image correlation effects.
    
    properties (Constant)
        name = "ModelConfigElectro"
        isFROSTEnabled = false
        isCOVISEnabled = false
        isMCLearningEnabled = false
        hasCriterialNoise = false
        hasStroopInterference = false
    end

    methods (Static)
        function vals = hebbianValues()
            vals = HebbianConst(1e-10, 300);
        end
    end
    
    methods
        function obj = ModelConfigElectro()
            obj@ModelConfig();
        end
        
        function [xs, ys, groups] = loadCoords(~)
            loaded_input = load('data/electro/coords.mat');
            xs = loaded_input.visualInput.x;
            ys = loaded_input.visualInput.y;
%             groups = loaded_input.visualInput.groups;
            groups = repmat(Category.NONE, length(xs), 1);
        end
        
        function config = doPreprocessing(~)
            return
        end
        
        function dispResults(config, neurons)
            dispResults@ModelConfig(config, neurons);
        end
    end
end

