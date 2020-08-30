classdef ModelConfigDualTask < ModelConfig
    % Corresponds to the Zeithamova, D. paper.
    % Observes dual task interference in category learning.
    
    properties (Constant)
        isFROSTEnabled = true
        isCOVISEnabled = true
        isMCLearningEnabled = false
        hasCriterialNoise = true
    end
    
    methods
        function obj = ModelConfigDualTask()
            obj@ModelConfig();
        end
    end
end

