classdef ModelConfigDualTask < ModelConfig
    % Corresponds to the Zeithamova, D. paper.
    % Observes dual task interference in category learning.
    
    methods
        function obj = ModelConfigDualTask()
            obj@ModelConfig(true, true);
        end
    end
end

