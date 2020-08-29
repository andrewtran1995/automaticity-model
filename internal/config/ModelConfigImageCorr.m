classdef ModelConfigImageCorr < ModelConfig
    % Corresponds to the Wallis, Jonathan D. paper.
    % Observes image correlation effects.
    
    methods
        function obj = ModelConfigImageCorr()
            obj@ModelConfig(false, false);
        end
    end
end

