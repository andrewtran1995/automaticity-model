classdef HebbianConst
    %HEBBIANCONSTANTS Class holding all Hebbian constants that can be
    %specified per type of neuron that engages in learning.
    % Strengthening occurs if "[voltage integral] > [NMDA]".
    % Weakening occurs if "[NMDA] - [voltage integral] > 0".
    
    properties (SetAccess = immutable)
        COEF (1,1) = 1e-8
        NMDA (1,1) = 400 % Upper threshold
    end
    
    methods
        % Set overrides if given.
        function obj = HebbianConst(coef, nmda)
            if nargin > 0
                obj.COEF = coef;
                obj.NMDA = nmda;
            end
        end
    end
end

