classdef HebbianConst < handle
    %HEBBIANCONSTANTS Class holding all Hebbian constants that can be
    %specified per type of neuron that engages in learning.
    % Strengthening occurs if "[voltage integral] > [NMDA]".
    % Weakening occurs if "[NMDA] - [voltage integral] - [AMPA] > 0".
    
    properties (SetAccess = immutable)
        COEF
        ANTI
        NMDA % Upper threshold
        AMPA % Lower threshold
    end
    
    methods
        function obj = HebbianConst(coef, anti, nmda, ampa)
            if nargin > 0
                obj.COEF = coef;
                obj.ANTI = anti;
                obj.NMDA = nmda;
                obj.AMPA = ampa;
            end
        end
    end
end

