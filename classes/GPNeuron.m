classdef GPNeuron
    %Global Pallidus
    %   Global pallidus. For FROST.
    
    properties
        W_OUT = 1
        out
        spikes = 0
        v
        activations
    end
    
    methods
        function obj = GP(n, trials)
            obj.out = zeros(n,1);
            obj.v = repmat(QIAF.rv,n,1);
            obj.activations = zeros(trials,1);
        end
        
        function obj = reset(obj)
            obj.spikes = 0;
            obj.v(:) = QIAF.rv;
            obj.out(:) = 0;
        end
    end
    
end

