classdef CNNeuron
    %Caudate Nucleus
    %   Caudate nucleus. For FROST.
    
    properties
        W_OUT = 1
        out
        spikes = 0
        v
        u
        activations
    end
    
    methods
        function obj = CN(n, trials)
            obj.out = zeros(n,1);
            obj.v = repmat(MSN.rv,n,1);
            obj.u = zeros(n,1);
            obj.activations = zeros(trials,1);
        end
        
        function obj = reset(obj)
            obj.spikes = 0;
            obj.v(:) = RSN.rv;
            obj.u(:) = 0;
            obj.out(:) = 0;
        end
    end
    
end

