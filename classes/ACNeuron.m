classdef ACNeuron
    %AC
    %   AC. For FROST.
    
    properties
        W_OUT = 1
        out
        spikes = 0
        v
        u
        rule_stim = 0.1
    end
    
    methods
        function obj = AC(n)
            obj.out = zeros(n,1);
            obj.v = repmat(RSN.rv,n,1);
            obj.u = zeros(n,1);
        end
        
        function obj = reset(obj)
            obj.spikes = 0;
            obj.v(:) = RSN.rv;
            obj.u(:) = 0;
            obj.out(:) = 0;
        end
    end
    
end

