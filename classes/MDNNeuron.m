classdef MDNNeuron
    %MDN
    %   MDN. For FROST.
    
    properties
        W_OUT
        out
        spikes = 0
        v
        u
    end
    
    methods
        function obj = MDN(W_OUT, n)
            obj.W_OUT = W_OUT;
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

