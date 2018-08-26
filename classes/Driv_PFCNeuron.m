classdef Driv_PFCNeuron
    %Driving Signal from PFC (FROST).
    %   Driving signal from the PFC.
    
    properties
        W_OUT
        out
        spikes = 0
        v
        u
        rule_stim = 0
    end
    
    methods
        function obj = Driv_PFC(W_OUT, n)
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

