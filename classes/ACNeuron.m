classdef ACNeuron < Neuron
    %AC
    %   AC. Input from rule stimulus (arbitrary value - constant). Input from PFC neuron. Output to PFC neuron. For FROST.
    
    properties
        W_OUT = 1
        rule_stim = 0.1
    end
    
    methods
        function obj = ACNeuron(n, TAU, LAMBDA)
            obj@Neuron(n, TAU, LAMBDA);
            obj.v = repmat(RSN.rv,n,1);
        end
        
        function obj = reset(obj)
            obj.spikes = 0;
            obj.v(:) = RSN.rv;
            obj.u(:) = 0;
            obj.out(:) = 0;
        end

        function obj = iterate(obj, PFC)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1)=((obj.v(i) + TAU*(RSN.k*(obj.v(i)-RSN.rv)*(obj.v(i)-RSN.vt)-obj.u(i)+ 10 + PFC.W_OUT_AC*PFC.out(i) + obj.rule_stim))/RSN.C);
            obj.u(i+1)=obj.u(i)+TAU*RSN.a*(RSN.b*(obj.v(i)-RSN.rv)-obj.u(i));
            if obj.v(i+1)>=RSN.vpeak
                obj.v(i)= RSN.vpeak;
                obj.v(i+1)= RSN.c;
                obj.u(i+1)= obj.u(i+1)+ RSN.d;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end
    end
    
end

