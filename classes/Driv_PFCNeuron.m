classdef Driv_PFCNeuron < Neuron
    %Driving Signal from PFC (FROST).
    %   Driving signal from the PFC. Input from rule stimulus (arbitrary value - constant). Output to CN Neuron.
    
    properties
        W_OUT
        rule_stim = 0
    end
    
    methods
        function obj = Driv_PFCNeuron(n, TAU, LAMBDA, W_OUT)
            obj@Neuron(n, TAU, LAMBDA);
            obj.W_OUT = W_OUT;
            obj.v = repmat(RSN.rv,n,1);
        end
        
        function obj = reset(obj)
            obj.spikes = 0;
            obj.v(:) = RSN.rv;
            obj.u(:) = 0;
            obj.out(:) = 0;
            obj.restartTime();
        end

        function obj = iterate(obj)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1)=((obj.v(i) + TAU*(RSN.k*(obj.v(i)-RSN.rv)*(obj.v(i)-RSN.vt)-obj.u(i) + RSN.E + obj.rule_stim))/RSN.C);
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

