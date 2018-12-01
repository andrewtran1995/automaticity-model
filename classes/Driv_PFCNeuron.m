classdef Driv_PFCNeuron < RSN
    %Driving Signal from PFC (FROST).
    %   Driving signal from the PFC. Input from rule stimulus (arbitrary value - constant). Output to CN Neuron.
    
    properties
        W_OUT
        rule_stim = 0
    end
    
    methods
        function obj = Driv_PFCNeuron(n, TAU, LAMBDA, W_OUT)
            obj@RSN(n, TAU, LAMBDA);
            obj.W_OUT = W_OUT;
            obj.v = repmat(obj.rv,n,1);
        end
        
        function obj = reset(obj)
            obj = reset@RSN(obj);
        end

        function obj = iterate(obj)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1)=((obj.v(i) + TAU*(obj.k*(obj.v(i)-obj.rv)*(obj.v(i)-obj.vt)-obj.u(i) + obj.E + obj.rule_stim))/obj.C);
            obj.u(i+1)=obj.u(i)+TAU*obj.a*(obj.b*(obj.v(i)-obj.rv)-obj.u(i));
            if obj.v(i+1)>=obj.vpeak
                obj.v(i)= obj.vpeak;
                obj.v(i+1)= obj.c;
                obj.u(i+1)= obj.u(i+1)+ obj.d;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end
    end
    
end

