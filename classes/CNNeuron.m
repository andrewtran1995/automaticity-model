classdef CNNeuron < MSN
    %Caudate Nucleus
    %   Caudate nucleus. Input from Driv_PFC neuron. Output to GP neuron. For FROST.
    
    properties
        W_OUT = 1
        activations
    end
    
    methods
        function obj = CNNeuron(n, TAU, LAMBDA, trials)
            obj@MSN(n, TAU, LAMBDA);
            obj.v = repmat(obj.rv,n,1);
            obj.activations = zeros(trials,1);
        end
        
        function obj = reset(obj)
            obj = reset@MSN(obj);
        end

        function obj = iterate(obj, Driv_PFC)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1)=((obj.v(i) + TAU*(obj.k*(obj.v(i)-obj.rv)*(obj.v(i)-obj.vt)- obj.u(i) + Driv_PFC.W_OUT*Driv_PFC.out(i) + obj.E ))/obj.C);
            obj.u(i+1)= obj.u(i)+TAU*obj.a*(obj.b*(obj.v(i)-obj.rv)-obj.u(i));
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

