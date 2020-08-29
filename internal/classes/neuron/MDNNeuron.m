classdef MDNNeuron < RSN
    %MDN
    %   MDN. Input from GP neuron. Input from PFC neuron. Output to PFC neuron. For FROST.
    
    properties
        W_OUT
    end
    
    methods
        function obj = MDNNeuron(W_OUT)
            obj@RSN();
            obj.W_OUT = W_OUT;
            obj.v = repmat(obj.rv,obj.n,1);
        end

        function obj = iterate(obj, PFC, GP)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1)=((obj.v(i) + TAU*(obj.k*(obj.v(i)-obj.rv)*(obj.v(i)-obj.vt)-obj.u(i)+ 10 + PFC.W_OUT_MDN*PFC.out(i) - GP.W_OUT*GP.out(i)))/obj.C);
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

