classdef MDNNeuron < RSN
    %MDN
    %   MDN. Input from GP neuron. Input from PFC neuron. Output to PFC neuron. For FROST.
    
    properties (Constant)
        W_OUT = 1
    end
    
    methods
        function obj = MDNNeuron()
            obj@RSN();
            obj.v = repmat(obj.rv,obj.n,1);
        end

        function obj = iterate(obj, PFC, GP)
            % Create local variables for readability
            i = obj.i;
            TAU = obj.TAU;
            
            obj.v(i+1)=(obj.v(i) + TAU*(obj.k*(obj.v(i)-obj.rv)*(obj.v(i)-obj.vt)-obj.u(i) + obj.E + PFC.W_OUT_MDN*PFC.out(i) - GP.W_OUT*GP.out(i)))/obj.C;
            obj.u(i+1) = obj.u(i) ...
                       + TAU * obj.a * (obj.b * (obj.v(i) - obj.rv) - obj.u(i));
            if obj.v(i+1) >= obj.vpeak
                obj.v(i:obj.i+1) = [obj.vpeak, obj.c];
                obj.u(i+1) = obj.u(i+1) + obj.d;
                obj.out(i:end) = obj.out(i:end) + obj.LAMBDA_PRECALC(1:obj.n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end
    end
    
end

