classdef Driv_PFCNeuron < RSN
    %Driving Signal from PFC (FROST).
    %   Driving signal from the PFC. Input from rule stimulus (arbitrary value - constant). Output to CN Neuron.
    
    properties
        W_OUT = 1
        rule_stim = 0
    end
    
    methods
        function obj = Driv_PFCNeuron()
            obj@RSN();
            obj.v = repmat(obj.rv,obj.n,1);
        end

        function obj = iterate(obj)
            % Create local variables for readability
            i = obj.i;
            TAU = obj.TAU;
            
            obj.v(i+1)=((obj.v(i) + TAU*(obj.k*(obj.v(i)-obj.rv)*(obj.v(i)-obj.vt)-obj.u(i) + obj.E + obj.rule_stim))/obj.C);
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

