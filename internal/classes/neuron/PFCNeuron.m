classdef PFCNeuron < RSN
    %Prefontal Cortex Neuron
    %   This class represents a neuron from the prefontal cortex (PFC).

    properties
        W_OUT = 9
        W_OUT_MDN
        W_OUT_AC = 1
        v_stim = 0
    end
    
    properties (Constant)
        V_SCALE = 1 % scaling factor for visual input into PFC neurons
        W_LI = 2 % lateral inhibition between PFC A / PFC B
        NOISE = getconstants().NOISE_PFC
    end
    
    methods
        function obj = PFCNeuron(W_OUT_MDN)
            obj@RSN();
            obj.W_OUT_MDN = W_OUT_MDN;
            obj.v = repmat(obj.rv,obj.n,1);
        end
        
        function obj = iterate_FROST(obj, PFC_OTHER, MDN, AC)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1) = obj.v(i) ...
                       + TAU*(obj.k*(obj.v(i)-obj.rv)*(obj.v(i)-obj.vt)-obj.u(i) + obj.E + MDN.W_OUT*MDN.out(i) + AC.W_OUT*AC.out(i) + obj.v_stim - obj.W_LI*PFC_OTHER.out(i))/obj.C ...
                       + normrnd(0, obj.NOISE);
            obj.u(i+1) = obj.u(i) + TAU*obj.a*(obj.b*(obj.v(i)-obj.rv)-obj.u(i));
            if obj.v(i+1) >= obj.vpeak
                obj.v(i:i+1) = [obj.vpeak, obj.c];
                obj.u(i+1) = obj.u(i+1) + obj.d;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end
        
        function obj = iterate(obj, PFC_OTHER)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1) = obj.v(i) ...
                       + TAU*(obj.k*(obj.v(i)-obj.rv)*(obj.v(i)-obj.vt)-obj.u(i) + obj.E + obj.v_stim - obj.W_LI*PFC_OTHER.out(i))/obj.C ...
                       + normrnd(0, obj.NOISE);
            obj.u(i+1) = obj.u(i) + TAU*obj.a*(obj.b*(obj.v(i)-obj.rv)-obj.u(i));
            if obj.v(i+1) >= obj.vpeak
                obj.v(i:i+1) = [obj.vpeak, obj.c];
                obj.u(i+1) = obj.u(i+1) + obj.d;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end
    end
    
end

