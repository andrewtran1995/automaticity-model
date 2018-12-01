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
    end
    
    methods
        function obj = PFCNeuron(n, TAU, LAMBDA, W_OUT_MDN)
            obj@RSN(n, TAU, LAMBDA);
            obj.W_OUT_MDN = W_OUT_MDN;
            obj.v = repmat(obj.rv,n,1);
        end
        
        function obj = reset(obj)
            obj = reset@RSN(obj);
        end
        
        function obj = iterate_FROST(obj, NOISE_PFC, PFC_OTHER, PMC, MDN, AC)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1) = (obj.v(i) + TAU*(obj.k*(obj.v(i)-obj.rv)*(obj.v(i)-obj.vt)-obj.u(i) + obj.E + MDN.W_OUT*MDN.out(i) + AC.W_OUT*AC.out(i) + obj.v_stim + PMC.W_OUT*PMC.out(i) - obj.W_LI*PFC_OTHER.out(i))/obj.C) + normrnd(0,NOISE_PFC);
            obj.u(i+1) = obj.u(i) + TAU*obj.a*(obj.b*(obj.v(i)-obj.rv)-obj.u(i));
            if obj.v(i+1) >= obj.vpeak
                obj.v(i) = obj.vpeak;
                obj.v(i+1) = obj.c;
                obj.u(i+1) = obj.u(i+1) + obj.d;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end
        
        function obj = iterate(obj, NOISE_PFC, PFC_OTHER, PMC)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1) = (obj.v(i) + TAU*(obj.k*(obj.v(i)-obj.rv)*(obj.v(i)-obj.vt)-obj.u(i) + obj.E + obj.v_stim + PMC.W_OUT*PMC.out(i) - obj.W_LI*PFC_OTHER.out(i))/obj.C) + normrnd(0,NOISE_PFC);
            obj.u(i+1) = obj.u(i) + TAU*obj.a*(obj.b*(obj.v(i)-obj.rv)-obj.u(i));
            if obj.v(i+1) >= obj.vpeak
                obj.v(i) = obj.vpeak;
                obj.v(i+1) = obj.c;
                obj.u(i+1) = obj.u(i+1) + obj.d;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end
    end
    
end

