classdef PFCNeuron
    %Prefontal Cortex Neuron
    %   This class represents a neuron from the prefontal cortex (PFC).
    
    properties
        W_OUT = 9
        W_OUT_MDN
        W_OUT_AC = 1
        out
        spikes = 0
        v
        u
        pos_volt
        v_stim = 0
    end
    
    properties (Constant)
        V_SCALE = 1 % scaling factor for visual input into PFC neurons
        W_LI = 2 % lateral inhibition between PFC A / PFC B
    end
    
    methods
        function obj = PFC(W_OUT_MDN, n)
            obj.W_OUT_MDN = W_OUT_MDN;
            obj.out = zeros(n,1);
            obj.v = repmat(RSN.rv,n,1);
            obj.u = zeros(n,1);
            obj.pos_volt = zeros(n,1);
        end
        
        function obj = reset(obj)
            obj.spikes = 0;
            obj.pos_volt(:) = 0;
            obj.v(:) = RSN.rv;
            obj.u(:) = 0;
            obj.out(:) = 0;
        end
        
        function obj = neuron_iterate_FROST(obj, i, TAU, LAMBDA_PRECALC, NOISE_PFC, PFC_OTHER, PMC, MDN, AC)
            obj.v(i+1) = (obj.v(i) + TAU*(RSN.k*(obj.v(i)-RSN.rv)*(obj.v(i)-RSN.vt)-obj.u(i) + RSN.E + MDN.W_OUT*MDN.out(i) + AC.W_OUT*AC.out(i) + obj.v_stim + (PMC.W_OUT*PMC.out(i)) - PFCNeuron.W_LI*PFC_OTHER.out(i))/RSN.C) + normrnd(0,NOISE_PFC);
            obj.u(i+1) = obj.u(i) + TAU*RSN.a*(RSN.b*(obj.v(i)-RSN.rv)-obj.u(i));
            if obj.v(i+1) >= RSN.vpeak
                obj.v(i) = RSN.vpeak;
                obj.v(i+1) = RSN.c;
                obj.u(i+1) = obj.u(i+1) + RSN.d;
                obj.out(i:n) = obj.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
            end
        end
        
        function obj = neuron_iterate(obj, i, TAU, LAMBDA_PRECALC, NOISE_PFC, PFC_OTHER, PMC)
            obj.v(i+1) = (obj.v(i) + TAU*(RSN.k*(obj.v(i)-RSN.rv)*(obj.v(i)-RSN.vt)-obj.u(i) + RSN.E + obj.v_stim + (PMC.W_OUT*PMC.out(i)) - PFCNeuron.W_LI*PFC_OTHER.out(i))/RSN.C) + normrnd(0,NOISE_PFC);
            obj.u(i+1) = obj.u(i) + TAU*RSN.a*(RSN.b*(obj.v(i)-RSN.rv)-obj.u(i));
            if obj.v(i+1) >= RSN.vpeak
                obj.v(i) = RSN.vpeak;
                obj.v(i+1) = RSN.c;
                obj.u(i+1) = obj.u(i+1) + RSN.d;
                obj.out(i:n) = obj.out(i:n) + LAMBDA_PRECALC(1:n-i+1);
            end
        end
    end
    
end

