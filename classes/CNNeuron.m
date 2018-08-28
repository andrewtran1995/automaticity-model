classdef CNNeuron < Neuron
    %Caudate Nucleus
    %   Caudate nucleus. Input from Driv_PFC neuron. Output to GP neuron. For FROST.
    
    properties
        W_OUT = 1
        activations
    end
    
    methods
        function obj = CNNeuron(n, TAU, LAMBDA, trials)
            obj@Neuron(n, TAU, LAMBDA);
            obj.v = repmat(MSN.rv,n,1);
            obj.activations = zeros(trials,1);
        end
        
        function obj = reset(obj)
            obj.spikes = 0;
            obj.v(:) = RSN.rv;
            obj.u(:) = 0;
            obj.out(:) = 0;
            obj.restartTime();
        end

        function obj = iterate(obj, Driv_PFC)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1)=((obj.v(i) + TAU*(MSN.k*(obj.v(i)-MSN.rv)*(obj.v(i)-MSN.vt)- obj.u(i) + Driv_PFC.W_OUT*Driv_PFC.out(i) + MSN.E ))/MSN.C);
            obj.u(i+1)= obj.u(i)+TAU*MSN.a*(MSN.b*(obj.v(i)-MSN.rv)-obj.u(i));
            if obj.v(i+1)>=MSN.vpeak
                obj.v(i)= MSN.vpeak;
                obj.v(i+1)= MSN.c;
                obj.u(i+1)= obj.u(i+1)+ MSN.d;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end
    end
    
end

