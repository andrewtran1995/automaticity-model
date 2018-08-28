classdef GPNeuron < Neuron
    %Global Pallidus
    %   Global pallidus. Input from CN Neuron. Output to MDN neuron. For FROST.
    
    properties
        W_OUT = 1
        activations
    end
    
    methods
        function obj = GPNeuron(n, TAU, LAMBDA, trials)
            obj@Neuron(n, TAU, LAMBDA);
            obj.v = repmat(QIAF.rv,n,1);
            obj.activations = zeros(trials,1);
        end
        
        function obj = reset(obj)
            obj.spikes = 0;
            obj.v(:) = QIAF.rv;
            obj.out(:) = 0;
            obj.restartTime();
        end

        function obj = iterate(obj, CN)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            
            dGP = (-1)*CN.W_OUT*CN.out(i) + QIAF.beta + QIAF.gamma*(obj.v(i)- QIAF.rv)*(obj.v(i)-QIAF.vt);               
            obj.v(i+1) = obj.v(i) + dGP;
            if (obj.v(i+1) >= QIAF.vpeak)
                obj.v(i) = QIAF.vpeak;
                obj.v(i+1) = QIAF.vreset;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end
    end
    
end

