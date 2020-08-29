classdef GPNeuron < QIAF
    %Global Pallidus
    %   Global pallidus. Input from CN Neuron. Output to MDN neuron. For FROST.
    
    properties
        W_OUT = 1
        activations
    end
    
    methods
        function obj = GPNeuron(trials)
            obj@QIAF();
            obj.v = repmat(obj.rv,obj.n,1);
            obj.activations = zeros(trials,1);
        end

        function iterate(obj, CN)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            
            dGP = (-1)*CN.W_OUT*CN.out(i) + obj.beta + obj.gamma*(obj.v(i)- obj.rv)*(obj.v(i)-obj.vt);               
            obj.v(i+1) = obj.v(i) + dGP;
            if (obj.v(i+1) >= obj.vpeak)
                obj.v(i) = obj.vpeak;
                obj.v(i+1) = obj.vreset;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-i+1);
            end
            
            % Increment time
            obj.i = obj.i + 1;
        end
    end
    
end

