classdef ACNeuron < RSN
    %AC
    %   AC. Input from rule stimulus (arbitrary value - constant). Input from PFC neuron. Output to PFC neuron. For FROST.
    
    properties
        W_OUT = 1
        rule_stim = 0.1
    end
    
    methods
        function obj = ACNeuron(n, TAU, LAMBDA)
            obj@RSN(n, TAU, LAMBDA);
            obj.v = repmat(obj.rv,n,1);
        end
        
        function obj = reset(obj)
            obj = reset@RSN(obj);
        end

        function obj = iterate(obj, PFC)
            % Create local variables for readability
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;
            
            obj.v(i+1)=((obj.v(i) + TAU*(obj.k*(obj.v(i)-obj.rv)*(obj.v(i)-obj.vt)-obj.u(i)+ 10 + PFC.W_OUT_AC*PFC.out(i) + obj.rule_stim))/obj.C);
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

