classdef PFCStroopNeuron < PFCNeuron
    %PFCSTROOPNEURON Summary of this class goes here
    %   Detailed explanation goes here 
    %{
    - Stimulus is independent of categorization units, based off of own
    constant.
    - MDN and stroop units project onto each other to form a reverberating
    feedback loop.
    - Individual stroop units connect to only their respective MDN units.
    However, both stroop units connect to both PFC units.
    
    Other notes
    Only enable MC learning if button switch enabled
    %}
    
    properties (Constant)
        STROOP_STIM = 250
        NOISE = getconstants().NOISE_PFC
    end
    
    methods
        function obj = PFCStroopNeuron(inputArg1,inputArg2)
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function obj = iterate(obj, MDN)
            % Create local variables for readability.
            i = obj.i;
            n = obj.n;
            TAU = obj.TAU;

            obj.v(i+1) = obj.v(i) ...
                + ( ...
                    TAU * obj.k * (obj.v(i) - obj.rv) * (obj.v(i) - obj.vt) ...
                    - obj.u(i) ...
                    + obj.E ...
                    + obj.STROOP_STIM ...
                    + 20 * MDN.out(i) ...
                )/obj.C ...
                + normrnd(0, obj.NOISE);
            obj.u(i+1) = obj.u(i) + TAU * obj.a * (obj.b * (obj.v(i) - obj.rv) - obj.u(i));
            if obj.v(i+1) >= obj.vpeak
                obj.v(i:i+1) = [obj.vpeak, obj.c];
                obj.u(i+1) = obj.u(i+1) + obj.d;
                obj.out(i:n) = obj.out(i:n) + obj.LAMBDA_PRECALC(1:n-1:1);
            end

            % Increment time.
            obj.i = obj.i + 1;
        end
    end
end
