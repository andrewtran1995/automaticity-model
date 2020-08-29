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
    
    properties
        Property1
    end
    
    methods
        function obj = PFCStroopNeuron(inputArg1,inputArg2)
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

