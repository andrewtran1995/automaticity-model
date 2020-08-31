classdef Category < uint8
    %AREA Enum that represent the categories a neuron or visual stimulus
    %area may refer to.
    enumeration
        NONE(0), A(1), B(2)
    end
    
    methods
        function ret = ID(obj)
            ret = uint8(obj);
        end
    end
end

