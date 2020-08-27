classdef Area < uint32
    %AREA Enum that maps areas to the IDs of the neurons associated to
    %those areas.
    enumeration
        A(1), B(2)
    end
    
    methods
        function ret = ID(obj)
            ret = uint32(obj);
        end
    end
end

