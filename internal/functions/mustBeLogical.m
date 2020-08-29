function mustBeLogical(val)
%MUSTBELOGICAL Errors if val is not logical.
    if not(islogical(val))
        error('Value must be logical');
    end
end

