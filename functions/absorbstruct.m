function [ retstruct ] = absorbstruct( parent, child )
%ABSORBSTRUCT Merge fields from child into parent into new struct
%   Assign fields in the parent struct from corresponding fields in the
%   child struct. Throws error if unexpected field in child struct appears.
    retstruct = parent;
    child_names = fieldnames(child);
    for i = 1:numel(child_names)
        if isfield(parent, child_names{i})
            retstruct.(child_names{i}) = child.(child_names{i});
        else
            error('Field in arg_struct not found in PARAMS: %s\n Verify that arg_struct is valid.', child_names{i});
        end
    end
end