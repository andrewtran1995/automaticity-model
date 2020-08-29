function [merged_struct] = mergestructs(struct_a, struct_b)
    if isempty(struct_a)
        merged_struct=struct_b;
        return
    end
    if isempty(struct_b)
        merged_struct=struct_a;
        return
    end
    %%insert struct a
    merged_struct=struct_a;
    %%insert struct b
    f = fieldnames(struct_b);
    for i = 1:length(f)
        merged_struct(1).(f{i}) = struct_b(1).(f{i});
    end
end

