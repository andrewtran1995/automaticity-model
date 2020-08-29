function [bounded_array] = bound_array(num, lower_bound, upper_bound)
    % Bound the input to the given range, inclusive
    bounded_array = min(max(num, lower_bound), upper_bound);
end