function [ vectorized, non_vec, is_equal ] = testVectorizeRBV()
%TESTVECTORIZERBV Summary of this function goes here
%   Detailed explanation goes here
    r_y = 40;
    r_x = 50;
    radius = 1;
    Visual.stim = 50;
    
    [Y, X] = meshgrid(1:100, 1:100);
    vectorized = exp( -(sqrt((r_y-Y)^2 + (r_x-X)^2))/radius ) * Visual.stim;
    
    non_vec = zeros(100);
    for y=1:100
        for x=1:100
            non_vec(y,x) = exp( -(sqrt((r_y-y)^2 + (r_x-x)^2))/radius ) * Visual.stim;
        end
    end
    
    is_equal = isequal(vectorized, non_vec);
    
    return
end

