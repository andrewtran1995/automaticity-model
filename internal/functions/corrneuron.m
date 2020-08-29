function [ corr_vec ] = corrneuron( neuron, accuracy )
%CORRNEURON Compute correlation vector from neuron and accuracy
%   Given a 2D matrix of neural convolution data and a 2D accuracy
%   matrix, compute corr_vec (correlation per session). Matrix is assumed
%   to have 4 columns (for 4 distinct FMRI sessions).
    corr_vec = [corr(neuron(:,1), accuracy(:,1)), ...
                corr(neuron(:,2), accuracy(:,2)), ...
                corr(neuron(:,3), accuracy(:,3)), ...
                corr(neuron(:,4), accuracy(:,4))];
end

