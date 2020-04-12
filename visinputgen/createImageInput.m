function [visualInput] = createImageInput()
%CREATEIMAGEINPUT Create visual input for the IMAGE configuration.
    TRIALS = 15000;

    %% Load the image vectors.
    loadedInput = load('datasets/imageVectors.mat');
    imageVectors = loadedInput.imageVectors;
    
    %% Create correlation matrix.
    numImages = size(imageVectors, 1);
    corrMatrix = ones(numImages, numImages);
    for i=1:numImages
        for j=1:numImages
            corr = corrcoef(imageVectors(i,:), imageVectors(j,:));
            corrMatrix(i,j) = abs(corr(1,2));
        end
    end
    
    %% Create x_coordinate values.
    xValues = zeros(TRIALS, 1);
    for i=1:length(xValues)
        if rand(1) <= 0.5
            xValues(i) = 0.8;
        else
            xValues(i) = corrMatrix(randi(numImages), randi(numImages));
        end
    end
    
    %% Return the visual input struct.
    visualInput = struct(...
        'x',      xValues, ...
        'y',      zeros(TRIALS, 1)+50, ...
        'groups', zeros(TRIALS, 1) ...
    );
end