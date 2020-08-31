function [coords] = createDualTaskCoords()
    % Create coordinates mapped from the methods used to generate gabors
    % from the referenced Zeithamova paper, translating spatial frequency
    % and orientation to X and Y coordinates, respectively.
    NUM_GROUP = 40;
    NUM_COORDS = NUM_GROUP*50;
    a_coords = createBaseCoords(NUM_GROUP, 280, sqrt(75), 125, sqrt(9000), Category.A);
    b_coords = createBaseCoords(NUM_GROUP, 320, sqrt(75), 125, sqrt(9000), Category.B);
    
    coords(NUM_COORDS,1) = Coord();
    % Alternate between category A stimulus and category B stimulus every
    % NUM_GROUP coordinates.
    for i=1:NUM_GROUP*2:NUM_COORDS
        coords(i:i+NUM_GROUP-1) = a_coords;
        coords(i+NUM_GROUP:i+2*NUM_GROUP-1) = b_coords;
    end
end

function [coords] = createBaseCoords(num_coords, sf_mean, sf_std, so_mean, so_std, category)
    % Create base coordinates given:
    % - Spatial frequency mean
    % - Spatial frequency standard deviation
    % - Spatial orientation mean
    % - Spatial orientation standard deviation
    coords(num_coords,1) = Coord();
    for i=1:num_coords
        coords(i) = Coord( ...
            (0.25 + normrnd(sf_mean, sqrt(sf_std))/50) * 14 - 38, ...
            normrnd(so_mean, sqrt(so_std)) * pi/500 * 57.3 * 0.6 + 20, ...
            category ...
        );
    end
end