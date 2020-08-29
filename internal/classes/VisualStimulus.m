classdef VisualStimulus
    %VISUALSTIMULUS A set of visuaul stimulus rules.
    
    properties
        coord = struct('x', 0, 'y', 0, 'group', 0)
        RULES
    end
    
    properties (Constant)
        STIM = 50
    end
    
    methods
        function obj = VisualStimulus()
            AREA = obj.AREA();
            obj.RULES = [ ...
                %% One-dimensional rules
                % The left half is A. The right half is B.
                struct('A_X', AREA.LOWER_HALF, 'A_Y', AREA.ALL,        'B_X', AREA.UPPER_HALF, 'B_Y', AREA.ALL,        'INVERSE', 2); ...
                % The right half is A. The left half is A.
                struct('A_X', AREA.UPPER_HALF, 'A_Y', AREA.ALL,        'B_X', AREA.LOWER_HALF, 'B_Y', AREA.ALL,        'INVERSE', 1); ...
                % The lower half is A. The upper half is B.
                struct('A_X', AREA.ALL,        'A_Y', AREA.LOWER_HALF, 'B_X', AREA.ALL,        'B_Y', AREA.UPPER_HALF, 'INVERSE', 4); ...
                % The upper half is A. The lower half is B.
                struct('A_X', AREA.ALL,        'A_Y', AREA.UPPER_HALF, 'B_X', AREA.ALL,        'B_Y', AREA.LOWER_HALF, 'INVERSE', 3); ...

                %% Disjunctive rules
                % The far left and far right are A. The inner section is B.
                struct('A_X', AREA.OUTER,      'A_Y', AREA.ALL,        'B_X', AREA.INNER,      'B_Y', AREA.ALL,        'INVERSE', 6); ...
                % The inner section is A. The far left and far right are B.
                struct('A_X', AREA.INNER,      'A_Y', AREA.ALL,        'B_X', AREA.OUTER,      'B_Y', AREA.ALL,        'INVERSE', 5); ...
                % The far lower and upper are A. The inner section is B.
                struct('A_X', AREA.ALL,        'A_Y', AREA.OUTER,      'B_X', AREA.ALL,        'B_Y', AREA.INNER,      'INVERSE', 8); ...
                % The inner section is A. The far lower and upper are B.
                struct('A_X', AREA.ALL,        'A_Y', AREA.INNER,      'B_X', AREA.ALL,        'B_Y', AREA.OUTER,      'INVERSE', 7) ...
            ];
        end
    end
    methods (Static)
        %{
        The AREA variable is a struct that contains convenient short-hand
        descriptions based off of:
        - STIMULUS_GRID_SIZE: The length of the grid where stimulus is given.
        - BORDER_SIZE: The width of the border around the stimulus grid.
        - GRID_SIZE: The overall length of the grid, expected to be the length
        of the stimulus grid + the widths of the border.

        The different areas that can be described by AREA are:
        - LOWER_HALF: The lower half of a dimension of the whole grid.
        - UPPER_HALF: The upper half of a dimension of the whole grid.
        - OUTER: The first 1/4 and last 1/4 of a dimension of the stimulus grid.
        - INNER: The area between the 1/4 and 3/4 marks of a dimension of the stimulus grid.
        - ALL: The full length of a dimension of the whole grid.
        %}
        function AREA = AREA()
            STIMULUS_GRID_SIZE = ModelConfig.STIMULUS_GRID_SIZE;
            GRID_SIZE = ModelConfig.GRID_SIZE;
            BORDER_SIZE = ModelConfig.BORDER_SIZE;
            
            AREA = struct('LOWER_HALF', 1:GRID_SIZE/2, 'UPPER_HALF', GRID_SIZE/2+1:GRID_SIZE, ...
    			  'OUTER', [1:STIMULUS_GRID_SIZE/4+BORDER_SIZE, STIMULUS_GRID_SIZE*3/4+BORDER_SIZE+1:GRID_SIZE], ...
    			  'INNER', STIMULUS_GRID_SIZE/4+BORDER_SIZE+1:STIMULUS_GRID_SIZE*3/4+BORDER_SIZE, ...
    			  'ALL', 1:GRID_SIZE);
        end
    end
end

