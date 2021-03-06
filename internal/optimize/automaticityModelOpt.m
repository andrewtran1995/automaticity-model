function [ opt_val ] = automaticityModelOpt( arg_vector )
%AUTOMATICITYMODELOPT Cost function for global optimization
%   Wrapper for automaticityModel with a function signature that adheres to
%   global optimization function requirements
    addpath(genpath('.'));
    GROUP_SIZE = 12;
    CONFIG = ModelConfigButtonSwitch();
    NAN_CEILING = 10;
    CN  = zeros(GROUP_SIZE, 4);
    MDN = zeros(GROUP_SIZE, 4);
    PMC = zeros(GROUP_SIZE, 4);
    acc = zeros(GROUP_SIZE, 4);
    
    parfor i = 1:GROUP_SIZE
        arg_struct = argvectortostruct(arg_vector, CONFIG);

        % Populate second function argument
        optional_parms = struct('VIS_INPUT_FROM_PARM', 0, ...
                                'visualinput', zeros(2));

        % Call automaticity model function
        [config, neurons] = automaticityModel_mex(arg_struct, optional_parms);
        retval = calc_FMRI_corr_data(config, neurons)
        
        
        % Assign values for correlation
        CN(i,:)  = retval(1,:);
        MDN(i,:) = retval(2,:);
        PMC(i,:) = retval(3,:);
        acc(i,:) = retval(4,:);
    end
    
    % Calculate SSE
    actual_corr = [corrneuron(CN, acc); ...
                   corrneuron(MDN, acc); ...
                   corrneuron(PMC, acc)];
    target = load('fmri/targetFMRICorrelations.mat');
    target_diff = actual_corr - [target.CN; target.MDN; target.PMC];
    target_diff(isnan(target_diff)) = NAN_CEILING;
    opt_val = sum(sum(target_diff.^2));
end