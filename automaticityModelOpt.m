function [ opt_val ] = automaticityModelOpt( arg_vector )
%AUTOMATICITYMODELOPT Cost function for global optimization
%   Wrapper for automaticityModel with a function signature that adheres to
%   global optimization function requirements
    GROUP_SIZE = 12;
    CONFIG = 'FMRI';
    NAN_CEILING = 10;
    PFC = zeros(GROUP_SIZE, 4);
    CN  = zeros(GROUP_SIZE, 4);
    GP  = zeros(GROUP_SIZE, 4);
    MDN = zeros(GROUP_SIZE, 4);
    PMC = zeros(GROUP_SIZE, 4);
    acc = zeros(GROUP_SIZE, 4);
    
    parfor i = 1:GROUP_SIZE
        arg_struct = argvectortostruct(arg_vector, CONFIG);

        % Populate second function argument
        optional_parms = struct('FMRI_META_GROUP_RUN', 1, ...
                                'VIS_INPUT_FROM_PARM', 0, ...
                                'visualinput', zeros(2));

        % Call automaticity model function
        [~, retval] = automaticityModelFast_mex(arg_struct, optional_parms);
        
        % Assign values for correlation
        PFC(i,:) = retval(1,:);
        CN(i,:)  = retval(2,:);
        GP(i,:)  = retval(3,:);
        MDN(i,:) = retval(4,:);
        PMC(i,:) = retval(5,:);
        acc(i,:) = retval(6,:);
    end
    
    % Calculate SSE
    actual_corr = [corrneuron(PFC, acc); ...
                   corrneuron(CN, acc); ...
                   corrneuron(GP, acc); ...
                   corrneuron(MDN, acc); ...
                   corrneuron(PMC, acc)];
    target = load('fmri/targetFMRICorrelations.mat');
    target_diff = actual_corr - [target.PFC; target.CN; target.GP; target.MDN; target.PMC];
    target_diff(isnan(target_diff)) = NAN_CEILING;
    opt_val = sum(sum(target_diff.^2));
end