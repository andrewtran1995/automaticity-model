loaded = load('fmri/particleswarm_target_17_12_14.mat');
arg_vector = loaded.x;
arg_struct = argvectortostruct(arg_vector, 'FMRI');
optional_parms = struct('FMRI_META_GROUP_RUN', 0, ...
                        'VIS_INPUT_FROM_PARM', 0, ...
                        'visualinput', zeros(2));
                    
automaticityModel(arg_struct, optional_parms);

% automaticityModel_mex(arg_struct, optional_parms);