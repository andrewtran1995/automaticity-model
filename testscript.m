addpath(genpath('.'));
optional_parms = struct('FMRI_META_GROUP_RUN', 0, ...
                        'VIS_INPUT_FROM_PARM', 0, ...
                        'visualinput', zeros(2));

automaticityModel_mex(getmodelparams(ModelConfigButtonSwitch()), optional_parms);