optional_parms = struct('FMRI_META_GROUP_RUN', 0, ...
                        'VIS_INPUT_FROM_PARM', 0, ...
                        'visualinput', zeros(2));
                    
% automaticityModel(getautoparams('FMRI'), optional_parms);
% automaticityModel_mex(getautoparams('FMRI'), optional_parms);
automaticityModel_mex(getautoparams('WALLIS'), optional_parms);