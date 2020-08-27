addpath('config');
optional_parms = struct('FMRI_META_GROUP_RUN', 0, ...
                        'VIS_INPUT_FROM_PARM', 0, ...
                        'visualinput', zeros(2));

% automaticityModel(getmodelparams('FMRI'), optional_parms);
% automaticityModel(getmodelparams(AutomaticityConfiguration.FMRI, optional_parms));
automaticityModel_mex(getmodelparams(AutomaticityConfiguration.FMRI), optional_parms);