function [retval] = calc_FMRI_corr_data(config, neurons)
    % Extract parameters.
    FMRI_META = config.meta.optimization;
    CN = neurons.CN;
    MDN = neurons.MDN;
    PMC = neurons.PMC;
    
    % Set parameter values of hrf.
    t1 = 1; n = 4; lamda = 2;
    % Define time axis.
    t = 1:LEARNING_TRIALS;
    % Create hrf.
    hrf = ((t-t1).^(n-1)).*exp(-(t-t1)/lamda)/((lamda^n)*factorial(n-1));
    % Get correlation matrix for return value.
    retval = [get_FMRI_corr_data(CN.activations, FMRI_META, hrf); ...
              get_FMRI_corr_data(MDN.activations, FMRI_META, hrf); ...
              get_FMRI_corr_data(PMC.activations, FMRI_META, hrf); ...
              mean(config.accuracy(FMRI_META.SES_1)), mean(config.accuracy(FMRI_META.SES_4)), mean(config.accuracy(FMRI_META.SES_10)), mean(config.accuracy(FMRI_META.SES_20))];
end

% Finds correlation between different neurons and accuracy
function [corr_vec] = get_FMRI_corr_data(activations, FMRI_META, hrf)
    % Compute convolution for each trial
    boldPMC = conv(activations, hrf');
    corr_vec = [mean(boldPMC(FMRI_META.SES_1)), ...
                mean(boldPMC(FMRI_META.SES_4)), ...
                mean(boldPMC(FMRI_META.SES_10)), ...
                mean(boldPMC(FMRI_META.SES_20))];
end