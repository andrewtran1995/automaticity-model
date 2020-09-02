function [params] = getconstants()
%GETCONSTANTS Get configuration-agnostic constants.
    % Initialize configuration-agnostic parameters used in optimization
    agn_names = {'HEB_CONSTS';'NMDA';'AMPA';'W_MAX';'NOISE_PFC';'NOISE_PMC';'NOISE_MC';'PMC_A_W_OUT';'PMC_B_W_OUT'};
    agn_vals  = {       1e-8;  450;    300;     10;        2.4;          2;         5;            1;            1};
    % FROST Params
    frost_names = {'PFC_A_W_OUT_MDN';'PFC_B_W_OUT_MDN';'DRIV_PFC_W_OUT';'MDN_A_W_OUT';'MDN_B_W_OUT'};
    frost_vals  = {                5;                5;               5;            5;            5};
    % COVIS Params
    covis_names = {'COVIS_DELTA_C';'COVIS_DELTA_E';'COVIS_PERSEV';'COVIS_LAMBDA'};
    covis_vals  = {         9.6545;         6.1809;        8.5273;        3.8482};
    params = cell2struct(vertcat( agn_vals,  frost_vals, covis_vals), ...  % Parameter values
                         vertcat(agn_names, frost_names, covis_names), ... % Parameter names
                         1);
end

