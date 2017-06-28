LOADED_PMC = load('pmcAlphaVectors_17_6_2.mat');
LOADED_STIM = load('vectors/stimVectors.mat');

% modelAlphaVectorMap is ridiculously huge (120+ MB), so don't bother
% saving it to a .mat file, but use it to calculate TRBoldMap instead
modelAlphaVectorMap = replaceNeuronOutput(LOADED_PMC.PMC_ALPHA, ...
                                          LOADED_PMC.LATENCY, ...
                                          LOADED_STIM.stimVectors);

% TRBoldMap is much more suitable to save to a .mat file
TRBoldMap = generateTRBoldForMap(modelAlphaVectorMap);