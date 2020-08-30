function [ TRBoldMap ] = generateTRBoldForMap( vectorsMap )
%GENERATETRBOLDFORMAP Generate TR Bold for each vector in map
%   Uses generatebold to create BOLD for each vector in map, then
%   samples at the start of every TR
    TR_LENGTH = 2000;
    BOLD_STEP = 200;
    TRBoldMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
    for subjNum = keys(vectorsMap)
        subj = vectorsMap(subjNum{1});
        newSubj = cell(20,1);
        for sessionNum = find(~cellfun(@isempty,subj))'
            stimVector = subj{sessionNum};
            BOLD = generatebold(stimVector);
            newSubj{sessionNum} = BOLD(1:BOLD_STEP:length(stimVector)/10);
        end
        TRBoldMap(subjNum{1}) = newSubj;
    end
end

