function saveSubjectMapToFiles( vectorsMap, directory )
%SAVESUBJECTMAPTOFILES Save map to subject{}_session{} format
% vectorsMap takes in a map object (i.e., from crosshairsVectorsMap.mat,
% responseVectorsMap.mat, or stimVectors.mat)
% directory specifies an output directory to save the files in
% Example call: saveSubjectMapToFiles(responseVectorsMap, 'response');
    for subjNum = keys(vectorsMap)
        subj = vectorsMap(subjNum{1});
        for sessionNum = find(~cellfun(@isempty,subj))'
            disp(sessionNum);
            str = sprintf('%s/subject%d_session%d', directory, subjNum{1}, sessionNum);
            csvwrite(str, subj{sessionNum});
        end
    end
end

