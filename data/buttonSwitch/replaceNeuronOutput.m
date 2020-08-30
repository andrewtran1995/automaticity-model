function [ modelAlphaVectorMap ] = replaceNeuronOutput(outputVectors, latency, stimVectorMap)
%REPLACENEURONOUTPUT Replace stimulus vectors w/ summed neural output
%   Replace stimulus vectors (binary vector) with summed PMC A and PMC B 
%   output from model, such that trials are aligned between data
    FMRI_META = struct('NUM_TRIALS', 11520, ...
                       'SES_1',      1:480, 'SES_4',    1681:2160, ...
                       'SES_10', 5161:5640, 'SES_20', 11041:11520);
    % Assuming there will be no more than 2000 TRs for a session
    MAX_OUTPUT_LENGTH = 4000000;
    TR_LENGTH = 2000;
    stimKeys = keys(stimVectorMap);
    modelAlphaVectorMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for i = 1:length(stimVectorMap)
        % Initialize subject cell array
        newSubject = cell(20,1);
        % Get handle on subject
        subject = stimVectorMap(stimKeys{i});
        for j = find(~cellfun(@isempty,subject))'
            session = subject{j};
            switch j
                case 1
                    outputRegion = outputVectors(FMRI_META.SES_1,:);
                case 4
                    outputRegion = outputVectors(FMRI_META.SES_4,:);
                case 10
                    outputRegion = outputVectors(FMRI_META.SES_10,:);
                case 20
                    outputRegion = outputVectors(FMRI_META.SES_20,:);
            end
            newSession = ones(MAX_OUTPUT_LENGTH,1)*-1;
            sessionNum = 1;
            for k = 1:TR_LENGTH:length(session)
                if session(k) == 1
                    newSession(k:k+TR_LENGTH-1) = [outputRegion(sessionNum, 1:latency(sessionNum))'; zeros(TR_LENGTH - latency(sessionNum), 1)];
                    sessionNum = sessionNum + 1;
                else
                    newSession(k:k+TR_LENGTH-1) = 0;
                end
            end
            % Verify newSession's buffer did not fill up, then trim it
            assert(any(newSession == -1));
            newSession(newSession == -1) = [];
            % Assign newSession to appropriate cell of newSubject
            newSubject{j} = newSession;
        end
        modelAlphaVectorMap(stimKeys{i}) = newSubject;
    end
end

