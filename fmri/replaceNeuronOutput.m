function [ outputVectorMap ] = replaceNeuronOutput(outputVectors, latency, stimVectorMap)
%REPLACENEURONOUTPUT Replace stimulus vectors w/ summed neural output
%   Replace stimulus vectors (binary vector) with summed PMC A and PMC B 
%   output from model, such that trials are aligned between data
    FMRI_META = struct('NUM_TRIALS', 11520, ...
                       'SES_1',      1:480, 'SES_4',    1681:2160, ...
                       'SES_10', 5161:5640, 'SES_20', 11041:11520);
    % Assuming there will be no more than 2000 TRs for a session
    MAX_OUTPUT_LENGTH = 4000000;
    stimKeys = keys(stimVectorMap);
    outputVectorMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for i = 1:length(stimVectorMap)
        % Initialize subject cell array
        newSubject = cell(20,1);
        % Get handle on subject
        subject = stimVectorMap(stimKeys{i});
        for j = find(~cellfun(@isempty,subject))'
            arr = (subject{j})';
            % Find indices of start of runs, prepending 1 (since there
            % is always a run at the start)
            idx = find([1,diff(arr)] ~= 0);
            % Find the length of each run by taking the difference
            % between elements, appending one more for length of last run
            len = [diff(idx), length(arr) - idx(end) + 1];
            % Only look at runs of 0
%             stimulusOff = arr(idx) == 0;
%             idx = idx(stimulusOff);
%             len = len(stimulusOff);
            % Determine session
            switch j
                case 1
                    session = FMRI_META.SES_1;
                case 4
                    session = FMRI_META.SES_4;
                case 10
                    session = FMRI_META.SES_10;
                case 20
                    session = FMRI_META.SES_20;
            end
            % Fill in runs of no stimulus (0s), and replace runs of
            % stimulus (1s) with outputVectors' data
            % Use value of -1 to indicate uninitialized value
            output = outputVectors(session,:);
            newSession = ones(MAX_OUTPUT_LENGTH,1)*-1;
            lastIdx = 1;
            sessionNum = 1;
            for k = 1:length(idx)
                if arr(idx(k)) == 0
                    newSession(lastIdx:lastIdx+len(k)) = 0;
                    lastIdx = lastIdx + len(k);
                else
                    newSession(lastIdx:lastIdx+latency(sessionNum)) = output(sessionNum, latency(sessionNum));
                    lastIdx = lastIdx + latency(sessionNum);
                    sessionNum = sessionNum + 1;
                end
            end
            % Verify newSession's buffer did not fill up, then trim it
            assert(any(newSession == -1));
            newSession(newSession == -1) = [];
            % Assign newSession to appropriate cell of newSubject
            newSubject{j} = newSession;
        end
        outputVectorMap(stimKeys{i}) = newSubject;
    end
end

