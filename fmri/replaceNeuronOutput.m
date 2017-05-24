function [ outputVectorMap ] = replaceNeuronOutput(outputVectors, stimVectorMap)
%REPLACENEURONOUTPUT Replace stimulus vectors w/ summed neural output
%   Replace stimulus vectors (binary vector) with summed PMC A and PMC B 
%   output from model, such that trials are aligned between data
    FMRI_META = struct('NUM_TRIALS', 11520, ...
                       'SES_1',      1:480, 'SES_4',    1681:2160, ...
                       'SES_10', 5161:5640, 'SES_20', 11041:11520);
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
            % Only look at runs of 1
            stimulusOn = arr(idx) == 1;
            idx = idx(stimulusOn);
            len = len(stimulusOn);
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
            % Replace runs with output
            output = outputVectors(session,:);
            newSession = zeros(length(arr),1);
            for k = 1:length(idx)
                % Discrepancy where outputVectors only has output for
                % 1000 time units max, whereas stimVectors has TRs of
                % 2000 time units
                newSession(idx(k):idx(k)+len(k)-1) = output(k,1:len(k));
            end
            newSubject{j} = newSession;
        end
        outputVectorMap(stimKeys{i}) = newSubject;
    end
end

