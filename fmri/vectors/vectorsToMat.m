nameRegex = '[0-9]+';

stimVectors = containers.Map('KeyType','double','ValueType','any');
formatSpec = '%f';
% sizeArr = [Inf 1];


files = dir('subject*');
for file = files'
    matchStr = regexp(file.name,nameRegex,'match');
    matchNum = str2double(matchStr);
    if isKey(stimVectors,matchNum(1))
        subject = stimVectors(matchNum(1));
    else
        stimVectors(matchNum(1)) = cell(20,1);
        subject = stimVectors(matchNum(1));
    end
    
    fileID = fopen(file.name, 'r');
    trVector = fscanf(fileID, formatSpec);
    fclose(fileID);
    
    % Convert TR vector to stim vector (0s and 1s, millisecond scale)
    trVector = floor(trVector*1000);
    vector = repelem(trVector,1,2000);
    for i = 1:size(vector,1)
        if vector(i) ~= 0
            vector(i,:) = [repelem(1,vector(i)), repelem(0,2000-vector(i))];
        end
    end
    subject{matchNum(2)} = reshape(vector',[],1);
    stimVectors(matchNum(1)) = subject;
end

% save('stimVectors','stimVectors');