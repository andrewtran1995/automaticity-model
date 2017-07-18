nameRegex = '[0-9]+';

stimVectors = containers.Map('KeyType','double','ValueType','any');
formatSpec = '%f';
% sizeArr = [Inf 1];


CROSSHAIRS = 'crosshairs';
STIMULUS = 'stimulus';
RESPONSE = 'response';

VECTOR_TYPE = RESPONSE;

files = dir(strcat(VECTOR_TYPE, '/subject*'));
for file = files'
    matchStr = regexp(file.name,nameRegex,'match');
    matchNum = str2double(matchStr);
    if isKey(stimVectors,matchNum(1))
        subject = stimVectors(matchNum(1));
    else
        stimVectors(matchNum(1)) = cell(20,1);
        subject = stimVectors(matchNum(1));
    end
    
    fileID = fopen(strcat(VECTOR_TYPE, '/', file.name), 'r');
    trVector = fscanf(fileID, formatSpec);
    fclose(fileID);
    
    % Convert TR vector to stim/crosshairs/response vector (0s and 1s, millisecond scale)
    trVector = floor(trVector*1000);
    if strcmp(VECTOR_TYPE, STIMULUS)
        vector = repelem(trVector,1,2000);
        for i = 1:size(vector,1)
            if vector(i) ~= 0
                vector(i,:) = [repelem(1,vector(i)), repelem(0,2000-vector(i))];
            end
        end
        subject{matchNum(2)} = reshape(vector',[],1);
    elseif strcmp(VECTOR_TYPE, CROSSHAIRS)
        vector = repelem(trVector,1,2000);
        for i = 1:size(vector,1)
            if vector(i) ~= 0
                vector(i,:) = [repelem(0, 1000), repelem(1,1000)];
            end
        end
        subject{matchNum(2)} = reshape(vector',[],1);
    elseif strcmp(VECTOR_TYPE, RESPONSE)
        vector = repelem(trVector,1,2000);
        for i = 1:size(vector,1)
            if vector(i) ~= 0
                vector(i,:) = repelem(1,2000);
            end
        end
        subject{matchNum(2)} = reshape(vector',[],1);
    end
    stimVectors(matchNum(1)) = subject;
end

% save('stimVectors','stimVectors');