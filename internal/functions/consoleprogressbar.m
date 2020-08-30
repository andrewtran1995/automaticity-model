function consoleprogressbar(str, iter, total)
    % Console progress bar
    % Resources used:
    % https://www.mathworks.com/matlabcentral/fileexchange/28067-text-progress-bar
    % https://www.mathworks.com/matlabcentral/fileexchange/30297-consoleprogressbar
    
    % Constants and variable initialization
    coder.extrinsic('num2str');
    coder.varsize('STR_CR');
    BAR_LENGTH = 40;
    persistent STR_CR_COUNT;
    if isempty(STR_CR_COUNT)
        STR_CR_COUNT = 0;
    end
    
    % Create progress bar
    percentage = iter/total;
    fullbar = repmat('.', [1, BAR_LENGTH]);
    completedbar = round(percentage * BAR_LENGTH);
    fullbar(1:completedbar) = '#';
    fullbar = ['[', fullbar, ']'];
    
    strout = [fullbar, ' ', num2str(iter), '/', num2str(total), ' ', str];
    
    % Print bar
    if iter == 1
        fprintf(strout);
    else
        fprintf([repmat('\b', 1, STR_CR_COUNT), strout]);
    end
    
    % Update carriage return
    STR_CR_COUNT = length(strout);
end