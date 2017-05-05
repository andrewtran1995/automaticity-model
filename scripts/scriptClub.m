subject_regex = '(?<=sub_)[0-9]+';
trial_regex = '[0-9]+(?=data)';

% Get 1D
cd('D:\Users\Andrew\Documents\Google Drive\UCSB\Computational Cognitive Neuroscience Research\fMRI Behavioral Data Folder\New Behavioral Data\1d');
listing_1d = dir('*.dat');
trials_1d = {};
for i=1:size(listing_1d)
    disp(listing_1d(i).name);
    subject = str2double(regexp(listing_1d(i).name, subject_regex, 'match'));
    trial = str2double(regexp(listing_1d(i).name, trial_regex, 'match'));
    data = importdata(listing_1d(i).name, ' ');
    trials_1d{subject, trial} = data;
end


% Get Disjunctive
cd('D:\Users\Andrew\Documents\Google Drive\UCSB\Computational Cognitive Neuroscience Research\fMRI Behavioral Data Folder\New Behavioral Data\Disjunctive');
listing_disj = dir('*.dat');
trials_disj = {};
for i=1:size(listing_disj)
    disp(listing_disj(i).name);
    subject = str2double(regexp(listing_disj(i).name, subject_regex, 'match'));
    trial = str2double(regexp(listing_disj(i).name, trial_regex, 'match'));
    data = importdata(listing_disj(i).name, ' ');
    trials_disj{subject, trial} = data;
end

% Accessing non-zero entries
% [row, col] = find(~cellfun(@isempty,trials_disj))