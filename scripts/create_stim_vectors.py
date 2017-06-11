"""
Iterate through different folders in order to create box car vectors.
"""
import numpy as np
import os
import re


# Directories
SOURCE_DIR = 'D:\\Users\\Andrew\\Google Drive\\UCSB\\Computational Cognitive Neuroscience Research\\fMRI Behavioral Data Folder\\rb_automaticity_no_study_folders'
TARGET_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'fmri/crosshairs')

# Regular expressions
STIMULUS = 'Stimulus'
CROSSHAIRS = 'Crosshairs'
RESPONSE = '(Correct)|(Incorrect)|(Too)'


def tr_lines(f):
    for line in f:
        if line.rstrip() and line[0:2] == 'TR':
            yield line


def save_crosshairs_vectors_for_analysis(subjects_dict):
    for subject in subjects_dict.keys():
        for session in subjects_dict[subject].keys():
            old_vector = subjects_dict[subject][session]
            vector = np.repeat(old_vector, 2)
            vector[vector != 0] = 1
            vector[np.insert(np.diff(vector) == 1, 0, False)] = 0
            vector = np.repeat(vector, 1000)
            np.savetxt(os.path.join(TARGET_DIR, "subject{}_session{}".format(subject, session)), vector)


def save_response_vectors_for_analysis(subjects_dict):
    for subject in subjects_dict.keys():
        for session in subjects_dict[subject].keys():
            vector = subjects_dict[subject][session]
            vector = np.repeat(vector, 2000)
            vector[vector != 0] = 1
            np.savetxt(os.path.join(TARGET_DIR, "subject{}_session{}".format(subject, session)), vector)


def get_tr_dict(event_type, save_results):
    subjects = dict()
    for dirname in (filepath for filepath in os.listdir(SOURCE_DIR) if
                    re.fullmatch('Sub[0-9]+_Ses[0-9]+', filepath) is not None):
        [subject, session] = re.findall('[0-9]+', dirname)
        print(subject + ' ' + session)
        # Get handles on files of interest
        expertise_data = None
        log_file = None
        tr_count = 0
        for filename in os.listdir(os.path.join(SOURCE_DIR, dirname)):
            if re.fullmatch('.*data\.dat', filename) is not None:
                expertise_data = np.loadtxt(os.path.join(SOURCE_DIR, dirname, filename))
            elif re.fullmatch('.*log\.dat', filename) is not None:
                log_file = filename
            elif re.fullmatch('block_length_[0-9]+\.dat', filename):
                tr_count = sum((int(line) for line in open(os.path.join(SOURCE_DIR, dirname, filename))))
            match = re.fullmatch('.*(data\.dat|log\.dat)', filename)
            if match is not None:
                if re.fullmatch('.*data\.dat', match.group()) is not None:
                    expertise_data = np.loadtxt(os.path.join(SOURCE_DIR, dirname, filename))
                else:
                    log_file = filename
        assert (expertise_data, log_file) is not (None, None)
        # If new subject, create a new entry in the subjects dict
        if subject not in subjects:
            subjects[subject] = dict()
        # Initialize array for the new session
        subjects[subject][session] = np.zeros(tr_count, dtype=float)
        with open(os.path.join(SOURCE_DIR, dirname, log_file), 'r') as f:
            for i, line in enumerate(tr_lines(f)):
                str_arr = line.split()
                # if str_arr[3] == 'Stimulus':
                if re.match(event_type, str_arr[3]):
                    if re.match(STIMULUS, str_arr[3]):
                        subjects[subject][session][i] = expertise_data[int(str_arr[-1]) - 1][-1]
                    else:
                        subjects[subject][session][i] = 2
                else:
                    subjects[subject][session][i] = 0
    print("Done!")
    if save_results:
        for subject in subjects.keys():
            for session in subjects[subject].keys():
                np.savetxt(os.path.join(TARGET_DIR, "subject{}_session{}".format(subject, session)), subjects[subject][session])
    return subjects

if __name__ == '__main__':
    get_tr_dict(CROSSHAIRS, True)
