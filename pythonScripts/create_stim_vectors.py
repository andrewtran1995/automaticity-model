"""
Iterate through different folders in order to create box car vectors.
"""
import numpy as np
import os
import re


# Directories
SOURCE_DIR = 'D:\\Users\\Andrew\\Google Drive\\UCSB\\Computational Cognitive Neuroscience Research\\fMRI Behavioral Data Folder\\rb_automaticity_no_study_folders'
TARGET_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'fmri_data')

# Regex


if __name__ == '__main__':
    for dirname in (filepath for filepath in os.listdir(SOURCE_DIR) if re.fullmatch('Sub[0-9]+_Ses[0-9]+', filepath) is not None):
        [subject, session] = re.findall('[0-9]+', dirname)
        print(subject + ' ' + session)
        expertise_data = None
        log_file = None
        for filename in os.listdir(os.path.join(SOURCE_DIR, dirname)):
            match = re.fullmatch('.*(data\.dat|log\.dat)', filename)
            if match is not None:
                if re.fullmatch('.*data\.dat', match.group()) is not None:
                    expertise_data = np.loadtxt(os.path.join(SOURCE_DIR, dirname, filename))
                else:
                    log_file = filename
        assert (expertise_data, log_file) is not (None, None)
        with open(os.path.join(SOURCE_DIR, dirname, log_file), 'r') as f:
            print("all is well")
