"""
Iterate through different folders in order to create box car vectors.
"""
import numpy as np
import os


# Directories
SOURCE_DIR = 'D:\\Users\\Andrew\\Google Drive\\UCSB\\Computational Cognitive Neuroscience Research\\fMRI Behavioral Data Folder\\rb_automaticity_no_study_folders'
TARGET_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'fmri_data')

# Regex


if __name__ == '__main__':
    for dirname in os.listdir(SOURCE_DIR):
        print(dirname)