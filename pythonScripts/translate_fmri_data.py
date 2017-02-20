"""
SCRIPT DESCRIPTION
Given the path to a directory containing .dat files and an (optional) path of a destination directory,
map any found .dat files to the new destination, standardizing file names and removing duplicate data sets.

REQUIREMENTS
* Python 2.7
* Numpy 1.12.0

REFERENCES
* Python List Comprehensions: https://docs.python.org/2/tutorial/datastructures.html#list-comprehensions
* Regular Expressions: https://docs.python.org/2/library/re.html#
* Regular Expressions Introduction: http://www.aivosto.com/vbtips/regex.html
"""
# Import required Python modules
import os
import sys
import shutil
import re  # regular expressions
import numpy as np  # really freakin' cool matrix features

# Define constants
DEFAULT_DST_DIR = os.path.join(os.path.dirname(__file__), "output")  # default destination is currentDirectory/output
DAT_SUFFIX = ".dat"
RE_SUBJECT_NUM = re.compile(r"(?<=sub_)[0-9]+(?=[_0][0-9])")
RE_SCAN_NUM = re.compile(r"[0-9](?=_data)")

# Define helper functions


def getDatFilesFromDirectory(dir_path):
    """
    Given a directory path, return a list of file objects non-recursively
    :param dir_path:
    :return: files
    """
    # Use a list comprehension to build a list of files if f meets a certain condition:
    # * is a file
    # * ends with correct suffix
    files = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f)) and f.endswith(DAT_SUFFIX)]
    # Return files
    return files


def standardize_name(filename):
    """
    Standardize a data file name such that subject number always takes up three characters, and underscores are used
    to separate subject number and trial number.
    :param filename:
    :return:
    """
    subject_num = re.search(RE_SUBJECT_NUM, filename).group()
    trial_num = re.search(RE_SCAN_NUM, filename).group()
    return "sub_{:03d}_{}_data.dat".format(int(subject_num), trial_num)


def flagged_name(standardized_name):
    return "dup_" + standardized_name


def load_dat_file(datname):
    return np.genfromtxt(datname, dtype=np.float)


def dat_files_identical(dat1, dat2):
    """
    Compares two dat files' contents and determines if they are equal.
    Unlike regular integers, equality for floats can not be done exactly, and must be treated
    with a level of error tolerance (i.e., "is it close enough?").
    :param dat1:
    :param dat2:
    :return: boolean indicating whether dat1 is the same matrix as dat2
    """
    return np.allclose(load_dat_file(dat1), load_dat_file(dat2))


# Following code is called if file is executed as such:
# python translate_fmri_data.py data_directory
if __name__ == '__main__':
    # Check if proper number of arguments (1 or 2 in addition to calling the script itself)
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Invalid number of arguments in sys.argv: {}".format(len(sys.argv)))
        print("Proper usage: python translate_fmri_data.py src_dir [dst_dir]")
        print("sys.argv: {}".format(sys.argv))
        sys.exit()
    # Get arguments
    src_dir = sys.argv[1]
    if len(sys.argv) is 3:
        dst_dir = sys.argv[2]
    else:
        dst_dir = DEFAULT_DST_DIR
    # Define dictionary to map new file names/paths to old file names/paths
    # to keep track of what "new" filenames are already taken
    filenames = dict()
    # Get all valid files
    for filename in getDatFilesFromDirectory(src_dir):
        new_name = standardize_name(filename)
        # If new_name is not "taken" by another file, insert it into the dictionary and point to the old filename
        if new_name not in filenames:
            filenames[new_name] = filename
        # Else, must resolve a name conflict between two files
        else:
            dat1 = filename
            dat2 = filenames[new_name]
            # If dat files are identical, the current filename is a duplicate of one we have already accounted for.
            # Therefore, do not insert it into the filenames dictionary and continue the loop.
            if dat_files_identical(dat1, dat2):
                continue
            # Otherwise, create a unique filename to indicate that this dat file is an anomaly.
            else:
                filenames[flagged_name(new_name)] = filename
    # Create dst_dir if it does not already exist
    if os.path.exists(dst_dir) is False:
        os.makedirs(dst_dir)
    # Create new files using dictionary
    for new_name, old_name in filenames.items():
        old_path = os.path.join(src_dir, old_name)
        new_path = os.path.join(dst_dir, new_name)
        shutil.copyfile(old_path, new_path)
