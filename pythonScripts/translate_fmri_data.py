"""
SCRIPT DESCRIPTION
Given the path to a directory containing .dat files and an (optional) path of a destination directory,
map any found .dat files to the new destination, standardizing file names and removing duplicate data sets.

REQUIREMENTS
* Python 2.7
* Numpy 1.12.0

REFERENCES
* Python List Comprehensions: https://docs.python.org/2/tutorial/datastructures.html#list-comprehensions
"""
# Import required Python modules
import os

# Define helper functions


def getDatFilesFromDirectory(dir_path):
    """
    Given a directory path, return a list of file objects non-recursively
    :param dir_name:
    :return: files
    """
    # Use a list comprehension to build a list of files if f meets a certain condition (isfile)
    files = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]
    # Return files
    return files


# Following code is called if file is executed as such:
# python translate_fmri_data.py data_directory
if __name__ == '__main__':
    pass