# Automaticity Model
Automaticity model for the Ashby lab.

## Requirements
* MATLAB R2016B
* MATLAB 9.1
* Curve Fitting Toolbox
* Global Optimization Toolbox
* MATLAB Coder
* MATLAB Compiler
* Optimization Toolbox
* Parallel Computing Toolbox
* Statistics and Machine Learning Toolbox

## Project Structure (at-a-glance)
* Directories
** codegen: Automatically created directory upon compiling automaticityModel.m
** datasets: MAT files containing sample datasets for use with the model
** figures
** fmri: ALl files/code relevant to the FMRI configuration and optimization
** reference: Old/experimental code kept for reference, but not used in the current code
** scripts: One-off scripts
** visinputgen: Functions to be used for generating visual input matrices
* .gitignore: Informs version control software what files/paths to ignore when updating code repositories
* absorbstruct.m