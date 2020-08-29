# Automaticity model
Automaticity model for the Ashby lab.

### Using the model
The model's main entry point is `automaticityModel.m`. For performance reasons, it is compatible with code generation.

### Requirements
* MATLAB R2020A
* MATLAB 9.1
* Curve Fitting Toolbox
* Global Optimization Toolbox
* MATLAB Coder
* MATLAB Compiler
* Optimization Toolbox
* Parallel Computing Toolbox
* Statistics and Machine Learning Toolbox

### Project structure
* Directories
	* **classes** - classes and constants used within the model
	* **codegen** - automatically-created directory upon executing `automaticityModel_script.m`
	* **datasets** - MAT files containing sample datasets for use with the model
	* **figures**
	* **fmri** - files relevant to FMRI configuration and optimization
	* **libraries** - externally-sourced dependencies
	* **reference** - old/experimental code kept for reference that is not executed normally
	* **scripts** - one-off scripts, such as those used to process many files of experimental data
	* **visinputgen** - functions used for generating visual input matrices
* `.gitignore` - informs git version control which files/paths to ignore when updating repositories
* `automaticityModel_script.m` - compiles `automaticityModel.m` into much more performant `mex` code
* `automaticityModel.m` - the automaticity model
* `automaticityModel.prj` - GUI-loadable representation of the code-generation
* `automaticitymodelOpt.m` - wrapper around `automaticityModel.m` for global parameter optimization
* `particleSwarmOpt.m` - script to run global parameter optimization using particle swarm optimization

### Referenced papers
* Wallis, Jonathan D., and Earl K. Miller. “From Rule to Response: Neuronal Processes in the Premotor and Prefrontal Cortex.” Journal of Neurophysiology, vol. 90, no. 3, 7 May 2003, pp. 1790–1806., doi:10.1152/jn.00086.2003.
* Helie, S., et al. “Evidence for Cortical Automaticity in Rule-Based Categorization.” Journal of Neuroscience, vol. 30, no. 42, 2010, pp. 14225–14234., doi:10.1523/jneurosci.2393-10.2010. Referenced when working with FMRI data and optimization.
* Zeithamova, Dagmar, and W. Todd Maddox. “Dual-Task Interference in Perceptual Category Learning.” Memory &amp; Cognition, vol. 34, no. 2, 2006, pp. 387–398., doi:10.3758/bf03193416.