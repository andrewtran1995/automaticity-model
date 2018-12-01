# Automaticity model
Automaticity model for the Ashby lab.

### Using the model
The model's main entry point is `automaticityModel.m`. For performance reasons, it has been written to be compatible with code generation.

### Requirements
* MATLAB R2018A
* MATLAB 9.1
* Curve Fitting Toolbox
* Global Optimization Toolbox
* MATLAB Coder
* MATLAB Compiler
* Optimization Toolbox
* Parallel Computing Toolbox
* Statistics and Machine Learning Toolbox

### Project structure (at-a-glance)
* Directories
	* codegen: Automatically created directory upon compiling automaticityModel.m
	* datasets: MAT files containing sample datasets for use with the model
	* figures
	* fmri: Files relevant to the FMRI configuration and optimization
	* reference: Old/experimental code kept for reference, but not used in the current code
	* scripts: One-off scripts
	* visinputgen: Functions to be used for generating visual input matrices
* `.gitignore`: Informs version control software what files/paths to ignore when updating code repositories
* absorbstruct.m

### Referenced papers
* Helie, S., et al. “Evidence for Cortical Automaticity in Rule-Based Categorization.” Journal of Neuroscience, vol. 30, no. 42, 2010, pp. 14225–14234., doi:10.1523/jneurosci.2393-10.2010. Referenced when working with FMRI data and optimization.
* Maddox, W. Todd, et al. “Response Time Distributions in Multidimensional Perceptual Categorization.” Perception &amp; Psychophysics, vol. 60, no. 4, June 1998, pp. 620–637., doi:10.3758/bf03206050.
* Wallis, Jonathan D., and Earl K. Miller. “From Rule to Response: Neuronal Processes in the Premotor and Prefrontal Cortex.” Journal of Neurophysiology, vol. 90, no. 3, 7 May 2003, pp. 1790–1806., doi:10.1152/jn.00086.2003.
