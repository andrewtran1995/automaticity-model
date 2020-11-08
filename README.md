# Model on automaticity of rule-guided behaviors
A biologically-detailed computational model of how rule-guided behaviors become automatic.
* [How the model works](#How-the-model-works)
* [Understanding the code](#Understanding-the-code)

## How the model works
The model is run by specifying a **model configuration** and passing it as an argument (along with other parameters, if desired) to the function **automaticityModel** (found in `automaticityModel.m`).

The **automaticityModel** runs through logic which can be thought of in the following stages:
* *Pre-processing* - The model initializes variables and objects.
	* Initialize configuration-specific and configuration-agnostic variables.
	* Create neuron objects.
* *Processing* - The model simulates a specified number of trials.
	* Set a visual stimulus and simulate its effects on the prefrontal cortex (PFC) and premotor cortex (PMC) neurons.
	* Simulate the effects of the visual stimulus on the set of neurons in the model over a discrete time interval.
	* Record data for the trial on which neurons exhibited the strongest reactions.
	* If applicable, apply Hebbian learning to modify local weights for a neuron.
* *Post-processing* - The model displays results or does other processing (e.g., calculating global parameter optimization values).

### Model configurations
A **model configuration** can affect how the model behaves, especially for the  purpose of replicating real-world phemonemon. This includes (but is not limited to):
* Affecting how many trials the model runs through.
* Determining for which trials learning is active.
* Deciding if the COVIS model of rule learning is enabled.
* Deciding if the FROST model of working memory maintenance is used.

Configurations are represented as (object-oriented programming) classes. Individual configurations inherit from `ModelConfig`, which contains logic shared between all configurations. These classes can be found in `internal/config`.

The following configurations are available, along with their in-code class name.
* **Electrophysiology**/`ModelConfigElectro` - Wallis, Jonathan D., and Earl K. Miller. “From Rule to Response: Neuronal Processes in the Premotor and Prefrontal Cortex.” Journal of Neurophysiology, vol. 90, no. 3, 7 May 2003, pp. 1790–1806., doi:10.1152/jn.00086.2003.
* **Dual-task interference**/`ModelConfigDualTask` - Zeithamova, Dagmar, and W. Todd Maddox. “Dual-Task Interference in Perceptual Category Learning.” Memory &amp; Cognition, vol. 34, no. 2, 2006, pp. 387–398., doi:10.3758/bf03193416.
* **Button-switch interference**/`ModelConfigButtonSwitch` - Helie, S., et al. “Evidence for Cortical Automaticity in Rule-Based Categorization.” Journal of Neuroscience, vol. 30, no. 42, 2010, pp. 14225–14234., doi:10.1523/jneurosci.2393-10.2010. Referenced when working with FMRI data and optimization.

### Sample usage
Run the model using the **button-switch interference** configuration.
```matlab
init;
compile;
optional_parms = struct( ...
    'FMRI_META_GROUP_RUN', 0, ...
	'VIS_INPUT_FROM_PARM', 0, ...
	'visualInput', zeros(2) ...
);
automaticityModel_mex( ...
    getmodelparams(ModelConfigButtonSwitch()), ...
	optional_parms ...
)
```

## Understanding the code

### Toolbox requirements
#### Required to run model
* MATLAB R2020A
* MATLAB 9.1
* Curve Fitting Toolbox
* MATLAB Coder
* MATLAB Compiler

#### Used in other parts of the project
* Global Optimization Toolbox
* Optimization Toolbox
* Parallel Computing Toolbox
* Statistics and Machine Learning Toolbox

### Project structure
* `/data` - Matrix files and related code, categorized by configuration.
* `/external` - External library dependencies.
* [`/internal`](/internal) - Code internal to this project. Most of the implementation can be found here.
* `.gitignore` - Version control file that specifies which files/directories should be ignored when updating this repository.
* `automaticityModel.m` - The function which runs an automaticity configuration through a series of trials. Can be thought of as the "main" function for this project.
* `compile.m` - Generates performant code for `automaticityModel.m`.
* `init.m` - Utility script to recursively load all folders from the current working directory. Useful for ensuring all classes/functions are in memory before doing work with the model.

### Code-generation and performance
Unfortunately, running `automaticityModel` by itself can be rather slow. In order to achieve acceptable performance, it is necessary to leverage [code generation](https://www.mathworks.com/help/coder/index.html) to compile `automaticityModel` into more performant MEX code.