# `/internal`
This directory contains the bulk of the code for the model.
* `classes` contains the majority of the (object-oriented programming) classes used in the model. These include objects that represent state (e.g., a neuron) or objects that group related constants or logic (e.g., Hebbian constants). Class inheritance is used extensively for the `Neuron` class and its sub-classes in order to represent the many variations of neurons used in the model.
* `config` contains classes and parameters that represent specific **model configurations**. Common behavior is captured in `ModelConfig`, whereas any configuration-specific behavior should be described in the appropriate subclasses. Certain parameters are defined in `config/params`, especially those that needed to be exposed for global parameter optimization.
* `functions` contains various functions that are used throughout the rest of the code.
* `optimize` contains functions/scripts used for global parameter optimization.