# **Hybrid Tuner for Hyperparameter Tuning**

## Installations
Install the following directory with:
'pip install -e .' 

## Dependencies
HybridTuner is usable with only Python libraries.

However, HybridTuner greatly benefits from the use of:
1. MATLAB
2. TOMLAB glcDIRECT

All of the results reported in the corresponding paper make use of the TOMLAB and MATLAB solvers described within. 

## Examples
We have included several examples in the Example directory.
To run an example file, modify 'matlab_path' to the location of your MATLAB installation, or remove it if you do not have MATLAB installed. 

To tune your own application, modify the myexec file to call your black-box function.

Required Format:
1. Parameters are read from myin into black-box function
2. Output is written to myout
3. Parameter lower and upper bounds are defined in the file myparams.json

### **Example outputs: when using MATLAB and solver = 8:**
Executing "python example.py myparams.json" in the ./examples/BanditDFO directory will return the following:

"Bandit DFO has completed"

"Best Solution = 1.0 found after 105 iterations"

![Image from Bandit DFO run](/examples/BanditDFO/banditResults.png)

## Papers
* B. Sauk and N.V. Sahindis. HybridTuner: Tuning with hybrid derivative-free optimization initialization strategies. 

## Support
* This work was conducted as part of the Institute for the Design of Advanced Energy Systems(IDAES) with funding from the Office of Fossil Energy, Cross-Cutting Research, U.S. Department of Energy. 
* This work used the Extreme Science and Engineering Discovery Environment (XSEDE),which is supported by National Science Foundation grant number ACI-1548562. Specifically, it used the Bridges system, which is supported by NSF award number ACI-1445606, at the Pittsburgh Super-computing Center (PSC). 
* We also gratefully acknowledge the support of the NVIDIA Corporation with the donation of the NVIDIA Tesla K40 GPU used for this research.
