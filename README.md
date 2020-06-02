# **Hybrid Tuner for Hyperparameter Tuning**

## Installations
Install the following directory with:
'pip install -e .' 

## Dependencies
HybridTuner is usable with only Python libraries.

However, HybridTuner greatly benefits from the use of:
1. Matlab
2. TOMLAB glcDIRECT

## Examples
We have included several examples in the Example directory.
To run an example file, modify 'matlab_path' to the location of your MATLAB installation, or remove it if you do not have MATLAB installed. 

To tune your own application, modify the myexec file to call your black-box function.

Required Format:
1. Parameters are read from myin into black-box function
2. Output is written to myout
3. Parameter lower and upper bounds are defined in the file myparams.json

### **Example outputs: when using MATLAB and solver = 8:**
Executing "python example.py myparams.json" in the BanditDFO directory will return the following:

"Bandit DFO has completed"

"Best Solution = 1.0 found after 102 iterations"

![Image from Bandit DFO run](/examples/BanditDFO/banditResults.png)

Hybrid DFO will return:

"Hybrid DFO has completed"

"Best Solution = 1.0 found after 28 iterations"

![Image from Hybrid DFO run](/examples/HybridDFO/hybridResults.png)

The single DFO solver SID-PSM will return:

"8 has completed"

"Best Solution = 1.0 found after 60 iterations"

![Image from HOPSPACK DFO run](/examples/SingleSolver/10Results.png)

## Papers
* B. Sauk and N.V. Sahindis. HybridTuner: Tuning with hybrid derivative-free optimization initialization strategies. *Submitted* to ACM: TACO, 2 June, 2020.

## Support
* This work was conducted as part of the Institute for the Design of Advanced Energy Systems(IDAES) with funding from the Office of Fossil Energy, Cross-Cutting Research, U.S. Department of Energy. 
* This work used the Extreme Science and Engineering Discovery Environment (XSEDE),which is supported by National Science Foundation grant number ACI-1548562. Specifically, it used the Bridges system, which is supported by NSF award number ACI-1445606, at the Pittsburgh Super-computing Center (PSC). 
* We also gratefully acknowledge the support of the NVIDIA Corporation with the donation of the NVIDIA Tesla K40 GPU used for this research.
