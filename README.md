# **Hybrid Tuner for Hyperparameter Tuning**

Install the following directory with:
'pip install -e .' 

We have included several examples in the Example directory.
In each example file, modify 'matlab_path' to the location of your MATLAB installation, or remove it if you do not have MATLAB installed.

To create your own example, modify the myexec file to call your black-box function!

Required Format:
1. Parameters are read from myin into black-box function
2. Output is written to myout

## **If you have MATLAB and attempt to run the example files with solver = 8:**

Bandit DFO will return:

"Bandit DFO has completed!"
"Best Solution = 1.0 found after 102 iterations!"

Hybrid DFO will return:

"Hybrid DFO has completed!"
"Best Solution = 1.0 found after 28 iterations!"

The single DFO solver SID-PSM will return:

"8 has completed!"
"Best Solution = 1.0 found after 60 iterations!"



