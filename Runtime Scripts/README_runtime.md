# Runtime Tests Readme
The runtime tests are run in two languages, Python 3 and Matlab. Example results have been included so that all code can be run as is.

# # Running the code:
The Matlab ILP tests must be run before the python code and the output saved to this directory. This can be done by running the runTime.m file in Matlab and saving the tests cell array. 

After the Matlab ILP tests have been run, the saved cell array is pre-processed in Python using the Load_matlab_runtime.py script. This will save out the pre-processed runtime tests to a dictionary.

Grobner basis and brute-force testing are then performed using hammingGrobner.py. This file can also take command line arguments that change the parameters of the testing as well as to create runtime plots.