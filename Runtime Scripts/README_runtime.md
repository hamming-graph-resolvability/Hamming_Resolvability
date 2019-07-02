# Runtime Tests Readme
The runtime tests are run in two languages, Python 3 and Matlab. Example results have been included so that all code can be run as is.

## Running the code
The Matlab ILP tests must be run before the python code and the output saved to this directory. This can be done by running the runTime.m file in Matlab and saving the tests cell array. 

After the Matlab ILP tests have been run, the saved cell array is pre-processed in Python using the Load_matlab_runtime.py script. This will save out the pre-processed runtime tests to a dictionary.

Grobner basis and brute-force testing are then performed using hammingGrobner.py. This file can also take command line arguments that change the parameters of the testing as well as to create runtime plots.

## Running Timing Tests
The grobnerTests.py file is used to experimentally compare Grobner basis and brute force appraoches to determing resolvability with respect to time. This code accepts a number of command line arguments. These include:

``--a`` The alphabet size to consider. This is assumed to be at most 10.

``--funcNum`` An integer representing the resolvability method to test, (1) for brute force (2) for Grobner basis (3) for a multiprocess version of the Grobner basis method.

``--data`` A file from which to extract test information. This file contains a single header line. Subsequent lines contain k, a, a set of vertices in H_{k,a} separated by commas, and True if the set is resolving or False otherwise all separated by tabs. Each line is a single test instance.

``--maxAK`` Only onsider test instances with a*k less than or equal to this value.

``--repeats`` The number of times to run a resolvability check for each test instance.

``--procs`` The number of processes to use.

``--combine`` If the combine flag is set, gather collected timing data into a single file for each resolavability check method.

``--test`` If the test flag is set, print the number of tests for which data has been collected.

For example, ``python grobnerTests.py --funcNum 1 --a 2 --data res_set_data_3.tsv --procs 4 --maxAK 25`` runs brute force checks on hypercubes (a=2) from the res_set_data_3.tsv file over 4 processes where 2*k is at most 25.


