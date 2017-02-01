# UFDiMuAnalysis
```
* Code used to process the Analyzed TTrees for the h2mu analysis at UF
    - Grabs info from the TTree to ...
    - make plots, categorize, synchronize with other groups, perform significance studies, resolution studies, etc
    - there are a variety of objects to help with
        % getting objects out of the TTree, plotting, bookkeeping, weighting, event selection, collection cleaning, parallelization, etc

* Create executables in the bin directory. Edit the makefile and set MAIN to the file you want to compile, e.g. MAIN=categorize.
    - This will compile categorize.cxx into categorize. Then run categorize with ./categorize.
    - bin/categorize.cxx is the most important executable
    - check it out to see how to use the code

* the makefile creates a shared library, so you need to update your LD_LIBRARY_PATH so that it can find it in the bin directory 
    - in bash: LD_LIBRARY_PATH=/path/to/UFDimuAnalysis/bin:${LD_LIBRARY_PATH}
    - in tcsh: setenv LD_LIBRARY_PATH /path/to/UFDimuAnalysis/bin:${LD_LIBRARY_PATH}

* bin/categorize.cxx produces the histograms and data vs MC stacks for the run1 or run2 categories depending on the inputs
    - ./categorize 1 -> run 1 categories
    - ./categorize 2 -> run 2 categories
    - see categorize.cxx for the other options

* create c++ executables in bin, python scripts in python
* the rest of the directories have objects that help with the h2mu analysis
    - Check out the ABOUT file in other folders to see how they help 
```
