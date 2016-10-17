# UFDiMuAnalysis
```
OKAY

Stage 1 Analyzer
------------------------------------------------------------------------
Latest version is the extra_mu_CMSSW_8_0_X branch
    -> https://github.com/acarnes/UfHMuMuCode/tree/extra_mu_CMSSW_8_0_X

Important scripts, besides the analyzer itself
    - Samples_v3.py
        +the database for all of the samples I'm using
        + can set files = [] if you don't want to copy the list of files from das to some local text file
            * then it will still run with crab, but not locally
    - make_crab_script.py
        + 'python make_crab_submit.py' to run
        + makes all of the crab submission scripts and analyzer python files needed for crab submission
        + tell it which json file and sample or collection of samples from Samples_v3.py to use by editing the file (obvious if you open it up)
        + takes templates for submission scripts and edits them appropriately
            * should change this: config.Data.outLFNDirBase = '/store/user/acarnes/'
            * and this:           config.Site.storageSite = 'T2_US_Florida'
            * in templates/crab_template.py to your own T2 storage info if you want the crab submission to work
        
    - crabsubmitall
        + 'source crabsubmitall' to submit all of the crab jobs that make_crab_script.py made
    - crabcleansubmit
        + 'source crabcleansubmit' to get rid of all of the submission scripts that make_crab_script.py made
    - crabstatusall
        + 'source crabstatusall' to check the status of all the crab jobs


Stage 2, plotting variables, limits, compare to fewz, etc (This repo)
------------------------------------------------------------------------    
Latest version is the more_muons branch
    -> https://github.com/acarnes/UFDimuAnalysis/tree/more_muons
    -> working at /home/acarnes/h2mumu/UFDimuAnalysis_v2 on the ihepa servers
    -> already made stage1 samples at /cms/data/store/user/t2/users/acarnes/h2mumu/samples/stage1
    
If you're not sure what some of the code does, it's usually commented reasonably well so check it out.

to compile a script edit the makefile and change MAIN to the name of the script without the .cxx extension

#### Directories #######################################################
    + bin
        - executables that you run to do stuff
            * make plots, compare to fewz, etc
        - write a .cxx file then compile it by running make, same process for compiling already available scripts
            * need to change MAIN = blahblah to the name of your script without the .cxx extension
            * then run with ./yourscript option1 option2 ...
        - will explain the important scripts later on in this documentation ...
    + lib
        - library of objects that you shouldn't need to edit that help with plotting, sample operations, etc
        - DataFormats.h needs to be the same DataFormats.h that was used with the analyzer to make the samples you are using
        - If the analyzer changes we need to update DataFormats.h and maybe VarSet.h, and Sample.cxx
            + DataFormats.h defines the data structures that we load the infro from the sample ttree into
            + VarSet.h is a datastructure that groups the different objects from DataFormats.h into a single collection
                * Samples.cxx uses VarSet.h like so, to get the recoMuon pt for mu 0 we would do Sample->vars.recoMuons[0].pt
                * so to access something you may want to plot or explore in your script it's usually like, vars.interesting_quantity
                * where interesting_quantity is some data structure from DataFormats.h
            + Sample.cxx allows us to access different variables, keeps track of information about the sample, what root file, xsec, etc
                * if the analyzer changes and DataFormats.h changes you may need to edit Sample->setBranchAddresses to find 
                  the new changes to the data structures 
                * sample->getEntry(int i) loads all of the info for that event into sample->vars
                * getWeight gets the weight for the event depending on the genWeight and the PU reweighting
                *
    + python
        - any python scripts are in here
        - currently use python to make root files and datacards for limit setting and to make the histograms from FEWZ output
    + selection
        - There are objects here that do the Muon and Event selection and the categorization for the analysis
        - Interfaces defined in .h files
        - The specific implementations for cuts and cateogory selections are in the .cxx files
            * I have different muon and event selections in SelectionCuts.cxx
            * and CategorySelectionRun1 is defined in CategorySelection.cxx
        - There are also tools to select good jets depending on our selection criteria
    + tools
        - General tools to help with different things: output event information to files, out event info to terminal,
          help make plots for pileup reweighting, make 4vecs from gen particle info, calculate dR, blah blah
          
#### Workflow ##########################################################
    * Make plots in every category for the run1 categories
        + code in categorize.cxx
        + run './categorize var_to_plot bins' in the bin directory 
            - var_to_plot goes from 0 to 9 tells the executable what variable to plot, 0 is dimuon mass for instance
              the numbers are mapped inside the .cxx file
            - bins is boolean, only matters when plotting the dimuon mass, 0 will plot mass from 50 to 200
              1 will plot mass for limit setting, mass from 110 to 160
            - so ./categorize 0 1 would plot dimu_mass in all run1 categories from 110 to 160
            - edit categorize.cxx and make if you want to change anything
            - will need to edit the make file and set MAIN to category then run 'make' to recompile the executable
            - right now 8_0_X MC has no HLT trigger info, so I just scale MC by about 91% to simulate trigger efficiency
            - in reality we need scale factors for pt and eta from muon pog, but this isn't implemented yet 
        + plots are saved in bin/rootfiles, some images automatically created in bin/imgs
        + some png files of the plots are automatically saved in bin/imgs
        + if you want to plot all the variables just run 'source run_categorize_all.sh'
        + there are stacks with the MC where the data is overlayed in the stacks folder in the root file 
            - with Data/MC ratio plots beneath, fit with straight lines
        + the individual histograms for each sample and category are also available in the histos folder in the root file
        + the net signal and net background and net data plots used for limit setting are in the net_histos folder in the root file
    * Limit setting
        + I made limit setting code for analytic shape fit and template 
        + code in python/limit_setting
        + check out WorkspaceAndDatacardMaker.py
        + intialize the object with the file and category name as shown in WorkspaceAndDatacardMaker.py
        + this limit setting code needs the net signal and net background dimu_mass histos so you have to make the plots in the category first 
            - go to bin folder and run ./categorize 0 1 
            - this will create the dimu_mass plots that you need for every run1 category
            - point maker object to the output file and tell it which category to make the .root and datacard for
            - then use the .root and datacard with higgs combine
            - systematics aren't accounted for yet
            - it's also set up to use the net mc instead of data for limit setting since we are blinded   
    * Compare FEWz to MC
        + make dimu_mass histograms (and others) from fewz output
            - python/fewz/parseFEWZtoTH1Fs.py, 'parseFEWZtoTH1Fs.py --infile=fewz_output_text_file_to_parse.txt'
            - root files created in python/fewz/fewz_predictions/
        + make dimu_mass histograms for MC and Data
            - bin/fewzCompare.cxx -> './fewzCompare var_to_plot useRecoMuCuts useRecoToPlotMu useRecoJetCuts useRecoToPlotJets cutNoGens'
            - var_to_plot is usually 0 to plot dimu_mass unless you want MC and Data plots for fewz categories for some other variable
            - the rest of the options are booleans
            - useRecoMuCuts: if true use reco values for muons for cuts and category selections otherwise use gen values
                % false if you want gen level comparisons
            - useRecoToPlotMu: if true plot the reco value not the gen value
                % false if you want gen level comparisons
            - useRecoJetCuts: if true use the reco values for the jets for cuts and category selections, should usually be true
                % fewz uses antikt jets so apples to apples is reco
            - useRecoToPlotJets: if true plot the reco value for jet vars and not the gen value in the category, usually set to true
            - cutNoGens: if true cut events where there is no gen muon pair from a Z or virutal gamma in the DY_MC
            - for gen level comparison between DYJets and fewz run ./fewzCompare 0 0 0 1 1 1
                # dimu_mass, cut on gen mu values, plot gen mu values, cut on reco jet values, plot reco jet values, cut events with no gen mu values
                # I apply muon isolation and tight muon id on reco to get rid of fakes, you can change this and recompile the executable if you want
        + overlay the dimu_mass spectrum for the different fewz categories and make comparison and ratio plots b/w DY_MC and FEWZ
            - code for this is bin/overlayFewz.cxx
            - have to run bin/fewzCompare first to create the dy_mc category plots
            - have to run python/fewz/parseFEWZtoTH1Fs.py to create the fewz category plots from the txt files
            - './overlayFewz dimu_mass_plots_created_by_fewzCompare.root lumi_of_plots make_log_plots' 
            - e.g. './overlayFewz rootfiles/11110_validate_dimu_mass_DY-FEWZ_MC_categories_3990.root 3990 0'
            - images automatically saved in python/fewz/img
            - root files saved in bin/rootfiles
```
