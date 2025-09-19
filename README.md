# NN_Classifier_Calibration_Stat_Code

# Xbb calibration using Zbb+Jet Event Topology

*  The repository contains a Binned template fit setup for the Xbb Calibration using Zbb+Jet event topolgy. The read-me given below is a basic documentation for the mechanism (running on Nominal datasets) and will be modififed as the analysis proceeds.
*  The fit mechanism is originally adapted from the [H--\>bb+j fitting software](https://gitlab.cern.ch/atlas-phys-exotics-dijetisr/xmlfit_boostedhbb).Some files, directory names and parts of the documentation outline are inherited from H--\>bb+j fitting.

*  The fit machinery uses `xmlAnaWSBuilder` as the baseline tool for fitting. The fit settings and links to inputs are given in the form of **xml** cards. The two main sources (for further reading and help) are:
    - [twiki](https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/XmlAnaWSBuilder)
    - @mstankai's [notes](https://docs.google.com/document/d/1NTpbmw9fV77SyHrFDUVxFCgFbC6EeOjuyI7TaNWNH74/edit?usp=sharing).


A `python` script is used to create the inputs to the fitter.
This script takes in a `json` file as input with links to input root files and various procedure flags.
The outputs are the following root files containing:
 - A data histogram to be given to the `xmlAnaWSBuilder` in the combination card;
 - Templates for all MC  with corresponding `RooWorkspace` to be given to the `xmlAnaWSBuilder` model cards;
 - All histograms for cross-checks (sanity-checks). 


### 0: Initial Set-Up (required only for the very first time) 

Cloning the repository:
```
git clone https://gitlab.cern.ch/shbansal/zjet_dxbb_templatefit.git
```
Getting `xmlAnaWSBuilder`:
```  
git clone https://gitlab.cern.ch/atlas-hgam-sw/xmlAnaWSBuilder.git
cd xmlAnaWSBuilder
mkdir -p ./Input/data/toy/
mkdir -p ./Input/data/histo/
```
Setting up soft links:
```
ln -s ../config/ config 
cd ..
ln -s xmlAnaWSBuilder/dtd dtd
```
Getting `quickFit` (fit framework for the xmlAnaWSBuilder output):
```
git clone https://gitlab.cern.ch/asciandr/quickFit.git
cd quickFit
git remote add upstream https://gitlab.cern.ch/atlas_higgs_combination/software/quickFit.git
ln -s ../xmlAnaWSBuilder/workspace workspace
```
-------------------------
### 1: Running workspace creation
The very formost step is creating the model workspace which is automatised and handled by a python script `modelMaker/simple_auto.py`
The configuration files are stored in the `json` folder (one example stored file is Xbb.json).

To set-up:
```
source setup_lxplus.sh
```
To run:
```
python modelMaker/simple_auto.py <path-to-config-json> <GeV-per-bin> <model-name> <tag>
```

Examples:
```
python modelMaker/simple_auto.py json/Xbb.json 2 Vbbfull test
```
(As `RooWorkspace` needs a `down` component for all systematics, the nominal histogram is given for one-sided systematics).

> **Note:** 
> It is advised to handle whatever functionalities possible in functions which can be added to `modelMaker/model_helper.py` script and called in `modelMaker/simple_auto.py` with `helpers.your_function()`, or write separate "helper" scripts that contain your functionalities (e.g. signal fitting is contained in `modelMaker/sigfit.py`)

##### Changing Fit Range
The fit range can be altered based on the flags set in the json file. Usually we are fitting from 50 to 150 GeV, but you can change these by adding:
```
"m_min" : "50.",
"m_max" : "200.",
"m_fitmin" : "50."
"m_fitmax" : "200."
```
to the `data` object in your json file.
  
##### Resonance Fitting
If the  `fit` flag in the json file is set for a template, parametric fitting will be applied to "smooth" the template. The fitting/histograming code is contained in `modelMaker/sigfit.py`, and the fitted functions and starting parameters are in `modelMaker\sigfit_init.py`. The fitting can be done in different binning schemes. This can be changed at the top of `modelMaker\simple_auto.py`

##### Ingnoring template in the RooWorkspace
Use the `ignoreWS` flag to not record a template in the workspace:
```
"ignoreWS" : true
``` 
this is useful if you need to use a  template object for some other reason but then omit it in the final output. The template histograms will still be recorded in the `crosscheck` output. By default this is `false` with the exception of systematic templates (see below).

-------------------------
### 2: Producing xml-cards

This step involves the creation of the the xml-cards which are the main basis and instructions for the fitting framework.
The code ran from **1** creates a `.json` file with the istructions to the xml-card-creator. One can find them at: `./genxml/models/<model-name>__<tag>.json`

Run command for xml-card creation:
```
python genxml/generate.py <list_of_json_names> --title <fit_title> -c <l/s> 
```

**Example run**
```
python genxml/generate.py Vbbfull__test --title Vbbfull_Xbb60 --tag test --poi 'mu_Zboson' --fr [50,150] --bins 50
```
Although here we tend to float the mu_Zboson (i.e the strength of the Z--\>bb process), this command should also be able to float the yield_Zboson, by simply replacing 'mu' with 'yield'

> **Note:** 
> Due to some modifications which are yet to be made in the modelMaker files, the mu_Zboson, mu_ttbar etc have to set by hand in the category level xml cards (for example in : config/Vbbj/Vbbfull_Xbb60_test/Vbbfull_category.xml): like for example, setting (in order to run the above command with 'mu_Zboson', 'mu_ttbar' etc) Subsequently here changes can be made to the NormFactor Name easily, if one wants to fix the strength, or constraint it with in a certain limit (stat. only)
```
<NormFactor Name="mu_Zboson[1,-50,50]"/>
```

-------------------------
### 3: Running the fit

Data and RooWorkspace with templates is then given to the xmlAnaWSBuilder to perform the fit.
To set up:
```
cd xmlAnaWSBuilder
source setup.sh
```
To produce a `RooWorkspace` file (running the fit):
```
./exe/XMLReader -x <path-to-top-level-card>  
```
**Example run (following the examples shown previously)** 
```
./exe/XMLReader -x ./config/Vbbj/Vbbfull_Xbb60_test/Vbbfull_Xbb60.xml
```
-------------------------
### 4: Performing a single best-fit with quickFit functionality:

After fitting with the `xmlAnaWSBuilder`, it produces a RooWorkspace which can then be manipulated using other, more sophisticated fitting frameworks.

Start by compiling the framework (Dont forget to do the initial make setup while cloning the quickFit repository in **0**):
```
cd ..
cd quickFit
source setup_lxplus.sh
```

For carrying out the single best-fit with `hesse+minos`, run the following command:
```
quickFit -f workspace/hbbj/Vbbfull_Xbb60/Vbbfull_Xbb60_model_test.root -d combData -p mu_Zboson=1_-50_50 -o output_Vbbfull_Xbb60.root --savefitresult 1 --hesse 1 --minos 1 --saveWS 1 --ssname quickfit
```
Additional documentation on `quickFit` can be found here: [git repo](https://gitlab.cern.ch/atlas_higgs_combination/software/quickFit) .

_Post-Fit plots can be made via the makePlot.C file stored in the Plotting directory, some paths to the input files needs to be changed there.__

-------------------------
### Background Modelling : Scripts to choose different functional forms optimsied on  Data sidebands and QCD dijet background:

The scripts in the directory OptFuncChoice_scripts contains two directories, one containing the macros for the optimsiation on data sidebands for different pT bins and inclusive in pT, and second for the optimsation on the full mass range of qcd dijet background. They can be simply run in root after making some adjustments to the input path for the root files, according your own inputs. 

P.S: The root files are here the in the form of unbinned dataset (as the optimsiation of the scripts are done via an unbinned Maximum Likelihood Fit)

-------------------------
### Reweighting of dijet simulated samples:

Due to the discrepancy between data and MC, the dijet simulated samples are reweighted w.r.t to data side bands, via a fit to data/MC distribution in the side band regions , the scripts to do this are collected in DijetReweight_Files directory, which contains two scripts, namely the first file (fittodatavMC.cxx) is to fit the data/MC (signal sideband) histogram via an analytical function and the second GetRewdijet.cxx is the macro to get the weight
The reweighted dijet MC is used for the evaluation of spurious signal

For any further questions while running the scripts, contact shubham.bansal@cern.ch
