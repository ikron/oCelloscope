# oCelloscope
This repository contains scripts to load, process, and analyse data from the oCelloscope machine.

## Usage
First thing to do is to export the data from the oCelloscope UniExplorer program using the "export csv" option. This csv file then contains the data output by the analysis module that was used. The format in which the data is exported requires that the file has to be parsed first before it can be read into R.

**Some important things to note**: UniExplorer uses the settings of the computer about field and decimal separators (and file encoding?). In our current lab computer these are decimal separator = ".", field separator = "," and encoding = "UTF-16". Use these as parameters for the R read.ocelloscope.data function. Note that if you had any files that were run using the old setup and computer, decimal separator is "," and field separator is ";" for these files.

### Tutorial

First load the some R functions for the oCelloscope. These are in /scripts/oCelloscope_functions.R

´source("oCelloscope_functions.R")´
