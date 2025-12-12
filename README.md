# oCelloscope
This repository contains scripts to load, process, and analyse data from the oCelloscope machine.

## Usage
First thing to do is to export the data from the oCelloscope UniExplorer program using the "export csv" option. This csv file then contains the data output by the analysis module that was used. The format in which the data is exported requires that the file has to be parsed first before it can be read into R.

**Some important things to note**: UniExplorer uses the settings of the computer about field and decimal separators (and file encoding?). In our current lab computer these are decimal separator = ".", field separator = "," and encoding = "UTF-16". Use these as parameters for the R read.ocelloscope.data function. Note that if you had any files that were run using the old setup and computer, decimal separator is "," and field separator is ";" for these files. It also used a different UniExplorer version, and it looks like some feature names have changed to the new one.

### Tutorial

First load the some R functions for the oCelloscope. These are in /scripts/oCelloscope_functions.R

`source("oCelloscope_functions.R")`

Load example data file: "Test1_scan_areas.csv" from /data. You have to set the correct number of repetions, for your data. We use often 100 repetitions. 

`testdata <- load.ocelloscope.data2("Test1_scan_areas.csv", sep = ";", dec = ",", repetitions = 100, features = c("EstimatedVolume", "Area", "NumberOfTips", "TipsPerArea", "TimestampInSeconds"))`

Note: this file was produced with the old computer and earlier UniExplorer version. So sep = ";" and dec = ",". For newer runs these should be sep = "," and dec = ".". Furthermore, in the old UniExplorer version feature names were different: Are occupied by the hyphae is called "Area" instead of "HyphaeArea" in the newer version. If loading a file that produced with the new setup you can leave the features parameter out as the default ones should work.

After loading the data, each well is a separate variable and different features are list items. For calculating growth rates, we'll convert the data into a long format

`longdata <- convert.ocelloscope.data.long(testdata)`

Then we can fir linear growth curves to the data, using the window based method

`gr.res <- ocelloscope.growth.rate(longdata)`

```> head(gr.res)
  sample        gr        R2      window tips.at.limit1
1     A1 13529.628 0.9377004     0:15000            898
2     A2  9371.063 0.9957154  5000:20000            263
3     A3  2190.374 0.9915097 15000:30000             56
4     B1 14749.068 0.9532670     0:15000            838
5     B2  6809.822 0.9970386 10000:25000            276
6     B3  4517.472 0.9666366 20000:35000             33 ```
