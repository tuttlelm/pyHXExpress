# [pyHXExpress](https://github.com/tuttlelm/pyHXExpress)
This is a python adaptation of [HX-Express](https://www.hxms.com/HXExpress/) [[1]](#1) to perform polymodal fits of Hydrogen-Deuterium Exchange Mass Spectrometry (HDX-MS) data. This package can be used for high throughput analysis of entire peptide pools and beyond. The pyHXExpress tools can be used in python scripts or used interactively within Jupyter Notebooks. The goal of this project is to allow users to easily interact with their data and subsequent outputs to gain statistically supported insights into the polymodality of the HDX data. 

## 01 - Overview

This is a general overview of how to use pyHXExpress. The easiest way to get started will be to try some of the provided Jupyter Notebooks that contain a variety of examples of processing and analyzing HDX-MS spectra. 


### Input Data

#### I. The Spectra

Input data is m/z vs Intensity for peptide/charge/timepoint/replicate. This can be in the format used by (1) [HX-Express](https://www.hxms.com/HXExpress/) or (2) as individual spectra as output by [HDExaminer](https://massspec.com/hdexaminer/). 


1. <b>HX-Express like tabular format</b>. All timepoints and replicates for a sample/peptide/charge are in one .xlsx or .csv file as adjacent columns. Some additional flexibility is possible here compared to the HX-Express requirements: time point labels do not need to be artifically incremented to reflect replicates as these will automatically be tallied and undeuterated (undeut) and total deuteration (TD) data can have replicates. Either the m/z or Intensity columns (but not both) must contain a timepoint label: undeut, TD, or time unit (e.g. 5 sec or 60 min). Supported units of time are sec, min, hour (or s, m, h). The other column label should be empty.

2. <b>HDExaminer peptide pool exported spectra</b> [[2]](#2). This output is typically in a folder called 'SpecExport' with subfolders matching the sample names. Within each sample folder are peptide folders that contain a csv file corresponding to every timepoint/charge/replicate.

A note on spectral quality: the best results will be obtained on clean data. This means the peaks in the spectra should fall on the expected m/z values and any contaminating peptides should be minimized. Because the purpose of this software is to detect polymodal data, it is necessarily assumed that any Intensities at the expected m/z values are real. A sign that the spectra for a petpide/charge contains contaminating peptides is when the UN/TD spectra are polymodal. This will be flagged in the output report.

#### II. The Metadata 

A metadata dataframe ('metadf') is used to control what peptide/charge/replicate data will be fit. This is a csv file that can be provided by the user or pyHXexpress can create it by inferring the information from the file names. 

1. To generate the metadf file For HX-Express type data, the filenames must be of the following format: SAMPLE_START-END-PEPTIDE-zCHARGE-usertext.xlsx or .csv. SAMPLE: a sample identifier, START: starting residue number of the peptide, END: last residue number, PEPTIDE: the entire peptide single letter amino acid sequence, CHARGE: the charge state. 

2. To generate the metadf file for HDExaminer exported spectra, the user need only include 'SAMPLE.fasta' files with the correct sequences (accounting for mutations) in the SpecExport directory where SAMPLE exactly matches the name of the sample subfolders in that directory. **It is highly recommended to remove any spaces in these file names. Underscores and dashes are your friends.


The minimum required columns in a metadf file are 'Index' (integer row identifier that becomes the 'data_id' in outputs),'file' (the filename with extension for tabular data or the peptide folder name for SpecExport data),'sample' (sample name),'start_seq' (begin residue number),'end_seq' (end residue number),'peptide_range' (begin-end),'charge' (integer charge state), and 'peptide' (correct single letter sequence). Any additional columns are allowed and will be ignored.

**If there are any peptide modifications, the user will need to manually create or modify the inferred metadf file to include a 'modification' column that contains any peptide modification present for each peptide. A phosphorylation modification (+HPO3) would be indicated by the composition: 'H:1 P:1 O:3'. Currently, the modification is assumed to be non-exchangable. Future versions will address this.


#### III. The User Configuration File
The default config.py file is required along with the hxex.py file (this is automatic if using pip install pyhxexpress). This file contains defaults for all of the potentially user-configurable settings that will specify input and output directories and numerous settings related to how the fits should be performed. Users should create their own user_config.py (any .py filename is allowed) and update the settings as appropriate. 

The example notebooks will go into greater detail about what each setting controls. Features include specifying the type of data (tabular or SpecExport), the min and max number of populations to be fit (*this can also be specified per dataset, see advanced examples), the p-value threshold to use to allow fitting of additional populations, a minimum population requirement, and much more.

Every effort has been made to allow the user to fully control how fits are performed, but the defaults in the config.py file are a good place to start for new users.

### References
<a id="1">[1]</a> 
Guttman M., Weis, D.D., Engen, J.R. & Lee, K. K. (2013).  Analysis of overlapped and noisy Hydrogen/Deuterium exchange data.  J. Amer. Soc. Mass Spectrom. 24(12), 1906-1912.
[DOI: 10.1007/s13361-013-0727-5](https://pubs.acs.org/doi/10.1007/s13361-013-0727-5).

<a id="2">[2]</a> From the [HDExaminer](https://massspec.com/hdexaminer/) documentation: "To export the raw spectra used for all calculated results in a Peptide Pool: switch to the Peptides View, then select a Peptide Pool or any peptide in that pool. Select “Tools”, then “Export”, then “Pool Spectra…” or right-click on any Peptide Pool or peptide and select “Export Pool Spectra…”. Select a folder. HDExaminer will save the raw spectrum data for each result corresponding to the selected Peptide Pool. The exported data will be saved to a subfolder called SpecExport, with another subfolder called Fragments in the case of middle down data. The individual spectra will be saved in csv format."

###### Last edited by LMT on 14April2024