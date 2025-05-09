# [pyHXExpress](https://github.com/tuttlelm/pyHXExpress)

This work has recently been published on BioRxiv [[1]](#1) https://www.biorxiv.org/content/10.1101/2025.03.13.643099v1 <br>
Please cite accordingly if you use this repository. See the ['paper' branch](https://github.com/tuttlelm/pyHXExpress/tree/paper) for the exact versions of scripts used and all outputs. 

<br><br>
*** This is under active development. Please report any comments, issues, requests *** 

This is a python adaptation of [HX-Express](https://www.hxms.com/HXExpress/) [[2,3]](#2) to perform multimodal fits of Hydrogen-Deuterium Exchange Mass Spectrometry (HDX-MS) data. This package can be used for high throughput analysis of entire peptide pools and beyond. The pyHXExpress tools can be used in python scripts or used interactively within Jupyter Notebooks. The goal of this project is to allow users to easily interact with their data and subsequent outputs to gain statistically supported insights into the modality of the HDX data. 

## How to use this tool

You'll need python 3.8+ and Jupyter Notebook. You can use pip install to install, which should install all dependencies. See [00_Installation.md](Documentation/00_Installation.md) for more details and alternatives to this approach.You may want to grab the latest pyhxexpress/*.py files to your installed location to have the most up-to-date features. 


    pip install pyhxexpress


Once installed, there are additional guides in the [Documentation](Documentation/01_Overview.md) to help you get started.

<br>

##### New to python and Jupyter?
If you are brand new to python and Jupyter Notebook, there are many guides out there to help you get started. I highly recommend a course put together by the RCSB that is avialble online: [Python Scripting for Biochemistry & Molecular Biology](https://pdb101.rcsb.org/train/training-events/python)

A common way to use Jupyter Notebook is as installed through Anaconda, where you then work in a web browser tab. This is fine. Personally, I strongly prefer working in [VS Code](https://code.visualstudio.com/docs/datascience/jupyter-notebooks) which is going to be more similar to an RStudio environment if you have any experience there.  

<br>

### References
<a id="1">[1]</a> 
Tuttle, L. M., Klevit, R. E., Guttman M. (2025).
A framework for automated multimodal HDX-MS analysis. 
bioRxiv 2025.03.13.643099; [DOI: 10.1101/2025.03.13.643099](https://www.biorxiv.org/content/10.1101/2025.03.13.643099v1).
<br><br>
<a id="2">[2]</a> 
Guttman M., Weis, D.D., Engen, J.R. & Lee, K. K. (2013).
Analysis of overlapped and noisy Hydrogen/Deuterium exchange data.  J. Amer. Soc. Mass Spectrom. 24(12), 1906-1912.
[DOI: 10.1007/s13361-013-0727-5](https://pubs.acs.org/doi/10.1007/s13361-013-0727-5).
<br><br>
<a id="3">[3]</a> 
Tuttle, L. M., E. I. James, F. Georgescauld, T. E. Wales, D. D. Weis, J. R. Engen, A. Nath, R. E. Klevit and M. Guttman (2025). 
"Rigorous Analysis of Multimodal HDX-MS Spectra." J Am Soc Mass Spectrom.
[DOI: 10.1021/jasms.4c00471](https://pubs.acs.org/doi/10.1021/jasms.4c00471).
<br><br>

###### Last edited by LMT on 18Apr2025
