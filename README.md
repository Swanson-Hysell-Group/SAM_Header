#.SAM file creation in Python

This repository contains Python code that takes paleomagnetic core orientation data from a .csv template and converts into the .sam format. This file format is used by the RAPID paleomag software (http://sourceforge.net/projects/paleomag/) that is used to generate paleomagnetic and rock magnetic data with 2G Enterprise superconducting rock magnetometers.

##How to use the code:

- Enter data into the spreadsheet template (sam_sample_template.xlsx or sam_sample_template.csv) and then save as a csv. 
- Run the python script mk_sam_file.py using command line specifying use of the .csv file you have saved:
```
~/$ python mk_sam_file.py site.csv
```
- The code should then generate a .sam header file as well as sample files for each sample in the site.

##Dependencies

This code utilizes Python functions from the PmagPy library (https://github.com/ltauxe/PmagPy). If the code is not working for you, download PmagPy and make sure you have added it to $PYTHONPATH.
