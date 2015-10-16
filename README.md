#.SAM file creation in Python

This repository contains Python code that takes paleomagnetic core orientation data from a .csv template and converts into the .sam format. This file format is used by the RAPID paleomag software (http://sourceforge.net/projects/paleomag/) that is used to generate paleomagnetic and rock magnetic data with 2G Enterprise superconducting rock magnetometers.

##How to use the code:

- Enter data into the spreadsheet template (sam_sample_template.xlsx or sam_sample_template.csv) and then save as a csv. 
- Run the python script mk_sam_file.py using command line specifying use of the .csv file you have saved:
```bash
~/$ python mk_sam_file.py site.csv [optional - output_path]
```
- The code should then generate a .sam header file as well as sample files for each sample in the site.

##Things to know:

- The code is currently set up so that if there are sun compass data those data are preferentially used for the sample orientations. If there are no sun compass data, the magnetic compass data are used.
- When there are coexisting magnetic and sun compass data, the difference between them (which is the local magnetic declination) is printed into the .csv file. If this calculated local magnetic declination is more than 5º away from the model IGRF field, a warning is printed to the terminal while the program is executing.
- Bedding orientation should be entered as strike and dip. The column 'correct_bedding_using_local_dec' takes either 'yes' or 'no'. If 'yes' the local calculated IGRF declination will be used to correct the bedding strike. If 'no', the bedding strike will be left uncorrected.

##Dependencies

The code requires the standard scientific python modules of numpy and scipy. Installing the Enthought Canopy python distribution (https://www.enthought.com/products/canopy/) is a way you can get quickly setup with python and the dependencies needed to run this code. Other necessary functions that were originally part of the PmagPy project (https://github.com/ltauxe/PmagPy/) that are dependencies for mk_sam_file.py have been collected in  mk_sam_utilities.py which is included in the repository.
