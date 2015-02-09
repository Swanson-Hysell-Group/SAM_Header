#.SAM creation system in python

###how to:

simply type data into the spreadsheet template provided and then save as a csv. Feed that csv into the python script mk_sam_file.py using command line and watch as the .sam file and all the sample files apear.

###reversing this process:

use the reverse_sam_file.py on the .sam file and sample files and you'll get a csv file of the format used by this system making instant spreadsheets from sam files. (currently not functional)

**dependencies** - this system relies on python functions from the PymagPy library by Lisa Tauxe which can be found here https://github.com/ltauxe/PmagPy, if files are not running check to make sure you've added PmagPy to $PYTHONPATH.
