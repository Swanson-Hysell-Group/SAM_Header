# .SAM file creation in Python

This repository contains Python code that takes paleomagnetic core orientation data from a .csv template and converts into the .sam format. This file format is used by the RAPID paleomag software (http://sourceforge.net/projects/paleomag/) to generate paleomagnetic and rock magnetic data with 2G Enterprise superconducting rock magnetometers.

## How to use the code:

- Download the [zip of this repository](https://github.com/Swanson-Hysell-Group/SAM_Header/archive/master.zip)
- Enter data into the spreadsheet template (sam_sample_template.xlsx or sam_sample_template.csv) and then save as a csv.
- Navigate into the folder with the program and template. Run the python script mk_sam_file.py using command line specifying use of the .csv file you have saved:
```bash
~/$ python mk_sam_file.py site.csv [optional - output_path]
```
- The code should then generate a .sam header file as well as sample files for each sample in the site.

## Site fields:

The first six rows are site information. These following values pertain to the site.

- *site_id* is required and is the name of the site (e.g. GB20-)
- *site_name* is optional text that provides additional site information (e.g lava flows in the Gooseberry Basalts)
- *site_lat* is the latitude of the site in decimal degrees (required)
- *site_lon* is the longitude of the site in decimal degrees (required)
- *site_elevation* is the elevation of the site in meters

## Required sample fields:

Rows 8 and onward are sample information. These following values pertain to the samples.

- *magnetic_core_strike* is the strike of the ''core plate'' which is the plane perpendicular to the core axis. This number is the trend of the core axis + 90º.
- *core_dip* is the dip of the ``core plate'' which is the plane perpendicular to the core axis. This number is the conjugate of the plunge of the core axis.
- *bedding_strike* is the strike of the bedding using the right-hand rule. This value will be used for tilt-correcting the data. If no tilt-correction is necessary enter 0 for bedding_strike and 0 for bedding_dip.
- *bedding_dip* is the dip of the bedding. This value will be used for tilt-correcting the data.
- *mass* is the mass of the specimen in grams. If you do not wish to enter mass data, enter 1.0 for each specimen.

## Optional sample fields:

- *correct_bedding_using_local_dec* This field should either be 'yes' or 'no'. If 'yes' ['yes' is the default if the field is left blank which is why this field is optional] the local calculated IGRF declination will be used to correct the bedding strike. If 'no', the bedding strike will be left uncorrected.
- *shadow_angle* is the angle read from a sun compass. The code processes these data using the convention of a counter-clockwise sun compass (the type used on a Pomeroy orienting fixture). If a clockwise sun compass is used instead (we use these in our lab for block sampling), then the data need to be transformed to be counter-clockwise upon entry.
- *GMT_offset* is the time difference between the local time and Greenwich Mean Time. What should be entered is the number of hours to SUBTRACT from local time to get GMT. For example, Ethiopia is 3 hours ahead of GMT so the value that should be entered is 3. In the summer months in Minnesota, the time is CDT which is 5 hours behind GMT so the value that should be entered is -5.
- *year*,	*month*,	*days*,	*hours*,	*minutes* are required date/time information if sun compass data are provided.

## Things to know:

- **Sun compass data are preferentially used:** The code is currently set up so that if there are sun compass data those data are preferentially used for the sample orientations. If there are no sun compass data, the magnetic compass data are used and they are corrected for the local magnetic declination calculated from the model IGRF field. Note that in both cases, the local magnetic declination value in the .sam file is set to be zero since the orientations are already corrected.
- **Inspect local magnetic declination vs. IGRF when running program:** When there are coexisting magnetic and sun compass data, the difference between them (which is the local magnetic declination) is printed into the .csv file. If this calculated local magnetic declination is more than 5º away from the model IGRF field, a warning is printed to the terminal while the program is executing. It is recommended to examine the modified .csv file after the code is executed to inspect these calculated local magnetic declination values. If the values are all over the place, it is likely that something is wrong related to data entry (such as GMT value or CW instead of CCW sun compass values).
- **Decide whether or not bedding orientations need to be corrected for the local magnetic declination:** Bedding orientation should be entered as strike and dip collected using the right-hand rule. The column 'correct_bedding_using_local_dec' takes either 'yes' or 'no'. If 'yes' the local calculated IGRF declination will be used to correct the bedding strike. If 'no', the bedding strike will be left uncorrected.
- **Recognize that the program assumes counter-clockwise sun compass data:** see explanation in the *shadow_angle* entry above.
- **Make sure that the GMT_offset is the hours to subtract from local time to get to GMT:** see explanation in the *shadow_angle* entry above.
- **To run in the RAPID software these files need to be a folder that corresponds to the site name** The RAPID software will not recognize the .sam if it is not in a folder with the same name.

## Dependencies

The code requires the standard scientific python modules of numpy, scipy and pandas. Installing the Enthought Canopy python distribution (https://www.enthought.com/products/canopy/) is a way you can get quickly setup with python and the dependencies needed to run this code. The Anaconda distribution (Python 2.7) is another way to get set-up that is relatively straight forward (https://www.continuum.io/downloads). Other necessary functions from the PmagPy project (https://github.com/PmagPy/PmagPy/) that are dependencies for mk_sam_file.py have been collected in  mk_sam_utilities.py which is included in the repository such that you don't need to download PmagPy for the program to run.

## Tips
If your directory structure follows the general format of ```./<site>/<template>.csv``` and you have multiple templates ready for conversion, you might find the following command line regex useful:

```find */*\.csv -exec python mk_sam_file.py '{}' \;```

This regex navigates into all subdirectories, generates .sam header files for any CSV files it finds, and outputs them into the same subdirectory. Note that you must add the SAM_Header folder to your path so that the ```mk_sam_file.py``` program is accessible outside of your current directory. To do this, open ```.profile``` or ```.bash_profile``` (located in your home directory) in a text editor and add the following line:

```export PATH=<absolute path to ‘SAM_Header’ folder>/SAM_Header/:./:$PATH```
