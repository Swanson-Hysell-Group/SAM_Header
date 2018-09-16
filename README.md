# .SAM file creation in Python

This repository contains Python code that takes paleomagnetic core orientation data from a .csv template and converts into the .sam format. This file format is used by the RAPID paleomag software (http://sourceforge.net/projects/paleomag/) that is used to generate paleomagnetic and rock magnetic data with 2G Enterprise superconducting rock magnetometers.

## How to use the code
- Download the [ZIP of this repository](https://github.com/Swanson-Hysell-Group/SAM_Header/archive/master.zip)
- Enter data into the spreadsheet template (sam_sample_template.xlsx or sam_sample_template.csv) and then save as a csv. 
- At the command line, navigate into the folder with the program and template. Run the python script `mk_sam_file.py`, specifying the name of the .csv file created from the template.

    ```bash
    ~/$ mk_sam_file.py site.csv -o [output path - optional]
    ```
- The code should then generate a .sam header file as well as sample files for each sample in the site.
- For quick reference to basic instructions and examples, use `mk_sam_file.py --help`

##### Additional setup (optional, but you probably want to do this!)

- To call the `mk_sam_file.py` script outside of this directory, add the following line to your `~/.profile` or `~/.bash_profile`

    ```bash
    export PATH="<path to this repository>/SAM_Header/:$PATH"
    ```

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
- *core_dip* is the dip of the ''core plate'' which is the plane perpendicular to the core axis. This number is the conjugate of the plunge of the core axis.
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

## Dependencies

The code requires the standard scientific python modules numpy and pandas. The Anaconda distribution
(Python 3.6) is a relatively straight-forward way to get set up quickly
(https://www.anaconda.com/download/). Other necessary functions from the PmagPy project
(https://github.com/PmagPy/PmagPy/) that are dependencies for `mk_sam_file.py` have been collected
in `mk_sam_utilities.py` which is included in the repository such that you don't need to download
PmagPy for the program to run.

## Other built-in features
### Working with multiple `<site>.csv` files
It is often the case that header files need to be generated for multiple sites, each with their own
.csv file. To run `mk_sam_file.py` on multiple files, simply run the script with a sequence of file
names:
```bash
~/$ mk_sam_file.py Z01.csv Z02.csv Z03.csv Z04.csv [...]
```
You can also pass glob patterns:
```bash
~/$ mk_sam_file.py Z*.csv # matches all files listed obove
~/$ mk_sam_file.py Z0[1-3].csv # excludes Z04.csv
```
To create header files for all .csv files in your current working directory:
```bash
~/$ mk_sam_file.py --all
```

### Automate output directory
By default, header files are written to your current working directory. As shown in [How to
use](#how-to-use), you can specify an output directory for header files using the option `-o` or
`--output-path`. To automatically write header files to a subdirectory with the same name as the
.csv file, use the `-ad` or `--auto-dirs` flag:
```bash
~/$ mk_sam_file.py Z01.csv Z02.csv -ad
# writes Z01 header files to ./Z01/, Z02 header files to ./Z02/
```
The output directory will be created if it does not yet exist.
