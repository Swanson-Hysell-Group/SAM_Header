#!/usr/bin/env python

import os
import sys
import argparse
import textwrap
import re
import math
import pandas as pd
from mk_sam_utilities import *
from datetime import datetime as dt
from functools import reduce


def main(csv_file, output_directory):
    """
    NAME
        mk_sam_file.py

    DESCRIPTION
        Takes formated CSV and creates and writes a .sam header file and a set
        of sample files for any number of samples.

    SYNTAX
        ~/$ python mk_sam_file.py site.csv [optional - output_directory]

    OUTPUT
        .sam and sample files

    """

    ###########################################################################
    #                              Read In Files                              #
    ###########################################################################

    # fetching comand line data
    # file_name = sys.argv[1]
    file_name = csv_file
    directory = os.path.split(file_name)[0]

    # try:
    #     output_directory = sys.argv[2]
    # except IndexError:
    #     output_directory = directory

    if output_directory is not None:
        output_directory = os.path.expanduser(output_directory)
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
    else:
        output_directory = directory

    print('Reading in file - ' + file_name)

    df_cols = ['sample_name', 'comment', 'strat_level',
               'magnetic_core_strike', 'core_dip', 'bedding_strike',
               'bedding_dip', 'correct_bedding_using_local_dec',
               'mass', 'runs', 'sun_core_strike', 'calculated_IGRF',
               'IGRF_local_dec', 'calculated_mag_dec', 'core_strike',
               'corrected_bedding_strike']
    sdf_cols = ['sample_name', 'shadow_angle', 'GMT_offset',
                'year', 'month', 'days', 'hours', 'minutes']

    # file read in
    hdf = pd.read_csv(file_name, header=0, index_col=0, nrows=5, usecols=[0, 1])
    df = pd.read_csv(file_name, header=6, index_col=0,
                     usecols=df_cols, dtype=object).transpose()
    sdf = pd.read_csv(file_name, header=6, index_col=0,
                      usecols=sdf_cols).transpose()

    # variable assignments
    samples = df.keys()
    attributes = ['core_strike', 'core_dip',
                  'bedding_strike', 'bedding_dip', 'mass']
    site_values = ['site_lat', 'site_long']
    time_types = ['year', 'month', 'days', 'hours', 'minutes']

    ###########################################################################
    #                         Find Calculated Values                          #
    ###########################################################################

    print('---------------------LOCAL MAGNETIC DECLINATION-----------------------')

    # calculate sun_core_strike for all samples
    for sample in samples:
        if (not sdf[sample].isnull().any()):
            # df[sample]['sun_core_strike'] = float('nan')
            # continue
            # else:
            time_values = []
            for i in range(len(time_types)):
                time_values.append(str(int(sdf[sample][time_types[i]])))
            assert (len(time_values[0]) == 4),\
                "must input full year for sun compass calculation (i.e. YYYY)"
            sundata = {}
            if (len(time_values[1]) == 1):
                time_values[1] = '0' + time_values[1]
            if (len(time_values[2]) == 1):
                time_values[2] = '0' + time_values[2]
            if (len(time_values[3]) == 1):
                time_values[3] = '0' + time_values[3]
            if (len(time_values[4]) == 1):
                time_values[4] = '0' + time_values[4]
            sundata['date'] = reduce(lambda x, y: x + ':' + y, time_values)
            sundata['lat'] = hdf['site_info']['site_lat']
            sundata['lon'] = hdf['site_info']['site_long']
            sundata['shadow_angle'] = sdf[sample]['shadow_angle']
            sundata['delta_u'] = sdf[sample]['GMT_offset']
            df[sample]['sun_core_strike'] = round(sundec(sundata), 1)

    # calculate IGRF
        if (sdf[sample]['GMT_offset':'month'].isnull()).any():
            raise ValueError("not enough data to calculate IGRF to correct "
                             "bedding please input at least GMT_offset, "
                             "year, month, day of measurement\n")
        else:
            if math.isnan(float(hdf['site_info']['site_elevation'])):
                hdf['site_info']['site_elevation'] = 0.0
            for time_type in time_types:
                if math.isnan(float(sdf[sample][time_type])):
                    sdf[sample][time_type] = 1
            date = to_year_fraction(dt(int(sdf[sample]['year']),
                                       int(sdf[sample]['month']),
                                       int(sdf[sample]['days']),
                                       int(sdf[sample]['hours']),
                                       int(sdf[sample]['minutes'])))
            df[sample]['calculated_IGRF'] = list(
                    igrf([date,
                          float(hdf['site_info']['site_elevation'])/1000,
                          float(hdf['site_info']['site_lat']),
                          float(hdf['site_info']['site_long'])]))
            if float(df[sample]['calculated_IGRF'][0]) > 180:
                df[sample]['IGRF_local_dec'] = df[sample]['calculated_IGRF'][0] - 360
            else:
                df[sample]['IGRF_local_dec'] = df[sample]['calculated_IGRF'][0]
            # print out the local IGRF
            print(hdf['site_info']['site_id'] + str(sample) + " has local IGRF declination of: ")
            print(df[sample]['IGRF_local_dec'])

        # calculate magnetic declination
        print('The local declination calculated through magnetic and sun compass comparison is:')
        if math.isnan(float(df[sample]['sun_core_strike'])) or \
                math.isnan(float(df[sample]['magnetic_core_strike'])):
            df[sample]['calculated_mag_dec'] = 'insufficient data'
            print('insufficient data')
        else:
            df[sample]['calculated_mag_dec'] = (float(df[sample]['sun_core_strike']) -
                                                float(df[sample]['magnetic_core_strike']))
            print(str(round(float(df[sample]['sun_core_strike']) -
                            float(df[sample]['magnetic_core_strike']), 2)))
            if abs(float(df[sample]['IGRF_local_dec']) -
                   float(df[sample]['calculated_mag_dec'])) > 5:
                print("WARNING: local IGRF declination & calculated magnetic "
                      "declination are more than 5 degree different")
        print('')
    print('')
    print('Site averages:')
    print('Average of local IGRF declination is: ' + str(df.transpose()['IGRF_local_dec'].mean()))
    # print('Average of calculated local declination: ' +
    #       str(df.transpose()['calculated_mag_dec'].mean()))
    print('')
    print('---------------------OUTPUT-----------------------')

    ###########################################################################
    #                         Create .SAM Header File                         #
    ###########################################################################

    # setting name
    sam_header = hdf['site_info']['site_name'] + '\r\n'

    # creating long lat and dec info
    for value in site_values:
        hdf['site_info'][value] = str(round(float(hdf['site_info'][value]), 1))
        if value == 'site_lat':
            sam_header += ' ' + hdf['site_info'][value]
        else:
            sam_header += ' '*(5-len(hdf['site_info'][value]) + 1) + hdf['site_info'][value]
    sam_header += ' '*(3) + '0.0'
    sam_header += '\r\n'

    # making writing sample info
    for sample in samples:
        sam_header += hdf['site_info']['site_id'] + str(sample) + '\r\n'

    # creating and writing file
    print('Writing file - ' + os.path.join(output_directory, hdf['site_info']['site_id'] + '.sam'))
    sam_file = open(os.path.join(output_directory,
                                 hdf['site_info']['site_id'] + '.sam'), 'w+')
    sam_file.write(sam_header)
    sam_file.close()

    ###########################################################################
    #                           Create Sample Files                           #
    ###########################################################################

    for sample in samples:
        # assign variables for easy refrence
        site_id = hdf['site_info']['site_id']
        if not math.isnan(df[sample]['runs']):
            runs = df[sample]['runs'].split(';')
        else:
            runs = []

        # decide which core_strike to use, default is sun_core_strike but if not supplied
        # magnetic_core_strike will be used
        if type(df[sample]['correct_bedding_using_local_dec']) == float and \
                math.isnan(df[sample]['correct_bedding_using_local_dec']):
            df[sample]['correct_bedding_using_local_dec'] = 'yes'
        if not math.isnan(df[sample]['IGRF_local_dec']):
            if math.isnan(df[sample]['sun_core_strike']):
                df[sample]['core_strike'] = (float(df[sample]['magnetic_core_strike']) +
                                             float(df[sample]['IGRF_local_dec']))
                df[sample]['comment'] = 'mag compass orientation (IGRF corrected)'
            else:
                df[sample]['core_strike'] = float(df[sample]['sun_core_strike'])
                df[sample]['comment'] = 'sun compass orientation'

        if ((df[sample]['correct_bedding_using_local_dec']) == 'yes' or
                (df[sample]['correct_bedding_using_local_dec']) == 'Yes' or
                (df[sample]['correct_bedding_using_local_dec']) == 'YES'):
            df[sample]['corrected_bedding_strike'] = (float(df[sample]['bedding_strike']) +
                                                      float(df[sample]['IGRF_local_dec']))

        comment = df[sample]['comment']

        # check for no comment
        if type(comment) == float and math.isnan(comment):
            comment = ''

        # insure input is valid
        assert (len(site_id) <= 5),\
            "Locality ID excedes 5 characters: refer too "\
            "http://cires.colorado.edu/people/jones.craig/PMag_Formats.html "\
            "(although that says that 4 is the limit)"
        assert (len(comment) <= 255),\
            "Sample comment excedes 255 characters: refer too "\
            "http://cires.colorado.edu/people/jones.craig/PMag_Formats.html"
        assert (len(str(sample)) <= 9),\
            "Sample name excedes 9 characters: refer too "\
            "http://cires.colorado.edu/people/jones.craig/PMag_Formats.html"

        # write sample name and comment for sample file
        new_file = site_id + ' ' + str(sample) + ' ' + comment + '\r\n'

        # start second line strat_level get's special treatment
        if (math.isnan(float(df[sample]['strat_level']))):
            df[sample]['strat_level'] = "     0"
        df[sample]['strat_level'] = str((df[sample]['strat_level']))
        assert (len(df[sample]['strat_level']) <= 6),\
            "Length of strat_level excedes 6 characters: refer too "\
            "http://cires.colorado.edu/people/jones.craig/PMag_Formats.html"
        new_file += ' ' + ' '*(6-len(df[sample]['strat_level'])) + df[sample]['strat_level']

        # write in sample attributes on the second line
        for attribute in attributes:
            # assert (str(df[sample][attribute]).isdigit()),\
            #         str( attribute) + ' is a requred numeric value'
            if attribute == 'bedding_strike' and \
                        ((df[sample]['correct_bedding_using_local_dec']) == 'yes' or
                         (df[sample]['correct_bedding_using_local_dec']) == 'Yes' or
                         (df[sample]['correct_bedding_using_local_dec']) == 'YES') and\
                    not math.isnan(df[sample]['corrected_bedding_strike']):
                attribute = 'corrected_bedding_strike'

            if type(df[sample][attribute]) == float and math.isnan(df[sample][attribute]):
                if attribute == 'mass':
                    df[sample][attribute] = '1.0'
                    print(
                        "no mass found for sample %s, setting to default = 1.0 g" % (sample))
                else:
                    df[sample][attribute] = ''
            else:
                df[sample][attribute] = str(round(float(df[sample][attribute]), 1))

            # attributes must follow standard sam format
            assert (len(df[sample][attribute]) <= 5),\
                "Length of " + attribute + \
                " excedes 5 characters: refer too " + \
                "http://cires.colorado.edu/people/jones.craig/PMag_Formats.html"

            new_file += ' ' + ' '*(5-len(df[sample][attribute])) + df[sample][attribute]

        new_file += '\r\n'

        # if there are previous sample runs write that to the bottem of the file
        for run in runs:
            new_file += run + '\r\n'

        # create and write sample file
        new_file = new_file.rstrip('\r\n') + '\r\n'
        print('Writing file - ' + os.path.join(output_directory, site_id + str(sample)))
        sample_file = open(os.path.join(output_directory, site_id + str(sample)), 'w+')
        sample_file.write(new_file)
        sample_file.close()

    ###########################################################################
    #                        Write New Values to .csv                         #
    ###########################################################################

    csv_file = open(file_name, 'rU')
    csv_str = ''

    for i in range(5):
        csv_str += csv_file.readline()

    comma_count = csv_file.readline().count(',')
    csv_str += 'site_elevation' + ',' + \
               str(hdf['site_info']['site_elevation']) + ','*(comma_count-1) + '\n'
    # elev_line = csv_file.readline().split(',')
    # elev_line[1] = str(hdf['site_info']['site_elevation'])
    # reduce(lambda x,y: x + ',' + y, elev_line)

    header = csv_file.readline()
    csv_str += header
    header = header.strip('\r\n').split(',')

    for sample in samples:
        line = csv_file.readline()
        items = line.split(',')
        for i in range(len(header)):
            if i == 0:
                continue
            elif header[i] == 'calculated_IGRF':
                if type(df[sample][header[i]]) == str or type(df[sample][header[i]]) == float:
                    items[i] = str(df[sample][header[i]])
                else:
                    items[i] = str(list(df[sample][header[i]])).replace(',', ';')
            elif header[i] in df[sample].keys():
                items[i] = str(df[sample][header[i]])
            elif header[i] in sdf[sample].keys():
                items[i] = str(sdf[sample][header[i]])
            else:
                raise KeyError('there is no item: ' + header[i])
        csv_str += reduce(lambda x, y: x + ',' + y, items) + '\r\n'

    print('Writing file - ' + os.path.join(output_directory, hdf['site_info']['site_id'] + '.csv'))
    new_csv_file = open(os.path.join(output_directory, hdf['site_info']['site_id'] + '.csv'), 'w+')
    new_csv_file.write(csv_str)
    new_csv_file.close()

    generate_inp_file(output_directory, df, hdf)


def fix_line_breaks(file_name):
    """ Reads in the file given as a command line argument and rewrites it both line
        break types '\ r' and '\ n' so that python will for sure register all lines
    """
    # fix line breaks between different OS and python's default
    try:
        csv_file = open(file_name, 'r')
        csv_str = csv_file.read()
    except UnicodeDecodeError:
        csv_file.close()
        # Both pandas and the python open built-in sometimes complain about
        # encoding of these csv files. Switching to the encoding below generally
        # does the trick...no idea what the underlying problem is though
        # <09-08-18, Luke Fairchild> #
        csv_file = open(file_name, 'r', encoding="ISO-8859-1")
        csv_str = csv_file.read()
    if csv_str.find('\r\n') != -1:
        fixed_lines = csv_str.replace('\r\n', '\n')
    else:
        fixed_lines = csv_str.replace('\r', '\n')
    new_csv_file = open(file_name, 'w')
    new_csv_file.write(fixed_lines)
    csv_file.close()


def generate_inp_file(od, df, hdf):
    """
    DESCRIPTION
        Uses sample and site DataFrames from mk_sam_file.main function to generate inp file

        @param: od - output directory
        @param: df - sample Dataframe
        @param: hdf - site DataFrame

    OUTPUT
        .inp file

    """

    inps = ""
    inps += "CIT\n"
    inps += "sam_path\tfield_magic_codes\tlocation\tnaming_convention\tnum_terminal_char\tdont_average_replicate_measurements\tpeak_AF\ttime_stamp\n"
    inps += (os.path.join('.', hdf['site_info']['site_id'] + '.sam')) + '\t'
    if all(df.T['comment'] == 'sun compass orientation'):
        inps += 'SO-SUN\t'
    elif all(df.T['comment'] == 'mag compass orientation (IGRF corrected)'):
        inps += 'SO-MAG\t'
    else:
        inps += 'SO-SM\t'

    inps += (hdf['site_info']['site_name'] if hdf['site_info']['site_name']
             != '' or hdf['site_info']['site_name'] is not None else 'unknown') + '\t'

    # for some files the value below was not a string and was not being
    # converted to one like it is now in the line below---this was one of the
    # infrequent errors I came across testing this code
    #  <09-08-18, Luke Fairchild> #
    first_sample_id = str(df.keys()[0])

    if '-' in first_sample_id or '-' in hdf['site_info']['site_id']:
        inps += '2\t'
        first_sample_id.replace('-', '')
    elif '.' in first_sample_id:
        inps += '3\t'
        first_sample_id.replace('.', '')
    elif hdf['site_info']['site_id'] == first_sample_id:
        inps += '5\t'
    else:
        inps += '4\t'

    inps += str(len(first_sample_id)) + '\t'
    inps += "True\t"
    inps += "None\t"
    inps += '0.0\n'

    print('Writing file - ' + os.path.join(od, hdf['site_info']['site_id'] + '.inp'))
    if od != '' and not os.path.exists(od):
        os.makedirs(od)
    inpf = open(os.path.join(od, hdf['site_info']['site_id'] + '.inp'), 'w+')
    inpf.write(inps)
    inpf.close()


if __name__ == "__main__":
    # fix_line_breaks()
    prog_desc = textwrap.dedent("""\
            Using a formatted CSV, creates and writes a .sam header file and a
            set of sample files. The program can run on multiple files provided
            either as an explicit sequence:
                $ mk_sam_file.py P1.csv P2.csv P3.csv
            or as a glob pattern:
                $ mk_sam_file.py P*.csv
            """)
    parser = argparse.ArgumentParser(prog="mk_sam_file.py", add_help=False,
                                     description=prog_desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('csv_file', nargs='*',
                        help=""".csv file(s). Required unless the --all option
                        is given""")
    config_opts = parser.add_argument_group(title="additional options")
    config_opts.add_argument('-h', '--help', action='help',
                             help='show this help message and exit')
    config_opts.add_argument('-a', '--all', action='store_true',
                             help="""Create .sam header files for all CSV files
                             in the current directory.""")
    # config_out = parser.add_mutually_exclusive_group()
    config_opts.add_argument('--dirout', dest='output_directory', metavar='',
                             help="""Output directory. If path does not exist,
                             it will be created. To automatically configure
                             output directories (based on csv name), use
                             --auto-dirs.""")
    config_opts.add_argument('--auto-dirs', action='store_true',
                             help="""Write contents to a directory with same
                             name as the csv file.""")
    args = vars(parser.parse_args())
    csv_file_list = args.pop('csv_file')

    # if no 'read' option specified...
    if len(csv_file_list) == 0 and not args['all']:
        parser.error("Nothing to read.\nYou must provide file name(s) unless "
                     "running with the --all option.")
    # handle csv_file/--all argument conflicts
    elif len(csv_file_list) != 0 and args['all']:
        parser.error("Ambiguous; use file name(s) or --all (not both)")
    # handle the --all option
    elif len(csv_file_list) == 0 and args['all']:
        # find all csv files in CWD
        csv_file_list = []
        filt = re.compile('.*\.csv')
        # search for csv files in CWD (filtering directories)
        for s in list(filter(os.path.isfile, os.listdir())):
            if filt.match(s):
                csv_file_list.append(s)
    if args['auto_dirs']:
        for fname in csv_file_list:
            main(fname, fname.replace('.csv', ''))
    else:
        for fname in csv_file_list:
            main(fname, args["output_directory"])
