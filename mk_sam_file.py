import pandas as pd
import math,sys,time
from pmag import dosundec as sundec
from ipmag import igrf
from datetime import datetime as dt
from to_year_fraction import *
from functools import reduce

def main():
    """
    NAME
        mk_sam_file.py

    DESCRIPTION
        Takes formated CSV and creates and writes a .sam header file and a set of sample files for any number of samples.

    SYNTAX
        ~/$ python mk_sam_file.py site.csv

    OUTPUT
        .sam and 

    """


    #################READ IN FILES####################

    if "-h" in sys.argv:
        print main.__doc__
        sys.exit()

    #fetching comand line data
    file_name = sys.argv[1]
    directory = reduce(lambda x,y: x+'/'+y, [sub for sub in file_name.split('/') if not sub.endswith('.csv')]) + '/'

    #file read in
    hdf = pd.read_csv(file_name,header=0,index_col=0,nrows=4,usecols=[0,1])
    df = pd.read_csv(file_name,header=5,index_col=1,usecols=[0,1,2,3,4,5,13,14,15,16,17,18,19,20,21]).transpose()
    sdf = pd.read_csv(file_name,header=5,index_col=0,usecols=[1,6,7,8,9,10,11,12]).transpose()

    #variable assignments
    samples = df.keys()
    attributes = ['core_strike','core_dip','bedding_strike','bedding_dip','mass']
    site_values = ['lat','long']
    time_types = ['year','month','days','hours','minutes']

    ##########Find Calculated Values##################

    #calculate sun_core_strike for all samples
    for sample in samples:
        time_values = []
        for i in range(len(time_types)):
             time_values.append(str(int(sdf[sample][time_types[i]])))
        assert (len(time_values[0]) == 4),"must input full year for sun compass calculation (i.e. YYYY)"
        sundata = {}
        if float('nan') in sdf[sample]:
            df[sample]['sun_core_strike'] = float('nan')
        else:
            if (len(time_values[1]) == 1):
                time_values[1] = '0' + time_values[1]
            if (len(time_values[2]) == 1):
                time_values[2] = '0' + time_values[2]
            if (len(time_values[3]) == 1):
                time_values[3] = '0' + time_values[3]
            if (len(time_values[4]) == 1):
                time_values[4] = '0' + time_values[4]
            sundata['date'] = reduce(lambda x,y: x + ':' + y, time_values)
            sundata['lat'] = hdf['site_info']['lat']
            sundata['lon'] = hdf['site_info']['long']
            sundata['shadow_angle'] = sdf[sample]['shadow_angle']
            sundata['delta_u'] = sdf[sample]['GMT_offset']
            df[sample]['sun_core_strike'] = round(sundec(sundata),1)

    #calculate IGRF
        if math.isnan(sdf[sample]['year']) or math.isnan(sdf[sample]['month']) or math.isnan(sdf[sample]['days']) or math.isnan(sdf[sample]['hours']) or math.isnan(sdf[sample]['minutes']):
            df[sample]['calculated_IGRF'] = 'insufficient data'
        else:
            if math.isnan(float(hdf['site_info']['elevation'])):
                hdf['site_info']['elevation'] = 0
            date = to_year_fraction(dt(int(sdf[sample]['year']),int(sdf[sample]['month']),int(sdf[sample]['days']),int(sdf[sample]['hours']),int(sdf[sample]['minutes'])))
            df[sample]['calculated_IGRF'] = igrf([date,float(hdf['site_info']['elevation']),float(hdf['site_info']['lat']),float(hdf['site_info']['long'])])
            df[sample]['IGRF_local_dec'] = df[sample]['calculated_IGRF'][0]

    #calculate magnetic declination
        if math.isnan(df[sample]['sun_core_strike']) or math.isnan(df[sample]['magnetic_core_strike']):
            df[sample]['calculated_mag_dec'] = 'insufficient data'
        else:
            df[sample]['calculated_mag_dec'] = df[sample]['sun_core_strike'] - df[sample]['magnetic_core_strike']

    if abs(float(df.transpose()['IGRF_local_dec'].mean()) - float(df.transpose()['calculated_mag_dec'].mean())) > 5:
        print('WARNING: local IGRF declination & calculated magnetic declination where ' + str(abs(round(float(df.transpose()['IGRF_local_dec'].mean()) - float(df.transpose()['calculated_mag_dec'].mean()),2))) + ' degrees different')

    ##########CREATE .SAM HEADER FILE##################

    #setting name
    sam_header = hdf['site_info']['name'] + '\r\n'

    #creating long lat and dec info
    for value in site_values:
        hdf['site_info'][value] = str(round(float(hdf['site_info'][value]),1))
        if value == 'lat':
            sam_header += ' ' + hdf['site_info'][value]
        else:
            sam_header += ' '*(5-len(hdf['site_info'][value]) + 1) + hdf['site_info'][value]
    sam_header += ' '*(5) + '0'
    sam_header += '\r\n'

    #making writing sample info
    for sample in samples:
        sam_header += df[sample]['site_id'] + sample + '\r\n'

    #creating and writing file
    sam_file = open(directory + df[sample]['site_id'] + '.sam', 'w')
    sam_file.write(sam_header)
    sam_file.close()

    ################Create Sample Files#################

    for sample in samples:

        #assign variables for easy refrence
        site_id = df[sample]['site_id']
        comment = df[sample]['comment']
        if type(df[sample]['runs']) == str:
            runs = df[sample]['runs'].split(';')

        #decide which core_strike to use, default is sun_core_strike but if not supplied 
        #magnetic_core_strike will be used
        if math.isnan(df[sample]['sun_core_strike']):
            df[sample]['core_strike'] = df[sample]['magnetic_core_strike'] + (df[sample]['calculated_IGRF'] - 360)
        else:
    ##########################change 'sun_core_strike'##########################################
            df[sample]['core_strike'] = df[sample]['sun_core_strike']

        #check for no comment
        if type(comment) == float and math.isnan(comment):
            comment = ''

        #insure input is valid
        assert (len(site_id) <= 4),'Locality ID excedes 4 characters: refer too http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
        assert (len(comment) <= 255),'Sample comment excedes 255 characters: refer too http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
        assert (len(sample) <= 9),'Sample name excedes 9 chaaracters: refer too http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
        
        #write sample name and comment for sample file
        new_file =  site_id + ' ' + sample + ' ' + comment + '\r\n'

        #start second line strat_level get's special treatment
        df[sample]['strat_level'] = str((df[sample]['strat_level']))
        assert (len(df[sample]['strat_level']) <= 6),'Length of strat_level excedes 6 characters: refer too http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
        new_file += ' ' + ' '*(6-len(df[sample]['strat_level'])) + df[sample]['strat_level']

        #write in sample attributes on the second line
        for attribute in attributes:

            assert (str(df[sample][attribute]).isdigit),str(attributes) + 'must all be numbers'

            if type(df[sample][attribute]) == float and math.isnan(df[sample][attribute]):
                df[sample][attribute] = ''
            else:
                df[sample][attribute] = str(round(float(df[sample][attribute]),1))


            #attributes must follow standard sam format
            assert (len(df[sample][attribute]) <= 5),'Length of ' + attribute + ' excedes 5 characters: refer too http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'

            new_file += ' ' + ' '*(5-len(df[sample][attribute])) + df[sample][attribute]

        new_file += '\r\n'
        
        #if there are previous sample runs write that to the bottem of the file
        for run in runs:
            new_file += run + '\r\n'
        
        #create and write sample file
        new_file = new_file.rstrip('\r\n') + '\r\n'
        sample_file = open(directory + site_id + sample, 'w')
        sample_file.write(new_file)
        sample_file.close()

    ################Write New Values to .csv###################

    csv_file = open(file_name)
    csv_str = ''

    for i in range(4):
        csv_str += csv_file.readline()

    elev_line = csv_file.readline().split(',')
    elev_line[1] = str(hdf['site_info']['elevation'])
    csv_str += reduce(lambda x,y: x + ',' + y, elev_line)

    header = csv_file.readline()
    csv_str += header
    header = header.strip('\n').split(',')

    for sample in samples:
        line = csv_file.readline()
        items = line.split(',')
        for i in range(len(header)):
            if i == 1:
                continue
            elif header[i] == 'calculated_IGRF':
                if type(df[sample][header[i]]) != str:
                    items[i] = str(list(df[sample][header[i]])).replace(',',';')
                else:
                    items[i] = str(df[sample][header[i]])
            elif header[i] in df[sample].keys():
                items[i] = str(df[sample][header[i]])
            elif header[i] in sdf[sample].keys():
                items[i] = str(sdf[sample][header[i]])
            else:
                raise KeyError('there is no item: ' + header[i])
        csv_str += reduce(lambda x,y: x + ',' + y, items) + '\n'

    new_csv_file = open(file_name,'w')
    new_csv_file.write(csv_str)
    new_csv_file.close()

main()
