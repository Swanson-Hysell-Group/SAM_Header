import pandas as pd
import sys

##Takes formated CSV and creates and writes a .sam header file and a set of sample files##
##for any number of samples.                                                            ##


#################READ IN FILES####################

#fetching comand line data
file_name = sys.argv[1]

#file read in
hdf = pd.read_csv(file_name,header=0,index_col=0,nrows=4)
df = pd.read_csv(file_name,header=5,index_col=0)

#variable assignments
samples = df.keys()
attributes = df[samples[0]].keys()[2:-1]
site_values = hdf['site_info'].keys()[1:]

##########CREATE .SAM HEADER FILE##################

#setting name
sam_header = hdf['site_info']['name'] + '\r\n'

#creating long lat and dec info
i=1
for value in site_values:
    sam_header += ' ' + hdf['site_info'][value] + ' '*i
    i += 1
sam_header = sam_header[0:-i+1]
sam_header += '\r\n'

#making writing sample info
for sample in samples:
    sam_header += df[sample]['local_id'] + sample + '\r\n'

#creating and writing file
sam_file = open(df[sample]['local_id'] + '.sam', 'w')
sam_file.write(sam_header)
sam_file.close()

################Create Sample Files#################

for sample in samples:

    #assign variables for easy refrence
    local_id = df[sample]['local_id']
    comment = df[sample]['comment']
    if type(df[sample]['runs']) == str:
        runs = df[sample]['runs'].split(';')

    #check for no comment
    if type(comment) != str:
        comment = ''

    #insure input is valid
    assert (len(local_id) <= 4),'Locality ID excedes 4 characters: refer too http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
    assert (len(comment) <= 255),'Sample comment excedes 255 characters: refer too http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
    assert (len(sample) <= 9),'Sample name excedes 9 chaaracters: refer too http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
    
    #write sample name and comment for sample file
    new_file =  local_id + ' '*(4-len(local_id)) + sample + ' '*(9-len(sample)) + comment + '\r\n '

    #write in sample attributes on the second line
    for attribute in attributes:
        if type(df[sample][attribute]) != str:
            df[sample][attribute] = ''
        if df[sample][attribute].isdigit():
            df[sample][attribute] = str(float(df[sample][attribute]))

        #attributes must follow standard sam format
        assert (len(df[sample][attribute]) <= 5),'Length of ' + attribute + ' excedes 5 characters: refer too http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'

        new_file += ' ' + df[sample][attribute] + ' '*(5-len(df[sample][attribute]))

    new_file += '\r\n'
    
    #if there are previous sample runs write that to the bottem of the file
    for run in runs:
        new_file += run + '\r\n'
    
    #create and write sample file
    sample_file = open(local_id + sample, 'w')
    sample_file.write(new_file)
    sample_file.close()
