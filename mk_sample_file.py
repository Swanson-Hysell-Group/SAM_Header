import pandas as pd
import sys

file_name = sys.argv[1]

df = pd.read_csv(file_name,header=0,index_col=0)
samples = df.keys()
attributes = df[samples[0]].keys()[2:-1]

for sample in samples:

    local_id = df[sample]['local_id']
    comment = df[sample]['comment']
    if type(df[sample]['runs']) == str:
        runs = df[sample]['runs'].split(';')

    if type(comment) != str:
        comment = ''

    assert (len(local_id) <= 4),'Locality ID excedes 4 characters: refer too http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
    assert (len(comment) <= 255),'Sample comment excedes 255 characters: refer too http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
    assert (len(sample) <= 9),'Sample name excedes 9 chaaracters: refer too http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
    

    new_file =  local_id + ' '*(4-len(local_id)) + sample + ' '*(9-len(sample)) + comment + '\r\n '

    for attribute in attributes:
        if type(df[sample][attribute]) != str:
            df[sample][attribute] = ''
        if df[sample][attribute].isdigit():
            df[sample][attribute] = str(float(df[sample][attribute]))

        assert (len(df[sample][attribute]) <= 5),'Length of ' + attribute + ' excedes 5 characters: refer too http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'

        new_file += ' ' + df[sample][attribute] + ' '*(5-len(df[sample][attribute]))

    new_file += '\r\n'
    
    for run in runs:
        new_file += run + '\r\n'
    
    filename = open(local_id + sample, 'w')
    filename.write(new_file)
    filename.close()
