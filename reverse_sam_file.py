import pandas as pd
import sys

sam_file = sys.argv[1]
sam_text = open(sam_file)
csv_file = ''

name = sam_text.readline().strip('\r\n')
sam_attributes = sam_text.readline().split()

line = 'mary had a little lamb'
samples = []
while line != '':
    line = sam_text.readline()
    samples.append(line.strip('\r\n'))
samples.remove('')

csv_file += ',site_info' + ','*(len(sys.argv)-3) + '\r\n'
csv_file += 'name,' + name + ','*(len(sys.argv)-3) + '\r\n'
csv_file += 'lat,' + sam_attributes[0] + ','*(len(sys.argv)-3) + '\r\n'
csv_file += 'long,' + sam_attributes[1] + ','*(len(sys.argv)-3) + '\r\n'
csv_file += 'mag_dec,' +  sam_attributes[2] + ','*(len(sys.argv)-3) + '\r\n'
csv_file += 'sample_name'

for sample in samples:
    if '.' in sample:
        sample = '.' + sample.split('.')[1]
    elif '-' in sample:
        sample = '-' + sample.split('-')[1]
    csv_file += ',' + sample

csv_file += '\r\nlocal_id'

for sample in samples:
    if '.' in sample:
        sample = sample.split('.')[0]
    elif '-' in sample:
        sample = sample.split('-')[0]
    csv_file += ',' + sample   

sample_dict = {'comment': [], 'strat_level': [],'core_strike': [],'core_dip': [],'bedding_strike': [],'bedding_dip': [],'mass': [], 'runs': []}

if len(sys.argv) > 2:
    for sample_file in sys.argv[2:]:
        text = open(sample_file)
        sample_name_comment = text.readline().strip('\r\n').split()
        attributes = text.readline().strip('\r\n').split()

        assert (sample_name_comment[0]+sample_name_comment[1] in samples),'Sample not in .sam file wrong same file used'

        if len(sample_name_comment) > 3:
            for word in sample_name_comment[3:]:
                comment += word + ' '
            sample_dict['comment'].append(comment)
        else:
            sample_dict['comment'].append('')

        sample_dict['strat_level'].append(attributes[0])
        sample_dict['core_strike'].append(attributes[1])
        sample_dict['core_dip'].append(attributes[2])
        sample_dict['bedding_strike'].append(attributes[3])
        sample_dict['bedding_dip'].append(attributes[4])
        sample_dict['mass'].append(attributes[5])

        line = 'mary had a little lamb'
        runs = ''
        while line != '':
            line = text.readline()
            line = line.replace('\r\n',';')
            runs += line
        sample_dict['runs'].append(runs)

keys = ['comment','strat_level','core_strike','core_dip','bedding_strike','bedding_dip','mass','runs']
csv_file += '\r\n' + reduce(lambda x,y: x + '\r\n' + y,[key + ',' + reduce(lambda x,y: x + ',' + y,sample_dict[key]) for key in keys])


csv_file_final = open(sam_file.replace('.sam','.csv'), 'w')
csv_file_final.write(csv_file)
csv_file_final.close()
