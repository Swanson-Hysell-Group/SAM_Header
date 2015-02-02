import pandas as pd
import sys

sam_file = sys.argv[1]
sam_text = open(sam_file)
csv_file = ''

name = sam_text.readline()
sam_attributes = sam_text.readline().split()

line = 'mary had a little lamb'
samples = []
while line != '':
    line = sam_text.readline()
    samples.append(line)

csv_file += ',site_info' + ','*(len(sys.argv)-3) + '\r\n'
csv_file += 'name/comment,' + name + ','*(len(sys.argv)-3) + '\r\n'
csv_file += 'lat,' + sam_attributes[0] + ','*(len(sys.argv)-3) + '\r\n'
csv_file += 'long,' + sam_attributes[1] + ','*(len(sys.argv)-3) + '\r\n'
csv_file += 'mag_dec,' +  sam_attributes[2] + ','*(len(sys.argv)-3) + '\r\n'
csv_file += 'sample_name,'

for sample in samples:
    if '.' in sample:
        sample = sample.split('.')[1]
    elif '-' in sample:
        sample = sample.split('-')[1]
    csv_file += sample + ','

csv_file += '\r\nlocal_id,'

for sample in samples:
    if '.' in sample:
        sample = sample.split('.')[0]
    elif '-' in sample:
        sample = sample.split('-')[0]
    csv_file += sample + ','    

csv_file += '\r\ncomment,'
comment = 'comment,'

if len(sys.argv) > 2:
    for sample_file in sys.argv[2:]:
        text = open(sample_file)
        sample_name_comment = text.readline().split
        attributes = text.readline().split()

        assert (sample_name in samples),'Sample not in .sam file wrong same file used'

        if len(sample_name_comment) > 3:
            csv_parts = csv_file.partition(comment)
            for word in sample_name_comment[3:]:
                comment += word + ' '
            csv_parts[1] += comment + ','
            csv_file = sum(csv_parts)
        else:
            csv_file += comment + ','

        line = 'mary had a little lamb'
        runs = ''
        while line != '':
            line = text.readline()
            line = line.replace(' \r\n',';')
            runs += line



csv_file_final = open(sam_file.replace('.sam','.csv'), 'w')
csv_file_final.write(csv_file)
csv_file_final.close()
