text = open('Z52/Z52.5a')
text.readline()
attributes = text.readline()
attributes = attributes.split()

line = 'mary had a little lamb'
runs = ''
while line != '':
    line = text.readline()
    line = line.replace('\n','') + line[line.index('\n'):]
    runs += line + ';'

print(runs)
