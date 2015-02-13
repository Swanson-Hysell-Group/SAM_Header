import sys

## Tests if two files have the same text characters                                     ##

args = sys.argv

assert (len(args) == 3),'Bad number of arguments please input 2 file names'

filename1 = args[1]
filename2 = args[2]

file1 = open(filename1,'r')
file2 = open(filename2,'r')

line1 = 'mary had'
line2 = 'a little lamb'
list1 = []
list2 = []
bools = []
false_lines = []
while line1 and line2:
    line1 = file1.readline()
    line2 = file2.readline()
    bools.append(line1 == line2)
    if line1 != line2:
        false_lines.append(len(bools))
    list1.append(line1)
    list2.append(line2)
#    print(repr(line1))
#    print(repr(line2))
#    print(line1 == line2)
print('all lines same: ' + str(all(bools)))
if all(bools) == False:
    print('check lines: ' + str(false_lines))
for line in false_lines:
    print(filename1 + ':\n' + list1[line-1])
    print(filename2 + ':\n' + list2[line-1])
