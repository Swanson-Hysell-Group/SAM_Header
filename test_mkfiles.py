import sys

args = sys.argv

assert (len(args) == 3),'Bad number of arguments please input 2 file names'

filename1 = args[1]
filename2 = args[2]

file1 = open(filename1,'r')
file2 = open(filename2,'r')

line1 = 'mary had'
line2 = 'a little lamb'
bools = []
false_lines = []
while line1 and line2:
    line1 = file1.readline()
    line2 = file2.readline()
    bools.append(line1 == line2)
    if line1 != line2:
        false_lines.append(len(bools))
#    print(repr(line1))
#    print(repr(line2))
#    print(line1 == line2)
print('all lines equal: ' + str(all(bools)))
if all(bools) == False:
    print('messed up lines: ' + str(false_lines))

