import sys

args = sys.argv

assert (len(args) == 3),'Bad number of arguments please input 2 file names'

filename1 = args[1]
filename2 = args[2]

file1 = open(filename1,'r')
file2 = open(filename2,'r')

line1 = 'mary had'
line2 = 'a little lamb'
while line1 and line2:
    line1 = file1.readline()
    line2 = file2.readline()
    print(repr(line1))
    print(repr(line2))
    print(line1 == line2)

