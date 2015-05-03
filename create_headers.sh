SITE=SLB
path=SLB
l="01 01b 02 02b 03 03b 04 04b 05 05b 06 06b 07 07b 08 08b 09 09b 10 10b 11 11b 12 12b 13 13b 14 14b 15 15b 16 16b 17 17b 18 18b 19 19b 20 20b 21 21b 22 22b 23 23b 24 24b 25 25b 26 26b 27 27b 28 28b 29 29b 30 30b 31 31b 32 32b 33 33b 34 34b 35 35b 36 36b 37 37b 38 38b 39 39b 40 40b"

for n in `echo $l`
do python mk_sam_file.py "$path"/"$SITE""$n".csv "$path"/"$SITE""$n"/
done
