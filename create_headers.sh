SITE=SLB
path=SLB
l="01 02 03 04 05 06 07 08 09 10 11 12 13 14 15"

for n in `echo $l`
do python mk_sam_file.py "$path"/"$SITE""$n".csv "$path"/"$SITE""$n"/
done
