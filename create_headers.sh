SITE=SLB
path=SLB
for n in 01 02 03 04 05 06 07 08 09 10
do python mk_sam_file.py "$path"/"$SITE""$n".csv "$path"/"$SITE""$n"/
done
