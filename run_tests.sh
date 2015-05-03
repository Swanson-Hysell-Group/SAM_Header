
SITE=SLB
path=SLB
l="01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"

for n in `echo $l`
    do for k in 1 2 3 4 5 6 7 8 9 10
        do python test_mkfiles.py "$path"/"$SITE""$n"/"$SITE""$n"."$k"a "$path"/"$SITE""$n"b/"$SITE""$n"."$k"b
    done
done

#for k in .sam 337.4b 337.7b 338.3b 338.5b 340.2b 340.3b 390.2b 390.3b 390.9b 391.5b 391.51b 391b 392.3b D1b D2b D3b D4b D6b D7b D8b D9b
#do python test_mkfiles.py tests/SI9/SI9-"$k" tests/SI9_original/SI9-"$k"
#done
