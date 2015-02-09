for k in sam 1a 2a 5a 6a 7a 8a
do python test_mkfiles.py Z52/Z52."$k" Z52_origional/Z52."$k"
done

for k in .sam 337.4b 337.7b 338.3b 338.5b 340.2b 340.3b 390.2b 390.3b 390.9b 391.5b 391.51b 391b 392.3b D1b D2b D3b D4b D6b D7b D8b D9b
do python test_mkfiles.py SI9/SI9-"$k" SI9_origional/SI9-"$k"
done
