#!/usr/bin/env sh

# Convert file containing only an markdown table to csv file
# Note: Extra Spaces are removed, so column headers will be concatenated
#
# usage: ./md2csv.sh experiment.md -> experiment.csv

in=$1
out=$(basename $1 .md).csv

#sed -e '/^[|]-/d' -e 's/^|//' -e 's/|$//' -e 's/,/;/g'  -e 's/|/,/g'  -e 's/ //g'  $in > $out
sed -e '/^[|]-/d' -e 's/^|//' -e 's/|$//' -e 's/,/;/g'  -e 's/|/,/g'  -e 's/ //g' -e 's/APOBEC/APO/g' -e 's/ActiveMotif//g' $in > $out

echo Required Headers:
echo 'SeqID', 'SampleID', 'Content (->[csv experiment]->)'
echo
echo Actual Headers:
head -1 $out
echo
echo wrote $out
