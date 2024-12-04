#!/bin/bash 

while IFS= read -r line
do
    #why do backticks work here and not '$()'
    prefix=`echo $line | pcregrep -o1 '(BT\\d+|scrBT\\d+)' -`
    filter="$(echo $line | pcregrep -o1 '(raw|filtered)' -)"
    filename="${prefix}_${filter}"
    filename+=".h5"
    echo $filename
    gcloud storage cp "${line}" "${filename}"
done < paths.txt
