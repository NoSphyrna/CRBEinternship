#!/bin/bash

indices="$1"
map="$2"

cat "$map" | while IFS= read -r line; do
	#echo "Processing line : $line"
	if [[ "$line" =~ Kinnex16S_Fwd_([0-9]{2})--Kinnex16S_Rev_([0-9]{2}),([A-Za-z0-9-]+)-ITS9-ITS4\n? ]]; then
		n1=${BASH_REMATCH[1]}
		n2=${BASH_REMATCH[2]}
		sample=${BASH_REMATCH[3]}
		adaptfwd=$(sed "${n1}q;d" "$indices")
		adaptrev=$(sed "${n2}q;d" "$indices" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev)
		echo -e ">$sample\n$adaptfwd...$adaptrev"
	fi
done
