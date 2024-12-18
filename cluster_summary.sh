#!/bin/bash

# cd-hit summary
# Example run: bash cluster_summary.sh file.clstr outfile_name


# Check max bin size with
MAX_SIZE=$(($(grep -E "^[0-9]{1,2}" ${1} |cut -f1|sort -un|tail -1) + 1))

echo $MAX_SIZE

# Use minimum this value as max bin size

N_CLUST=$(( $(grep -c "^>" ${1}) - 1))

TOTAL_SEQ=0
BIN_1=0
BIN_2=0
BIN_3=0
BIN_4=0
BIN_5=0
BIN_6=0
BIN_7=0
BIN_8=0
BIN_9=0
BIN_10=0
BIN_11=0

echo -n "" > ${2}
echo -e "Bin_size_1\tBin_size_2\tBin_size_3\tBin_size_4\tBin_size_5\tBin_size_6\tBin_size_7\tBin_size_8\tBin_size_9\tBin_size_10\tBin_size_11\tTotal_sequences" >> ${2}

for ((i=0;i<=N_CLUST;i++)); do

	CLST1=">Cluster ${i}"
	CLST2=">Cluster $((${i}+1))"

	CLST_SIZE=$(awk -v clst1="^$CLST1$" -v clst2="^$CLST2$" '$0 ~ clst1{flag=1;next} $0 ~ clst2{flag=0}flag' $1 | wc -l)

	TOTAL_SEQ=$(($TOTAL_SEQ+1))

	if [[ "$CLST_SIZE" == 1 ]]; then

		BIN_1=$(($BIN_1+1))

	elif [[ "$CLST_SIZE" == 2 ]]; then

		BIN_2=$(($BIN_2+1))

	elif [[ "$CLST_SIZE" == 3 ]]; then

		BIN_3=$(($BIN_3+1))

	elif [[ "$CLST_SIZE" == 4 ]]; then

		BIN_4=$(($BIN_4+1))

	elif [[ "$CLST_SIZE" == 5 ]]; then

		BIN_5=$(($BIN_5+1))

	elif [[ "$CLST_SIZE" == 6 ]]; then

		BIN_6=$(($BIN_6+1))

	elif [[ "$CLST_SIZE" == 7 ]]; then

		BIN_7=$(($BIN_7+1))

	elif [[ "$CLST_SIZE" == 8 ]]; then

                BIN_8=$(($BIN_8+1))

	elif [[ "$CLST_SIZE" == 9 ]]; then

                BIN_9=$(($BIN_9+1))

	elif [[ "$CLST_SIZE" == 10 ]]; then

                BIN_10=$(($BIN_10+1))

	elif [[ "$CLST_SIZE" == 11 ]]; then

                BIN_11=$(($BIN_11+1))

	else

		echo "Sequences at ${CLST1} out of range!\n\nThey are not counted in the defined clustering bins"

	fi

done

echo -e "$BIN_1\t$BIN_2\t$BIN_3\t$BIN_4\t$BIN_5\t$BIN_6\t$BIN_7\t$BIN_8\t$BIN_9\t$BIN_10\t$BIN_11\t$TOTAL_SEQ" >> ${2}
