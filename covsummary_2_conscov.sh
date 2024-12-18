#!/bin/bash

# Generate a consensus-wise coverage file from the output of buildSummary.pl (RepeatMasker script in RepeatMasker/utils/): copy count, bp coverage, consensus id and length
# output file contains 5 tab-separated columns: consensus original id, element type, consensus length, copy count, coverage in bp

# To run this script you need to run beforehand RepeatMasker and the util script buildSummary.pl to obtain the coverage summary of all TE families
# seqkit must be in $PATH
# Required arguments: coverage summary file, library fasta file, name of the output file
# Example run: bash covsummary_2_conscov.sh aedes_manuallib_RM_covsummary.tsv aedes_manuallib_edit.fa outfile_name


echo -n "" > $3

echo -e "consensus_id\telement_type\tconsensus_lenght\tcopy_count\tcoverage_bp\tcoverage_perc" >> $3

awk '/Repeat Stats/{flag=1;next}/By Sequence/{flag=0}flag' $1 |  grep -v "([A-T]*)n\|-rich"|tac | tail -n +3 | tac | sed '1,5d' > tmpfile

for SEQID in $(awk '{print$1}' tmpfile); do

	TETYPE=$(seqkit grep -r -p ${SEQID}# $2| seqkit fx2tab --length --name | cut -f1 | cut -d"#" -f2)
	CONSLEN=$(seqkit grep -r -p ${SEQID}# $2| seqkit fx2tab --length --name | cut -f2)
	COUNT=$(grep "${SEQID}\s" tmpfile| awk '{print$2}')
	COVBP=$(grep "${SEQID}\s" tmpfile| awk '{print$3}')
	COVPERC=$(grep "${SEQID}\s" tmpfile| awk '{print$4}' | sed 's/%//g')

echo -e "$SEQID\t$TETYPE\t$CONSLEN\t$COUNT\t$COVBP\t$COVPERC" >> $3

done

rm tmpfile

