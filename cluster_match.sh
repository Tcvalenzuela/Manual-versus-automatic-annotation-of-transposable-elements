#!/bin/bash

# cd-hit summary
# Example run: bash cluster_match.sh file.clstr outfile_name


# Check max bin size with
MAX_SIZE=$(($(grep -E "^[0-9]{1,2}" ${1} |cut -f1|sort -un|tail -1) + 1))

echo $MAX_SIZE

# define function to dynamically match REF_ORDER with the good order and increase the corresponding variables
# use with either "NOMATCH", "GOOD", or "BAD": match_order_and_count "$REF_ORDER" "NOMATCH"

match_order_and_count () {
	local reforder="$1"
	local matchtype="$2"
	declare -A order_matchtype=()
	order_matchtype[$1]=$(echo "${reforder}_${matchtype}")
	ORDER_MATCHTYPE=${order_matchtype[$1]}
	declare -A total_matchtype=()
	total_matchtype[$1]=$(echo "TOTAL_${matchtype}")
	TOTAL_MATCHTYPE=${total_matchtype[$1]}
	((${ORDER_MATCHTYPE}++))
	((${TOTAL_MATCHTYPE}++))

}


# initialize variables
N_CLUST=$(( $(grep -c "^>" ${1}) - 1))

TOTAL_SEQ=0

TOTAL_NOMATCH=0
TOTAL_GOOD=0
TOTAL_BAD=0

DNA_NOMATCH=0
DNA_GOOD=0
DNA_BAD=0
MITE_NOMATCH=0
MITE_GOOD=0
MITE_BAD=0
RC_NOMATCH=0
RC_GOOD=0
RC_BAD=0
ClassII_NOMATCH=0
ClassII_GOOD=0
ClassII_BAD=0
LINE_NOMATCH=0
LINE_GOOD=0
LINE_BAD=0
Penelope_NOMATCH=0
Penelope_GOOD=0
Penelope_BAD=0
SINE_NOMATCH=0
SINE_GOOD=0
SINE_BAD=0
LTR_NOMATCH=0
LTR_GOOD=0
LTR_BAD=0
ClassI_NOMATCH=0
ClassI_GOOD=0
ClassI_BAD=0


echo -n "" > ${2}
echo -e "DNA_nomatch\tDNA_goodmatch\tDNA_badmatch\tMITE_nomatch\tMITE_goodmatch\tMITE_badmatch\tRC_nomatch\tRC_goodmatch\tRC_badmatch\tClassII_Others_nomatch\tClassII_Others_goodmatch\tClassII_Others_badmatch\tLINE_nomatch\tLINE_goodmatch\tLINE_badmatch\tPenelope_nomatch\tPenelope_goodmatch\tPenelope_badmatch\tSINE_nomatch\tSINE_goodmatch\tSINE_badmatch\tLTR_nomatch\tLTR_goodmatch\tLTR_badmatch\tClassI_Others_nomatch\tClassI_Others_goodmatch\tClassI_Others_badmatch\tTotal_nomatch\tTotal_goodmatch\tTotal_badmatch\tTotal" >> ${2}


for ((i=0;i<=N_CLUST;i++)); do

	# count the consensus in for the total count
	TOTAL_SEQ=$(($TOTAL_SEQ+1))
	CLST1=">Cluster ${i}"
	CLST2=">Cluster $((${i}+1))"

	# number of clustering sequences including the ref consensus
	CLST_SIZE=$(awk -v clst1="^$CLST1$" -v clst2="^$CLST2$" '$0 ~ clst1{flag=1;next} $0 ~ clst2{flag=0}flag' $1 | wc -l)

	# order of the ref consensus
	REF_ORDER=$(awk -v clst1="^$CLST1$" -v clst2="^$CLST2$" '$0 ~ clst1{flag=1;next} $0 ~ clst2{flag=0}flag' $1| grep ^0| cut -d"#" -f2| cut -d"/" -f1 |  cut -d"." -f1)

	if [[ "$CLST_SIZE" == 1 ]]; then # if no matches for this consensus (only the ref consensus forms the cluster), hence

		match_order_and_count "$REF_ORDER" "NOMATCH" # increase <ORDER>_NOMATCH and TOTAL_NOMATCH variables

	elif [[ "$CLST_SIZE" > 1 ]]; then # if 1 or more matches for this consensus (at least one sequence other than the ref consensus is in the cluster), check whether at least 1 good match is in the cluster

		# tmp list of all matching orders excluding the ref consensus
		awk -v clst1="^$CLST1$" -v clst2="^$CLST2$" '$0 ~ clst1{flag=1;next} $0 ~ clst2{flag=0}flag' $1 | grep -v ^0| cut -d"#" -f2| cut -d"/" -f1| cut -d"." -f1 | sort -u > matching_orders_tmp

		#Checking for some loose matches:
		# if a ClassI consensus corresponds to another ClassI element (either way), the match is good
		if ([[ "$REF_ORDER" == "ClassI" ]] && grep -qw "ClassI\|LINE\|Penelope\|SINE\|LTR" matching_orders_tmp) || ([[ "$REF_ORDER" =~ ^(ClassI|LINE|Penelope|SINE|LTR)$ ]] && grep -wq "ClassI" matching_orders_tmp); then

			match_order_and_count "$REF_ORDER" "GOOD"

		# if a ClassII consensus corresponds to another ClassII element (either way), the match is good
		elif ([[ "$REF_ORDER" == "ClassII" ]] && grep -qw "ClassII\|DNA\|MITE\|RC" matching_orders_tmp) || ([[ "$REF_ORDER" =~ ^(ClassII|DNA|MITE|RC)$ ]] && grep -wq "ClassII" matching_orders_tmp); then

			match_order_and_count "$REF_ORDER" "GOOD"

		# if a Penelope consensus corresponds to a LINE consensus (either way), the match is good (no Penelope is detected by MCHelper)
		elif [[ "$REF_ORDER" =~ ^(Penelope|LINE)$ ]] && grep -qw "LINE\|Penelope" matching_orders_tmp; then

			match_order_and_count "$REF_ORDER" "GOOD"

		# in all other cases, a match is good if there is exact correspondance, hence
		elif grep -qw $REF_ORDER matching_orders_tmp; then # if the consensus has at least one match with the same exact order,

			match_order_and_count "$REF_ORDER" "GOOD" # increase <ORDER>_GOOD and TOTAL_GOOD variables

		else # otherwise, all other cases are no good matches, hence

			match_order_and_count "$REF_ORDER" "BAD" # increase <ORDER>_BAD and TOTAL_BAD variables

		fi

	fi

done

echo -e "$DNA_NOMATCH\t$DNA_GOOD\t$DNA_BAD\t$MITE_NOMATCH\t$MITE_GOOD\t$MITE_BAD\t$RC_NOMATCH\t$RC_GOOD\t$RC_BAD\t$ClassII_NOMATCH\t$ClassII_GOOD\t$ClassII_BAD\t$LINE_NOMATCH\t$LINE_GOOD\t$LINE_BAD\t$Penelope_NOMATCH\t$Penelope_GOOD\t$Penelope_BAD\t$SINE_NOMATCH\t$SINE_GOOD\t$SINE_BAD\t$LTR_NOMATCH\t$LTR_GOOD\t$LTR_BAD\t$ClassI_NOMATCH\t$ClassI_GOOD\t$ClassI_BAD\t$TOTAL_NOMATCH\t$TOTAL_GOOD\t$TOTAL_BAD\t$TOTAL_SEQ" >> ${2}
