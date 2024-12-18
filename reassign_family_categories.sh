#!/bin/bash

# Define list of perfect, good, present, fragmented, absent consensi in one library (query) compared to the other (reference)
# run this script after running the reciprocal masking and the get_family_summary_paper.sh script from Flynn et al. (2020). Supply the reference TE library as argument (the one used with the -lib option in the reciprocal masking), a suffix for reference library, a suffix for species name

# example run: bash reassign_family_categories.sh reference_TElib.fa manualref species


grep -Fvw -f good.perfect.families present.all.families > present.only.families
grep -Fvw -f present.all.families repbase.families > fragment.families 

for fam in $(cat fragment.families); do awk -v fam="$fam" '$1==fam && $3>0.30 {print $1}' repbase.families.covsumm80 >> fragment_80_30.families; done # query families with at least 80% similarity over more than 30% length of the reference sequence

grep -Fwv -f fragment_80_30.families fragment.families > absent1.families
sed 's/$/#/g' repbase.families | sed 's/^/>/g'> frag_pres.fam

for seqid in $(cat fragment_80_30.families); do grep -m1 $seqid'#' $1 | cut -d'#' -f2 >> fragment_80_30_classif; done
for seqid in $(cat absent1.families); do grep -m1 $seqid'#' $1 | cut -d'#' -f2 >> absent1_classif; done
grep ">" $1 | grep -Fv -f frag_pres.fam | sed 's/>//g' | cut -d'#' -f1 > absent2.families
grep ">" $1 | grep -Fv -f frag_pres.fam | sed 's/>//g' | cut -d'#' -f2 > absent2_classif
for seqid in $(cat present.only.families); do grep -m1 $seqid'#' $1 | cut -d'#' -f2 >> present.only_classif; done
for seqid in $(cat good.families); do grep -m1 $seqid'#' $1 | cut -d'#' -f2 >> good_classif; done
for seqid in $(cat perfect.families); do grep -m1 $seqid'#' $1 | cut -d'#' -f2 >> perfect_classif; done


paste absent1.families absent1_classif > absent1
paste absent2.families absent2_classif > absent2
cat absent1 absent2 > ${3}_${2}_absent.tsv
paste fragment_80_30.families fragment_80_30_classif > ${3}_${2}_fragment_80_30.tsv
paste present.only.families present.only_classif > ${3}_${2}_present.tsv
paste good.families good_classif > ${3}_${2}_good.tsv
paste perfect.families perfect_classif > ${3}_${2}_perfect.tsv

category_summary() {
echo -n > tmpcol
catcount=$($1 | wc -l ) # first arg is the category tsv file
range=$(seq 1 $catcount)
for i in $range ; do printf "${2}\n"; done >> tmpcol # second arg is the category type (Perfect, Good, Absent, etc)
paste $1 tmpcol > $3 # third arg is the output filename 
}

category_summary ${3}_${2}_absent.tsv Absent abstsv
category_summary ${3}_${2}_fragment_80_30.tsv Fragment fragtsv
category_summary ${3}_${2}_present.tsv Present prestsv
category_summary ${3}_${2}_good.tsv Good goodtsv
category_summary ${3}_${2}_perfect.tsv Perfect perftsv

echo -e "consensus_id\telement_type\toverlap_level" > ${3}_${2}_categories_summary.tsv
cat abstsv fragtsv prestsv goodtsv perftsv >> ${3}_${2}_categories_summary.tsv
