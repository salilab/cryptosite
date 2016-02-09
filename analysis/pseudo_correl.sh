#!/bin/bash
# in Qi diff, NaN signifies delta Q>0.99 and is changed to 999.0 => those residues are set so P=0.5
# p refers to a state/substate above 0, n refers to a state/substate below 0
# also, if there does not exist an i in state p then P(j is n | i is p)=0.5

analysis_dir="/netapp/sali/peterc/Undrugabble/AllosModRuns/analysis"
transval=0.0
LIST=list
BDIR=${analysis_dir}/bin

${BDIR}/get_p.sh
NLANDSCAPE=`awk '(NF>0){a+=1}END{print a}' ${LIST}`

if test -e tempcij881; then rm tempcij881; fi
for landscape in `cat ${LIST}`; do
    for dir in `ls -1d ${landscape}/*`; do
	if test -d ${dir}; then
	if (test ! -e ${dir}/pm.pdb.B10010001.pdb); then continue; fi
	sres=1
	eres=`awk '(NR==1){print NF}' ${dir}/qidiffoftsc.dat`

	awk '{print $1/'${NLANDSCAPE}'}' ${dir}/p.dat >>tempcijp
	awk '{for (a='${sres}'; a<='${eres}' ;a++) {if($a!="NaN"){printf $a" "}else{printf 999.0" "}}}{printf "\n"}' ${dir}/qidiffoftsc.dat >>tempcij881
	fi
    done
done

paste -d" " tempcijp tempcij881 | ${BDIR}/randomize_list_2.pl | tail -n70000 >tempcij880

awk '{for (a=2; a<=NF ;a++) {print $a}}' tempcij880 >tempcij88a

ZNEW=`awk '{a+=$1}END{print a}' tempcij880`
awk '{print $1/'${ZNEW}'}' tempcij880 >tempcijp
NPOINTS=`awk 'END{print NR}' tempcijp`
numres=`echo "${eres} - ${sres} + 1" | bc -l`

echo tempcij88a >datlist
echo $NPOINTS $numres >>datlist
echo $transval >>datlist
echo tempcijp >>datlist

${BDIR}/correl_ij
awk '{printf $1" "}(NR%'${numres}'==0){printf "\n"}' pjisn_iisn.out >tempcij991
awk '{printf $1" "}(NR%'${numres}'==0){printf "\n"}' pjisn_iisp.out >tempcij992
#take logs odd ratio but set minimum P to 0.01
paste pjisn_iisn.out pjisn_iisp.out | awk '{a=$1;b=$2}(a<0.01){a=0.01}(b<0.01){b=0.01}{printf log(a/b)" "}(NR%'${numres}'==0){printf "\n"}' >pjisn_iisn_o_pjisn_iisp.out
paste pjisn_iisn.out pjisn_iisp.out | awk '{a=1-$1;b=1-$2}(a<0.01){a=0.01}(b<0.01){b=0.01}{printf log(b/a)" "}(NR%'${numres}'==0){printf "\n"}' >pjisp_iisp_o_pjisp_iisn.out
mv tempcij991 pjisn_iisn.out
mv tempcij992 pjisn_iisp.out

rm tempcij88[a0-9] datlist tempcijp 

