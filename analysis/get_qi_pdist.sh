#!/bin/bash
#Script to create probability distributions of Qi_diff metric
#list is a set of paths for each landscape used to calculate the distribution

analysis_dir="/netapp/sali/peterc/Undrugabble/AllosModRuns/analysis"
LIST=list
BDIR=${analysis_dir}/bin
NLANDSCAPE=`awk '(NF>0){a+=1}END{print a}' ${LIST}`

PWD0=`pwd`

${BDIR}/get_p.sh

if test -e tempqip47; then rm tempqip47; fi
if test -e tempqip_p; then rm tempqip_p; fi
for landscape in `cat ${LIST}`; do
    for dir in `ls -1d ${landscape}/*`; do
        if test -d ${dir}; then
	    if (test ! -e ${dir}/pm.pdb.B10010001.pdb); then continue; fi
	    NRES=`awk '(NR==1){print NF}' ${dir}/qidiffoftsc.dat`
	    
	    awk '{print $1/'${NLANDSCAPE}'}' ${dir}/p.dat >tempqip_p
	    paste -d" " tempqip_p ${dir}/qidiffoftsc.dat | awk '{printf $1" "}{for (a=2;a<='${NRES}'+1;a++) {printf $a" "}}{printf "\n"}' >>tempqip47
	fi
    done
done


MMa=10.0
TOTPNTS=`awk 'END{print NR}' tempqip47`
MINPNTS=`echo "${TOTPNTS}/1000" |bc -l` #0.1%
for r in `echo ${NRES} | awk '{for(a=1;a<=$1;a++){print a}}'`; do
    rp1=$((${r}+1))
    java -classpath ${BDIR}/ gen_Pbinweight tempqip47 $rp1,1 m${MMa},${MMa} 100 | awk '($3>'${MINPNTS}'){print $0}' >qidiff${r}sc.dat
done

rm tempqip47 tempqip_p
