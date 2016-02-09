#!/bin/bash
#Script generates a rolling average of local Qi ($FIL_qi) as a function of global Q ($FIL_q)
#Predicts relative folding stability of regions in the protein under the assumption that folding stability is proportional to the rate of local unfolding
#list is a set of paths for each landscape used to calculate the distribution

analysis_dir="/netapp/sali/peterc/Undrugabble/AllosModRuns/analysis"
LIST=list
BDIR=${analysis_dir}/bin
NLANDSCAPE=`awk '(NF>0){a+=1}END{print a}' ${LIST}`
rcut=11

PWD0=`pwd`

${BDIR}/get_p.sh

#if test -e tempqi47; then rm tempqi47; fi
for landscape in `cat ${LIST}`; do
    for dir in `ls -1d ${landscape}/*`; do
        if test -d ${dir}; then
	    #if (test ! -e ${dir}/pm.pdb.B10010001.pdb); then continue; fi
	    FIL_qi=`ls -1 ${dir}/qioft_*_${rcut}sc.dat | head -n1`
	    FIL_q=`ls -1 ${dir}/qalloft_*_${rcut}sc.dat | head -n1`
	    NRES=`awk '(NR==1){print NF}' $FIL_qi`
	    awk '{print $1/'${NLANDSCAPE}'}' ${dir}/p.dat >temp_p
	    awk '{print $1}' $FIL_q >temp_q
	    paste -d" " temp_p temp_q $FIL_qi >>tempqi47
	fi
    done
done

TOTPNTS=`awk '(NF>1){a+=1}END{print a}' tempqi47`
MINPNTS=`echo "${TOTPNTS}/1000" |bc -l` #0.1%
for r in `echo ${NRES} | awk '{for(a=1;a<=$1;a++){print a}}'`; do
    rp2=$((${r}+2))
    java -classpath ${BDIR}/ gen_Abinweight tempqi47 2,${rp2},1 0,1 25 | awk '($4>'${MINPNTS}'){print $0}' >propq${r}.dat
done
    
rm tempqi47 temp_p temp_q

