#!/bin/bash

QEXEC=${1}/bin/getqiavg_sc
PWD0=`pwd`
rcut=11
LIST="list"

for landscape in `cat ${LIST}`; do
    for dir in `ls -1d ${landscape}/*`; do
    if test -d ${dir}; then
    cd $dir

    #if (test ! -e pm.pdb.B10010001.pdb); then continue; fi

    ls -1 pm.pdb.B[1-8]*pdb >templlist
    TEMP1=`awk '(NR==1){print "pm_"$1}' list`
    TEMP2=`awk 'BEGIN{a="none"}(NR==2){a="pm_"$1}END{print a}' list`
    NFILES=`awk 'END{print NR}' templlist`
    NRES=`awk 'BEGIN{a=0}($1=="ATOM"&&$3=="CA"){a+=1}END{print a}' ${TEMP1}`

    mkdir qi_${TEMP2}_${rcut} qi_${TEMP1}_${rcut}

    # Qi_${TEMP1}
    echo ${TEMP1} >targlist
    echo ${NRES} >>targlist
    echo $rcut >>targlist
    echo $NFILES >>targlist
    awk '{print $1}' templlist >>targlist
    $QEXEC
    mv qi_*dat qi_${TEMP1}_${rcut}

    # Qi_${TEMP2}
    if test $TEMP2 != "none"; then
    echo ${TEMP2} >targlist
    echo ${NRES} >>targlist
    echo $rcut >>targlist
    echo $NFILES >>targlist
    awk '{print $1}' templlist >>targlist
    $QEXEC
    mv qi_*dat qi_${TEMP2}_${rcut}
    fi

    # make Qi(t)
    if test -e qi_${TEMP1}_${rcut}/qi_sc1.dat;    then ls -1 qi_${TEMP1}_${rcut}/qi_sc[1-9].dat >templlist; fi
    if test -e qi_${TEMP1}_${rcut}/qi_sc10.dat;   then ls -1 qi_${TEMP1}_${rcut}/qi_sc[1-9][0-9].dat >>templlist; fi
    if test -e qi_${TEMP1}_${rcut}/qi_sc100.dat;  then ls -1 qi_${TEMP1}_${rcut}/qi_sc[1-9][0-9][0-9].dat >>templlist; fi
    if test -e qi_${TEMP1}_${rcut}/qi_sc1000.dat; then ls -1 qi_${TEMP1}_${rcut}/qi_sc[1-9][0-9][0-9][0-9].dat >>templlist; fi
    if test -e qioft_${TEMP1}_${rcut}sc.dat; then rm qioft_${TEMP1}_${rcut}sc.dat; fi
    for f in `cat templlist`; do
        awk '{printf $2" "}END{printf "\n"}' $f >>qioft_${TEMP1}_${rcut}sc.dat
    done
    awk '{q=0}{for(a=1;a<=NF;a++){q+=$a;if(a==NF){printf "%9.4f\n",q/NF}}}' qioft_${TEMP1}_${rcut}sc.dat >qalloft_${TEMP1}_${rcut}sc.dat

    if test $TEMP2 != "none"; then
    if test -e qi_${TEMP2}_${rcut}/qi_sc1.dat;    then ls -1 qi_${TEMP2}_${rcut}/qi_sc[1-9].dat >templlist; fi
    if test -e qi_${TEMP2}_${rcut}/qi_sc10.dat;   then ls -1 qi_${TEMP2}_${rcut}/qi_sc[1-9][0-9].dat >>templlist; fi
    if test -e qi_${TEMP2}_${rcut}/qi_sc100.dat;  then ls -1 qi_${TEMP2}_${rcut}/qi_sc[1-9][0-9][0-9].dat >>templlist; fi
    if test -e qi_${TEMP2}_${rcut}/qi_sc1000.dat; then ls -1 qi_${TEMP2}_${rcut}/qi_sc[1-9][0-9][0-9][0-9].dat >>templlist; fi
    if test -e qioft_${TEMP2}_${rcut}sc.dat; then rm qioft_${TEMP2}_${rcut}sc.dat; fi
    for f in `cat templlist`; do
        awk '{printf $2" "}END{printf "\n"}' $f >>qioft_${TEMP2}_${rcut}sc.dat
    done
    awk '{q=0}{for(a=1;a<=NF;a++){q+=$a;if(a==NF){printf "%9.4f\n",q/NF}}}' qioft_${TEMP2}_${rcut}sc.dat >qalloft_${TEMP2}_${rcut}sc.dat
    fi

    # Qi between allosteric states
    if test $TEMP2 != "none"; then
    echo ${TEMP1} >targlist
    echo ${NRES} >>targlist
    echo $rcut >>targlist
    echo 1 >>targlist
    echo ${TEMP2} >>targlist
    $QEXEC
    mv qi_sc1.dat qiassc.dat

    # Qidiff
    NLINES=`awk 'END{print NR}' qioft_${TEMP1}_${rcut}sc.dat`
    awk '{printf "%9.4f ",1.0-$2}END{printf "\n"}' qiassc.dat | awk '{for (a=1; a<='${NLINES}' ;a++) {print $0}}' >tempqidiffsc

    paste qioft_${TEMP2}_${rcut}sc.dat qioft_${TEMP1}_${rcut}sc.dat tempqidiffsc | \
        awk '{for (a=1; a<='${NRES}' ;a++) { {b='${NRES}'+a; c='${NRES}'*2+a} if($c>0.01){printf "%9.4f ",($b-$a)/$c} else {printf "NaN "} }}{printf "\n"}' >qidiffoftsc.dat
    fi

    rm -rf templlist qi_${TEMP2}_${rcut} qi_${TEMP1}_${rcut} qi_avg.dat qi_scavg.dat tempqidiffsc targlist

    cd ${PWD0}
    fi
    done
done
