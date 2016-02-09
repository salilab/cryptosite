#!/bin/bash

MAXVEL=0.0010
LIST=list
PWD0=`pwd`

for landscape in `cat ${LIST}`; do
    for dir in `ls -1d ${landscape}/*`; do
    if test -d ${dir}; then
    cd $dir

    if (test ! -e pm.pdb.B10010001.pdb); then 
	echo echo ${dir} "**** WARNING: This directory does not contain a simulation trajectory. ****"
	continue
    fi

    TEMP1=`awk '(NR==1){print "pm_"$1}' ${LIST}`
    QFIL=`ls -1 qalloft_${TEMP1}_*sc.dat | head -n1`

    if test -e ${QFIL}; then
	PSUM=`awk '{a+=$1}END{print a}' p.dat | awk '{printf "%9.4f\n",$1}'`
	QAVG=`awk '{a+=$1; n+=1}END{printf "%9.4f\n",a/n}' ${QFIL}`
	MAXV=`awk '(NR==1){a=-1000000}(a<$2&&($2>=0||$2<0)){a=$2}END{print a}' velocity.dat`

	if test `echo "${QAVG}<0.5" |bc -l` -eq 1 -o\
		`echo "${MAXV}>${MAXVEL}" |bc -l` -eq 1; then 
	    echo ${dir} $QAVG $MAXV $PSUM "WARNING: this simulation trajectory likely has errors and should not be used"
	else
	    echo ${dir} $QAVG $MAXV $PSUM "OK"
	fi
    fi

    cd $PWD0
    fi
    done
done
