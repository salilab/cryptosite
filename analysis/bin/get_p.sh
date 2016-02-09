#!/bin/bash

LIST="list"
PWD0=`pwd`

getstd()
{
XAVG=`awk '(NR==1){a=0;n=0}{a+=$'${2}'; n+=1}END{print a/n}' $1`
awk '(NR==1){n=0;a=0}{a+=($'${2}'-"'${XAVG}'")^2;n+=1}END{print sqrt(a/n)}' $1
}

if test -e ${PWD0}/p.log; then rm ${PWD0}/p.log; fi
for landscape in `cat ${LIST}`; do
    for dir in `ls -1d ${landscape}/*`; do
	if test -d ${dir}; then
	    if (test ! -e ${dir}/energy.dat); then continue; fi
	    cat ${dir}/energy.dat >>${PWD0}/tempgp77
	fi
    done

    emin=`awk '(NR==1){a=10000000}(a>$2&&($2>=0||$2<0)){a=$2}END{print a}' tempgp77`
    RT=`getstd tempgp77 2`
    zpart1=`awk '{a+=exp(-1*($2-'${emin}')/'${RT}')}END{print a}' tempgp77`
    rm tempgp77

    for dir in `ls -1d ${landscape}/*`; do
	if test -d ${dir}; then
	cd $dir
    
	if (test ! -e energy.dat); then continue; fi

	awk '{printf "%15.13f\n", exp(-1*($2-'${emin}')/'${RT}')/'${zpart1}'}' energy.dat >p.dat
	psum=(`awk '{a+=$1}END{printf "%9.5f\n",a}' p.dat`)
	echo $dir ${psum[0]} >>${PWD0}/p.log
	
	cd $PWD0
	fi
    done
done

