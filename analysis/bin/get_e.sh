#!/bin/bash
#input: a list of all run directories

LIST="list"
PWD0=`pwd`

for landscape in `cat ${LIST}`; do
    for dir in `ls -1d ${landscape}/*`; do
    if test -d ${dir}; then
    cd $dir

    if (test ! -e pm.pdb.B10010001.pdb); then continue; fi
    NMOV=`ls -1 pm.pdb.B[1-8]*pdb | awk 'END{print NR}'`

    grep -B2 Step pm.pdb.D00000001 | grep -v "#" | grep -v "-" >temp55
    tail -n1 pm.pdb.D00000001 >>temp55

    ENDMOV=`awk 'END{print NR}' temp55`
    STARTMOV=$((${ENDMOV} - ${NMOV} +1))

    awk 'BEGIN{ctr=1}(NR>='${STARTMOV}'&&NR<='${ENDMOV}'){print ctr,$2,$3,$4,$5,$6;ctr+=1}' temp55 >energy.dat

    rm temp55

    cd $PWD0
    fi
    done
done

