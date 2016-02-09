#!/bin/bash

LIST=$1
PWD0=`pwd`
PDIR=`awk '{print "dirname "$1}' $LIST | sh | head -n1`
cd $PDIR

awk '{print "basename "$1}' ${PWD0}/$LIST | sh >input

for s in `cat input`; do
awk '($1=="ATOM"||$1=="HETATM"){print $0}' $s |\
sed "s/CD  ILE/CD1 ILE/g" | sed "s/OT1/O  /g" | sed "s/OT2/OXT/g" >tempi994
mv tempi994 $s
done

perl /netapp/sali/pweinkam/amber/1_run_HBPLUS.pl

