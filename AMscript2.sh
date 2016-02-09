#!/bin/bash                         
#                                     
#$ -S /bin/bash                        
#$ -o /netapp/sali/leon/cryptosite/src_multichain/sge                          
#$ -e /netapp/sali/leon/cryptosite/src_multichain/sge                        
#$ -l arch=linux-x64
#$ -l scratch=2G
#$ -l mem_free=2G

##export CURDIR=`pwd`
##export TMPDIR="/scratch/peterc/$JOB_ID/$SGE_TAKS_ID"
##mkdir -p $TMPDIR

##cp %s.pdb %s.features $TMPDIR
##cd $TMPDIR
##cp /netapp/sali/peterc/cryptosite/src_v1/PREDICTER.py /netapp/sali/peterc/cryptosite/src_v1/PolyS*.pkl .

module load sali-libraries
export PYTHONPATH="/netapp/sali/peterc/lib64/python"

python /netapp/sali/leon/cryptosite/src_multichain/PREDICTER.py XXX 

##cp *.pol.pred *.pol.pred.pdb $CURDIR
rm -rf $TMPDIR
date
