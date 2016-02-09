#!/bin/bash                         
#                                     
#$ -S /bin/bash                        
#$ -o /netapp/sali/leon/cryptosite/src_multichain/sge                          
#$ -e /netapp/sali/leon/cryptosite/src_multichain/sge                        
#$ -l arch=linux-x64
#$ -l scratch=2G
#$ -l mem_free=6G
#$ -t 1-1



export MYTMP=`mktemp -d`
cd $MYTMP

TT='TNi261'
ARRAY=(/scrapp/$TT/pred_dECALCrAS1000/*)
DRIN=${ARRAY[$SGE_TASK_ID - 1]}
DIN=`echo ${DRIN} | cut -d '/' -f 5`
DROUT=`echo ${DRIN} | cut -d '/' -f 4,5`
DRT="./$DROUT"
mkdir -p $DRT
PDB=`echo $DIN | cut -d \. -f 1`



echo $ARRAY
echo $DRIN
echo $DROUT
echo $DRT
echo $PDB



cp -r /netapp/sali/leon/cryptosite/src_multichain/analysis/* .
cp $DRIN/pm.pdb.B299*1.pdb $DRIN/$PDB.pdb .

wc *.pdb

cp /netapp/sali/leon/cryptosite/src_multichain/AM_BMI.py .
cp /netapp/sali/leon/cryptosite/src_multichain/CHASA.py .
cp /netapp/sali/leon/cryptosite/src_multichain/soap_clean.py .
cp /netapp/sali/leon/cryptosite/src_multichain/gatherer.py .


/diva1/home/modeller/modpy.sh python soap_clean.py $SGE_TASK_ID
##python soap_clean.py $SGE_TASK_ID

cat SnapList.txt


sleep 20

## - pockets
cp /netapp/sali/leon/cryptosite/src_multichain/pocket_parser.py .

python pocket_parser.py $SGE_TASK_ID


mkdir -p /scrapp/AM/$TT/$DROUT

cp pockets.out /scrapp/AM/$TT/$DROUT

sleep 20



## - AM features
/diva1/home/modeller/modpy.sh python AM_BMI.py $file $SGE_TASK_ID
##python AM_BMI.py $file $SGE_TASK_ID


cp am_features.out /scrapp/AM/$TT/$DROUT


## - copy energy.dat
ls -l

pwd

cp $DRIN/* $DRT/

echo "./pred_dECALCrAS1000" > list
bash /netapp/sali/peterc/cryptosite/src_v1/analysis/allosmod_analysis.sh

cat check_runs.out
cp *.dat /scrapp/AM/$TT/$DROUT

## - gatherer
python gatherer.py


ls -l

#rm -rf $MYTMP
date


