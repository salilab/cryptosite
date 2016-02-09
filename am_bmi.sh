#!/bin/bash 
# 
#$ -S /bin/bash 
#$ -o /scrapp/peterc/Re-analysis/ 
#$ -e /scrapp/peterc/Re-analysis/
##$ -cwd 
#$ -r y 
#$ -j y 
#$ -l netappsali=3G  
#$ -l diva1=1G,database=1G,scratch=2G 
#$ -l h_rt=120:00:00 
#$ -l mem_free=6G 
#$ -p 0 
#$ -R yes
#$ -t 1-150




export TMPDIR=/scratch
export MYTMP=`mktemp -d`
cd $MYTMP

pwd


TT="tt12"
ARRAY=(/scrapp/peterc/Re-analysis/$TT/pred_dECALCrAS1000/*)
DRIN=${ARRAY[$SGE_TASK_ID - 1]}
DIN=`echo ${DRIN} | cut -d '/' -f 7`
DROUT=`echo ${DRIN} | cut -d '/' -f 6,7`
DRT="./$DROUT"
PDB=`echo $DIN | cut -d \. -f 1`

echo $DRIN
echo $DIN
echo $DROUT
echo $DRT
echo $PDB



cp -r /netapp/sali/peterc/cryptosite/src_v2/analysis/for_peter/* .
cp $DRIN/pm.pdb.B???????1.pdb $DRIN/$PDB.pdb .

wc *.pdb

cp /scrapp/peterc/Re-analysis/AM_BMI.py .
cp /netapp/sali/peterc/cryptosite/src/CHASA.py .
cp /scrapp/peterc/Re-analysis/soap_clean.py .

/diva1/home/modeller/modpy.sh python soap_clean.py $SGE_TASK_ID

cat SnapList.txt


sleep 20

## - pockets
cp /scrapp/peterc/Re-analysis/pocket_parser.py .

python pocket_parser.py $SGE_TASK_ID

cp pockets.out /scrapp/peterc/Re-analysis/$TT/$DROUT

sleep 20


## - AM featres
/diva1/home/modeller/modpy.sh python AM_BMI.py $file $SGE_TASK_ID

cp am_features.out /scrapp/peterc/Re-analysis/$TT/$DROUT


##rm -rf TMPDIR
rm -rf $MYTMP
date

