#!/bin/bash 
# 
#$ -S /bin/bash 
#$ -o /scrapp/peterc/AM
#$ -e /scrapp/peterc/AM
##$ -cwd 
#$ -r y 
#$ -j y 
#$ -l netappsali=3G  
#$ -l arch=linux-x64
#$ -l netappsali=3G 
#$ -l diva1=1G,database=1G,scratch=2G 
#$ -l h_rt=120:00:00 
#$ -l mem_free=2G 
#$ -p 0 
#$ -R yes
#$ -t 1-500


export TMPDIR="/scratch/peterc/$JOB_ID/$SGE_TAKS_ID"
mkdir -p $TMPDIR
cd $TMPDIR
pwd

TT="tt10"
ARRAY=(/scrapp/peterc/AM/$TT/pred_dECALCrAS1000/*)
DRIN=${ARRAY[$SGE_TASK_ID - 1]}
DIN=`echo ${DRIN} | cut -d '/' -f 7`
DROUT=`echo ${DRIN} | cut -d '/' -f 6,7`
DRT="./$DROUT"
PDB=`echo $DIN | cut -d \. -f 1`

echo
echo $ARRAY
echo $DIN
echo $DRIN
echo $DROUT
echo $DRT
echo $PDB
echo /scrapp/peterc/AM/$TT/frust/$DIN
echo

##mkdir -p scrapp/peterc/AM/$TT/frust/$DIN/
##cp $DRIN/pm.pdb.B1*1.pdb $DRIN/pm.pdb.B2*1.pdb $DRIN/pm.pdb.B3*1.pdb scrapp/peterc/AM/$TT/frust/$DIN/ 
##cp -r /scrapp/peterc/AM/$TT/frust/$DIN/pm.pdb.*.pdb .
cp -r /netapp/sali/peterc/Undrugabble/AllosModRuns/analysis/for_peter/* .
cp $DRIN/pm.pdb.B1*1.pdb $DRIN/pm.pdb.B2*1.pdb $DRIN/pm.pdb.B3*1.pdb .

cp /netapp/sali/peterc/Undrugabble/Server/SiteCrypt/AM_BMI.py .
cp /netapp/sali/peterc/Undrugabble/Server/SiteCrypt/CHASA.py .


for file in `ls pm.pdb.B*1.pdb | sed -n '1~4p'`
do
  /diva1/home/modeller/modpy.sh python AM_BMI.py $file
done


mkdir -p /scrapp/peterc/AM/$TT/frust/$DIN
cp *.feat pm.pdb.B*.pdb /scrapp/peterc/AM/$TT/frust/$DIN

date


