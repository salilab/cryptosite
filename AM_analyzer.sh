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


ARRAY=(/scrapp/tt1/pred_dECALCrAS1000/*)

export TMPDIR="/scratch/peterc/$JOB_ID/$SGE_TASK_ID"
mkdir -p $TMPDIR

cd $TMPDIR

DRIN=${ARRAY[$SGE_TASK_ID - 1]}
DROUT=`echo ${DRIN} | cut -d '/' -f 4,5`
DRT="./$DROUT"
mkdir -p $DRT

echo
echo $DRIN
echo $DROUT
echo $DRT
echo
#cp /scrapp/tt1/pred_dECALCrAS1000/1ALV.pdb_9/pm.pdb.B1*1.pdb ./pred_dECALCrAS1000/1ALV.pdb_9/ 
#cp /scrapp/tt1/pred_dECALCrAS1000/1ALV.pdb_9/pm.pdb.B2*1.pdb ./pred_dECALCrAS1000/1ALV.pdb_9/
#cp /scrapp/tt1/pred_dECALCrAS1000/1ALV.pdb_9/pm.pdb.B3*1.pdb ./pred_dECALCrAS1000/1ALV.pdb_9/
#cp /scrapp/tt1/pred_dECALCrAS1000/1ALV.pdb_9/pm.pdb.D00000001 ./pred_dECALCrAS1000/1ALV.pdb_9/

cp $DRIN/* $DRT/

pwd
ls -lt

echo "./pred_dECALCrAS1000" > list
bash /netapp/sali/peterc/Undrugabble/AllosModRuns/analysis/allosmod_analysis.sh

cat check_runs.out

cp -r $DRT /scrapp/peterc/AM/tt1/pred_dECALCrAS1000
cp check_runs.out /scrapp/peterc/AM/tt1/$DROUT


# pockets
ls $DRT/pm.pdb.B1*1.pdb $DRT/pm.pdb.B2*1.pdb $DRT/pm.pdb.B3*1.pdb > SnapList.txt

cp ~/Undrugabble/AllosModRuns/pocket_parser.py .

python pocket_parser.py

cp pockets.out /scrapp/peterc/AM/tt1/$DROUT 

date


