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

export amberpath='/netapp/sali/pweinkam/utils/AMBER/amber11/bin'
export AMBERHOME=/netapp/sali/pweinkam/utils/AMBER/amber11/

export PATH=/usr/sbin:/usr/local/bin:/sbin:/usr/local/lib:/netapp/sali/pweinkam/bin:/bin:/sbin:/usr/bin:/usr/include:/usr/sbin:/usr/local/bin:/usr/local/sbin:/netopt/mpi/mpich2/bin:/diva1/home/modeller/SVN/bin:/netapp/sali/pweinkam/utils/mmtsb_toolset/perl:/netapp/sali/pweinkam/utils/mmtsb_toolset/bin:/netapp/sali/pweinkam/utils/AMBER/amber11/bin:.:$PATH

TT="tt1"
ARRAY=(/scrapp/$TT/pred_dECALCrAS1000/*)
DRIN=${ARRAY[$SGE_TASK_ID - 1]}
DIN=`echo ${DRIN} | cut -d '/' -f 5`
DROUT=`echo ${DRIN} | cut -d '/' -f 4,5`
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


cp -r /scrapp/peterc/AM/$TT/frust/$DIN/* .

perl 3_amber_minimize.pl

rename pdb-sequence pdb-A-sequence *.pdb-sequence
perl 4_energy_decomose.pl

for file in `ls pm.pdb.B*1.pdb | sed -n '1~4p'`
do
  ./get_delEres.sh ${file}-A_statistics.out ${file}-A.pdb > ${file}.frust
done

cp *.frust /scrapp/peterc/AM/$TT/frust/$DIN

date

