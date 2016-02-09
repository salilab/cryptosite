#!/bin/bash 
# 
#$ -S /bin/bash 
#$ -o /scrapp/peterc/AM
#$ -e /scrapp/peterc/AM
##$ -cwd 
#$ -r y 
#$ -j y 
#$ -l netappsali=3G  
#$ -l arch=linux-x86 
#$ -l diva1=1G,database=1G,scratch=2G 
#$ -l h_rt=120:00:00 
#$ -l mem_free=1G 
#$ -p 0 
#$ -R yes
#$ -t 1-500


export TMPDIR="/scratch/peterc/$JOB_ID/$SGE_TASK_ID"
mkdir -p $TMPDIR
cd $TMPDIR
pwd


export amberpath='/netapp/sali/pweinkam/utils/AMBER/amber11/bin'
export AMBERHOME=/netapp/sali/pweinkam/utils/AMBER/amber11/

export PATH=/usr/sbin:/usr/local/bin:/sbin:/usr/local/lib:/netapp/sali/pweinkam/bin:/bin:/sbin:/usr/bin:/usr/include:/usr/sbin:/usr/local/bin:/usr/local/sbin:/netopt/mpi/mpich2/bin:/diva1/home/modeller/SVN/bin:/netapp/sali/pweinkam/utils/mmtsb_toolset/perl:/netapp/sali/pweinkam/utils/mmtsb_toolset/bin:/netapp/sali/pweinkam/utils/AMBER/amber11/bin:.:$PATH

TT='tt2'
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


cp -r /netapp/sali/peterc/Undrugabble/AllosModRuns/analysis/for_peter/* .
cp $DRIN/pm.pdb.B1*1.pdb $DRIN/pm.pdb.B2*1.pdb $DRIN/pm.pdb.B3*1.pdb .

#echo "${PDB}.pdb A" > input
ls pm.pdb.B*1.pdb | awk '{print $1,"A"}' > input

head input


# clean PDB (ATOM, TER, HETATM lines only)
for file in pm.pdb.B*1.pdb
do
  python cleaner_frust.py ${file}
done

perl 0_pdb_check.pl
perl 1_run_HBPLUS.pl

ls -lt

rename hb2 pdb-A.hb2 *.hb2
rename nb2 pdb-A.nb2 *.nb2

perl 2_Sulfate.pl

ls -lt

for file in pm.pdb.B*1.pdb
do
  cp ${file} ${file}-A.pdb
done

mkdir -p /scrapp/peterc/AM/$TT/frust/$DIN
cp -r * /scrapp/peterc/AM/$TT/frust/$DIN


date

