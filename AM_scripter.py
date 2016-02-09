import glob
import sys


#dirs = glob.glob('/scrapp/peterc/AM/tt3/output/*/pred_*/*/model_run.py')
dirs = glob.glob('/scrapp/peterc/AM/output/*/pred_*/*/model_run.py')
#                                      /scrapp/peterc/AllosModRuns/Test/output/input/pred_dECALCrAS1000/1BNC.pdb_0
#                                  tt1/output/1A8IA/pred_dE0.1rAS1000/1A8I.pdb_0/

outt = open('blaa.sh','w')
for d in dirs:
	#d = dirs[0]
	print d
	path = d.rsplit('/',1)[0]
	predfile = path.split('/')[-1]
	print predfile
	pdb,runnum = predfile.split('.')[0], predfile.split('_')[-1]
	print pdb,runnum
	outt.write('qsub am_%s_%s.sh\n' % (pdb,runnum))	

	text = """
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
#$ -t 1-5


ulimit -v 6000000

export TMPDIR="/scratch/peterc/$JOB_ID/$SGE_TASK_ID"
mkdir -p $TMPDIR
cd $TMPDIR
pwd

echo "here" > REAMDE


cp -r %s .
##cp -r /netapp/sali/peterc/Undrugabble/AllosModRuns/analysis .
pwd > list

cd %s
##cp /netapp/sali/peterc/Undrugabble/AllosModRuns/pocket_parser.py .
##cp /netapp/sali/peterc/Undrugabble/AllosModRuns/rename.py .
cp /netapp/sali/peterc/Undrugabble/Server/SiteCrypt/allosmod.py_fold .
mv allosmod.py_fold allosmod.py

""" % (path, predfile)
	
	
	text +="""
# - run at different temps
let "temp = $SGE_TASK_ID * 100 + 200"
echo $temp
/diva1/home/modeller/modSVN model_run.py $temp
#python rename.py $temp

""" 

	text +="""

mkdir -p %s
cp * %s

date
#rm -r $TMPDIR

""" % (path.rsplit('/',1)[0]+'_${temp}/'+predfile, path.rsplit('/',1)[0]+'_${temp}/'+predfile)


	out = open('am_%s_%s.sh' % (pdb,runnum),'w')
	out.write(text)
	out.close()
	
outt.close()



