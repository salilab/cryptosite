import glob
import os
import sys


dirs = glob.glob('/scrapp/peterc/AllosModRuns/%s/output/input/pred_*_*00/*/model_run.py' % sys.argv[-1])
#                                      output_11_21/1EXMA_1HA3B_MAU/pred_dE0.1rAS1000/1EXMA.pdb_0/model_run.py

'''
dirs = []
for fil in glob.glob("output_*/*/pred_dE0.1rAS1000/*/model_run.py"):
        f = fil.rsplit('/',1)[0]
        s = len([i for i in os.listdir(f) if i[-4:]=='.pdb' and 'pm.pdb.B1' in i])
        p = fil.split('/')[1]
        if s<500:
                if int(f.rsplit('_',1)[-1]) not in range(1,7): continue
                #print 'am_'+f.split('/')[-1].replace('.pdb','')
		dirs.append('/scrapp/peterc/AllosModRuns/'+fil)
'''

temp = 200

outt = open('blaa3.sh','w')
for d in dirs:
	#d = dirs[0]
	path = d.rsplit('/',1)[0]
	predfile = path.split('/')[-1]
	pdb,runnum = predfile.split('.')[0], predfile.split('_')[-1]
	if int(runnum) not in range(0,10): continue
	
	##/scrapp/peterc/AllosModRuns/output_41-50/2AKAA_1YV3A_BIT/pred_dE0.1rAS1000_300/2AKAA.pdb_3/
	print pdb,path,path.split('/')[-2].split('_')[-1], path.split('/')[-1].split('_')[-1]
	outt.write('qsub am_%s_%s_%s_1.sh\n' % (pdb,path.split('/')[-2].split('_')[-1], path.split('/')[-1].split('_')[-1]))	

	text1 = """
#!/bin/bash 
# 
#$ -S /bin/bash 
#$ -o /scrapp/peterc/AllosModRuns/%s
#$ -e /scrapp/peterc/AllosModRuns/%s
##$ -cwd 
#$ -r y 
#$ -j y 
#$ -l netappsali=3G  
#$ -l arch=linux-x86 
#$ -l netappsali=3G 
#$ -l diva1=1G,database=1G,scratch=2G 
#$ -l h_rt=120:00:00 
#$ -l mem_free=1G 
#$ -p 0 
#$ -R yes


ulimit -v 6000000

export TMPDIR="/scratch/peterc/$JOB_ID/"
mkdir -p $TMPDIR
cd $TMPDIR
pwd

export amberpath='/netapp/sali/pweinkam/utils/AMBER/amber11/bin'
export AMBERHOME=/netapp/sali/pweinkam/utils/AMBER/amber11/

export PATH=/usr/sbin:/usr/local/bin:/sbin:/usr/local/lib:/netapp/sali/pweinkam/bin:/bin:/sbin:/usr/bin:/usr/include:/usr/sbin:/usr/local/bin:/usr/local/sbin:/netopt/mpi/mpich2/bin:/diva1/home/modeller/SVN/bin:/netapp/sali/pweinkam/utils/mmtsb_toolset/perl:/netapp/sali/pweinkam/utils/mmtsb_toolset/bin:/netapp/sali/pweinkam/utils/AMBER/amber11/bin:.:$PATH

export PDB='%s'

cp -r /netapp/sali/peterc/Undrugabble/AllosModRuns/analysis/for_peter/* .
cp %s/pm.pdb.B1*01.pdb .

#echo "${PDB}.pdb A" > input
ls pm.pdb.B1*01.pdb | awk '{print $1,"A"}' > input

head input


# clean PDB (ATOM, TER, HETATM lines only)
for file in pm.pdb.B1*01.pdb
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

for file in pm.pdb.B1*01.pdb
do
  cp ${file} ${file}-A.pdb
done

mkdir -p %s
cp -r * %s


date

""" % (sys.argv[-1],sys.argv[-1],pdb, path, '/'.join(path.split('/')[:4]+[sys.argv[-1]+'/FrstOut']+path.split('/')[5:]), '/'.join(path.split('/')[:4]+[sys.argv[-1]+'/FrstOut']+path.split('/')[5:]) )


	text = """
#!/bin/bash 
# 
#$ -S /bin/bash 
#$ -o /scrapp/peterc/AllosModRuns/%s
#$ -e /scrapp/peterc/AllosModRuns/%s
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
##$ -t 1-6


ulimit -v 6000000

export TMPDIR="/scratch/peterc/$JOB_ID/"
mkdir -p $TMPDIR
cd $TMPDIR
pwd

export amberpath='/netapp/sali/pweinkam/utils/AMBER/amber11/bin'
export AMBERHOME=/netapp/sali/pweinkam/utils/AMBER/amber11/

export PATH=/usr/sbin:/usr/local/bin:/sbin:/usr/local/lib:/netapp/sali/pweinkam/bin:/bin:/sbin:/usr/bin:/usr/include:/usr/sbin:/usr/local/bin:/usr/local/sbin:/netopt/mpi/mpich2/bin:/diva1/home/modeller/SVN/bin:/netapp/sali/pweinkam/utils/mmtsb_toolset/perl:/netapp/sali/pweinkam/utils/mmtsb_toolset/bin:/netapp/sali/pweinkam/utils/AMBER/amber11/bin:.:$PATH

cp -r %s/* .

perl 3_amber_minimize.pl

rename pdb-sequence pdb-A-sequence *.pdb-sequence
perl 4_energy_decomose.pl

for file in pm.pdb.B1*01.pdb
do
  ./get_delEres.sh ${file}-A_statistics.out ${file}-A.pdb > ${file}.frust
done

cp *.frust %s

date

""" % (sys.argv[-1],sys.argv[-1],'/'.join(path.split('/')[:4]+[sys.argv[-1]+'/FrstOut']+path.split('/')[5:]), '/'.join(path.split('/')[:4]+[sys.argv[-1]+'/FrstOut']+path.split('/')[5:]))

	#/scrapp/peterc/AllosModRuns/output_41-50/2AKAA_1YV3A_BIT/pred_dE0.1rAS1000_300/2AKAA.pdb_3/	
	out = open('am_%s_%s_%s_1.sh' % (pdb,path.split('/')[-2].split('_')[-1], path.split('/')[-1].split('_')[-1]),'w')
	out.write(text1)
	out.close()

	out1 = open('am_%s_%s_%s_2.sh' % (pdb,path.split('/')[-2].split('_')[-1], path.split('/')[-1].split('_')[-1]),'w')
        out1.write(text)
        out1.close()
	
outt.close()



