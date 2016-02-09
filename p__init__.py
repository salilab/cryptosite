import saliweb.backend
import logging, os
from random import randint,choice
import string, re



class Job(saliweb.backend.Job):

    runnercls = saliweb.backend.SaliSGERunner

    def run(self):

        ### - read in the chain and sequence
        paramsFile = open('param.txt','r')
        jobName, email, chainid = [i.strip() for i in paramsFile.readlines()]
        paramsFile.close()

        if 'stage.out' in os.listdir('.'): 
            data = open('stage.out')
            D = data.readlines()
            data.close()
            stage = D[0].strip()

            ### - set random number
            outran= open('random.out')
            rfil = outran.read().strip()
            outran.close()


            if stage=="pre-AllosMod": 
                
                os.system('zip -r %s.zip XXX' % rfil)
                ### - submit to AllosMod
                r = saliweb.backend.SaliWebServiceRunner('http://modbase.compbio.ucsf.edu/allosmod/job',
                                 ['name=%s' % rfil, 'jobemail=peterc@salilab.org', 'zip=@%s.zip' % rfil])
                out = open('stage.out','w')
                out.write('AllosMod')
                out.close()
                return r

            elif stage=='AllosMod':
                
                self.logger.info("Post-processing AllosMod results -bmi")
                
                ### - setup the SGE script #1
                script = """
export MYTMP=`mktemp -d`
cd $MYTMP

TT=%s
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



cp -r /netapp/sali/peter/cryptosite/src_multichain/analysis/for_peter/* .
cp $DRIN/pm.pdb.B*.pdb $DRIN/$PDB.pdb .

wc *.pdb

cp /netapp/sali/peter/cryptosite/src_multichain/AM_BMI.py .
cp /netapp/sali/peter/cryptosite/src_multichain/CHASA.py .
cp /netapp/sali/peter/cryptosite/src_multichain/soap_clean.py .


/diva1/home/modeller/modpy.sh python soap_clean.py $SGE_TASK_ID
##python soap_clean.py $SGE_TASK_ID

cat SnapList.txt


sleep 20

## - pockets
cp /netapp/sali/peter/cryptosite/src_multichain/pocket_parser.py .

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

cp *.dat /scrapp/AM/$TT/$DROUT

cat check_runs.out


ls -l

#rm -rf $MYTMP
date
""" % (rfil,)
            
                r = self.runnercls(script)
                r.set_sge_options('-l arch=linux-x64 -l scratch=2G -l mem_free=6G -t 1-50')

                self.logger.info("Calculated pockets for AllosMod results")

                out = open('stage.out','w')
                out.write('AllosMod-bmi')
                out.close()

                return r


            elif stage=='AllosMod-bmi':

                self.logger.info("Applying SVM model")

                ### - setup the SGE script #1
                script = """
##export CURDIR=`pwd`
##export TMPDIR="/scratch/peterc/$JOB_ID/$SGE_TAKS_ID"
##mkdir -p $TMPDIR

##cp %s.pdb %s.features $TMPDIR
##cd $TMPDIR
##cp /netapp/sali/peterc/cryptosite/src_v1/PREDICTER.py /netapp/sali/peterc/cryptosite/src_v1/PolyS*.pkl .

module load sali-libraries
export PYTHONPATH="/netapp/sali/peterc/lib64/python"

python PREDICTER.py %s

##cp *.pol.pred *.pol.pred.pdb $CURDIR
rm -rf $TMPDIR
date""" % ('XXX'+chainid,'XXX'+chainid,'XXX'+chainid)

                r = self.runnercls(script)
                r.set_sge_options('-l arch=linux-x64 -l scratch=2G -l mem_free=2G')

                self.logger.info("Prediction DONE!")

                out = open('stage.out','w')
                out.write('DONE')
                out.close()

                return r


        else: 
            ### - set logging file
            self.logger.setLevel(logging.INFO)
            self.logger.info("Beginning preprocess() for job %s " %self.name)


            ### - set random number
            rint = randint(100,999)
            rfil = ''.join([choice(string.letters) for nl in range(3)])+'%3i' % (rint,)
            outran= open('random.out','w')
            outran.write(rfil)
            outran.close()


            ### - read in the chain and sequence
            #paramsFile = open('param.txt','r')
            #jobName, email, chainid = [i.strip() for i in paramsFile.readlines()]
            #paramsFile.close()

            ### - setup the SGE script #1
            script = """

##export TMPDIR="/scratch/peterc/$JOB_ID"
##mkdir -p $TMPDIR
ls -lt
##cp input.pdb input.seq param.txt $TMPDIR

##cd $TMPDIR
pwd

##cp -r /netapp/sali/peterc/cryptosite/src_v1 $TMPDIR
cp -r /netapp/sali/peterc/cryptosite/src_multichain/* .

##cd ${TMPDIR}/src
##cp ../input.pdb ../input.seq ../param.txt .
pwd
ls -lt

export LD_LIBRARY_PATH="/netapp/sali/peterc/Undrugabble/Software/concavity_distr/bin/x86_64"
echo "working on"

module load modeller

##/diva1/home/modeller/modpy.sh python mainer.py 
python mainer.py %s

module unload modeller

echo "pre-AllosMod" > stage.out

##mkdir -p /scrapp/cryptosite/%s
##cp -r * /scrapp/cryptosite/%s

date

""" % ('XXX'+chainid, rfil, rfil)
        
            r = self.runnercls(script)
            r.set_sge_options('-l arch=linux-x64 -l diva1=1G -l scratch=2G -l mem_free=2G')

	    self.logger.info("Calculated bioinformatics features for job: %s" % rfil)
            self.logger.info("\n\nSubmitting to AllosMod")

            return r

    def postprocess(self, results=None):
        ### reschedule run
        data = open('stage.out')
        D = data.readlines()
        data.close()
        stage = D[0].strip()

        outran= open('random.out')
        rfil = outran.read().strip()
        outran.close()

        paramsFile = open('param.txt','r')
        jobName, email, chainid = [i.strip() for i in paramsFile.readlines()]
        jobName = jobName.replace(' ','_')
        paramsFile.close()


        if stage=="pre-AllosMod": self.reschedule_run()
        elif stage=="AllosMod": self.reschedule_run() 
        elif stage=="AllosMod-bmi":

            #if stage=="pre-AllosMod": 

            ### - gather the AM data
            self.logger.info("Gathering AllosMod results")
            os.system('python /netapp/sali/peterc/cryptosite/src_multichain/gatherer.py %s %s' % (rfil,chainid))

            ### - run SVM
            ## TO DELETE
            data = open('stage.out','w')
            data.write('AllosMod-bmi')
            data.close()
            ##
            self.reschedule_run()
        elif stage=="DONE":
            ### - make chimera session file
            self.logger.info("Completing the job: writing a Chimera session file.")
            #paramsFile = open('param.txt','r')
            #jobName, email, chainid = [i.strip() for i in paramsFile.readlines()]
            #jobName = jobName.replace(' ','_')
            #paramsFile.close()

            os.system('cp XXX%s.pol.pred cryptosite.pol.pred' % chainid)
            os.system('cp XXX%s.pol.pred.pdb cryptosite.pol.pred.pdb' % chainid)
            
            data = open('script.chimerax')
            chimeraSession = data.read()
            data.close()

            (filepath, jobdir) = os.path.split(self.url)
            pdb_temp = "/cryptosite.pol.pred.pdb?"
            pdb_url = re.sub('\?', pdb_temp, jobdir)
            pdburl = filepath + "/" + pdb_url
            urlAddress = self.url

            ftr_temp = "/cryptosite.pol.pred?"
            ftr_url = re.sub('\?', ftr_temp, jobdir)
            ftrurl = filepath + "/" + ftr_url

            chimeraSession += """
open_files("%s", "%s")
]]>
</py_cmd>
</ChimeraPuppet>""" % (pdburl,ftrurl)# jobName, jobName)

            out = open('cryptosite.chimerax', 'w')
            out.write(chimeraSession)
            out.close()
            self.logger.info("Job completed! Results available at: %s" % urlAddress)

            os.system('zip chimera.zip cryptosite.chimerax')

        else: 
            self.logger.info("\n\nEXITING: no stage output!")
            pass


def get_web_service(config_file):
    db = saliweb.backend.Database(Job)
    config = saliweb.backend.Config(config_file)
    return saliweb.backend.WebService(config, db)

