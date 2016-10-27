from SiteCrypt import PATH2MUSCLE, PATH2BIOPYTHON

import os, sys, warnings, subprocess
sys.path.append(PATH2BIOPYTHON)
from Bio import PDB
from Bio.Blast import NCBIXML
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import numpy
import re
import subprocess
from modeller import *
from modeller.automodel import *
from operator import itemgetter


parser = PDBParser()


def get_pdb_seq(pdb, chain):
	'''
	Read a PDB[pdb] input and outputs a sequence of a desired chain(s)[chain].
	'''
	#TODO: implement for a set of chains

        pdbh = open(pdb)
        structure = parser.get_structure("protein", pdb)
        model = structure[0]
        chain = model[chain]

        ppb = PPBuilder()
        print chain
        ss = ppb.build_peptides(chain)
        print ss
        seq = ''.join([str(i.get_sequence()) for i in ss])

        #print [str(i.get_sequence()) for i in ss]

        return seq



def muscleAlign(qSeq, sSeq, pdb, chain):
	'''
	Align two sequences and produce Modeller-compatible output
	'''

	# --- write sequence file
	output = open('sequences.seq','w')
	output.write( '>%s%s\n' % (pdb,chain) )
	output.write( qSeq+'\n' )
	output.write( '>%s%spdb\n' % (pdb,chain) )
	output.write( sSeq )
	output.close()

	# --- align using Muscle
        cmd = [PATH2MUSCLE+"muscle", "-in", "sequences.seq", "-out", "alignment.ali"]
        print cmd
        prc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        prc.wait()

        data = open('alignment.ali')
        D = data.read().split('>')
        data.close()

        strc = D[1].split('\n')
        strcid,strcsq = strc[0], ''.join(strc[1:])
        seq = D[2].split('\n')
        seqid,seqsq = seq[0], ''.join(seq[1:])

	# --- clean
	os.system("rm alignment.ali sequences.seq")

	return (strcsq, seqsq)


def get_gaps(chainso,Ls):

	data = open('alignment.pir')
	D = data.read().split('>P1;')
	data.close()

	strc = D[1].split('\n')
	strcid,strcsq = strc[0], ''.join(strc[2:])[:-1]
	seq = D[2].split('\n')
	seqid,seqsq = seq[0], ''.join(seq[2:])[:-1]

	# --- find gaps using regular expression matching
	gaps = []
	L=0
        seqsq = seqsq.split('/')
	strcsq = strcsq.split('/')
	chains=map(chr, range(65, 65+len(strcsq)))
	for strcsq_index, strcsq_part in enumerate(strcsq):
		gap_reg = re.compile('\W-*\W')
		iterator = gap_reg.finditer(strcsq_part)
		for match in iterator:
			if len(chains)>1:
				gaps.append(','.join([str(match.span()[0]+L+1)+':'+chains[strcsq_index], str(match.span()[1]+L)+':'+chains[strcsq_index]]))
			else:
				gaps.append(','.join([str(match.span()[0]+L+1)+':', str(match.span()[1]+L)+':']))
                #L+=int(Ls[chainso[strcsq_index]][1])
                L+=len(seqsq[strcsq_index])

	return gaps




def concatenate_models(pdb,chains):

        openr={}

        for chain in chains:
                op=open('%s%s_mdl.pdb' % (pdb,chain))
                openr[chain]=op.readlines()
                op.close()


        w=open('XXX_.pdb', 'w')
        for chain in chains:
                for line in openr[chain]:
                        if line[:4]=='ATOM' or line[:3]=='TER':
                                w.write(line)

        w.close()




def build_model(pdb,chains,chainLs):
	'''
	Build model using Modeller, treating residues in the structure as rigid, and
	loop modeling for the rest.
	'''

	# --- get gaps
        gaps = get_gaps(chains,chainLs)
	out = open('gaps_%s.txt' % pdb,'w')
        out.write(pdb+'_mdl\t'+str(gaps))
        out.close()

	# --- set up modeling
        env = environ()

        env.io.atom_files_directory = ['.', '../atom_files']
        #env.io.hetatm = True
	#env.io.convert_modres = True	

        class MyLoop(loopmodel):
                def select_loop_atoms(self):
                        gaps = get_gaps(chains,chainLs)
                        return selection(self.residue_range(i.split(',')[0], i.split(',')[1]) for i in gaps)
		def special_restraints(self, aln):
                        rsr = self.restraints
                        wholeSel = selection(self) - self.select_loop_atoms()
                        r = rigid_body(wholeSel)
                        rsr.rigid_bodies.append(r)

        if len(gaps)>0:

                a = MyLoop(env,
                   alnfile  = 'alignment.pir',      # alignment filename
                   knowns   = pdb,                  # codes of the templates
                   sequence = pdb.lower()+'_'+'X',# code of the target
                   loop_assess_methods=assess.normalized_dope) # assess each loop with DOPE

                a.starting_model= 1                 # index of the first model
                a.ending_model  = 1                 # index of the last model

                a.loop.starting_model = 1           # First loop model
                a.loop.ending_model   = 10          # Last loop model

                a.make()                            # do modeling and loop refinement
	
		ok_models = [i for i in a.loop.outputs if i['failure'] is None]
	
		# -- select best model
		zscores = sorted([(i['name'], i['Normalized DOPE score']) for i in ok_models], key=itemgetter(1))
		bestModel = zscores[0][0]
		os.system('mv %s %s_mdl.pdb' % (bestModel, pdb))
		os.system('rm %s_X.*' % (pdb.lower(),))

        else:
                a = automodel(env, alnfile='alignment.pir',
                      knowns=pdb, sequence=pdb.lower()+'_X',
                      assess_methods=(assess.DOPE,))
                a.very_fast()
                a.starting_model = 1
                a.ending_model = 1
                a.make(exit_stage=2)

                os.system('mv %s_X.ini %s_mdl.pdb' % (pdb.lower(), pdb))
		os.system('rm %s_X.*' % (pdb.lower(),))
		#os.system('mv %s.pdb /netapp/sali/peterc/Undrugabble/PDB/' % (pdb+chains[-1],))

	if len(chains)==1:
        	mdl  = model(env, file='XXX_mdl')
	        mdl.rename_segments(segment_ids='A')
        	mdl.write(file='%s_mdl.pdb' % (pdb))

        # Read an alignment for the transfer
        #aln = alignment(env, file='alignment.pir', align_codes=('XXX', 'xxx_X'))
        # Read the template and target models:
        #mdl2 = model(env, file='XXX')
        #mdl  = model(env, file='XXX_mdl')
        # Transfer the residue and chain ids and write out the new MODEL:
        #mdl.res_num_from(mdl2, aln)
        #mdl.write(file='%s_.pdb' % (pdb))


        #mdl  = model(env)

        #for chain in chains:
        #        mdl.read(file='XXX_', model_format='PDB', model_segment=('FIRST:%s' % chain, 'LAST:%s' % chain))
	#mdl.rename_segments(segment_ids='A')
        #mdl.write(file='%s_mdl.pdb' % (pdb))

        #concatenate_models(pdb,chains)

	#os.system('mv %s_.pdb %s_mdl.pdb' % (pdb,pdb))
	
	#RESMAP={}
	#mdl  = model(env, file='XXX_mdl')
	#mdl2 = model(env, file='XXX_')
	#for indexr,res in enumerate(mdl.residues):
	#	res2=mdl2.residues[indexr]
	#	RESMAP[(res.name, res.num, res.chain.name)]=[res2.name, res2.num, res2.chain.name]
		

	#return RESMAP
