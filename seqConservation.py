from SiteCrypt import PATH2BLAST, PATH2UNIPROT, PATH2UCLUST, PATH2BIOPYTHON

import os, sys, warnings, subprocess
sys.path.append(PATH2BIOPYTHON)
from Bio.Blast import NCBIXML
from Bio.PDB.PDBParser import PDBParser
from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder
import numpy



def run_blast(pdb):
	'''
	Run blast of a given sequence
	(only works on netapp disc where Uniprot is located)
	'''

        cmd = ["-query","test.seq","-db",PATH2UNIPROT,"-evalue","0.00001","-out",pdb+".blast","-outfmt","5","-num_alignments","10000"]
        cmd = [PATH2BLAST+"blastp"] + cmd
        #os.system(cmd)
	
	print cmd
	print "Running BLAST ..."

        prc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        prc.wait()

	print "BLASTing finished!"




def parse_blast(blastOut, pdb, qseq, evalue=0.00001):
	'''
	Parse XML Blast outputs.
	Parameters:
	  - evalue: set an alignment cutoff
	  - qseq:   query sequence as a string
	'''

        blast_record = NCBIXML.read(open(blastOut))
        A = numpy.zeros((len(qseq),21))
        q = ['A','C','D','E','F','G','H','I','K','L','M','N','P','R','S','T','V','Y','W','Q','-']

        out = open(pdb+'.ali','w')
        out.write('>%sq\n' % pdb)
        out.write(qseq+'\n')

        Seqs = {pdb+'q':qseq}
        for alignment in blast_record.alignments:

                for hsp in alignment.hsps:
                        if float(hsp.expect) > evalue: continue
                        sseq = alignment.title.split('|')[-2]
                        if sseq not in Seqs:
                                Seqs[sseq] = 1
                                out.write('>'+sseq+'1'+'\n')
                                out.write(hsp.sbjct+'\n')
                        else:
                                Seqs[sseq] += 1
                                out.write('>'+sseq+str(Seqs[sseq])+'\n')
                                out.write(hsp.sbjct+'\n')
        out.close()
        # --- get A matrix
        Clusters = ucluster(pdb+'.ali')
        M = len(Clusters)
        Meff = 0.
        Seqs = {}
        for alignment in blast_record.alignments:

                for hsp in alignment.hsps:
                        if float(hsp.expect) > evalue: continue
                        sseq = alignment.title.split('|')[-2]
			qstart = hsp.query_start-1
                        qi = qstart
                        hquery = hsp.query
                        hsbjct = hsp.sbjct
                        if sseq not in Seqs:
                                Seqs[sseq] = 1
                                sseq += '1'
                        else:
                                Seqs[sseq] += 1
                                sseq += str(Seqs[sseq])
                        ma = 1./Clusters[sseq]
                        
                        for i in xrange(len(hquery)):
                                if hquery[i] != '-':
                                        if hquery[i]!=qseq[qi]:
                                                print 'Some Error in alignment!'
                                                print hquery[i],qseq[qi],hquery[i]==qseq[qi]
                                                print peter
                                        if hsbjct[i] in q: A[qi,q.index(hsbjct[i])] += ma #1.
                                        else: pass
                                        qi += 1
                                else: pass
                        Meff += ma
        for i,a in enumerate(qseq):
                try:
                        A[i,q.index(a)] += 1./Clusters[pdb+'q'] #1.
                        Meff += 1./Clusters[pdb+'q']
                except ValueError: pass

        # --- re-weight the A matrix, correct for lambda factor
        lmbd = Meff
        out = open(pdb+'.sqc','w')
        for i in xrange(len(A)):
                Si = 0
                Fa = sum(A[i])
                for j in xrange(len(q)):
                        si = (1./(Meff+lmbd)) * ( (lmbd/len(q)) + A[i,j] ) #A[i,j]/Fa 
                        A[i,j] = si
                        if si>0.: Si -= si*numpy.log(si)
                out.write('\t'.join([str(i+1), qseq[i], str(Si)]) + '\n')
        out.close()



def ucluster(ali, cutoff=0.8):
	'''
	Cluster homologs.
	Parameters:
	  - cutoff: a distance cutoff to form sequence clusters.
	'''

        cmd = [PATH2UCLUST + "usearch4.0.43", "--cluster", ali, "--uc", "results.uc", "--id", str(cutoff), "--usersort"]
        prc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        prc.wait()

        data = open('results.uc')
        D = data.readlines()
        data.close()

        Clusters = {}
        for d in D:
                if 'S'==d[0] or 'H'==d[0]:
                        d = d.strip().split('\t')
                        sbj,qry = d[8],d[9]
                        if sbj not in Clusters: Clusters[sbj] = int(d[1])
                        if qry not in Clusters: Clusters[qry] = int(d[1])
        Temp = {}
        for c in Clusters:
                i = Clusters[c]
                if i not in Temp: Temp[i] = 1
                else: Temp[i] += 1
        for c in Clusters.keys(): Clusters[c] = Temp[Clusters[c]]
        return Clusters

