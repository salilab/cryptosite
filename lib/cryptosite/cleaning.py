from __future__ import print_function, absolute_import
import os
import re
import subprocess
from operator import itemgetter

class InvalidResiduesError(Exception):
    pass

def get_pdb_seq(pdb, chain):
    '''
    Read a PDB[pdb] input and outputs a sequence of a desired chain[chain].
    '''
    import modeller

    e = modeller.environ()
    m = modeller.model(e, file=pdb,
                       model_segment=('FIRST:'+chain, 'LAST:'+chain))
    # BLAST will get confused if the PDB file contains non-standard ATOM
    # records (e.g. HIE rather than HIS)
    _detect_invalid_residue_types(m)
    return "".join(r.code for r in m.residues)

def _detect_invalid_residue_types(m):
    """Raise an error if the Modeller model contains any invalid restypes."""
    bad_res_types = {}
    for r in m.residues:
        if r.code == '.': # BLK residue
            bad_res_types[r.pdb_name] = None
    if bad_res_types:
        raise InvalidResiduesError(
               "Your PDB file contains the following invalid residue "
               "types: %s. Please either remove them or convert them to the "
               "closest equivalent standard PDB residue type."
               % ", ".join(bad_res_types.keys()))

def muscleAlign(qSeq, sSeq, pdb, chain):
    '''
    Align two sequences and produce Modeller-compatible output
    '''

    # todo: put these files in a temporary directory
    # --- write sequence file
    with open('sequences.seq', 'w') as output:
        output.write('>%s%s\n' % (pdb,chain))
        output.write(qSeq+'\n')
        output.write('>%s%spdb\n' % (pdb,chain))
        output.write(sSeq)

    # --- align using Muscle
    cmd = ["muscle", "-in", "sequences.seq", "-out", "alignment.ali"]
    print(cmd)
    subprocess.check_call(cmd)

    with open('alignment.ali') as data:
        D = data.read().split('>')

    strc = D[1].split('\n')
    strcid, strcsq = strc[0], ''.join(strc[1:])
    seq = D[2].split('\n')
    seqid, seqsq = seq[0], ''.join(seq[1:])

    # --- clean
    os.unlink('alignment.ali')
    os.unlink('sequences.seq')

    return (strcsq, seqsq)


def get_gaps(alnfile):
    """Get a list of all gaps in the given alignment file."""

    data = open(alnfile)
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
    chains = [chr(x) for x in range(65, 65+len(strcsq))]
    gap_reg = re.compile('\W-*\W')
    for strcsq_index, strcsq_part in enumerate(strcsq):
        iterator = gap_reg.finditer(strcsq_part)
        for match in iterator:
            if len(chains)>1:
                gaps.append(','.join([str(match.span()[0]+L+1)+':'+chains[strcsq_index], str(match.span()[1]+L)+':'+chains[strcsq_index]]))
            else:
                gaps.append(','.join([str(match.span()[0]+L+1)+':', str(match.span()[1]+L)+':']))
        L+=len(seqsq[strcsq_index])

    return gaps

def build_model(pdb,chains):
    '''
    Build model using Modeller, treating residues in the structure as rigid, and
    loop modeling for the rest.
    '''
    import modeller
    from modeller.automodel import loopmodel, automodel, assess

    # --- get gaps
    gaps = get_gaps('alignment.pir')
    out = open('gaps_%s.txt' % pdb,'w')
    out.write(pdb+'_mdl\t'+str(gaps))
    out.close()

    # --- set up modeling
    env = modeller.environ()

    env.io.atom_files_directory = ['.', '../atom_files']

    class MyLoop(loopmodel):
        def select_loop_atoms(self):
            gaps = get_gaps('alignment.pir')
            return modeller.selection(self.residue_range(
                              i.split(',')[0], i.split(',')[1]) for i in gaps)
        def special_restraints(self, aln):
            rsr = self.restraints
            wholeSel = modeller.selection(self) - self.select_loop_atoms()
            r = modeller.rigid_body(wholeSel)
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

    if len(chains)==1:
        mdl  = modeller.model(env, file='XXX_mdl')
        mdl.rename_segments(segment_ids='A')
        mdl.write(file='%s_mdl.pdb' % (pdb))
