import pandas as pd
import numpy as np
from ete3 import Tree
import subprocess
import os

raxmlCommand = 'raxmlHPC-PTHREADS-SSE3'

glot = Tree('asjp/glottolog.tre')

for l in glot.get_leaves():
    a, b, c = l.name.split('.')
    l.name = a+'.'+b.upper()+'.'+c


ccData = pd.read_csv('asjp/asjp18Clustered.csv')

ccData['longname'] = ['.'.join(x).replace('-', '_')
                      for x in ccData[['wls_fam', 'wls_gen',
                                       'doculect']].values]

ccData = ccData[ccData.longname.isin(glot.get_leaf_names())]


concepts = ccData.concept.unique()


with open('asjp/world_cc.phy') as f:
    raw = np.array([x.strip() for x in f.readlines()[1:]])
    ccMtx = pd.DataFrame([list(x.split()[1]) for x in raw],
                         index=[x.split()[0] for x in raw])


with open('asjp/world_sc.phy') as f:
    raw = np.array([x.strip() for x in f.readlines()[1:]])
    scMtx = pd.DataFrame([list(x.split()[1]) for x in raw],
                         index=[x.split()[0] for x in raw])


asjp = ccData[['longname', 'glot_fam']].drop_duplicates()

asjp.index = asjp.longname.values

db = pd.read_csv('asjp/dataset.tab', sep='\t')

db['longname'] = ['.'.join(x).replace('-', '_')
                  for x in db[['wls_fam', 'wls_gen', 'names']].values]

db = db[db.longname.isin(asjp.index)]
db = db[db['pop'] > 0]

asjp = asjp[asjp.index.isin(db.longname.values)]

ccData = ccData[ccData.longname.isin(db.longname.values)]

ccMtx = ccMtx[ccMtx.index.isin(db.longname.values)]
scMtx = scMtx[scMtx.index.isin(db.longname.values)]

taxa = ccMtx.index

fam10 = np.array(pd.value_counts(asjp.glot_fam)[
    pd.value_counts(asjp.glot_fam) >= 10].index)


ccMtx.loc['outgroup'] = '0'
scMtx.loc['outgroup'] = '0'

for fam in fam10:
    print(fam)
    fTaxa = ccData[ccData.glot_fam == fam].longname.unique()
    outgroup = 'outgroup'
    fTaxaO = np.r_[[outgroup], fTaxa]
    fCcMtx = ccMtx.loc[fTaxaO].copy()
    fScMtx = scMtx.loc[fTaxaO].copy()
    pad = max(map(len, fTaxaO))+2
    x = len(fCcMtx.columns)
    y = len(fScMtx.columns)
    with open('families/'+fam+'.part.txt', 'w') as f:
        f.write('BINX, cc=1-'+str(x)+'\n')
        f.write('BINX, sc='+str(x+1)+'-'+str(x+y)+'\n')
    fCcScMtx = pd.concat([fCcMtx, fScMtx], axis=1)
    with open('families/'+fam+'.cc_sc.phy', 'w') as f:
        f.write(' '.join(map(str, np.shape(fCcScMtx)))+'\n')
        for l in fCcScMtx.index:
            f.write(l.ljust(pad))
            rw = np.array(fCcScMtx.loc[l].values, str)
            rw[rw == '-1'] = '-'
            f.write(''.join(rw)+'\n')
    cmd = raxmlCommand + " -T 20 -g " + fam
    cmd += ".tre -c 4 -m BINGAMMAX -s " + fam
    cmd += ".cc_sc.phy -q " + fam + ".part.txt -n " + fam + "_rooted -p 12345"
    p = subprocess.Popen(cmd,
                         shell=True,
                         cwd=os.getcwd()+'/families/')
    os.waitpid(p.pid, 0)

    fTree = Tree('families/RAxML_bestTree.'+fam+'_rooted')
    fTree.resolve_polytomy()
    fTree.set_outgroup(fTree & outgroup)
    (fTree & outgroup).delete(preserve_branch_length=True)
    p = subprocess.Popen("rm -f RAxML*",
                         shell=True,
                         cwd=os.getcwd()+'/families/')
    os.waitpid(p.pid, 0)
    p = subprocess.Popen("rm -f " + fam + ".cc_sc.phy* " + fam + ".part.txt",
                         shell=True,
                         cwd=os.getcwd()+'/families/')
    os.waitpid(p.pid, 0)
    fTree.get_children()[0].write(
        outfile="outgroupRooted/"+fam+".outgroupRooted.tre", format=1)
