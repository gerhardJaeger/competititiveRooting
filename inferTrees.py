import pandas as pd
import numpy as np
from ete3 import Tree
import subprocess
import os


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


try:
    os.stat('families')
except os.error:
    os.mkdir('families')

for fam in fam10:
    print(fam)
    p = subprocess.Popen('rm -f *'+fam+'*',
                         shell=True,
                         cwd=os.getcwd()+'/families/')
    os.waitpid(p.pid, 0)
    fTaxa = ccData[ccData.glot_fam == fam].longname.unique()
    fCcMtx = ccMtx.loc[fTaxa].copy()
    fScMtx = scMtx.loc[fTaxa].copy()
    pad = max(map(len, fTaxa))+2
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
    fGlot = glot.copy()
    fGlot.prune([fGlot & l for l in fCcMtx.index])
    with open('families/glot'+fam+'.tre', 'w') as f:
        f.write(fGlot.write(format=9))
#    cmd = 'raxmlHPC-PTHREADS-SSE3 -c 4 -m BINGAMMAX -s '
    cmd = 'raxml -c 4 -m BINGAMMAX -s '
    cmd += fam
    cmd += '.cc_sc.phy -q '
    cmd += fam
    cmd += '.part.txt -n '+fam+' -g glot'+fam+'.tre -T 20 -p 12345 > /dev/null'
    p = subprocess.Popen(cmd,
                         shell=True,
                         cwd=os.getcwd()+'/families/')
    os.waitpid(p.pid, 0)
    fTree = Tree('families/RAxML_bestTree.'+fam)
    p = subprocess.Popen('rm -f *'+fam+'*',
                         shell=True,
                         cwd=os.getcwd()+'/families/')
    os.waitpid(p.pid, 0)
    fTree.write(format=1, outfile='families/'+fam+'.tre')
