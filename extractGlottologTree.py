from ete3 import Tree
import re
import pandas as pd
import numpy as np


with open('tree_glottolog_newick.txt') as f:
    raw = f.readlines()

trees = []

for i, ln in enumerate(raw):
    ln = ln.strip()
    ln = re.sub(r"\'[A-Z][^[]*\[", "[", ln)
    ln = re.sub(r"\][^']*\'", "]", ln)
    ln = re.sub(r"\[|\]", "", ln)
    ln = ln.replace(":1", "")
    trees.append(Tree(ln, format=1))


glot = Tree()
for t in trees:
    glot.add_child(t)


nonLeaves = [nd.name for nd in glot.traverse()
             if nd.name != '' and not nd.is_leaf()]

for i, nm in enumerate(nonLeaves):
    if i % 100 == 0:
        print(i)
    nd = glot & nm
    nd.name = ''
    nd.add_child(name=nm)

asjp = pd.read_csv('asjp/dataset.tab', sep='\t')

asjpLanguages = pd.read_csv('asjp/languages.csv')

glottocodes = np.unique([x for x in glot.get_leaf_names()
                         if x in asjpLanguages.Glottocode.values])

glot.prune([glot & l for l in glottocodes])

asjpLanguages = asjpLanguages[asjpLanguages.Glottocode.isin(glottocodes)]

asjp = asjp[asjp.names.isin(asjpLanguages.ID)]

asjp['longname'] = ['.'.join(x).replace('-', '_')
                    for x in asjp[['wls_fam', 'wls_gen', 'names']].values]

asjpLanguages = pd.merge(asjpLanguages, asjp[['names', 'longname']],
                         left_on='ID', right_on='names')


for l in asjpLanguages.longname.unique():
    gc = asjpLanguages[asjpLanguages.longname == l].Glottocode.values[0]
    (glot & gc).add_child(name=l)

for nd in glot.traverse():
    if len(nd.get_children()) == 1:
        nd.delete()

glot.write(outfile='asjp/glottolog.tre', format=9)
