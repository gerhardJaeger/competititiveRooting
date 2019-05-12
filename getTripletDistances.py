from ete3 import Tree
import os
import subprocess
import tempfile
import numpy as np
import pandas as pd

families = np.unique([x.split('.')[0]
                      for x in os.listdir('families')])

glot = Tree('asjp/glottolog.tre')

for l in glot.get_leaves():
    a, b, c = l.name.split('.')


def gtd(tree, glot):
    treeFile = tempfile.mkstemp()[1]
    glotFile = tempfile.mkstemp()[1]
    tree.write(outfile=treeFile, format=9)
    glot.write(outfile=glotFile, format=9)
    cmdString1 = 'tqDist-1.0.1/bin/triplet_dist -v  '+treeFile+' '+glotFile
    p1 = subprocess.check_output(cmdString1, shell=True)
    cmdString2 = 'tqDist-1.0.1/bin/triplet_dist -v  '+glotFile+' '+glotFile
    p2 = subprocess.check_output(cmdString2, shell=True)
    p = subprocess.Popen('rm -f ' + treeFile, shell=True)
    os.waitpid(p.pid, 0)
    p = subprocess.Popen('rm -f ' + glotFile, shell=True)
    os.waitpid(p.pid, 0)
    r1 = list(map(float, p1.split()))
    r2 = list(map(float, p2.split()))
    nTriplets = r1[1]
    nUnresolved = r2[-2]
    nResolved = nTriplets - nUnresolved
    if nResolved == 0:
        return np.nan
    nAgree = r1[4]
    return 1. - nAgree/nResolved


results = []
for fam in families:
    print(fam)
    fmGlot = glot.copy()
    tree = Tree('madRooted/' + fam + '.madRooted.tre')
    fmGlot.prune(tree.get_leaf_names())
    fmResults = [len(tree)]
    fmResults.append(gtd(Tree('madRooted/' + fam +
                              '.madRooted.tre'), fmGlot))
    fmResults.append(gtd(Tree('midpointRooted/' + fam +
                              '.midpointRooted.tre'), fmGlot))
    fmResults.append(gtd(Tree('outgroupRooted/' + fam +
                              '.outgroupRooted.tre'), fmGlot))
    fmResults.append(gtd(Tree('yuleRooted/' + fam +
                              '.yuleRooted.tre'), fmGlot))
    results.append(fmResults)

results = pd.DataFrame(results,
                       index=families,
                       columns=['size', 'mad', 'midpoint', 'outgroup', 'yule'])

results.to_csv('gtdEvaluation.csv')
