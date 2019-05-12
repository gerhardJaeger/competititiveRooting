from ete3 import Tree
import numpy as np
import os



def midpointRooting(tree):
    rooted = tree.copy()
    rooted.set_outgroup(rooted.get_leaves()[0])
    rooted.set_outgroup(rooted.get_midpoint_outgroup())
    return rooted


families = np.unique([x.split('.')[0]
                      for x in os.listdir('families')])


try:
    os.stat('midpointRooted')
except os.error:
    os.mkdir('midpointRooted')


for fam in families:
    print(fam)
    tree = Tree('families/'+fam+'.tre')
    rooted = midpointRooting(tree)
    rooted.write(outfile='midpointRooted/'+fam+'.midpointRooted.tre',
                 format=1)
