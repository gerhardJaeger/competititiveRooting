from ete3 import Tree
import numpy as np
import pandas as pd
import os


def yuleRooting(tree):
    t = tree.copy()
    taxa = np.array(t.get_leaf_names())
    t.set_outgroup(t & taxa[-1])
    t.set_outgroup(t & taxa[0])
    nodes = [nd.copy() for nd in t.traverse()
             if not nd.is_root()]

    def yuleLL(tree):
        return sum([-np.log(len(nd)-1) for nd in tree.traverse()
                    if not nd.is_leaf() and not nd.is_root()])
    loglikelihoods = []
    for nd in nodes:
        t.set_outgroup(t & taxa[-1])
        t.set_outgroup(t & taxa[0])
        ndTaxa = set(nd.get_leaf_names())
        if ndTaxa != set([taxa[0]]) & ndTaxa != set(taxa[1:]):
            if len(ndTaxa) == 1:
                ndNode = t & list(ndTaxa)[0]
            else:
                ndNode = t.get_common_ancestor([t & x for x in ndTaxa])
            t.set_outgroup(ndNode)
        loglikelihoods.append(yuleLL(t))
    loglikelihoods = pd.Series(loglikelihoods)
    nd = nodes[loglikelihoods.idxmax()]
    ndTaxa = set(nd.get_leaf_names())
    t.set_outgroup(t & taxa[-1])
    t.set_outgroup(t & taxa[0])
    if ndTaxa != set([taxa[0]]) & ndTaxa != set(taxa[1:]):
        if len(ndTaxa) == 1:
            ndNode = t & list(ndTaxa)[0]
        else:
            ndNode = t.get_common_ancestor([t & x for x in ndTaxa])
        t.set_outgroup(ndNode)
    return t


families = np.array([x.split('.')[0]
                     for x in os.listdir('families')])


try:
    os.stat('yuleRooted')
except os.error:
    os.mkdir('yuleRooted')


for fam in families:
    print(fam)
    tree = Tree('families/'+fam+'.tre')
    rooted = yuleRooting(tree)
    rooted.write(outfile='yuleRooted/'+fam+'.yuleRooted.tre',
                 format=1)
