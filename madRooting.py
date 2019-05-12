from ete3 import Tree
import numpy as np
import os
import tempfile
import subprocess


def mad(tree):
    treeFile = tempfile.mkstemp()[1]
    tree.write(outfile=treeFile, format=1)
    cmdString = 'python3 mad/mad.py ' + treeFile + ' > /dev/null'
    p = subprocess.Popen(cmdString, shell=True)
    os.waitpid(p.pid, 0)
    with open(treeFile + '.rooted') as f:
        rerooted = Tree(f.readlines()[0])
    p = subprocess.Popen('rm -f ' + treeFile + '*', shell=True)
    os.waitpid(p.pid, 0)
    return rerooted


families = np.unique([x.split('.')[0]
                      for x in os.listdir('families')])


try:
    os.stat('madRooted')
except os.error:
    os.mkdir('madRooted')


for fam in families:
    print(fam)
    tree = Tree('families/'+fam+'.tre')
    rooted = mad(tree)
    rooted.write(outfile='madRooted/'+fam+'.madRooted.tre',
                 format=1)
