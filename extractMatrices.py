import pandas as pd
import numpy as np


def writePhy(mtx, fn):
    pad = max(map(len, mtx.index))+5
    with open(fn, 'w') as f:
        f.write(' '.join(map(str, np.shape(mtx)))+'\n')
        for l in mtx.index:
            f.write(l.ljust(pad))
            f.write(''.join(map(str, mtx.loc[l].values))+'\n')


data = pd.read_csv('asjp/asjp18Clustered.csv',
                   na_filter=False, dtype=str)

data['longname'] = ['.'.join(x).replace('-', '_') for x in
                    data[['wls_fam', 'wls_gen', 'doculect']].values]

taxa = data.longname.unique()

concepts = data.concept.unique()


ccMtx = pd.DataFrame(index=taxa)
for c in concepts:
    cData = data[data.concept == c]
    cMtx = pd.crosstab(cData.longname, cData.cClass)
    cMtx[cMtx > 1] = 1
    cMtx = cMtx.reindex(taxa, fill_value='-')
    ccMtx = pd.concat([ccMtx, cMtx], axis=1)


writePhy(ccMtx, 'asjp/world_cc.phy')

sounds = np.unique(np.concatenate(list(map(list, data.simplified.values))))

scMtx = pd.DataFrame(index=taxa)
for c in concepts:
    print(c)
    cData = data[data.concept == c]
    cTaxa = cData.longname.unique()
    cWords = pd.Series([''.join(cData[cData.longname == l].simplified.values)
                        for l in cTaxa],
                       index=cTaxa)
    cMtx = pd.DataFrame([[int(s in cWords[l]) for s in sounds]
                         for l in cTaxa],
                        index=cTaxa,
                        columns=[c+':'+s
                                 for s in sounds]).reindex(taxa,
                                                           fill_value='-')
    scMtx = pd.concat([scMtx, cMtx], axis=1)


writePhy(scMtx, 'asjp/world_sc.phy')
