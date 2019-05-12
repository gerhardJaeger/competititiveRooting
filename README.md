# competititiveRooting
### comparing various rooting algorithms for phylogenetic language trees

To replicate my analysis, you need python3 to be installed. I assume the executable `python` to start python3. Also, an executable for RAxML (https://cme.h-its.org/exelixis/web/software/raxml/index.html) needs to be in your search path.

The file `asjp/asjp18Clustered.csv`was obtained by applying the procedure described in [Jäger, G., 2018. Global-scale phylogenetic linguistic inference from lexical resources, Scientific Data 5, 180189, 2018, doi.org/10.1038/sdata.2018.189.](https://www.nature.com/articles/sdata2018189) to the data from [ASJP v. 18](https://asjp.clld.org/). 

To install the required python packages:

`pip install -r requirements.txt`

To obtain the constraint tree from Glottolog:

`wget https://cdstar.shh.mpg.de/bitstreams/EAEA0-9478-C22F-4AAF-0/tree_glottolog_newick.txt`

To obtain the ASJP data:

`cd asjp
wget https://cdstar.shh.mpg.de/bitstreams/EAEA0-E32A-2C2D-B777-0/asjp_dataset.cldf.zip
unzip asjp_dataset.cldf.zip
wget https://cdstar.shh.mpg.de/bitstreams/EAEA0-E32A-2C2D-B777-0/asjp_dataset.tab.zip
unzip unzip asjp_dataset.tab.zip
cd ..`

You need to install the program tqDist:

`wget http://users-cs.au.dk/cstorm/software/tqdist/files/tqDist-1.0.1.zip
unzip tqDist-1.0.1.zip
cd tqDist-1.0.1`

Follow the "Build procedure" described in `README`. I got error messages when I ran "make", but adding the lines

   `#include <string.h>
  #include <stdio.h>`

to the headers of the files `all_pairs_quartet_distance.cpp, pairs_triplet_distance.cpp, test_triplet.cpp, pairs_quartet_distance.cpp, quartet_dist.cpp, triplet_dist.cpp` in the directory `tqDist` solved the problem.

`cd ..`

Finally, you need the script `mad.py` by the authors of Tria, F. D. K., Landan, G. and Dagan, T. Nat. Ecol. Evol. 1, 0193 (2017):

`wget https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip
unzip mad2-2.zip`

The actual analysis is carried out by the commands

`python extractGlottologTree.py
python inferTrees.py
python madRooting.py
python outgroupRooting.py
python midpointRooting.py
python yuleRooting.py
python getTripletDistances.py`

