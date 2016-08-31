Amy Campbell, 2016
----------

=======
staNMF
=======
Python 2.7 implementation of Siqi Wu et al.'s 2016 stability NMF(staNMF) for
selecting the number of NMF clusters(K) that results in the most stable output.

=======
Package Contents:
=======

staNMF.py
----------
staNMF.py class includes the necessary methods to perform stability
NMF on a user-specified .csv datase (see sourcecode/staNMF_Example for precise
usage instructions)

staNMF_Example.py
----------
Example of staNMF demonstrated on Wu et al.'s 2016
drosophila spatial expression data between K=15 and K=30; Generates
sample factorizations, calculates instability index, and plots instability
against K

data/WuExampleExpression.csv
----------
sample dataset

=======
Installation*:
=======

pip install staNMF

*Requires SPAMs package(version 2.5), available from Julien Mairal et al. at
http://spams-devel.gforge.inria.fr/downloads.html
