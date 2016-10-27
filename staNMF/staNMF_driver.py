import staNMF
import sys
import argparse
'''
staNMF Driver Script

Function:
This script is not required, but offers one way to use the three primary user
functions of staNMF in parallel. It can be called from the command or from a
shell script using:

$ python staNMF_driver.py <k1> <k2> <reps1> <reps2> <folder> <filename> <function>

:param: k1 (int) - lower bound(inclusive) for K calculations for this call

:param: k2 (int) - upper bound(inclusive) for K calculations in this call; note
that if you are breaking one K's calculation into multiple calls, k1 and k2 can
be the same integer

:param: reps1 (int) - lower bound(inclusive) for the NMF factorizations
generated in this call

:param: reps2 (int) -upper bound (not inclusive) for the NMF factorizations
generated in this call --NOTE: if you are calling instability or plot,
range(reps1, reps2) should include all the existing factorizations you wish to
input into the instability or plot function performed in this call

:param: folder (str) - folder ID for NMF factorizations that are either to be
generated(by runNMF) or are to be used by instabilty/plot

:param: filename (str) - .csv file containing the full dataset

;param: function (str) - the staNMF function to be run in this call; cannot be
more than one and must be 'runNMF', 'instability', or 'plot'

'''
parser = argparse.ArgumentParser()

parser.add_argument("k1", type=int)
parser.add_argument("k2", type=int)
parser.add_argument("reps1", type=int)
parser.add_argument("reps2", type=int)
parser.add_argument("folder", type=str)
parser.add_argument("filename", type=str)
parser.add_argument("function", choices=['runNMF', 'instability', 'plot'])
args = parser.parse_args()

if args.function == 'runNMF':
    staNMFobj = staNMF.staNMF(filename=args.filename, K1=args.k1, K2=args.k2,
                              replicates=(args.reps1, args.reps2),
                              folderID=args.folder, parallel=True)
    staNMFobj.runNMF()

elif args.function == 'instability':
    staNMFobj = staNMF.staNMF(filename=args.filename, K1=args.k1,
                              K2=args.k2, replicates=(args.reps1, args.reps2),
                              folderID=args.folder, NMF_finished=True,
                              parallel=True)
    staNMFobj.instability()

elif args.function == 'plot':
    staNMFobj = staNMF.staNMF(filename="PanCanMatrix.csv", K1=k1, K2=k2,
                              replicates1=reps1, replicates2=reps2,
                              folderID=folder, NMF_finished=True,
                              parallel=True)
    staNMFobj.plot(done=True)
