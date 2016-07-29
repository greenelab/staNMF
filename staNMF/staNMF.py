###########################
# Required Pacakges
##########################
import argparse
from timeit import default_timer as timer
import math
import spams
import csv
import numpy as np
import random
import os
import sys
import warnings
import matplotlib.pyplot as plt


class staNMF:
    '''
    (C) 2016 Amy E Campbell

    Python 2.7 implementation of Siqi Wu's 03/2016 Stability NMF (staNMF)

    Solves non-negative matrix factorization over a range of principal patterns
    with randomly sampled initial NMF "guesses" for an n x m matrix using
    the SPArse Modeling Software (J. Mairal, F. Bach, J. Ponce and G. Sapiro,
    2010)

    USAGE:
    Can be called and used by a script, imported as "staNMF"
    (see requirements.txt for required packages)

    To read in Wu et al's 2016 drosophila spatial expressoin data, set
    filename = 'example'

    INPUT:
    :param: filename (str, required) : Specify full file extension to a ".csv"
    containing a table with columns and rows labeled, or set filename='example'

    :param: folderID (str, optional with default ""): allows user to specify
    a unique (to the user's working directory) identifier for the 'staNMFDicts'
    folder that the runNMF() method creates

    :param: K1 (int, optional with default 15): lowest # PP's (K) tested

    :param: K2 (int, optional with default 30): highest # PP's (K) tested

    :param: weighted (bool, optional, default False): performs weighting step
    to account for multiple columns of the same name
    (True if filename='example')

    :param: seed (int, optional with default 123): sets numpy random seed

    :param: replicates(int, optional with default 100): Number of bootstrapped
    repetitions of NMF to be performed on each value of K, for use in stability
    analysis

    :param: finished (bool, optional with default False): True if runNMF has
    been completed for the dataset and k range for which you wish to calculate
    instability. To surpass NMF step if fileID file already contains
    factorization solutions for X in your range (K1, K2], set to True
    '''

    def __init__(self, filename, folderID="", K1=15, K2=30,
                 weighted=False, seed=123, replicates=100, finished=False):
        warnings.filterwarnings("ignore")
        self.K1 = K1
        self.K2 = K2
        self.weighted = weighted
        self.seed = seed
        self.guess = np.array([])
        self.guessdict = {}
        self.replicates = replicates
        self.X = []
        if filename == 'example':
            self.fn = "updatedexpressionpatterns.csv"
            self.weighted = True
        else:
            self.fn = filename
        self.folderID = folderID
        self.rowidmatrix = []
        self.finished = finished
        self.instabilitydict = {}
        self.load()
        self.instabilityarray = []

    def initialguess(self, X, K, i):

        '''
        Randomly samples K columns from X; sets input matrix guess to be a
        fortran array, and sets 'guesslist', a list of the column indices
        sampled from X for this replication of NMF

        Arguments:
        :param X :full matrix
        :param K : Number of columns to select at random

        Usage:
        Called by runNMF
        '''

        random.seed(self.seed)
        indexlist = random.sample(np.arange(1, len(X[1])), K)
        self.guess = np.asfortranarray(X[:, indexlist])
        self.guessdict[i] = indexlist

    def load(self):
        '''
        Loads full matrix from .csv file to

        Usage:
        Called by constructor

        '''
        if not self.finished:

            try:
                csvfile = open(self.fn, "r")
            except:
                print("File ''" + str(self.fn) + "' not found.")
                sys.exit()

            csvreader = csv.reader(csvfile)
            colnames = csvreader.next()[1:]
            workingmatrix = []
            for row in csvreader:
                line = row[1:]
                for w in range(len(line)):
                    line[w] = float(line[w])
                workingmatrix.append(line)
                self.rowidmatrix.append(row[0])

            if self.weighted:
                for g in range(len(colnames)):
                    if "." in str(colnames[g]):
                        index = colnames[g].find(".")
                        string = colnames[g]
                        colnames[g] = string[:index]

                colUnique = np.unique(colnames)
                colNum = np.zeros(len(colUnique))
                weight = np.zeros(len(colnames))

                for i in range(len(colUnique)):
                    colNum[i] = colnames.count(colUnique[i])
                    weight[i] = 1/(colNum[i])
                for j in range(len(workingmatrix)):
                    for i in range(len(colnames)):
                        num = float(workingmatrix[j][i])
                        workingmatrix[j][i] = math.sqrt(float(weight[i]) * num)
            X1 = (np.array(workingmatrix).astype(float))
            self.X = np.asfortranarray(X1)
            print(type(self.X))

    def runNMF(self):
        '''
        Iterates through range of integers between the K1 and K2 provided (By
        default, K1=15 and K2=30)
        '''
        self.finished = False
        numPatterns = np.arange(self.K1, self.K2)
        for k in range(len(numPatterns)):
            K = numPatterns[k]
            path = str("./staNMFDicts" + str(self.folderID) + "/K=" + str(K) +
                       "/")
            try:
                os.makedirs(path)
            except OSError:
                if not (os.path.isdir(path)):
                    raise
            m, n = np.shape(self.X)

            print("Working on " + str(K) + "...\n")

            for l in range(self.replicates):
                print(l)
                self.initialguess(self.X, K, l)
                Dsolution = spams.trainDL(
                    # Matrix
                    self.X,
                    # Initial guess as provided by initialguess()
                    D=self.guess,
                    # All available CPUs/cores
                    numThreads=-1,
                    # minibatch size
                    batchsize=min(1024, n),
                    # Number of columns in solution
                    K=K,
                    lambda1=0,
                    # Number of iterations to go into this round of NMF
                    iter=500,
                    # Specify optimization problem to solve
                    mode=2,
                    # Specify convex set
                    modeD=0,
                    # Positivity constraint on coefficients
                    posAlpha=True,
                    # Positivity constraint on solution
                    posD=True,
                    # Limited information printed about progress
                    verbose=False,
                    gamma1=0)

                # write solution to a csv file in the staNMFDicts/k=K/ folder
                outputfilename = "factorization_" + str(l) + ".csv"
                outputfilepath = os.path.join(path, outputfilename)

                outputfile = open(outputfilepath, "w")
                writer = csv.writer(outputfile)
                for i in range(len(self.rowidmatrix)):
                    dlist = Dsolution[i].tolist()
                    dlist.insert(0, str(self.rowidmatrix[i]).replace("'", ""))

                    writer.writerow(dlist)
                outputfile.close()

            indexoutputstring = "selectedcolumns" + str(K) + ".txt"
            indexoutputpath = os.path.join(path, indexoutputstring)
            indexoutputfile = open(indexoutputpath, "w")
            for m in sorted(self.guessdict):
                indexoutputfile.write(str(m) + '\t' + str(self.guessdict[m]) +
                                      '\n')
            indexoutputfile.close()

            self.finished = True

    def amariMaxError(self, correlation):
        '''
        Computes what Wu et al. (2016) described as a 'amari-type error'
        based on average distance between factorization solutions

        Return:
        Amari distance distM

        Arguments:
        :param: correlation: k by k matrix of pearson correlations

        Usage: Called by instability()
        '''

        n, m = correlation.shape
        maxCol = np.absolute(correlation).max(0)
        colTemp = np.mean((1-maxCol))
        maxRow = np.absolute(correlation).max(1)
        rowTemp = np.mean((1-maxRow))
        distM = (rowTemp + colTemp)/(2)

        return distM

    def findcorrelation(self, A, B, k):
        '''
        Construct k by k matrix of Pearson product-moment correlation
        coefficients for every combination of two columns in A and B

        :param: A : first NMF solution matrix
        :param: B : second NMF solution matrix, of same dimensions as A
        :param: k : number of columns in each matrix A and B

        Return: numpy array of dimensions k by k, where array[a][b] is the
        correlation between column 'a' of X and column 'b'

        Usage:
        Called by instability()

        '''
        corrmatrix = []

        for a in range(k):
            for b in range(k):
                c = np.corrcoef(A[:, a], B[:, b])
                corrmatrix.append(c[0][1])

        return np.asarray(corrmatrix).reshape(k, k)

    def instability(self, k1=0, k2=0):
        '''
        Performs instability calculation for NMF factorizations for each K
        within the range entered

        Arguments:

        :param: k1 (int, optional, default self.K1): lower bound of K to
        plot against stability

        :param: k2 (int, optional, default self.K2): upper bound of K to
        plot against instability

        Return:
        "instability.csv" containing instability index
        for each K between and including k1 and k2; updates
        self.instabilitydict (required for makeplot())
        '''
        if k1 == 0:
            k1 = self.K1
        if k2 == 0:
            k2 = self.K2

        numReplicates = self.replicates

        if self.finished is False:
            print("staNMF Error: runNMF is not complete\n")
        else:
            numPatterns = np.arange(k1, k2+1)

            modelK = numPatterns[0]
            path = str("./staNMFDicts" + str(self.folderID) + "/K=" +
                       str(modelK)+"/")
            inputfilename = "factorization_0.csv"
            inputfilepath = os.path.join(path, inputfilename)
            inputfile = open(inputfilepath, "rb")
            reader = csv.reader(inputfile, delimiter=',')
            matrix1 = np.array(list(reader))
            firstmatrix = matrix1[:, 1:]

            inputfile.close()
            d = np.size(firstmatrix, 0)
            for k in numPatterns:
                print("Calculating instability for " + str(k))
                path = str("./staNMFDicts" + str(self.folderID) + "/K=" +
                           str(k)+"/")
                Dhat = np.zeros((numReplicates, d, k))

                for L in range(numReplicates):
                    inputfilename = "factorization_" + str(L) + ".csv"
                    inputfilepath = os.path.join(path, inputfilename)
                    inputfile = open(inputfilepath, "rb")
                    reader = csv.reader(inputfile, delimiter=',')
                    inputmatrix = np.array(list(reader))
                    replicatematrix = inputmatrix[:, 1:]
                    Dhat[L] = replicatematrix
                    inputfile.close()

                distMat = np.zeros(shape=(numReplicates, numReplicates))

                for i in range(numReplicates):
                    for j in range(i, numReplicates):
                        x = Dhat[i]
                        y = Dhat[j]

                        CORR = self.findcorrelation(x, y, k)
                        distMat[i][j] = self.amariMaxError(CORR)
                        distMat[j][i] = distMat[i][j]

                self.instabilitydict[k] = (np.sum(distMat) / (numReplicates *
                                           (numReplicates-1)))
                outputfile = open("instability.csv", "w")
                for i in sorted(self.instabilitydict):
                    outputfile.write(str(i) + "," +
                                     str(self.instabilitydict[i]) + "\n")

    def plot(self, dataset_title="Drosophila Spatial Expression Data", xmax=0,
             xmin=-1, ymin=0, ymax=0, xlab="K", ylab="Instability Index"):
        '''
        Plots instability results for all K's between and including K1 and K2
        with K on the X axis and instability on the Y axis

        :param: dataset_title (str, optional, default "Drosophila
        Expression Data")

        :param: ymax (float, optional,  default
        largest Y + (largestY/# of points)

        :param: xmax (float, optional, default K2+1)

        :param: xlab (string, default "K") x-axis label

        :param: ylab (string, default "Instability Index") y-axis label

        Return:
        saves plot as <dataset_title>.png
        '''
        kArray = []
        for i in sorted(self.instabilitydict):
            kArray.append(i)
            self.instabilityarray.append(self.instabilitydict[i])
        if xmax == 0:
            xmax = self.K2 + 1
        if xmin == -1:
            xmin = self.K1
        ymin = 0
        ymax = max(self.instabilityarray) + (max(self.instabilityarray) /
                                             len(self.instabilityarray))
        plt.plot(kArray, self.instabilityarray)
        plt.axis([xmin, xmax, ymin, ymax])
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(str('Stability NMF Results: Principal Patterns vs.\
                       Instability in ' + dataset_title))
        plotname = str(dataset_title + ".png")
        plt.savefig(plotname)
