
#!/usr/bin/env python

###########################
# Required Pacakges
##########################
import math
import random
import os
import shutil
import sys
import warnings
import argparse
import collections
import csv
from timeit import default_timer as timer

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import spams


class staNMF:
    '''
    2016 Amy E Campbell

    Python 2.7 implementation of Siqi Wu's 03/2016 Stability NMF (staNMF)

    Solves non-negative matrix factorization over a range of principal patterns
    (PPs) with randomly sampled initial NMF "guesses" for an n x m matrix using
    the SPArse Modeling Software (J. Mairal, F. Bach, J. Ponce and G. Sapiro,
    2010)

    USAGE:
    Can be called and used by a script, imported as "staNMF"
    (see requirements.txt for required packages)

    To read in Wu et al's 2016 drosophila spatial expression data, set:
    filename = 'example'

    INPUT:
    :param: filename (str, required) : Specify full file extension to a ".csv"
    containing a table with columns and rows labeled, or set filename='example'

    :param: folderID (str, optional with default ""): allows user to specify
    a unique (to the user's working directory) identifier for the 'staNMFDicts'
    folder that the runNMF() method creates

    :param: K1 (int, optional with default 15): lowest # PP's (K) tested

    :param: K2 (int, optional with default 30): highest # PP's (K) tested

    :param: sample_weights (bool or list, optional, default False):
    performs weighting step on full data matrix to account for multiple columns
    of the same name (filename='example'defaults to special case of weighting
    (delimited by ".")) If sample_weights is a list of custom weights, it must
    be equal in length to the number of columns in the matrix, and this list of
    weights will be applied across columns

    :param: seed (int, optional with default 123): sets numpy random seed

    :param: replicates(int or tuple of ints of length 2, optional with default
    int 100):
    Specifies which/ how many bootstrapped repetitions to be performed on each
    value of K, for use in stability analysis; if a list of length 2 is given,
    self.replicates is set to a list of ints between the first and second
    elements of this tuple. If it is set to an integer, self.replicates is set
    to a list of ints between 0 and the given integer

    Number of bootstrapped
    repetitions of NMF to be performed on each value of K, for use in stability
    analysis

    :param: NMF_finished (bool, optional with default False): True if runNMF
    has been completed for the dataset and k range for which you wish to
    calculate instability. To surpass NMF step if fileID file already contains
    factorization solutions for X in your range [K1, K2], set to True

    :param: parallel (bool, optional with default False): True if NMF is to be
    run in parallel such that the instability calculation should write a file
    for each K containing its instability index

    '''

    def __init__(self, filename, folderID="", K1=15, K2=30,
                 sample_weights=False, seed=123, replicates=100,
                 NMF_finished=False, parallel=False):
        warnings.filterwarnings("ignore")
        self.K1 = K1
        self.K2 = K2
        self.sample_weights = sample_weights
        self.seed = seed
        self.guess = np.array([])
        self.guessdict = {}
        self.parallel = parallel
        if isinstance(replicates, int):
            self.replicates = range(replicates)
        elif isinstance(replicates, tuple):
            start, stop = replicates
            self.replicates = range(replicates[0], replicates[1])
        self.X = []
        if filename == 'example':
            self.fn = os.path.join("data", "WuExampleExpression.csv")
            self.sample_weights = True
        else:
            self.fn = filename
        self.folderID = folderID
        self.rowidmatrix = []
        self.NMF_finished = NMF_finished
        self.instabilitydict = {}
        self.load_data()
        self.instabilityarray = []
        self.stability_finished = False
        random.seed(self.seed)

    def initialguess(self, X, K, i):

        '''
        Randomly samples K columns from X; sets input matrix guess to be a
        fortran array, and sets 'guesslist', a list of the column indices
        sampled from X for this replication of NMF;

        Arguments:
        :param: X(numpy array, required): full matrix
        :param: K(int, required): Number of columns to select at random to be
        used as the 'initial guess' for the K PPs in the current simulation
        of NMF
        :param: i(int, required): Key at which indexlist will be stored(current
        replicate of NMF)

        Usage:
        Called by runNMF
        '''

        indexlist = random.sample(np.arange(1, X.shape[1]), K)
        self.guess = np.asfortranarray(X[:, indexlist])
        self.guessdict[i] = indexlist

    def load_data(self):
        '''
        Loads full data matrix from .csv file into numpy array; if the
        self.sample_weights variable is True, weights column names by their
        number of replicate occurances

        Usage:
        Called by constructor

        '''
        if not self.NMF_finished:
            csvfile = open(self.fn, "r")
            workingmatrix = pd.read_csv(csvfile, index_col=0)
            self.rowidmatrix = workingmatrix.index.values
            colnames = workingmatrix.columns.values

            if self.sample_weights is not False:
                if isinstance(self.sample_weights, list):
                    if len(self.sample_weights) != len(colnames):
                        raise ValueError("sample_weights length must equal the"
                                         " number of columns.")
                    else:
                        weight = self.sample_weights
                else:
                    # Special formatting case for Wu et al. expression data
                    if self.fn == os.path.join("data",
                                               "WuExampleExpression.csv"):
                        colnames = [(str(x).split('.'))[0] for x in colnames]

                    colUnique = np.unique(colnames)
                    colNum = np.zeros(len(colUnique))
                    weight = np.zeros(len(colnames))

                    for i in range(len(colUnique)):
                        colNum[i] = list(colnames).count(colUnique[i])
                        weight[i] = 1/(colNum[i])

                workingmatrix = workingmatrix.apply(lambda x: weight * x,
                                                    axis=1)
                workingmatrix = workingmatrix.applymap(lambda x: math.sqrt(x))

            X1 = (np.array(workingmatrix).astype(float))
            self.X = np.asfortranarray(X1)

    def runNMF(self, **kwargs):
        '''
        Iterates through range of integers between the K1 and K2 provided (By
        default, K1=15 and K2=30), runs NMF using SPAMS package; outputs
        NMF matrix files (.csv form) and updates self.guessdict containing the
        columns selected for the initial guess input(as calculated by
        staNMF.initialguess())

        Usage: Called by user (ex: '$ instance.runNMF()')

        Arguments: Optional **kwargs allows user to update spams.trainDL()
        parameters

        Return: None

        Output:
        (k2-k1) folders, each containing files for every replicate
        (labeled factorization_<replicate>.csv) , and each containing
        a 'selectedcolumns.txt' file, which prints 'self.guessdict', a
        dictionary with keys <factorzation #>, values <columns selected>

        '''

        self.NMF_finished = False
        numPatterns = np.arange(self.K1, self.K2+1)
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

            param = {"numThreads": -1,
                     # minibatch size
                     "batchsize": min(1024, n),
                     # Number of columns in solution
                     "K": K,
                     "lambda1": 0,
                     # Number of iterations to go into this round of NMF
                     "iter": 500,
                     # Specify optimization problem to solve
                     "mode": 2,
                     # Specify convex set
                     "modeD": 0,
                     # Positivity constraint on coefficients
                     "posAlpha": True,
                     # Positivity constraint on solution
                     "posD": True,
                     # Limited information about progress
                     "verbose": False,
                     "gamma1": 0}

            for p in param:
                if p not in kwargs:
                    kwargs[p] = param[p]

            for l in self.replicates:
                self.initialguess(self.X, K, l)
                Dsolution = spams.trainDL(
                    # Matrix
                    self.X,
                    # Initial guess as provided by initialguess()
                    D=self.guess,
                    **kwargs)

                # write solution to a csv file in the staNMFDicts/k=K/ folder
                outputfilename = "factorization_" + str(l) + ".csv"
                outputfilepath = os.path.join(path, outputfilename)

                Dsolution1 = pd.DataFrame(Dsolution, index=self.rowidmatrix)
                Dsolution1.to_csv(outputfilepath, header=None)

            indexoutputstring = "selectedcolumns" + str(K) + ".csv"
            indexoutputpath = os.path.join(path, indexoutputstring)

            with open(indexoutputpath, "w") as indexoutputfile:
                for m in sorted(self.guessdict):
                    indexoutputfile.write(str(m) + '\t' +
                                          str(self.guessdict[m]) + '\n')

            self.NMF_finished = True

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

        numReplicates = len(self.replicates)

        if self.NMF_finished is False:
            ("staNMF Error: runNMF is not complete\n")
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
                ("Calculating instability for " + str(k))
                path = str("./staNMFDicts" + str(self.folderID) + "/K=" +
                           str(k)+"/")
                Dhat = np.zeros((numReplicates, d, k))

                for replicate in range(numReplicates):
                    inputfilename = "factorization_" + str(replicate) + ".csv"
                    inputfilepath = os.path.join(path, inputfilename)
                    with open(inputfilepath, "rb") as inputfile:
                        matrix1 = pd.read_csv(inputfile, header=None)
                        inputmatrix = matrix1.drop(0, axis=1)
                        inputmatrix.columns = np.arange(0, matrix1.shape[1]-1)

                    Dhat[replicate] = inputmatrix

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

                if self.parallel:
                    outputfile = open(str(path + "instability.csv"), "w")
                    outputfile.write("\n{},{}".format(
                            k, self.instabilitydict[k]))
                    output.close()
        if not self.parallel:
            outputfile = open("instability.csv", "w")
            outputwriter = csv.writer(outputfile)
            for i in sorted(self.instabilitydict):
                outputwriter.writerow([i], [self.instabilitydict[i]])

    def get_instability(self):
        '''
        Retrieves instability values calculated in this instance of staNMF

        Returns:
        dictionary with keys K, values instability index

        Usage: Called by user (not required for output of 'instablity.csv', but
        returns usable python dictionary of these calculations)
        '''
        if self.stability_finished:
            return self.instabilitydict
        else:
            print("Instability has not yet been calculated for your NMF"
                  "results. Use staNMF.instability() to continue.")

    def plot(self, dataset_title="Drosophila Spatial Expression Data", xmax=0,
             xmin=-1, ymin=0, ymax=0, xlab="K", ylab="Instability Index"):
        '''
        Plots instability results for all K's between and including K1 and K2
        with K on the X axis and instability on the Y axis

        Arguments:

        :param: dataset_title (str, optional, default "Drosophila
        Expression Data")

        :param: ymax (float, optional,  default
        largest Y + (largest Y/ # of points)

        :param: xmax (float, optional, default K2+1)

        :param: xlab (string, default "K") x-axis label

        :param: ylab (string, default "Instability Index") y-axis label

        Returns: None, saves plot as <dataset_title>.png

        Usage: Called by user to generate plot
        '''
        kArray = []

        if self.parallel:
            for K in range(x1, x2+1):
                kpath = "./staNMFDicts{}/K={}/instability.csv".format(
                                                      self.folderID, K)
                df = pd.read_csv(kpath)
                kArray.append(int(df.columns[0]))
                self.instabilityarray.append(float(df.columns[1]))
        else:
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
        plt.axes.titlesize = 'smaller'
        plt.title(str('Stability NMF Results: Principal Patterns vs.'
                      'Instability in ' + dataset_title))
        plotname = str(dataset_title + ".png")
        plt.savefig(plotname)

    def ClearDirectory(self, k_list):
        '''
        A storage-saving option that clears the entire directory of each K
        requested, including the instability.csv file in each folder

        :param: k_list (list, required) list of K's to delete corresponding
        directories of

        NOTE: this should only be used after stability has been calculated for
        each K you wish to delete.
        '''
        for K in k_list:
            path = "./staNMFDicts{}/K={}/".format(self.folderID, K)
            shutil.rmtree(path)
