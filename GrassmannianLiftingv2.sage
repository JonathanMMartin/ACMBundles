import numpy as np
import sage.libs.lrcalc.lrcalc as lrcalc

#   TO DO:
#       Correct for twisting issue
#           Keep track of twists?
#           Only apply twists when multiplying, and the undo them?
#       Improve documentation
#           Add docstrings for all classes and methods
#       Rename classes
#           "GLInvariantBundle" instead of "VectorBundle"?
#           There must be a better name then "DirectSumOfVectorBundles"?
#       Should we use tuples or lists for the partitions?
#       Fix ImodIsquared function to created the untwisted version of the ImI2 bundle, also, rename the function 
#


#   Example of how to use this code:
#      VB = VectorBundle((2,2),(1,1,1))
#      L = liftingBundle(VB)
#      L.cohomology(2)
#
#   This will return the 2nd cohomology of VB tensor VB^V tensor I/I^2
#   L.cohomology(n) will return the nth such cohomology
#



class VectorBundle:
    """

    """
    Qpartition = ()
    SVpartition = ()
    
    _cohomologyRank = -1
    _cohomology = ()
    

    def __init__(self, Qpartition: tuple, SVpartition: tuple):
        if not(all(Qpartition[i] >= Qpartition[i+1] for i in range(len(Qpartition)-1))):
            raise Exception("Qpartition must be weakly decreasing")
        if not(all(SVpartition[i] >= SVpartition[i+1] for i in range(len(SVpartition)-1))):
            raise Exception("SVpartition must be weakly decreasing")
        m = min(Qpartition[-1], SVpartition[-1])
        self.Qpartition = tuple([x-m for x in Qpartition])
        self.SVpartition = tuple([x-m for x in SVpartition])
        self._cohomologyRank = None
        self._cohomology = ()

    def __repr__(self):
        return "𝛴^{}Q ⨂ 𝛴^{}S^V".format(self.Qpartition, self.SVpartition)

    def __eq__(self,other):
        return (self.Qpartition == other.Qpartition) and (self.SVpartition == other.SVpartition)

    def __mul__(self,other):
        Qprod = lrcalc.mult(self.Qpartition,other.Qpartition,min(len(self.Qpartition),len(other.Qpartition)))
        SVprod = lrcalc.mult(self.SVpartition,other.SVpartition,min(len(self.SVpartition),len(other.SVpartition)))

        bundles = []
        for q in list(Qprod.keys()):
            correctedq = list(q)
            if len(correctedq) < len(self.Qpartition):
                correctedq.extend([0]*(len(self.Qpartition)-len(correctedq)))
            for sv in list(SVprod.keys()):
                correctedsv = list(sv)
                if len(correctedsv) < len(self.SVpartition):
                    correctedsv.extend([0]*(len(self.SVpartition)-len(correctedsv)))
                bundles.append((VectorBundle(correctedq,correctedsv), Qprod[q]*SVprod[sv]))
        return DirectSumOfVectorBundles(bundles)


    def _computeCohomology(self):
        if self._cohomologyRank is not None:
            return
        alpha = self.Qpartition + self.SVpartition
        alphaRho = [alpha[i]+len(alpha)-i for i in range(len(alpha))]
        if len(alphaRho) != len(set(alphaRho)):
            self._cohomologyRank = -1
            self._cohomology = 0
        else:
            self._cohomologyRank = Permutation(np.argsort(list(reversed(alphaRho)))+1).length()
            alphaRho.sort(reverse=True)
            self._cohomology = tuple([alphaRho[i]-len(alphaRho)+i for i in range(len(alphaRho))])
        
    
    def Qrank(self):
        return len(self.Qpartition)

    def SVrank(self):
        return len(self.SVpartition)

    def grassmannianRank(self):
        return self.Qrank(), self.Qrank()+self.SVrank()

    def cohomologyRank(self):
        if self._cohomologyRank is None:
            self._computeCohomology()
        return self._cohomologyRank

    def cohomology(self):
        if self._cohomologyRank is None:
            self._computeCohomology()
        return self._cohomology
            
    def dual(self):
        return VectorBundle(tuple([-1*x for x in reversed(self.Qpartition)]),
                            tuple([-1*x for x in reversed(self.SVpartition)]))


def _mulhelper(x: tuple):
    return (VectorBundle(x[0][0],x[0][1]),x[1])


class DirectSumOfVectorBundles:
    """

    """
    Bundles = []

    def __init__(self, Bundles: list):
        if len(Bundles) < 1:
            raise Exception("DirectSumOfVectorBundles needs at least one bundle")
        # TO DO: Check that all entries of Bundles are (VectorBundle, integer) pairs
        self.Bundles = Bundles

    def __repr__(self):
        s = "({})".format(str(self.Bundles[0][0]))
        for _ in range(1,self.Bundles[0][1]):
            s = "{} ⨁ ({})".format(s, str(self.Bundles[0][0]))
        
        for i in range(1,len(self.Bundles)):
            for _ in range(self.Bundles[i][1]):
                s = "{} ⨁ ({})".format(s, str(self.Bundles[i][0]))
        return s

    def __mul__(self,other):
        # TO DO: introduce a helper function, there is too much indentation here!
        helperDict = {}
        for VB1 in self.Bundles:
            for VB2 in other.Bundles:
                prod = VB1[0]*VB2[0]
                for vb in prod.Bundles:
                    k = (vb[0].Qpartition, vb[0].SVpartition)
                    if helperDict.get(k) is None:
                        helperDict[k] = vb[1]
                    else:
                        helperDict[k] += vb[1]
        allBun = map(_mulhelper, helperDict.items())
        return DirectSumOfVectorBundles(list(allBun))

    def cohomology(self, n:int):
        """
        Returns the nth cohomology of this direct sum of VectorBundles
        :param n: an integer
        :return: A list of lists of integers
        """
        coho = []

        for bun in self.Bundles:
            if bun[0].cohomologyRank() == n:
                coho.append((bun[0].cohomology(),bun[1]))
        if len(coho) == 0:
            return 0
        return coho

    def prettyPrint(self):
        return str(self.Bundles)
        

def zeroFunctor(rank: int):
    """
    Returns a partition of all 0 of length equal to rank
    :param rank: an integer
    :return: A list of integers
    """
    if rank < 0:
        raise Exception("Rank must be non-negative")
    return [0]*rank


def wedgeFunctor(n: int, rank: int):
    """
    Returns a partition corresponding to the nth wedge product length equal to rank
    :param n: an integer
    :param rank: an integer
    :return: A list of integers
    """
    if (rank < 0) or (n < 1):
        raise Exception("Rank must be non-negative and n must be positive")
    if rank < n:
        return []
    p = [1]*n
    p.extend([0]*(rank-n))
    return p


def symetricFunctor(n: int, rank: int):
    """
    Returns a partition corresponding to the nth symetric functor of length equal to rank
    :param n: an integer
    :param rank: an integer
    :return: A list of integers
    """
    if (rank < 0) or (n < 1):
        raise Exception("Rank must be non-negative and n must be positive")
    if rank == 0:
        return []
    p = [0]*rank
    p[0]=n
    return p

def ImodIsquared(k: int, n: int):
    """
    Returns I/I^2 for the Grassmannian G(k,n)
    :param k: a positive integer
    :param n: a positive integer strictly greater than k
    :return: A VectorBundle
    """
    if (k < 0) or (n <= k):
        raise Exception("Preconditions not met, must have 0 < k < n")

    if (k != 2) or (n < 4):
        raise NotImplementedError("Currently only the case where k = 2 and n >=4 is supported")

    svpartition = [1]*(n-2)
    svpartition[0] = 2
    svpartition[1] = 2
    return VectorBundle([0,0],svpartition)


def liftingBundle(bundle: VectorBundle):
    k,n = bundle.grassmannianRank()
    ImI2 = ImodIsquared(k,n)

    return (bundle*(bundle.dual()))*DirectSumOfVectorBundles([(ImI2,1)])

def liftable(bundle: VectorBundle):
    return liftingBundle(bundle).cohomology(2) == 0