#   TO DO:
#       Improve documentation
#           Add docstrings for all classes and methods
#       Keep track of twists using seperate twist parameter
#       Spead this across multiply files, that is put non-class function into a seperate file and then import this one into it.

import numpy as np
import sage.libs.lrcalc.lrcalc as lrcalc

class IrreducibleGLInvariantBundle:
    """

    """
    # Note we are using tuples here to help with multiplication of GLInvariantBundles later, don't swtch them to lists!!
    _Qpartition = ()
    _SVpartition = ()
    
    _cohomologyRank = -1
    _cohomology = ()
    

    def __init__(self, Qpartition: tuple, SVpartition: tuple, twist:int = 0):
        if not(all(Qpartition[i] >= Qpartition[i+1] for i in range(len(Qpartition)-1))):
            raise Exception("Qpartition must be weakly decreasing")
        if not(all(SVpartition[i] >= SVpartition[i+1] for i in range(len(SVpartition)-1))):
            raise Exception("SVpartition must be weakly decreasing")
        self._Qpartition = Qpartition
        self._SVpartition = SVpartition
        self._cohomologyRank = -1
        self._cohomology = ()

    def __repr__(self):
        return "𝛴^{}Q ⨂ 𝛴^{}S^V".format(self.qPart(), self.svPart())

    def __eq__(self,other):
        return (self.qPart() == other.qPart()) and (self.svPart() == other.svPart())

    def __mul__(self,other):
        # TO DO: Double check that these rank calculations make sence
        minQRank = min(len(self.qPart()),len(other.qPart()))
        minSVRank = min(len(self.svPart()),len(other.svPart()))

        selfShift = min(self.qPart()[-1],self.svPart()[-1])
        otherShift = min(other.qPart()[-1],other.svPart()[-1])

        Qprod = lrcalc.mult([x-selfShift for x in self.qPart()],[x-otherShift for x in other.qPart()],minQRank)
        SVprod = lrcalc.mult([x-selfShift for x in self.svPart()],[x-otherShift for x in other.svPart()],minSVRank)

        bundles = []
        for q in list(Qprod.keys()):
            correctedq = list(q)
            if len(correctedq) < minQRank:
                correctedq.extend([0]*(minQRank-len(correctedq)))
            correctedq = tuple([x + selfShift + otherShift for x in correctedq])
            for sv in list(SVprod.keys()):
                correctedsv = list(sv)
                if len(correctedsv) < minSVRank:
                    correctedsv.extend([0]*(minSVRank-len(correctedsv)))
                correctedsv = tuple([x + selfShift + otherShift for x in correctedsv])
                bundles.append((IrreducibleGLInvariantBundle(correctedq,correctedsv), Qprod[q]*SVprod[sv]))
        return GLInvariantBundle(bundles)

    def qPart(self):
        return self._Qpartition

    def svPart(self):
        return self._SVpartition

    def cohomologyRank(self):
        if self._cohomologyRank == -1:
            self._computeCohomology()
        return self._cohomologyRank

    def cohomology(self):
        if self._cohomologyRank == -1:
            self._computeCohomology()
        return self._cohomology

    def qRank(self):
        return len(self._Qpartition)

    def svRank(self):
        return len(self._SVpartition)

    def grassRank(self):
        return self.qRank(), self.qRank() + self.svRank()

    def dual(self):
        return IrreducibleGLInvariantBundle(
                    [-1*x for x in reversed(self.qPart())],
                    [-1*x for x in reversed(self.svPart())])

    def twist(self,n:int = 1):
        return IrreducibleGLInvariantBundle(
                    [x+n for x in self.qPart()],
                    [x+n for x in self.svPart()])

    def _computeCohomology(self):
        if self._cohomologyRank != -1:
            return
        alpha = self.qPart() + self.svPart()
        alphaRho = [alpha[i]+len(alpha)-i for i in range(len(alpha))]
        if len(alphaRho) != len(set(alphaRho)):
            self._cohomologyRank = None
            self._cohomology = 0
        else:
            self._cohomologyRank = Permutation(np.argsort(list(reversed(alphaRho)))+1).length()
            alphaRho.sort(reverse=True)
            self._cohomology = tuple([alphaRho[i]-len(alphaRho)+i for i in range(len(alphaRho))])



class GLInvariantBundle:
    """

    """
    _bundles = []

    def __init__(self, Bundles: list):
        if len(Bundles) < 1:
            raise Exception("GLInvariantBundle needs at least one bundle")
        # TO DO: Check that all entries of Bundles are (IrreducibleGLInvariantBundle, integer) pairs
        self._bundles = Bundles

    def __repr__(self):
        s = "({})".format(str(self.bundles()[0][0]))
        for _ in range(1,self.bundles()[0][1]):
            s = "{} ⨁ ({})".format(s, str(self.bundles()[0][0]))
        
        for i in range(1,len(self.bundles())):
            for _ in range(self.bundles()[i][1]):
                s = "{} ⨁ ({})".format(s, str(self.bundles()[i][0]))
        return s

    def __mul__(self,other):
        # TO DO: introduce a helper function, there is too much indentation here!
        helperDict = {}
        for (VB1, mult1) in self.bundles():
            for (VB2, mult2) in other.bundles():
                prod = VB1*VB2
                for (vb, mult3) in prod.bundles():
                    k = (vb.qPart(), vb.svPart())
                    if helperDict.get(k) is None:
                        helperDict[k] = mult1*mult2*mult3
                    else:
                        helperDict[k] += mult1*mult2*mult3
        allBun = map(lambda x: (IrreducibleGLInvariantBundle(x[0][0],x[0][1]),x[1]), helperDict.items())
        return GLInvariantBundle(list(allBun))

    def bundles(self):
        return self._bundles

    def twist(self, n: int = 1):
        return list(map(lambda x: (x[0].twist(),x[1]), self.bundles()))

    def cohomology(self, n:int):
        """
        Returns the nth cohomology of this direct sum of VectorBundles
        :param n: an integer
        :return: A list of lists of integers
        """
        coho = []

        for bun in self.bundles():
            if bun[0].cohomologyRank() == n:
                coho.append((bun[0].cohomology(),bun[1]))
        if len(coho) == 0:
            return 0
        return coho


def zeroFunctor(rank: int):
    """
    Returns a partition of all 0 of length equal to rank
    :param rank: an integer
    :return: A tuple of integers
    """
    if rank < 0:
        raise Exception("Rank must be non-negative")
    return tuple([0]*rank)


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
    return tuple(p)


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
    return tuple(p)

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

    svpartition = [0]*(n-2)
    svpartition[0] = 1
    svpartition[1] = 1
    return IrreducibleGLInvariantBundle((-1,-1),tuple(svpartition))


def liftingBundle(bundle: IrreducibleGLInvariantBundle):
    k,n = bundle.grassRank()
    ImI2 = ImodIsquared(k,n)

    return (bundle*(bundle.dual()))*GLInvariantBundle([(ImI2,1)])

def liftable(bundle: IrreducibleGLInvariantBundle):
    return liftingBundle(bundle).cohomology(2) == 0


