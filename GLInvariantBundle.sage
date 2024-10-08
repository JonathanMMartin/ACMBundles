# Python Imports
from numpy import argsort
from itertools import combinations_with_replacement
import pickle

# Sage Imports
import sage.libs.lrcalc.lrcalc as lrcalc


# TODO: Split this across multiple files
# TODO: Improve error handling (try and except) for larger methods like multiplication
# TODO: When multiplication of IrreducibleGLInvariantBundles results in a list of exactly one bundle with multiplicity 1, should
#       it return an IrreducibleGLInvariantBundle? Or would that cause some annoying problems?

class IrreducibleGLInvariantBundle:
    """

    """
    _Qpartition = ()
    _SVpartition = ()
    
    _cohomologyRank = -1
    _cohomology = ()
    

    def __init__(self, Qpartition: tuple, SVpartition: tuple):
        try:
            assert all(Qpartition[i] >= Qpartition[i+1] for i in range(len(Qpartition)-1)), "Qpartition must be weakly decreasing"
            assert all(SVpartition[i] >= SVpartition[i+1] for i in range(len(SVpartition)-1)), "SVpartition must be weakly decreasing"
            self._Qpartition = Qpartition
            self._SVpartition = SVpartition    
            self._cohomologyRank = -1
            self._cohomology = ()
        except AssertionError as err:
            print("Error:", err)

    def __repr__(self):
        if self.noQ() and self.noSV():
            return "𝒪_G({},{})".format(self.qRank(), self.qRank()+self.svRank())
        elif self.noQ():
            return "𝛴^{}S^V".format(self.svPart())
        elif self.noSV():
            return "𝛴^{}Q".format(self.qPart())
        return "𝛴^{}Q ⨂ 𝛴^{}S^V".format(self.qPart(), self.svPart())

    def __eq__(self,other):
        return (self.qPart() == other.qPart()) and (self.svPart() == other.svPart())

    def __mul__(self,other):
        # TODO: Add documentation to this method
        if type(other) not in [IrreducibleGLInvariantBundle, GLInvariantBundle]:
            raise Exception("Can only multiply a GLInvariantBundle by another GLInvariantBundle")
        if type(other) == GLInvariantBundle:
            return other*(GLInvariantBundle([(self, 1)]))

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

    def qRank(self):
        return len(self._Qpartition)

    def svRank(self):
        return len(self._SVpartition)

    def grassRank(self):
        return self.qRank(), self.qRank() + self.svRank()

    def noQ(self):
        return self.qPart()[0] == 0 and self.qPart()[-1] == 0

    def noSV(self):
        return self.svPart()[0] == 0 and self.svPart()[-1] == 0

    def cohomologyRank(self):
        if self._cohomologyRank == -1:
            self._computeCohomology()
        return self._cohomologyRank

    def cohomology(self, n:int = -1):
        if self._cohomologyRank == -1:
            self._computeCohomology()
        if n == -1:
            return self.cohomology(self.cohomologyRank())
        elif n != self.cohomologyRank():
            return 0
        return self._cohomology

    def dual(self):
        return IrreducibleGLInvariantBundle(
                    [-1*x for x in reversed(self.qPart())],
                    [-1*x for x in reversed(self.svPart())])

    def shift(self,n:int = 1):
        return IrreducibleGLInvariantBundle(
                    [x+n for x in self.qPart()],
                    [x+n for x in self.svPart()])

    def _computeCohomology(self):
        """
        Computes the non-zero cohomology using BBW based on the conventions of Costa et al.
        :return: Either 0 or a tuple of length n where G(k,n) is the underlying Grassmannian
        """
        if self._cohomologyRank != -1:
            return
        alpha = self.qPart() + self.svPart()
        alphaRho = [alpha[i]+len(alpha)-i for i in range(len(alpha))]
        if len(alphaRho) != len(set(alphaRho)):
            self._cohomologyRank = None
            self._cohomology = 0
        else:
            self._cohomologyRank = Permutation(argsort(list(reversed(alphaRho)))+1).length()
            alphaRho.sort(reverse=True)
            self._cohomology = tuple([alphaRho[i]-len(alphaRho)+i for i in range(len(alphaRho))])



class GLInvariantBundle:
    """

    """
    _bundles = []

    def __init__(self, Bundles: list):
        try:
            assert len(Bundles) > 0, "GLInvariantBundle needs at least one bundle"
            # TODO: Check that all entries of Bundles are (IrreducibleGLInvariantBundle, integer) pairs
            self._bundles = Bundles
        except AssertionError as err:
            print("Error: ", err)
        
    def __repr__(self):
        # TODO: Find a better way to print this, is there an easy way to use superscripts for multiplicities?
        s = "({})".format(str(self.bundles()[0][0]))
        for _ in range(1,self.bundles()[0][1]):
            s = "{} ⨁ ({})".format(s, str(self.bundles()[0][0]))
        
        for i in range(1,len(self.bundles())):
            for _ in range(self.bundles()[i][1]):
                s = "{} ⨁ ({})".format(s, str(self.bundles()[i][0]))
        return s

    def __mul__(self,other):
        if type(other) not in [IrreducibleGLInvariantBundle, GLInvariantBundle]:
            raise Exception("Can only multiply a GLInvariantBundle by another GLInvariantBundle")
        if type(other) == IrreducibleGLInvariantBundle:
            return self*(GLInvariantBundle([(other, 1)]))

        # TODO: introduce a helper function/method, there is too much indentation here!
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

    def shift(self, n: int = 1):
        return GLInvariantBundle(list(map(lambda x: (x[0].shift(n),x[1]), self.bundles())))

    def cohomology(self, n:int):
        """
        Returns the nth cohomology of this GLInvariantBundle
        :param n: an integer
        :return: Either 0 or a list of (list of integers, integer) pairs
        """
        coho = []

        for bun in self.bundles():
            if bun[0].cohomologyRank() == n:
                coho.append((bun[0].cohomology(),bun[1]))
        if len(coho) == 0:
            return 0
        return coho


# Here we define helper functions to give us special bundles
def trivialBundle(k:int, n: int):
    """
    Returns the trivial bundle on G(k,n)
    :param k: an integer
    :param n: an integer
    :return: An IrreducibleGLInvariantBundle
    """
    try:
        assert (0 < k) and (k < n), "Preconditions not met, must have 0 < k < n"
        return IrreducibleGLInvariantBundle(tuple([0]*k), tuple([0]*(n-k)))
    except AssertionError as err:
        print("Error: ", err)

def wedgeBundle(k: int, n: int, m: int, Q: bool = True):
    """
    Returns the mth wedge product of Q (or S^V) on G(k,n)
    :param k: an integer
    :param n: an integer
    :param m: an integer
    :param Q: a boolean
    :return: An IrreducibleGLInvariantBundle
    """
    try:
        assert (0 < k) and (k < n), "Preconditions not met, must have 0 < k < n"
        if Q:
            assert m <=k, "For the wedge product of Q, m must be at most the dimension of Q (which is k)"
            p = [1]*m
            p.extend([0]*(k-m))
            return IrreducibleGLInvariantBundle(tuple(p), tuple([0]*(n-k)))
        assert m <= n-k, "For the wedge product of S^V, m must be at most the dimenstion of S^V (which is n-k)"
        p = [1]*m
        p.extend([0]*(n-k-m))
        return IrreducibleGLInvariantBundle(tuple([0]*k), tuple(p))
    except AssertionError as err:
        print("Error: ", err)

def symetricBundle(k: int, n: int, m: int, Q: bool = True):
    """
    Returns the mth symetric of Q (or S^V) on G(k,n)
    :param k: an integer
    :param n: an integer
    :param m: an integer
    :param Q: a boolean
    :return: An IrreducibleGLInvariantBundle
    """
    try:
        assert (0 < k) and (k < n), "Preconditions not met, must have 0 < k < n"
        if Q:
            p = [0]*k
            p[0] = m
            return IrreducibleGLInvariantBundle(tuple(p), tuple([0]*(n-k)))
        p = [0]*(n-k)
        p[0] = m
        return IrreducibleGLInvariantBundle(tuple([0]*k), tuple(p))
    except AssertionError as err:
        print("Error: ", err)

def tangentBundle(k: int, n: int):
    """
    Returns the cotangent bundle of G(n,k), which is Q ⨂ S = 𝛴^(1,0) ⨂ 𝛴^(0,-1)S^V 
    :param k: an integer
    :param n: an integer
    :return: An IrreducibleGLInvariantBundle
    """
    qPart = [0]*k
    qPart[0]=1
    
    svPart = [0]*(n-k)
    svPart[-1] = -1
    return IrreducibleGLInvariantBundle(tuple(qPart),tuple(svPart))

def cotangentBundle(k: int, n: int):
    """
    Returns the cotangent bundle of G(n,k), which is Q^V ⨂ S^V = 𝛴^(0,-1) ⨂ 𝛴^(1,0)S^V 
    :param k: an integer
    :param n: an integer
    :return: An IrreducibleGLInvariantBundle
    """
    return tangentBundle(k,n).dual()

def endomorphism(bundle):
    """
    Returns End(bundle) = bundle ⨂ bundle^V
    :param bundle: Either a IrreducibleGLInvariantBundle or a GLInvariantBundle
    :return: A GLInvariantBundle
    """
    return bundle*(bundle.dual())

def ImodIsquared(k: int, n: int):
    """
    Returns I/I^2 for the Grassmannian G(k,n)
    :param k: a positive integer
    :param n: a positive integer strictly greater than k
    :return: A VectorBundle
    """
    try:
        assert (0 < k) and (k < n), "Preconditions not met, must have 0 < k < n"
        assert (k == 2) and (4 <= n), "Currently only the case where k = 2 and n >=4 is supported"
        svpartition = [0]*(n-2)
        svpartition[0] = 1
        svpartition[1] = 1
        return IrreducibleGLInvariantBundle((-1,-1),tuple(svpartition))
    except AssertionError as err:
        print("Error: ", err)

# Here we define functions that help us compute the obstruction space of a vector bundle
def liftingBundle(bundle: IrreducibleGLInvariantBundle):
    k,n = bundle.grassRank()
    ImI2 = ImodIsquared(k,n)
    return endomorphism(bundle)*ImI2

def liftable(bundle: IrreducibleGLInvariantBundle):
    l = liftingBundle(bundle)
    for (bun, mult) in l.bundles():
        if bun.cohomologyRank() == 2:
            return False
    return True

# Here we define the function that will do the main search work for us
# In order to use this code the minimum requirement is to understand the parameters of the following function
def exhuastiveObstructionSearch(end: int, k: int = 2, n: int = 5, skip: bool = False, maxNew: int = -1, checkpoint: int = 100, logging:bool = False):
    """
    Checks all IrreducibleGLInvariantBundle s of the form 𝛴^α Q ⨂ 𝛴^β S^V on G(k,n) where
    α, β are partitions where all numbers are at most end to see if any of them lift
    :param end: an integer, determines how high we are willing to check
    :param k: an integer, defines the Grassmannian we are working over
    :param n: an integer, defines the Grassmannian we are working over
    :param skip: a boolean, determines if we should skip bundles of the form 𝛴^(n,n) Q and 𝛴^(n,n,n) S^V
    :param checkpoint: an integer, how often to report on progress and save the results to a .pickle file
    :param logging: a boolean, whether or not to send output to a .txt file to record the results
    :return: a list of non-linebundle IrreducibleGLInvariantBundles that have 0 obstruction space
    """
    if end < 0:
        raise Exception("Currently we must have end >= 0")
    
    liftResults = _getObstructionResults(k,n)

    # TODO: Add "infinite" search mode
    
    possiblitycounter = 0
    checkcounter = 0
    
    VB = trivialBundle(k,n)
    # TODO: clean this up, it is confusing the way that we are manipulating these tuples, it is correct, but there must be a cleaner way
    for qPart in combinations_with_replacement(range(end+1),k):
        if skip and qPart[0] == qPart[-1]:
            continue
        for svPart in combinations_with_replacement(range(end+1),n-k):            
            if skip and svPart[0] == svPart[-1]:
                continue
            
            # Note that a bundle has 0 cohomology iff any shifted version of it has 0 second cohomology
            # Therefore it suffices to "normalize" the bundle and check if that version has 0 obstruction space
            m = min(qPart[0],svPart[0])
            correctedq = tuple([x-m for x in qPart[::-1]])
            correctedsv = tuple([x-m for x in svPart[::-1]])
            key = (correctedq, correctedsv)

            # If liftResults already has a value for this, then we have already check if an equivalent
            # bundles has zero obstruction space, so it is not needed to check this one.
            possiblitycounter += 1
            if liftResults.get(key) is not None:
                continue
            
            VB = IrreducibleGLInvariantBundle(correctedq,correctedsv)
            liftResults[key] = liftable(VB)
            if liftResults[key]:
                print("{} has zero obstruction space!!!".format(VB))
            
            checkcounter+=1
            if checkcounter % checkpoint == 0:
                print("We have checked {} new bundles so far, the last one checked was: {}".format(checkcounter,VB))
                _saveObstructionResults(liftResults, k, n)
            if checkcounter == maxNew:
                break
        if checkcounter == maxNew:
            break
    if checkcounter == 0:
        print("There were {} bundles in the given search space. Of them no new results were found".format(possiblitycounter))
    elif checkcounter == maxNew:
        print("{} new bundles were checked. The last one was {}".format(checkcounter, VB))
    else:
        print("There were {} bundles in the given search space. Of them {} were new bundles, the last one was: {}".format(possiblitycounter, checkcounter, VB))
    _saveObstructionResults(liftResults, k, n)

    if logging:
        _outputSearchResults(end, k, n, skip, maxNew)

# Here are some helper functions to aid with logging the results.
def bundlesWith0ObstructionSpace(k: int, n: int):
    results = _getObstructionResults(k,n)
    bundles = []
    for (bun, value) in results.items():
        if value:
            bundles.append(bun)
    return bundles

def _outputSearchResults(end: int, k: int, n: int, skip: bool, maxNew: int):
    bundles = bundlesWith0ObstructionSpace(k,n)
    file = open("g{}{}obstructionSearch.txt".format(k,n),'a')
    file.write("The results from using the parameters: end={}, k={}, n={}, skip={}, maxNew={} are: \n".format(end, k, n, skip, maxNew))
    file.write("Results: {}\n\n".format(bundles))
    file.close()

def _getObstructionResults(k: int, n: int):
    try:
        with open("g{}{}ObstructionResults.pickle".format(k,n), "rb") as f:
            liftResults = pickle.load(f)
    except:
        print("No results found for G({},{})".format(k,n))
        liftResults = {}
    return liftResults

def _saveObstructionResults(results: dict, k: int, n: int):
    try:
        with open("g{}{}ObstructionResults.pickle".format(k,n), "wb") as f:
            pickle.dump(results,f)
    except:
        print("There was an error while saving the results. Please verify the integrety of the results file")
