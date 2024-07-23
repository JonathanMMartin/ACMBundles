from numpy import argsort
import sage.libs.lrcalc.lrcalc as lrcalc

class IrreducibleGLInvariantBundle:
    """

    """
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
        if all(x == 0 for x in [self.qPart()[0], self.qPart()[-1], self.svPart()[0], self.svPart()[-1]]):
            return "𝒪_G({},{})".format(self.qRank(), self.qRank()+self.svRank())
        elif self.qPart()[0] == 0 and self.qPart()[-1]==0:
            return "𝛴^{}S^V".format(self.svPart())
        elif self.svPart()[0] == 0 and self.svPart()[-1]==0:
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
        if len(Bundles) < 1:
            raise Exception("GLInvariantBundle needs at least one bundle")
        # TODO: Check that all entries of Bundles are (IrreducibleGLInvariantBundle, integer) pairs
        self._bundles = Bundles

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
    if (k < 0) or (n <= k):
        raise Exception("Preconditions not met, must have 0 < k < n")
    return IrreducibleGLInvariantBundle(tuple([0]*k), tuple([0]*(n-k)))


def wedgeBundle(k: int, n: int, m: int, Q: bool = True):
    """
    Returns the mth wedge product of Q (or S^V) on G(k,n)
    :param k: an integer
    :param n: an integer
    :param m: an integer
    :param Q: a boolean
    :return: An IrreducibleGLInvariantBundle
    """
    if (k < 0) or (n <= k):
        raise Exception("Preconditions not met, must have 0 < k < n")
    if Q:
        if k < m:
            raise Exception("m must be at most the dimension of Q")
        p = [1]*m
        p.extend([0]*(k-m))
        return IrreducibleGLInvariantBundle(tuple(p), tuple([0]*(n-k)))
    
    if m > n-k:
        raise Exception("m must be at most the dimension of S^V")
    p = [1]*m
    p.extend([0]*(n-k-m))
    return IrreducibleGLInvariantBundle(tuple([0]*k), tuple(p))

def symetricBundle(k: int, n: int, m: int, Q: bool = True):
    """
    Returns the mth symetric of Q (or S^V) on G(k,n)
    :param k: an integer
    :param n: an integer
    :param m: an integer
    :param Q: a boolean
    :return: An IrreducibleGLInvariantBundle
    """
    if (k < 0) or (n <= k):
        raise Exception("Preconditions not met, must have 0 < k < n")
    if Q:
        p = [0]*k
        p[0] = m
        return IrreducibleGLInvariantBundle(tuple(p), tuple([0]*(n-k)))
    p = [0]*(n-k)
    p[0] = m
    return IrreducibleGLInvariantBundle(tuple([0]*k), tuple(p))

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
    return endomorphism(bundle)*ImI2

def liftable(bundle: IrreducibleGLInvariantBundle):
    l = liftingBundle(bundle)
    for (bun, mult) in l.bundles():
        if bun.cohomologyRank() == 2:
            return False
    return True

def exhuastiveObstructionSearch(end: int, start: int = 0, k: int = 2, n: int = 5, checkpoint: int = 10):
    """
    Checks all IrreducibleGLInvariantBundle s of the form 𝛴^α Q ⨂ 𝛴^β S^V on G(k,n) where
    α, β are partitions where all numbers are at most end to see if any of them lift
    :param k: an integer
    :param n: an integer
    :param start: an integer
    :param end: an integer
    :param checkpoint: an integer
    :return: 
    """
    if k != 2 or n != 5:
        raise Exception("Currently only the case where k = 2 and n = 5 is supported")
    if end < 0:
        raise Exception("Currently we must have end >= 0")
    if end < start:
        raise Exception("Preconditions not met must have start <= end")
    
    check = 0

    # TODO: Implement an iterator to help make this cleaner and easier to generatlize
    # TODO: Add "infinite" search mode
    for a in range(start,end+1):
        for b in range(a+1):
            for c in range(end+1):
                for d in range(c+1):
                    for e in range(d+1):
                        VB = IrreducibleGLInvariantBundle((a,b),(c,d,e))
                        if liftable(VB):
                            print(str(VB) + " has zero obstruction space!!!")
                        check+=1
                        if check % checkpoint == 0:
                            print("We have checked " + str(check) + " Bundles, the last one checked was: " + str(VB))