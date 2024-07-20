import numpy as np
import sage.libs.lrcalc.lrcalc as lrcalc


class SchurFunctor:
    """
    The class of SchurFunctors
    """
    partition = []
    rank = 0

    def __init__(self, partition: list):
        #if not(all(type(x) is Integer for x in partition)):
        #    raise Exception("SchurFunctor only accepts partitions of integers")
        if not(all(partition[i] >= partition[i+1] for i in range(len(partition)-1))):
            raise Exception("SchurFunctor only accepts partitions, that is list must be weakly decreasing")
        if len(partition)>0 and partition[-1] < 0:
            raise Exception("SchurFunctor only accepts partitions of non-negative integers")
        
        self.partition = partition
        self.rank = len(partition)

    def __repr__(self):
        return "𝛴^{}".format(self.partition)

def zeroFunctor(rank: int):
    """
    Returns the SchurFunctor corresponding to the 0 partition of the desired rank
    :param rank: an integer
    :return: A SchurFunctor
    """
    if rank < 0:
        raise Exception("Rank must be non-negative")
    return SchurFunctor([0]*rank)


def wedgeFunctor(n: int, rank: int):
    """
    Returns the SchurFunctor corresponding to the nth wedge product of a vector bundle of the desired rank
    :param n: an integer
    :param rank: an integer
    :return: A SchurFunctor
    """
    if (rank < 0) or (n < 1):
        raise Exception("Rank must be non-negative and n must be positive")
    if rank < n:
        return SchurFunctor([])
    p = [1]*n
    p.extend([0]*(rank-n))
    return SchurFunctor(p)


def symetricFunctor(n: int, rank: int):
    """
    Returns the SchurFunctor corresponding to the nth symetric functor of a vector bundle of the desired rank
    :param n: an integer
    :param rank: an integer
    :return: A SchurFunctor
    """
    if (rank < 0) or (n < 1):
        raise Exception("Rank must be non-negative and n must be positive")
    if rank == 0:
        return SchurFunctor([])
    p = [0]*rank
    p[0]=n
    return SchurFunctor(p)
    

def productOfSchurFunctors(SF1: SchurFunctor, SF2: SchurFunctor):
    """
    Compute the product of two SchurFunctors. Note the resulting list order is not canoncial
    :param SF2: a SchurFunctor
    :param SF2: a SchurFunctor
    :return: A list of SchurFunctors whose direct sum is the product of SF1 and SF2
    """
    # Question: should the rank actually be SF1.rank?
    r = min(SF1.rank,SF2.rank)

    # Leverages the lrcalc package to do the computations for us
    sfdict = lrcalc.mult(SF1.partition, SF2.partition, r)
    functors = list(sfdict.keys())
    prod = []
    for fun in functors:
        # We want to preserve the notion of rank in our SchurFunctors
        # Here we exetend the partitions by appending 0 until its length is equal to the desired rank
        f = list(fun)
        if len(f) < r:
            f.extend([0]*(r-len(f)))
        # There is a notion of multiplicity in the product of Schur Functors, so we append a given functor
        # multiple times to match the multiplicity 
        for _ in range(sfdict[fun]):
            prod.append(SchurFunctor(f))
    return prod


def normalizeListPair(lst1: list, lst2: list):
    """
    Normalize a pair of lists integers, so that they are both non-negative but the difference between elements is the same
    Precondition: Both lists are weakly decreasing
    :param lst1: a list of integers
    :param lst2: a list of integers
    :return: a pair of lists (m1, m2) that are the normalized forms of the original lists
    """
    m = min(lst1[-1], lst2[-1])
    return ([x-m for x in lst1], [x-m for x in lst2])


class VectorBundle:
    """
    The class of vector bundles of the form S^α Q ⊗ S^β S^v, where α and β are partitions
    """
    def __init__(self, qFunctor: SchurFunctor, sVFunctor: SchurFunctor):
        self.QFunctor = qFunctor
        self.SVFunctor = sVFunctor
        self._normalize()

    def _normalize(self):
        qpart, svpart = normalizeListPair(self.QFunctor.partition, self.SVFunctor.partition)
        self.QFunctor = SchurFunctor(qpart)
        self.SVFunctor = SchurFunctor(svpart)

    def __repr__(self):
        return "{}Q ⊗ {}S^V".format(str(self.QFunctor), str(self.SVFunctor))

    def grassmannianRank(self):
        return self.QFunctor.rank, self.QFunctor.rank + self.SVFunctor.rank

    def dual(self):
        Qdual = [-1*x for x in reversed(self.QFunctor.partition)]
        SVdual = [-1*x for x in reversed(self.SVFunctor.partition)]
        Qdual, SVdual = normalizeListPair(Qdual, SVdual)
        return VectorBundle(SchurFunctor(Qdual), SchurFunctor(SVdual))


def productOfVectorBundles(VB1: VectorBundle, VB2: VectorBundle):
    """
    Compute the product of two VectorBundles. Note the resulting list order is not canoncial
    :param VB1: a VectorBundle
    :param VB2: a VectorBundle
    :return: A list of VectorBundles whose direct sum is the product of VB1 and VB2
    """
    Qproduct = productOfSchurFunctors(VB1.QFunctor, VB2.QFunctor)
    SVproduct = productOfSchurFunctors(VB1.SVFunctor, VB2.SVFunctor)
    
    prodlist = []
    for q in Qproduct:
        for sv in SVproduct:
            prodlist.append(VectorBundle(q,sv))
    return prodlist


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
    return VectorBundle(SchurFunctor([0,0]),SchurFunctor(svpartition))


def transpositions_needed(arr: list):
    """
    :param arr: a list of integers
    :return: The number of transpositions need to put arr in order
    """
    return Permutation(np.argsort(arr)+1).length()



def needs_two_transpositions(arr: list):
    """
    Checks if a list needs exactly two transpositions to be sorted in decreasing order
    :param arr: a list of integers
    :return: True iff the list needs exactly two transpositions False otherwise
    """
    return transpositions_needed(list(reversed(arr))) == 2


def h2(bundle: VectorBundle):
    """
    Computes the second cohomology of the given Vector Bundle on Gr(k,n) using Borel-Bott-Weil
    where k = bundle.QFunctor.rank, and n = k + bundle.SVFunctor.rank
    :param bundle: A vector bundle
    :return: 0 if there is no second cohomology, otherwise return a weakly decreasing list of integers, λ
             such that H^2(G(k,n), bundle) = 𝛴^λ V^*
    """
    alpha = []
    alpha.extend(bundle.QFunctor.partition)
    alpha.extend(bundle.SVFunctor.partition)

    alphaRho = [alpha[i]+len(alpha)-i for i in range(len(alpha))]
    if len(alphaRho) != len(set(alphaRho)):
        # This checks if the list contains any repetitions
        # By BBW all cohomogies are 0 if there are repetitions in alphaRho
        return 0
    elif not needs_two_transpositions(alphaRho):
        # By BBW H2 is non-zero if alphaRho needs exactly two transpositions to be sorted
        return 0

    alphaRho.sort(reverse=True)
    return [alphaRho[i]-len(alphaRho)+i for i in range(len(alphaRho))]


def _checkForLiftingHelper(bundle: VectorBundle):
    k,n = bundle.grassmannianRank()
    ImI2 = ImodIsquared(k,n)
    
    bundleTensorBundleDual = productOfVectorBundles(bundle, bundle.dual())
    
    completeTensorProduct = []
    for vb in bundleTensorBundleDual:
        completeTensorProduct.extend(productOfVectorBundles(vb, ImI2))
    
    return completeTensorProduct



def checkForLifting(bundle: VectorBundle):
    """
    Computes the second cohomology of bundle ⊗ bundle^V ⊗ I/I^2, where I/I^2 is from the underlying Grassmannian
    :param bundle: A VectorBundle
    :return: Return list consisting of 0s and of a weakly decreasing list of integers, λ_i
             such that H^2(G(k,n), bundle ⊗ bundle^V ⊗ I/I^2) = ⨁ 𝛴^λ_i V^*
    """
    completeTensorProduct = _checkForLiftingHelper(bundle)
    secondCohomology = []
    for vb in completeTensorProduct:
        secondCohomology.append(h2(vb))

    return secondCohomology


def lifts(bundle: VectorBundle):
    """
    Returns True iff the second cohomology is 0, returns False otherwise
    :param bundle: A VectorBundle
    :return: Boolean
    """
    completeTensorProduct = _checkForLiftingHelper(bundle)
    for vb in completeTensorProduct:
        if h2(vb) != 0:
            return False
    return True




