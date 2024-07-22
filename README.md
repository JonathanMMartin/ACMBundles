# ACMBundles
This project is to aid in the computation of cohomologies of vector bundles on Grassmannians.

# Example 
VB = IrreducibleGLInvariantBundle((2,0),(1,1,1))

L = liftingBundle(VB)

L.cohomology(2)

This will output the second cohomology of End(𝛴^(2,0)Q ⨂ 𝛴^(1,1,1)S^V) ⨂ I/I^2, for G(2,5)
