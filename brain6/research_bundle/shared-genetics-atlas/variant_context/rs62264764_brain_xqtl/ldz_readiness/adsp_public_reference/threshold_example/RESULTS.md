# Why zero-filling thresholded LD is not a valid reconstruction

This is a deliberately constructed three-variant example, **not an ADSP observation**. Begin with an actual positive-definite correlation matrix whose off-diagonal values are 0.81, 0.60 and 0.19. Its smallest eigenvalue is 0.0769. If a pairwise file suppresses the 0.19 correlation under a |R| ≥ 0.2 rule and the missing value is filled with zero, the resulting matrix has smallest eigenvalue **−0.00802** and is not positive semidefinite. It cannot be the exact correlation matrix of real random variables.

This demonstrates a possible failure, not a claim that every thresholded matrix is indefinite or that all absent records are known to be below 0.2. Missing variants, low counts or other processing can also cause absence. Projecting a filled matrix to positive semidefinite form would alter values and would not reconstruct the missing correlations without additional information. The full actual matrix and study matching remain needed for valid multi-signal modeling.

[Code](verify.R) and [results](results.tsv) reproduce the example. See the [reference assessment](../REFERENCE_ASSESSMENT.md) for the publisher's threshold and access distinctions.
