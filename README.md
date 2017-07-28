# redsvd

A fork of the original repo from google code (find it [here](https://code.google.com/archive/p/redsvd/)).

## Building

*Note*: If you want to build the tests and samples, just uncomment out the commented lines in `wscript` and remove the very last line.

## Changes to the original

When using sparse matrices and the SVD option, the program spits out some additional files:

- the `US` matrix
- the `VS` matrix
- the all vs. all cosine dissimilarity matrix of either `US` or `VS`, whichever has fewer rows
- the `US` vs `VS` cosine dissimilarity matrix if `US` has fewer rows than `VS`, or `VS` vs `US` otherwise.

Also, it builds against the Eigen library that is included in this repository to make it an all in one package.
