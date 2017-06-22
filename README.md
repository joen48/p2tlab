# p2tlab

Solve sets of polynomial equations using higher-order tensor decompositions. 
A Matlab implementation of the algorithms presented in Sets of polynomial equations, decompositions of higher-order tensors and multidimensional harmonic retrieval: connections and algorithms, J. Vanderstukken, promoter: L. De Lathauwer. Master's thesis, KULeuven, 2017.

The algorithms are mostly wrappers in Matlab 9.1 (R2016b) that construct a third-order tensor containing the null space of the so-called Macaulay matrix and forward the computation of its decomposition to [Tensorlab, version 3.0][1]. [PNLA][2] is used to construct the Macaulay matrix.

## Getting started

Once you have fetched P2Tlab containing the programs, Tensorlab and PNLA to the current directory, add them to the Matlab search path.
```matlab
addpath('p2tlab', 'tensorlab', 'PNLA_MATLAB_OCTAVE')
```

Appendix H of the text highlights the use of the algorithms.

[1]: http://www.tensorlab.net/
[2]: https://github.com/kbatseli/PNLA_MATLAB_OCTAVE/
