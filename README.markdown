# CLODE - Interface to the GSL ODE solver.

## Usage

`git clone` the repository.

Ensure that quicklisp can find `CLODE`. (e.g. `git clone` into your `quicklisp/local-projects`)

Then in an instance of Common Lisp: 

```CL-USER> (ql:quickload :clode)
To load "clode":
  Load 1 ASDF system:
    clode
; Loading "clode"

(:CLODE)
CL-USER> (clode:with-ode-system (vanderpol u v v (- (* μ v (- 1 (* u u))) u) μ 10.0)
           (clode:integrate vanderpol '(1.0 0.0) :from-to '(0 . 100) :steps 100))
1.000000 -1.456875 -11.547242
2.000000 -1.956084 0.069065
3.000000 -1.884810 0.073642
.
.
.
98.000000 -1.913287 0.071734
99.000000 -1.839039 0.076957
100.000000 -1.758888 0.083643
0
CL-USER> 
```

The second command reproduces the example of the 'Van Der Pol' oscillator from the GSL documentation.

Additional arguments can be given to `INTEGRATE` in order to make use of a subset of GSL stepping functions (those which do not require jacobians), and to dump the results to a stream (e.g. a file). Please take a look at the keyword arguments for `INTEGRATE`.

## Installation

CLODE makes use of GSL library. In particular, `libgslcblas` is required. The source code can be obtained [[here][https://www.gnu.org/software/gsl/doc/html/intro.html#obtaining-gsl]. Alternatively, install `libgsl` with your package manager.

## Author

* Lewis Grozinger (l.grozinger2@ncl.ac.uk)

## Copyright

Copyright (c) 2019 Lewis Grozinger (l.grozinger2@ncl.ac.uk)
