# CSNumerics

<b>Portable numerical algorithms in C#</b>

Copyright (c) 2012-1014 Anders Gustafsson, Cureos AB. Licensed under the General Public License version 3.<br />

## Usage

*CSNumerics* is a Portable Class Library (PCL) of various numerical algorithms written in C#. It currently consists of a C# implementation of Michael J.D. Powell's *LINCOA* 
algorithm for minimizing a nonlinear objective function subject to linear constraints.

The *CSNumerics* assembly is immediately applicable with the following targets:

* .NET Framework 4 and later
* Windows Store apps (Windows 8 and later)
* Windows Phone 8
* Silverlight 5
* Xamarin.iOS
* Xamarin.Android

### LINCOA

The C# implementation is relatively faithful to the original Fortran 77 implementation. It should be noted however that the  
indexing of the variables and constraints arrays in the *public* interface is zero-based, i.e. for a problem with 3 variables, x[0],
x[1] and x[2] should be employed in the objective function evaluations. Furthermore, the signature of the objective function
evaluation method is different from the original Fortran implementation.

*LINCOA* solves the following optimization task:

> minimize F(x)<br />
> subject to A.x <= B

where *x* is the *n*-dimensional variable vector, *A* is the *m x n*-dimensional linear constraints matrix and *B* is the *m*-dimensional constraint vector.

To invoke the C# *LINCOA* algorithm, first implement a method for computing the objective function *F* with the following signature:

    double objective(int n, double[] x, bool isXFeasible)

where `n` is the number of variables, `x` is the variable array and `isXFeasible` is a flag indicating whether current `x` is feasible with respect to the linear constraints. 
The method should return the calculated objective function value.

The linear constraint matrix *A* is in the C# implementation represented by the multidimensional array `a`, with a size equal to or larger than `[m, n]`. The constraint vector *B* is represented by the array `b` with size at least `m`.

To minimize the objective function subject to bounds, call the static _Lincoa.FindMinimum_ method:

    Lincoa.Result Lincoa.FindMinimum(Func<int, double[], bool, double> objective, int n, int npt, int m,
            double[,] a, double[] b, double[] x, double rhobeg, double rhoend, int iprint, int maxfun, 
			TextWriter logger)

where `x` on input is the initial variable array, `n` is the number of variables in `x`, `npt` is the number
of interpolation conditions (recommended value `2 * n + 1`), `m` is the number of linear constraints specified in detail by `a` and `b`, `rhobeg` and `rhoend` 
are initial and final values of a trust region radius, `iprint` (0..3) specifies the level of output, `maxfun` is the maximum allowed number of function evaluations, 
and _logger_ is a text writer to where *LINCOA*'s log will be output. 

The method returns the result as an `Lincoa.Result` object containing:

    Lincoa.Status Status   /* Optimization status upon exit from LINCOA */
	int           Evals    /* Total number of function evaluations */
	double        F        /* Optimal value of the objective function
	double[]      X        /* Optimal values of the optimization variables */

The `FindMinimum` method implements default values as follows: `rhobeg = 1.0`, `rhoend = 1.0e-6`, `iprint = 1`, `maxfun = 10000` and `logger = null`. 

#### Fortran 77 README excerpt

> /.../ LINCOA is attached. Its purpose is to seek the least value of a function F of several variables subject to general linear inequality constraints on the variables, when derivatives of F are not available. The name LINCOA denotes LINearly Constrained Optimization Algorithm. F is specified by the user through a subroutine called CALFUN. The algorithm is intended to change the variables to values that are close to a local constrained minimum of F. The user, however, should assume responsibility for finding out if the calculations are adequate. It may be helpful to employ several starting points in the space of the variables and to try different values of the parameters NPT and RHOEND. I intend to write a paper that explains briefly the main features of the software.

> LINCOA is not suitable for very large numbers of variables because no attention is given to any sparsity. A few calculations with 1000 variables, however, have been run successfully overnight, and the performance of LINCOA is satisfactory usually for small numbers of variables. Several calculations of the objective function may be required at points that do not satisfy the linear constraints, especially if an equality constraint is expressed as two inequalities.

> /.../

> In addition to providing CALFUN, the linear constraints, and an initial vector of variables, the user has to set the values of RHOBEG, RHOEND and NPT. After scaling the individual variables, so that the magnitudes of their expected changes are similar, RHOBEG is the initial steplength for changes to the variables, a reasonable choice being the mesh size of a coarse grid search. Further, RHOEND should be suitable for a search on a very fine grid. Typically, the final vector of variables is within distance 10*RHOEND of a local minimum. The parameter NPT specifies the number of interpolation conditions on each quadratic model, the value NPT=2*N+1 being recommended for a start, where N is the number of variables.

> /.../ I hope that the time and effort I have spent on developing the package will be helpful to much and to many applications.

> December 6th, 2013                    M.J.D. Powell (mjdp@cam.ac.uk)

## Links

* [Wikipedia article on LINCOA](http://en.wikipedia.org/wiki/LINCOA)
* [Original Fortran 77 implementation of LINCOA](http://plato.asu.edu/ftp/lincoa.zip)

## Revision history

* February 3, 2014: Initial document.
