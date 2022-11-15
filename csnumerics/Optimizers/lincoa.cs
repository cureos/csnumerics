/*
 *  Copyright (c) 2012-2022, Cureos AB.
 *  All rights reserved.
 *  http://www.cureos.com
 *
 *	This file is part of CSNumerics.
 *
 *  CSNumerics is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as 
 *  published by the Free Software Foundation, either version 3 of the 
 *  License, or (at your option) any later version.
 *
 *  CSNumerics is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with CSNumerics.  If not, see <http://www.gnu.org/licenses/>.
 *  
 *  Remarks:
 * 
 *  The original Fortran 77 version of this code was developed by 
 *  Michael Powell (mjdp@cam.ac.uk) and can be downloaded from this location:
 *  http://plato.asu.edu/ftp/lincoa.zip
 */

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Cureos.Numerics.Optimizers
{
    #region DELEGATES

    /// <summary>
    /// Delegate for the LINCOA objective function formulation.
    /// </summary>
    /// <param name="n">Number of variables.</param>
    /// <param name="x">Variable array.</param>
    /// <param name="constraintsSatisfied">true if <paramref name="x"/> satisfies all constraints within workign accuract, false otherwise.</param>
    /// <returns>Value of the objective function at <paramref name="x"/>.</returns>
    public delegate double LincoaObjectiveFunctionDelegate(int n, double[] x, bool constraintsSatisfied);

    #endregion

    /// <summary>
    /// LINCOA, Derivative-free optimizer of nonlinear objective function with linear constraints.
    /// </summary>
    public class Lincoa : IQuadraticModelOptimizer
    {
// ReSharper disable InconsistentNaming

        #region DELEGATES

        private delegate void LincoaCalfunDelegate(int n, double[] x, ref double f);

        #endregion

        #region INNER TYPES

        private class LincoaCalfunAdapter
        {
            #region FIELDS

            private readonly LincoaObjectiveFunctionDelegate _objective;

            #endregion

            #region CONSTRUCTORS

            internal LincoaCalfunAdapter(LincoaObjectiveFunctionDelegate objective)
            {
                _objective = objective;
            }

            #endregion

            #region METHODS

            internal void CALFUN(int n, double[] xx, ref double f)
            {
                var x = new double[n];
                Array.Copy(xx, 1, x, 0, n);
                f = _objective(n, x, f > 0.0);
            }

            #endregion
        }

        #endregion

        #region FIELDS

        private const double ZERO = 0.0;
        private const double HALF = 0.5;
        private const double ONE = 1.0;
        private const double TENTH = 0.1;
        private const double TINY = 1.0E-60;
        private const double CTEST = 0.01;

        private static readonly string LF = Environment.NewLine;

// ReSharper disable ConvertToConstant.Local
        private static readonly string LINCOA_10 = "Return from LINCOA because N is less than 2.";
        private static readonly string LINCOA_20 = "Return from LINCOA because NPT is not in the required interval.";
        private static readonly string LINCOA_30 = "Return from LINCOA because MAXFUN is less than NPT+1.";
        private static readonly string LINCOA_50 = "Return from LINCOA because the gradient of a constraint is zero.";
        private static readonly string LINCOA_70 = "LINCOA has made the initial X feasible by increasing part(s) of B.";

        private static readonly string LINCOB_230 = "Return from LINCOA because CALFUN has been called MAXFUN times.";
        private static readonly string LINCOB_250 = "Return from LINCOA because rounding errors prevent reasonable changes to X.";
        private static readonly string LINCOB_260 = "Function number {0,6:D}    F = {1,18:E10}    The corresponding X is:" + LF + "{2}";
        private static readonly string LINCOB_320 = "Return from LINCOA because the denominator of the updating formula is zero.";
        private static readonly string LINCOB_570 = LF;
        private static readonly string LINCOB_580 = "New RHO = {0,11:E4}     Number of function values = {1,6:D}";
        private static readonly string LINCOB_590 = "Least value of F = {0,23:E15}         The corresponding X is:" + LF + "{1}";
        private static readonly string LINCOB_620 = "At the return from LINCOA     Number of function values = {0,6:D}";

        private static readonly string PRELIM_140 = LINCOB_260;
// ReSharper restore ConvertToConstant.Local

        #endregion

        #region FIELDS

        private readonly LincoaObjectiveFunctionDelegate _objective;

        private readonly double[,] _a;
        private readonly double[] _b;

        private readonly int _n;
        private readonly int _m;

        private double _rhobeg;
        private double _rhoend;

        private int _npt;
        private int _maxfun;
        private int _iprint;

        private TextWriter _logger;

        #endregion

        #region CONSTRUCTORS

        /// <summary>
        /// Initializes an instance of the LINCOA optimizer.
        /// </summary>
        /// <param name="n">Number of variables.</param>
        /// <param name="m">Number of linear constraints.</param>
        /// <param name="objective">Objective function subject to minimization.</param>
        /// <param name="a">Linear constraints matrix.</param>
        /// <param name="b">Linear constraints vector.</param>
        public Lincoa(int n, int m, LincoaObjectiveFunctionDelegate objective, double[,] a, double[] b)
        {
            _n = n;
            _m = m;
            _objective = objective;

            _a = new double[1 + n, 1 + m];
            for (var j = 0; j < m; ++j)
                for (var i = 0; i < n; ++i)
                    _a[1 + i, 1 + j] = a[j, i];

            _b = new double[1 + m];
            Array.Copy(b, 0, _b, 1, m);

            _npt = 2 * _n + 1;
            _rhobeg = 1.0;
            _rhoend = 1.0e-6;
            _maxfun = 10000;
            _iprint = 0;
            _logger = null;
        }

        /// <summary>
        /// Initializes an instance of the LINCOA optimizer, assuming that number of variables and constraints are given by the size of <paramref name="a"/>.
        /// </summary>
        /// <param name="objective">Objective function subject to minimization.</param>
        /// <param name="a">Linear constraints matrix.</param>
        /// <param name="b">Linear constraints vector.</param>
        public Lincoa(LincoaObjectiveFunctionDelegate objective, double[,] a, double[] b) 
            : this(a.GetLength(1), a.GetLength(0), objective, a, b)
        {
        }

        #endregion

        #region PROPERTIES

        /// <summary>
        /// Gets or sets the number of interpolation conditions.
        /// </summary>
        public int InterpolationConditions
        {
            get { return _npt; }
            set { _npt = value; }
        }

        /// <summary>
        /// Gets or sets the start value of the trust region radius.
        /// </summary>
        public double TrustRegionRadiusStart
        {
            get { return _rhobeg; }
            set { _rhobeg = value; }
        }

        /// <summary>
        /// Gets or sets the final value of the trust region radius.
        /// </summary>
        public double TrustRegionRadiusEnd
        {
            get { return _rhoend; }
            set { _rhoend = value; }
        }

        /// <summary>
        /// Gets or sets the number of maximum function calls.
        /// </summary>
        public int MaximumFunctionCalls
        {
            get { return _maxfun; }
            set { _maxfun = value; }
        }

        /// <summary>
        /// Gets or sets the print level to the logger.
        /// </summary>
        public int PrintLevel {
            get { return _iprint; }
            set { _iprint = value; }
        }

        /// <summary>
        /// Gets or sets the logger to which LINCOA log information should be sent.
        /// </summary>
        public TextWriter Logger
        {
            get { return _logger; }
            set { _logger = value; }
        }


        #endregion

        #region METHODS

        /// <summary>
        /// Find a local minimum of provided objective function satisfying the provided linear constraints.
        /// </summary>
        /// <param name="x0">Initial variable array.</param>
        /// <returns>Summary of the optimization result.</returns>
        public OptimizationSummary FindMinimum(double[] x0)
        {
            var x = new double[1 + _n];
            Array.Copy(x0, 0, x, 1, _n);

            double f;
            int nf;
            var calfun = new LincoaCalfunDelegate(new LincoaCalfunAdapter(_objective).CALFUN);
            var status = LINCOA(calfun, _n, _npt, _m, _a, _b, x, _rhobeg, _rhoend, _iprint, _maxfun, out f, out nf, _logger);

            var xopt = new double[_n];
            Array.Copy(x, 1, xopt, 0, _n);

            return new OptimizationSummary(status, nf, xopt, f);
        }

// ReSharper disable SuggestUseVarKeywordEvident
        private static OptimizationStatus LINCOA(LincoaCalfunDelegate calfun, int n, int npt, int m, double[,] a, double[] b,
            double[] x, double rhobeg, double rhoend, int iprint, int maxfun, out double f, out int nf, TextWriter logger)
        {
            f = Double.MaxValue;
            nf = 0;
//
//     This subroutine seeks the least value of a function of many variables,
//       subject to general linear inequality constraints, by a trust region
//       method that forms quadratic models by interpolation. Usually there
//       is much freedom in each new model after satisfying the interpolation
//       conditions, which is taken up by minimizing the Frobenius norm of
//       the change to the second derivative matrix of the model. One new
//       function value is calculated on each iteration, usually at a point
//       where the current model predicts a reduction in the least value so
//       far of the objective function subject to the linear constraints.
//       Alternatively, a new vector of variables may be chosen to replace
//       an interpolation point that may be too far away for reliability, and
//       then the new point does not have to satisfy the linear constraints.
//       The arguments of the subroutine are as follows.
//
//     N must be set to the number of variables and must be at least two.
//     NPT must be set to the number of interpolation conditions, which is
//       required to be in the interval [N+2,(N+1)(N+2)/2]. Typical choices
//       of the author are NPT=N+6 and NPT=2*N+1. Larger values tend to be
//       highly inefficent when the number of variables is substantial, due
//       to the amount of work and extra difficulty of adjusting more points.
//     M must be set to the number of linear inequality constraints.
//     A is a matrix whose columns are the constraint gradients, which are
//       required to be nonzero.
//     IA is the first dimension of the array A, which must be at least N.
//     B is the vector of right hand sides of the constraints, the J-th
//       constraint being that the scalar product of A(.,J) with X(.) is at
//       most B(J). The initial vector X(.) is made feasible by increasing
//       the value of B(J) if necessary.
//     X is the vector of variables. Initial values of X(1),X(2),...,X(N)
//       must be supplied. If they do not satisfy the constraints, then B
//       is increased as mentioned above. X contains on return the variables
//       that have given the least calculated F subject to the constraints.
//     RHOBEG and RHOEND must be set to the initial and final values of a
//       trust region radius, so both must be positive with RHOEND<=RHOBEG.
//       Typically, RHOBEG should be about one tenth of the greatest expected
//       change to a variable, and RHOEND should indicate the accuracy that
//       is required in the final values of the variables.
//     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
//       amount of printing. Specifically, there is no output if IPRINT=0 and
//       there is output only at the return if IPRINT=1. Otherwise, the best
//       feasible vector of variables so far and the corresponding value of
//       the objective function are printed whenever RHO is reduced, where
//       RHO is the current lower bound on the trust region radius. Further,
//       each new value of F with its variables are output if IPRINT=3.
//     MAXFUN must be set to an upper bound on the number of calls of CALFUN,
//       its value being at least NPT+1.
//     W is an array used for working space. Its length must be at least
//       M*(2+N) + NPT*(4+N+NPT) + N*(9+3*N) + MAX [ M+3*N, 2*M+N, 2*NPT ].
//       On return, W(1) is set to the final value of F, and W(2) is set to
//       the total number of function evaluations plus 0.5.
//
//     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
//       F to the value of the objective function for the variables X(1),
//       X(2),...,X(N). The value of the argument F is positive when CALFUN
//       is called if and only if the current X satisfies the constraints
//       to working accuracy.
//
//     Check that N, NPT and MAXFUN are acceptable.
//
            double smallx = 1.0E-6 * rhoend;
            int np = n + 1;
            if (n <= 1)
            {
                PRINT(logger, LINCOA_10);
                return OptimizationStatus.N_TooSmall;
            }
            if (npt < n + 2 || npt > ((n + 2) * np) / 2)
            {
                PRINT(logger, LINCOA_20);
                return OptimizationStatus.NPT_OutOfRange;
            }
            if (maxfun <= npt)
            {
                PRINT(logger, LINCOA_30);
                return OptimizationStatus.MAXFUN_NotLargerThan_NPT;
            }
//
//     Normalize the constraints, and copy the resultant constraint matrix
//       and right hand sides into working space, after increasing the right
//       hand sides if necessary so that the starting point is feasible.
//
            double[,] amat = new double[1 + n, 1 + m];
            double[] bnorm = new double[1 + m]; // B in LINCOB
            int iflag = 0;
            if (m > 0)
            {
                for (int j = 1; j <= m; ++j)
                {
                    double sum = ZERO;
                    double temp = ZERO;
                    for (int i = 1; i <= n; ++i)
                    {
                        sum += a[i, j] * x[i];
                        temp += a[i, j] * a[i, j];
                    }
                    if (temp == ZERO)
                    {
                        PRINT(logger, LINCOA_50);
                        return OptimizationStatus.ConstraintGradientIsZero;
                    }
                    temp = Math.Sqrt(temp);
                    if (sum - b[j] > smallx * temp) iflag = 1;
                    bnorm[j] = Math.Max(b[j], sum) / temp;
                    for (int i = 1; i <= n; ++i)
                    {
                        amat[i, j] = a[i, j] / temp;
                    }
                }
            }
            if (iflag == 1)
            {
                if (iprint > 0) PRINT(logger, LINCOA_70);
            }
            return LINCOB(calfun, n, npt, m, amat, bnorm, x, rhobeg, rhoend, iprint, maxfun, out f, out nf, logger);
        }

        private static OptimizationStatus LINCOB(LincoaCalfunDelegate calfun, int n, int npt, int m, double[,] amat, double[] b,
            double[] x, double rhobeg, double rhoend, int iprint, int maxfun, out double f, out int nf, TextWriter logger)
        {
            OptimizationStatus? status = null;
//
//     The arguments N, NPT, M, X, RHOBEG, RHOEND, IPRINT and MAXFUN are
//       identical to the corresponding arguments in SUBROUTINE LINCOA.
//     AMAT is a matrix whose columns are the constraint gradients, scaled
//       so that they have unit length.
//     B contains on entry the right hand sides of the constraints, scaled
//       as above, but later B is modified for variables relative to XBASE.
//     XBASE holds a shift of origin that should reduce the contributions
//       from rounding errors to values of the model and Lagrange functions.
//     XPT contains the interpolation point coordinates relative to XBASE.
//     FVAL holds the values of F at the interpolation points.
//     XSAV holds the best feasible vector of variables so far, without any
//       shift of origin.
//     XOPT is set to XSAV-XBASE, which is the displacement from XBASE of
//       the feasible vector of variables that provides the least calculated
//       F so far, this vector being the current trust region centre.
//     GOPT holds the gradient of the quadratic model at XSAV = XBASE+XOPT.
//     HQ holds the explicit second derivatives of the quadratic model.
//     PQ contains the parameters of the implicit second derivatives of the
//       quadratic model.
//     BMAT holds the last N columns of the big inverse matrix H.
//     ZMAT holds the factorization of the leading NPT by NPT submatrix
//       of H, this factorization being ZMAT times Diag(DZ) times ZMAT^T,
//       where the elements of DZ are plus or minus one, as specified by IDZ.
//     NDIM is the first dimension of BMAT and has the value NPT+N.
//     STEP is employed for trial steps from XOPT. It is also used for working
//       space when XBASE is shifted and in PRELIM.
//     SP is reserved for the scalar products XOPT^T XPT(K,.), K=1,2,...,NPT,
//       followed by STEP^T XPT(K,.), K=1,2,...,NPT.
//     XNEW is the displacement from XBASE of the vector of variables for
//       the current calculation of F, except that SUBROUTINE TRSTEP uses it
//       for working space.
//     IACT is an integer array for the indices of the active constraints.
//     RESCON holds useful information about the constraint residuals. Every
//       nonnegative RESCON(J) is the residual of the J-th constraint at the
//       current trust region centre. Otherwise, if RESCON(J) is negative, the
//       J-th constraint holds as a strict inequality at the trust region
//       centre, its residual being at least |RESCON(J)|; further, the value
//       of |RESCON(J)| is at least the current trust region radius DELTA.
//     QFAC is the orthogonal part of the QR factorization of the matrix of
//       active constraint gradients, these gradients being ordered in
//       accordance with IACT. When NACT is less than N, columns are added
//       to QFAC to complete an N by N orthogonal matrix, which is important
//       for keeping calculated steps sufficiently close to the boundaries
//       of the active constraints.
//     RFAC is the upper triangular part of this QR factorization, beginning
//       with the first diagonal element, followed by the two elements in the
//       upper triangular part of the second column and so on.
//     PQW is used for working space, mainly for storing second derivative
//       coefficients of quadratic functions. Its length is NPT+N.
//     The array W is also used for working space. The required number of
//       elements, namely MAX[M+3*N,2*M+N,2*NPT], is set in LINCOA.
//
//     Set some constants.
//
            int np = n + 1;
            int nh = (n * np) / 2;
            int nptm = npt - np;
            int iamat = Math.Max(m + 3 * n, Math.Max(2 * m + n, 2 * npt)) + 1;
            int ndim = npt + n;
//
//     Partition the working space array, so that different parts of it can be
//     treated separately by the subroutine that performs the main calculation.
//
            double[] xbase = new double[1 + n];
            double[,] xpt = new double[1 + npt, 1 + n];
            double[] fval = new double[1 + npt];
            double[] xsav = new double[1 + n];
            double[] xopt = new double[1 + n];
            double[] gopt = new double[1 + n];
            double[] hq = new double[1 + (n * np) / 2];
            double[] pq = new double[1 + npt];
            double[,] bmat = new double[1 + ndim, 1 + n];
            double[,] zmat = new double[1 + npt, 1 + nptm];
            double[] step = new double[1 + n];
            double[] sp = new double[1 + npt + npt];
            double[] xnew = new double[1 + n];
            int[] iact = new int[1 + n];
            double[] rescon = new double[1 + m];
            double[,] qfac = new double[1 + n, 1 + n];
            double[] rfac = new double[1 + (n * np) / 2];
            double[] pqw = new double[1 + npt + n];
            double[] w = new double[iamat];
//
//     Set the elements of XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT,
//       ZMAT and SP for the first iteration. An important feature is that,
//       if the interpolation point XPT(K,.) is not feasible, where K is any
//       integer from [1,NPT], then a change is made to XPT(K,.) if necessary
//       so that the constraint violation is at least 0.2*RHOBEG. Also KOPT
//       is set so that XPT(KOPT,.) is the initial trust region centre.
//
            int kopt, idz;
            PRELIM(calfun, n, npt, m, amat, b, x, rhobeg, iprint, xbase, xpt, fval, xsav, xopt, gopt, out kopt, hq, pq,
                bmat, zmat, out idz, ndim, sp, rescon, logger);
//
//     Begin the iterative procedure.
//
            nf = npt;
            double fopt = fval[kopt];
            double rho = rhobeg;
            double delta = rho;
            int ifeas = 0;
            int nact = 0;
            int itest = 3;

            int knew, nvala, nvalb;
            double fsave, xoptsq;

            LINCOB_10:

            knew = 0;
            nvala = 0;
            nvalb = 0;
//
//     Shift XBASE if XOPT may be too far from XBASE. First make the changes
//       to BMAT that do not depend on ZMAT.
//
            LINCOB_20:

            fsave = fopt;
            xoptsq = ZERO;
            for (int i = 1; i <= n; ++i)
                xoptsq += xopt[i] * xopt[i];
            if (xoptsq >= 1.0E4 * delta * delta)
            {
                double qoptsq = 0.25 * xoptsq;
                for (int k = 1; k <= npt; ++k)
                {
                    double sum = ZERO;
                    for (int I = 1; I <= n; ++I)
                        sum += xpt[k, I] * xopt[I];
                    sum -= HALF * xoptsq;
                    w[npt + k] = sum;
                    sp[k] = ZERO;
                    for (int i = 1; i <= n; ++i)
                    {
                        xpt[k, i] -= HALF * xopt[i];
                        step[i] = bmat[k, i];
                        w[i] = sum * xpt[k, i] + qoptsq * xopt[i];
                        int ip = npt + i;
                        for (int j = 1; j <= i; ++j)
                            bmat[ip, j] += step[i] * w[j] + w[i] * step[j];
                    }
                }
//
//     Then the revisions of BMAT that depend on ZMAT are calculated.
//
                for (int k = 1; k <= nptm; ++k)
                {
                    double sumz = ZERO;
                    for (int i = 1; i <= npt; ++i)
                    {
                        sumz += zmat[i, k];
                        w[i] = w[npt + i] * zmat[i, k];
                    }
                    for (int j = 1; j <= n; ++j)
                    {
                        double sum = qoptsq * sumz * xopt[j];
                        for (int i = 1; i <= npt; ++i)
                            sum += w[i] * xpt[i, j];
                        step[j] = sum;
                        if (k < idz) sum = -sum;
                        for (int i = 1; i <= npt; ++i)
                            bmat[i, j] += sum * zmat[i, k];
                    }
                    for (int i = 1; i <= n; ++i)
                    {
                        int ip = i + npt;
                        double temp = step[i];
                        if (k < idz) temp = -temp;
                        for (int j = 1; j <= i; ++j)
                            bmat[ip, j] += temp * step[j];
                    }
                }
//
//     Update the right hand sides of the constraints.
//
                if (m > 0)
                {
                    for (int j = 1; j <= m; ++j)
                    {
                        double temp = ZERO;
                        for (int i = 1; i <= n; ++i)
                            temp += amat[i, j] * xopt[i];
                        b[j] -= temp;
                    }
                }
//
//     The following instructions complete the shift of XBASE, including the
//       changes to the parameters of the quadratic model.
//
                for (int ih = 0, j = 1; j <= n; ++j)
                {
                    w[j] = ZERO;
                    for (int k = 1; k <= npt; ++k)
                    {
                        w[j] += pq[k] * xpt[k, j];
                        xpt[k, j] -= HALF * xopt[j];
                    }
                    for (int i = 1; i <= j; ++i)
                    {
                        ih++;
                        hq[ih] += w[i] * xopt[j] + xopt[i] * w[j];
                        bmat[npt + i, j] = bmat[npt + j, i];
                    }
                }
                for (int j = 1; j <= n; ++j)
                {
                    xbase[j] += xopt[j];
                    xopt[j] = ZERO;
                    xpt[kopt, j] = ZERO;
                }
            }
//
//     In the case KNEW=0, generate the next trust region step by calling
//       TRSTEP, where SNORM is the current trust region radius initially.
//       The final value of SNORM is the length of the calculated step,
//       except that SNORM is zero on return if the projected gradient is
//       unsuitable for starting the conjugate gradient iterations.
//
            f = ZERO;
            double vquad = ZERO;
            double snorm = ZERO;
            double delsav = delta;
            int ksave = knew;
            if (knew == 0)
            {
                snorm = delta;
                for (int i = 1; i <= n; ++i)
                    xnew[i] = gopt[i];
                TRSTEP(n, npt, m, amat, xpt, hq, pq, ref nact, iact, rescon, qfac, rfac, ref snorm, step, xnew);
//
//     A trust region step is applied whenever its length, namely SNORM, is at
//       least HALF*DELTA. It is also applied if its length is at least 0.1999
//       times DELTA and if a line search of TRSTEP has caused a change to the
//       active set. Otherwise there is a branch below to label 530 or 560.
//
                double temp = HALF * delta;
                if (xnew[1] >= HALF) temp = 0.1999 * delta;
                if (snorm <= temp)
                {
                    delta *= HALF;
                    if (delta <= 1.4 * rho) delta = rho;
                    ++nvala;
                    ++nvalb;
                    temp = snorm / rho;
                    if (delsav > rho) temp = ONE;
                    if (temp >= HALF) nvala = 0;
                    if (temp >= TENTH) nvalb = 0;
                    if (delsav > rho) goto LINCOB_530;
                    if (nvala < 5 && nvalb < 3) goto LINCOB_530;
                    if (snorm > ZERO) ksave = -1;
                    goto LINCOB_560;
                }
                nvala = 0;
                nvalb = 0;
//
//     Alternatively, KNEW is positive. Then the model step is calculated
//       within a trust region of radius DEL, after setting the gradient at
//       XBASE and the second derivative parameters of the KNEW-th Lagrange
//       function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively.
//
            }
            else
            {
                double del = Math.Max(TENTH * delta, rho);
                for (int i = 1; i <= n; ++i)
                    w[i] = bmat[knew, i];
                for (int k = 1; k <= npt; ++k)
                    pqw[k] = ZERO;
                for (int j = 1; j <= nptm; ++j)
                {
                    double temp = zmat[knew, j];
                    if (j < idz) temp = -temp;
                    for (int k = 1; k <= npt; k++)
                        pqw[k] += temp * zmat[k, j];
                }
                QMSTEP(n, npt, m, amat, xpt, xopt, nact, iact, rescon, qfac, kopt, knew, del, step, w, pqw, out ifeas);
            }
//
//     Set VQUAD to the change to the quadratic model when the move STEP is
//       made from XOPT. If STEP is a trust region step, then VQUAD should be
//       negative. If it is nonnegative due to rounding errors in this case,
//       there is a branch to label 530 to try to improve the model.
//
            for (int ih = 0, j = 1; j <= n; ++j)
            {
                vquad = vquad + step[j] * gopt[j];
                for (int i = 1; i <= j; ++i)
                {
                    ++ih;
                    double temp = step[i] * step[j];
                    if (i == j) temp *= HALF;
                    vquad = vquad + temp * hq[ih];
                }
            }
            for (int k = 1; k <= npt; ++k)
            {
                double temp = ZERO;
                for (int j = 1; j <= n; ++j)
                {
                    temp = temp + xpt[k, j] * step[j];
                    sp[npt + k] = temp;
                }
                vquad = vquad + HALF * pq[k] * temp * temp;
            }
            if (ksave == 0 && vquad >= ZERO) goto LINCOB_530;
//
//     Calculate the next value of the objective function. The difference
//       between the actual new value of F and the value predicted by the
//       model is recorded in DIFF.
//
            LINCOB_220:

            ++nf;
            if (nf > maxfun)
            {
                --nf;
                if (iprint > 0) PRINT(logger, LINCOB_230);
                status = OptimizationStatus.MAXFUN_Reached;
                goto LINCOB_600;
            }
            double xdiff = ZERO;
            for (int i = 1; i <= n; ++i)
            {
                xnew[i] = xopt[i] + step[i];
                x[i] = xbase[i] + xnew[i];
                xdiff += Math.Pow(x[i] - xsav[i], 2.0);
            }
            xdiff = Math.Sqrt(xdiff);
            if (ksave == -1) xdiff = rho;
            if (xdiff <= TENTH * rho || xdiff >= delta + delta)
            {
                ifeas = 0;
                if (iprint > 0) PRINT(logger, LINCOB_250);
                status = OptimizationStatus.X_RoundingErrorsPreventUpdate;
                goto LINCOB_600;
            }
            if (ksave <= 0) ifeas = 1;
            f = ifeas;
            calfun(n, x, ref f);
            if (iprint == 3)
                PRINT(logger, LINCOB_260, nf, f, FORMAT("  ", "15:E6", x, 1, n));
            if (ksave == -1) goto LINCOB_600;
            double diff = f - fopt - vquad;
//
//     If X is feasible, then set DFFALT to the difference between the new
//       value of F and the value predicted by the alternative model.
//
            double dffalt = ZERO;
            if (ifeas == 1 && itest < 3)
            {
                for (int k = 1; k <= npt; ++k)
                {
                    pqw[k] = ZERO;
                    w[k] = fval[k] - fval[kopt];
                }
                for (int j = 1; j <= nptm; ++j)
                {
                    double sum = ZERO;
                    for (int i = 1; i <= npt; ++i)
                        sum += w[i] * zmat[i, j];
                    if (j < idz) sum = -sum;
                    for (int k = 1; k <= npt; ++k)
                        pqw[k] = pqw[k] + sum * zmat[k, j];
                }
                double vqalt = ZERO;
                for (int k = 1; k <= npt; ++k)
                {
                    double sum = ZERO;
                    for (int j = 1; j <= n; ++j)
                        sum += bmat[k, j] * step[j];
                    vqalt = vqalt + sum * w[k];
                    vqalt += pqw[k] * sp[npt + k] * (HALF * sp[npt + k] + sp[k]);
                }
                dffalt = f - fopt - vqalt;
            }
            if (itest == 3)
            {
                dffalt = diff;
                itest = 0;
            }
//
//     Pick the next value of DELTA after a trust region step.
//
            double ratio = ZERO;
            if (ksave == 0)
            {
                ratio = (f - fopt) / vquad;
                if (ratio <= TENTH)
                {
                    delta *= HALF;
                }
                else if (ratio <= 0.7)
                {
                    delta = Math.Max(HALF * delta, snorm);
                }
                else
                {
                    double temp = Math.Sqrt(2.0) * delta;
                    delta = Math.Max(HALF * delta, snorm + snorm);
                    delta = Math.Min(delta, temp);
                }
                if (delta <= 1.4 * rho) delta = rho;
            }
//
//     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
//       can be moved. If STEP is a trust region step, then KNEW is zero at
//       present, but a positive value is picked by subroutine UPDATE.
//
            UPDATE(n, npt, xpt, bmat, zmat, idz, sp, step, kopt, ref knew);
            if (knew == 0)
            {
                if (iprint > 0) PRINT(logger, LINCOB_320);
                status = OptimizationStatus.UpdatingFormulaDenominatorZero;
                goto LINCOB_600;
            }
//
//     If ITEST is increased to 3, then the next quadratic model is the
//       one whose second derivative matrix is least subject to the new
//       interpolation conditions. Otherwise the new model is constructed
//       by the symmetric Broyden method in the usual way.
//
            if (ifeas == 1)
            {
                ++itest;
                if (Math.Abs(dffalt) >= TENTH * Math.Abs(diff)) itest = 0;
            }
//
//     Update the second derivatives of the model by the symmetric Broyden
//       method, using PQW for the second derivative parameters of the new
//       KNEW-th Lagrange function. The contribution from the old parameter
//       PQ(KNEW) is included in the second derivative matrix HQ. W is used
//       later for the gradient of the new KNEW-th Lagrange function.       
//
            if (itest < 3)
            {
                for (int k = 1; k <= npt; ++k)
                    pqw[k] = ZERO;
                for (int j = 1; j <= nptm; ++j)
                {
                    double temp = zmat[knew, j];
                    if (temp != ZERO)
                    {
                        if (j < idz) temp = -temp;
                        for (int k = 1; k <= npt; ++k)
                            pqw[k] += temp * zmat[k, j];
                    }
                }
                for (int ih = 0, i = 1; i <= n; ++i)
                {
                    w[i] = bmat[knew, i];
                    double temp = pq[knew] * xpt[knew, i];
                    for (int j = 1; j <= i; ++j)
                    {
                        ++ih;
                        hq[ih] += temp * xpt[knew, j];
                    }
                }
                pq[knew] = ZERO;
                for (int k = 1; k <= npt; ++k)
                    pq[k] += diff * pqw[k];
            }
//
//     Include the new interpolation point with the corresponding updates of
//       SP. Also make the changes of the symmetric Broyden method to GOPT at
//       the old XOPT if ITEST is less than 3.
//
            fval[knew] = f;
            sp[knew] = sp[kopt] + sp[npt + kopt];
            double ssq = ZERO;
            for (int i = 1; i <= n; ++i)
            {
                xpt[knew, i] = xnew[i];
                ssq += step[i] * step[i];
            }
            sp[npt + knew] = sp[npt + kopt] + ssq;
            if (itest < 3)
            {
                for (int k = 1; k <= npt; ++k)
                {
                    double temp = pqw[k] * sp[k];
                    for (int i = 1; i <= n; ++i)
                        w[i] += temp * xpt[k, i];
                }
                for (int i = 1; i <= n; ++i)
                    gopt[i] += +diff * w[i];
            }
//
//     Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the
//       least calculated value so far with a feasible vector of variables.
//
            if (f < fopt && ifeas == 1)
            {
                fopt = f;
                for (int j = 1; j <= n; ++j)
                {
                    xsav[j] = x[j];
                    xopt[j] = xnew[j];
                }
                kopt = knew;
                snorm = Math.Sqrt(ssq);
                for (int j = 1; j <= m; ++j)
                {
                    if (rescon[j] >= delta + snorm)
                    {
                        rescon[j] = snorm - rescon[j];
                    }
                    else
                    {
                        rescon[j] += +snorm;
                        if (rescon[j] + delta > ZERO)
                        {
                            double temp = b[j];
                            for (int I = 1; I <= n; ++I)
                                temp -= xopt[I] * amat[I, j];
                            temp = Math.Max(temp, ZERO);
                            if (temp >= delta) temp = -temp;
                            rescon[j] = temp;
                        }
                    }
                }
                for (int k = 1; k <= npt; ++k)
                    sp[k] = sp[k] + sp[npt + k];
//
//     Also revise GOPT when symmetric Broyden updating is applied.
//
                if (itest < 3)
                {
                    for (int ih = 0, j = 1; j <= n; ++j)
                    {
                        for (int i = 1; i <= j; ++i)
                        {
                            ++ih;
                            if (i < j) gopt[j] += hq[ih] * step[i];
                            gopt[i] += hq[ih] * step[j];
                        }
                    }
                    for (int k = 1; k <= npt; ++k)
                    {
                        double temp = pq[k] * sp[npt + k];
                        for (int i = 1; i <= n; ++i)
                            gopt[i] += temp * xpt[k, i];
                    }
                }
            }
//
//     Replace the current model by the least Frobenius norm interpolant if
//       this interpolant gives substantial reductions in the predictions
//       of values of F at feasible points.
//
            if (itest == 3)
            {
                for (int k = 1; k <= npt; ++k)
                {
                    pq[k] = ZERO;
                    w[k] = fval[k] - fval[kopt];
                }
                for (int j = 1; j <= nptm; ++j)
                {
                    double sum = ZERO;
                    for (int i = 1; i <= npt; ++i)
                        sum += w[i] * zmat[i, j];
                    if (j < idz) sum = -sum;
                    for (int k = 1; k <= npt; ++k)
                        pq[k] = pq[k] + sum * zmat[k, j];
                }
                for (int j = 1; j <= n; ++j)
                {
                    gopt[j] = ZERO;
                    for (int i = 1; i <= npt; ++i)
                        gopt[j] += w[i] * bmat[i, j];
                }
                for (int k = 1; k <= npt; ++k)
                {
                    double temp = pq[k] * sp[k];
                    for (int i = 1; i <= n; ++i)
                        gopt[i] += temp * xpt[k, i];
                }
                for (int ih = 1; ih <= nh; ++ih)
                    hq[ih] = ZERO;
            }
//
//     If a trust region step has provided a sufficient decrease in F, then
//       branch for another trust region calculation. Every iteration that
//       takes a model step is followed by an attempt to take a trust region
//       step.
//
            knew = 0;
            if (ksave > 0) goto LINCOB_20;
            if (ratio >= TENTH) goto LINCOB_20;
//
//     Alternatively, find out if the interpolation points are close enough
//       to the best point so far.
//
            LINCOB_530:

            double distsq = Math.Max(delta * delta, 4.0 * rho * rho);
            for (int k = 1; k <= npt; ++k)
            {
                double sum = ZERO;
                for (int j = 1; j <= n; ++j)
                    sum += Math.Pow(xpt[k, j] - xopt[j], 2.0);
                if (sum > distsq)
                {
                    knew = k;
                    distsq = sum;
                }
            }
//
//     If KNEW is positive, then branch back for the next iteration, which
//       will generate a "model step". Otherwise, if the current iteration
//       has reduced F, or if DELTA was above its lower bound when the last
//       trust region step was calculated, then try a "trust region" step
//       instead.
//
            if (knew > 0) goto LINCOB_20;
            knew = 0;
            if (fopt < fsave) goto LINCOB_20;
            if (delsav > rho) goto LINCOB_20;
//
//     The calculations with the current value of RHO are complete.
//       Pick the next value of RHO.
//
            LINCOB_560:

            if (rho > rhoend)
            {
                delta = HALF * rho;
                if (rho > 250.0 * rhoend)
                {
                    rho *= TENTH;
                }
                else if (rho <= 16.0 * rhoend)
                {
                    rho = rhoend;
                }
                else
                {
                    rho = Math.Sqrt(rho * rhoend);
                }
                delta = Math.Max(delta, rho);
                if (iprint >= 2)
                {
                    if (iprint >= 3) PRINT(logger, LINCOB_570);
                    PRINT(logger, LINCOB_580, rho, nf);
                    PRINT(logger, LINCOB_590, fopt, FORMAT("  ", "15:E6", xbase.Zip(xopt, (xb, xo) => xb + xo), 1, n));
                }
                goto LINCOB_10;
            }
//
//     Return from the calculation, after branching to label 220 for another
//       Newton-Raphson step if it has not been tried before.
//
            if (ksave == -1) goto LINCOB_220;

            LINCOB_600:

            if (fopt <= f || ifeas == 0)
            {
                for (int i = 1; i <= n; ++i)
                    x[i] = xsav[i];
                f = fopt;
            }
            if (iprint >= 1)
            {
                PRINT(logger, LINCOB_620, nf);
                PRINT(logger, LINCOB_590, f, FORMAT("  ", "15:E6", x, 1, n));
            }

            return status.GetValueOrDefault(OptimizationStatus.Normal);
        }

        private static double GETACT(int n, int m, double[,] amat, ref int nact, int[] iact, double[,] qfac,
            double[] rfac, double snorm, double[] resnew, double[] resact, double[] g, double[] dw)
        {
            double[] vlam = new double[1 + n];
            double[] w = new double[1 + n];
//
//     N, M, AMAT, B, NACT, IACT, QFAC and RFAC are the same as the terms
//       with these names in SUBROUTINE LINCOB. The current values must be
//       set on entry. NACT, IACT, QFAC and RFAC are kept up to date when
//       GETACT changes the current active set.
//     SNORM, RESNEW, RESACT, G and DW are the same as the terms with these
//       names in SUBROUTINE TRSTEP. The elements of RESNEW and RESACT are
//       also kept up to date.
//     VLAM and W are used for working space, the vector VLAM being reserved
//       for the Lagrange multipliers of the calculation. Their lengths must
//       be at least N.
//     The main purpose of GETACT is to pick the current active set. It is
//       defined by the property that the projection of -G into the space
//       orthogonal to the active constraint normals is as large as possible,
//       subject to this projected steepest descent direction moving no closer
//       to the boundary of every constraint whose current residual is at most
//       0.2*SNORM. On return, the settings in NACT, IACT, QFAC and RFAC are
//       all appropriate to this choice of active set.
//     Occasionally this projected direction is zero, and then the final value
//       of W(1) is set to zero. Otherwise, the direction itself is returned
//       in DW, and W(1) is set to the square of the length of the direction.
//
//     Set some constants and a temporary VLAM.
//
            double violmx = ZERO;
            double tdel = 0.2 * snorm;
            double ddsav = ZERO;
            for (int i = 1; i <= n; ++i)
            {
                ddsav += g[i] * g[i];
                vlam[i] = ZERO;
            }
            ddsav *= 2.0;
//
//     Set the initial QFAC to the identity matrix in the case NACT=0.
//
            if (nact == 0)
            {
                for (int i = 1; i <= n; ++i)
                {
                    for (int j = 1; j <= n; ++j)
                        qfac[i, j] = ZERO;
                    qfac[i, i] = ONE;
                }
                goto GETACT_100;
            }
//
//     Remove any constraints from the initial active set whose residuals
//       exceed TDEL.
//
            int iflag = 1;
            int ic = nact;

            GETACT_40:

            if (resact[ic] > tdel) goto GETACT_800;

            GETACT_50:

            --ic;
            if (ic > 0) goto GETACT_40;
//
//     Remove any constraints from the initial active set whose Lagrange
//       multipliers are nonnegative, and set the surviving multipliers.
//
            iflag = 2;

            GETACT_60:

            if (nact == 0) goto GETACT_100;
            ic = nact;

            GETACT_70:

            double temp = ZERO;
            for (int i = 1; i <= n; ++i)
                temp += qfac[i, ic] * g[i];
            int idiag = (ic * ic + ic) / 2;
            if (ic < nact)
            {
                int jw = idiag + ic;
                for (int j = ic + 1; j <= nact; ++j)
                {
                    temp -= rfac[jw] * vlam[j];
                    jw += j;
                }
            }
            if (temp >= ZERO) goto GETACT_800;
            vlam[ic] = temp / rfac[idiag];
            --ic;
            if (ic > 0) goto GETACT_70;
//
//     Set the new search direction D. Terminate if the 2-norm of D is zero
//       or does not decrease, or if NACT=N holds. The situation NACT=N
//       occurs for sufficiently large SNORM if the origin is in the convex
//       hull of the constraint gradients.
//
            GETACT_100:

            if (nact == n) goto GETACT_290;
            for (int j = nact + 1; j <= n; ++j)
            {
                w[j] = ZERO;
                for (int i = 1; i <= n; ++i)
                    w[j] = w[j] + qfac[i, j] * g[i];
            }
            double dd = ZERO;
            for (int i = 1; i <= n; ++i)
            {
                dw[i] = ZERO;
                for (int j = nact + 1; j <= n; ++j)
                    dw[i] -= w[j] * qfac[i, j];
                dd += dw[i] * dw[i];
            }
            if (dd >= ddsav) goto GETACT_290;
            if (dd == ZERO) goto GETACT_300;
            ddsav = dd;
            double dnorm = Math.Sqrt(dd);
//
//     Pick the next integer L or terminate, a positive value of L being
//       the index of the most violated constraint. The purpose of CTOL
//       below is to estimate whether a positive value of VIOLMX may be
//       due to computer rounding errors.
//
            int l = 0;
            violmx = ZERO;
            double ctol = ZERO;
            if (m > 0)
            {
                double test = dnorm / snorm;
                for (int j = 1; j <= m; ++j)
                {
                    if (resnew[j] > ZERO && resnew[j] <= tdel)
                    {
                        double sum = ZERO;
                        for (int i = 1; i <= n; ++i)
                            sum += amat[i, j] * dw[i];
                        if (sum > test * resnew[j])
                        {
                            if (sum > violmx)
                            {
                                l = j;
                                violmx = sum;
                            }
                        }
                    }
                }
                temp = 0.01 * dnorm;
                if (violmx > ZERO && violmx < temp)
                {
                    if (nact > 0)
                    {
                        for (int k = 1; k <= nact; ++k)
                        {
                            int j = iact[k];
                            double sum = ZERO;
                            for (int i = 1; i <= n; ++i)
                                sum += dw[i] * amat[i, j];
                            ctol = Math.Max(ctol, Math.Abs(sum));
                        }
                    }
                }
            }
            w[1] = ONE;
            if (l == 0) goto GETACT_300;
            if (violmx <= 10.0 * ctol) goto GETACT_300;
//
//     Apply Givens rotations to the last (N-NACT) columns of QFAC so that
//       the first (NACT+1) columns of QFAC are the ones required for the
//       addition of the L-th constraint, and add the appropriate column
//       to RFAC.
//
            int nactp = nact + 1;
            idiag = (nactp * nactp - nactp) / 2;
            double rdiag = ZERO;
            for (int j = n; j >= 1; --j)
            {
                double sprod = ZERO;
                for (int i = 1; i <= n; ++i)
                    sprod = sprod + qfac[i, j] * amat[i, l];
                if (j <= nact)
                {
                    rfac[idiag + j] = sprod;
                }
                else
                {
                    if (Math.Abs(rdiag) <= 1.0E-20 * Math.Abs(sprod))
                    {
                        rdiag = sprod;
                    }
                    else
                    {
                        temp = Math.Sqrt(sprod * sprod + rdiag * rdiag);
                        double cosv = sprod / temp;
                        double sinv = rdiag / temp;
                        rdiag = temp;
                        for (int i = 1; i <= n; ++i)
                        {
                            temp = cosv * qfac[i, j] + sinv * qfac[i, j + 1];
                            qfac[i, j + 1] = -sinv * qfac[i, j] + cosv * qfac[i, j + 1];
                            qfac[i, j] = temp;
                        }
                    }
                }
            }

            if (rdiag < ZERO)
            {
                for (int i = 1; i <= n; ++i)
                    qfac[i, nactp] = -qfac[i, nactp];
            }
            rfac[idiag + nactp] = Math.Abs(rdiag);
            nact = nactp;
            iact[nact] = l;
            resact[nact] = resnew[l];
            vlam[nact] = ZERO;
            resnew[l] = ZERO;
//
//     Set the components of the vector VMU in W.
//
            GETACT_220:

            w[nact] = ONE / Math.Pow(rfac[(nact * nact + nact) / 2], 2.0);
            if (nact > 1)
            {
                for (int i = nact - 1; i >= 1; --i)
                {
                    idiag = (i * i + i) / 2;
                    int jw = idiag + i;
                    double sum = ZERO;
                    for (int j = i + 1; j <= nact; ++j)
                    {
                        sum -= rfac[jw] * w[j];
                        jw += +j;
                    }
                    w[i] = sum / rfac[idiag];
                }
            }
//
//     Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
//
            double vmult = violmx;
            ic = 0;
            {
                int j = 1;
                while (j < nact)
                {
                    if (vlam[j] >= vmult * w[j])
                    {
                        ic = j;
                        vmult = vlam[j] / w[j];
                    }
                    ++j;
                }
            }
            for (int j = 1; j <= nact; ++j)
                vlam[j] = vlam[j] - vmult * w[j];
            if (ic > 0) vlam[ic] = ZERO;
            violmx = Math.Max(violmx - vmult, ZERO);
            if (ic == 0) violmx = ZERO;
//
//     Reduce the active set if necessary, so that all components of the
//       new VLAM are negative, with resetting of the residuals of the
//       constraints that become inactive.
//
            iflag = 3;
            ic = nact;

            GETACT_270:

            if (vlam[ic] < ZERO) goto GETACT_280;
            resnew[iact[ic]] = Math.Max(resact[ic], TINY);
            goto GETACT_800;

            GETACT_280:

            --ic;
            if (ic > 0) goto GETACT_270;
//
//     Calculate the next VMU if VIOLMX is positive. Return if NACT=N holds,
//       as then the active constraints imply D=0. Otherwise, go to label
//       100, to calculate the new D and to test for termination.
//
            if (violmx > ZERO) goto GETACT_220;
            if (nact < n) goto GETACT_100;

            GETACT_290:

            dd = ZERO;

            GETACT_300:

            return dd;
//
//     These instructions rearrange the active constraints so that the new
//       value of IACT(NACT) is the old value of IACT(IC). A sequence of
//       Givens rotations is applied to the current QFAC and RFAC. Then NACT
//       is reduced by one.
//
            GETACT_800:

            resnew[iact[ic]] = Math.Max(resact[ic], TINY);
            int jc = ic;
            while (jc < nact)
            {
                int jcp = jc + 1;
                idiag = jc * jcp / 2;
                int jw = idiag + jcp;
                temp = Math.Sqrt(rfac[jw - 1] * rfac[jw - 1] + rfac[jw] * rfac[jw]);
                double cval = rfac[jw] / temp;
                double sval = rfac[jw - 1] / temp;
                rfac[jw - 1] = sval * rfac[idiag];
                rfac[jw] = cval * rfac[idiag];
                rfac[idiag] = temp;
                if (jcp < nact)
                {
                    for (int j = jcp + 1; j <= nact; ++j)
                    {
                        temp = sval * rfac[jw + jc] + cval * rfac[jw + jcp];
                        rfac[jw + jcp] = cval * rfac[jw + jc] - sval * rfac[jw + jcp];
                        rfac[jw + jc] = temp;
                        jw += j;
                    }
                }
                int jdiag = idiag - jc;
                for (int i = 1; i <= n; ++i)
                {
                    if (i < jc)
                    {
                        temp = rfac[idiag + i];
                        rfac[idiag + i] = rfac[jdiag + i];
                        rfac[jdiag + i] = temp;
                    }
                    temp = sval * qfac[i, jc] + cval * qfac[i, jcp];
                    qfac[i, jcp] = cval * qfac[i, jc] - sval * qfac[i, jcp];
                    qfac[i, jc] = temp;
                }
                iact[jc] = iact[jcp];
                resact[jc] = resact[jcp];
                vlam[jc] = vlam[jcp];
                jc = jcp;
            }
            --nact;
            switch (iflag)
            {
                case 1:
                    goto GETACT_50;
                case 2:
                    goto GETACT_60;
                case 3:
                    goto GETACT_280;
                default:
                    throw new InvalidOperationException("Invalid IFLAG value");
            }
        }

        private static void PRELIM(LincoaCalfunDelegate calfun, int n, int npt, int m, double[,] amat, double[] b, double[] x,
            double rhobeg, int iprint, double[] xbase, double[,] xpt, double[] fval, double[] xsav, double[] xopt,
            double[] gopt, out int kopt, double[] hq, double[] pq, double[,] bmat, double[,] zmat, out int idz, int ndim,
            double[] sp, double[] rescon, TextWriter logger)
        {
            double[] step = new double[1 + n];
            double[] w = new double[1 + npt + n];
//
//     The arguments N, NPT, M, AMAT, B, X, RHOBEG, IPRINT, XBASE, XPT, FVAL,
//       XSAV, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SP and RESCON are the
//       same as the corresponding arguments in SUBROUTINE LINCOB.
//     KOPT is set to the integer such that XPT(KOPT,.) is the initial trust
//       region centre.
//     IDZ is going to be set to one, so that every element of Diag(DZ) is
//       one in the product ZMAT times Diag(DZ) times ZMAT^T, which is the
//       factorization of the leading NPT by NPT submatrix of H.
//     STEP, PQW and W are used for working space, the arrays STEP and PQW
//       being taken from LINCOB. The length of W must be at least N+NPT.
//
//     SUBROUTINE PRELIM provides the elements of XBASE, XPT, BMAT and ZMAT
//       for the first iteration, an important feature being that, if any of
//       of the columns of XPT is an infeasible point, then the largest of
//       the constraint violations there is at least 0.2*RHOBEG. It also sets
//       the initial elements of FVAL, XOPT, GOPT, HQ, PQ, SP and RESCON.
//
//     Set some constants.
//
            int nptm = npt - n - 1;
            double rhosq = rhobeg * rhobeg;
            double recip = ONE / rhosq;
            double reciq = Math.Sqrt(HALF) / rhosq;
            double test = 0.2 * rhobeg;
            kopt = 0;
            idz = 1;
            const int kbase = 1;
//
//     Set the initial elements of XPT, BMAT, SP and ZMAT to zero. 
//
            for (int j = 1; j <= n; ++j)
            {
                xbase[j] = x[j];
                for (int k = 1; k <= npt; ++k)
                    xpt[k, j] = ZERO;
                for (int i = 1; i <= ndim; ++i)
                    bmat[i, j] = ZERO;
            }
            for (int k = 1; k <= npt; ++k)
            {
                sp[k] = ZERO;
                for (int j = 1; j <= npt - n - 1; ++j)
                    zmat[k, j] = ZERO;
            }
//
//     Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT],
//       but they may be altered later to make a constraint violation
//       sufficiently large. The initial nonzero elements of BMAT and of
//       the first min[N,NPT-N-1] columns of ZMAT are set also.
//
            for (int j = 1; j <= n; ++j)
            {
                xpt[j + 1, j] = rhobeg;
                if (j < npt - n)
                {
                    int jp = n + j + 1;
                    xpt[jp, j] = -rhobeg;
                    bmat[j + 1, j] = HALF / rhobeg;
                    bmat[jp, j] = -HALF / rhobeg;
                    zmat[1, j] = -reciq - reciq;
                    zmat[j + 1, j] = reciq;
                    zmat[jp, j] = reciq;
                }
                else
                {
                    bmat[1, j] = -ONE / rhobeg;
                    bmat[j + 1, j] = ONE / rhobeg;
                    bmat[npt + j, j] = -HALF * rhosq;
                }
            }
//
//     Set the remaining initial nonzero elements of XPT and ZMAT when the
//       number of interpolation points exceeds 2*N+1.
//
            if (npt > 2 * n + 1)
            {
                for (int k = n + 1; k <= npt - n - 1; ++k)
                {
                    int itemp = (k - 1) / n;
                    int ipt = k - itemp * n;
                    int jpt = ipt + itemp;
                    if (jpt > n) jpt -= n;
                    xpt[n + k + 1, ipt] = rhobeg;
                    xpt[n + k + 1, jpt] = rhobeg;
                    zmat[1, k] = recip;
                    zmat[ipt + 1, k] = -recip;
                    zmat[jpt + 1, k] = -recip;
                    zmat[n + k + 1, k] = recip;
                }
            }
//
//     Update the constraint right hand sides to allow for the shift XBASE.
//
            if (m > 0)
            {
                for (int j = 1; j <= m; ++j)
                {
                    double temp = ZERO;
                    for (int i = 1; i <= n; ++i)
                        temp += amat[i, j] * xbase[i];
                    b[j] -= temp;
                }
            }
//
//     Go through the initial points, shifting every infeasible point if
//       necessary so that its constraint violation is at least 0.2*RHOBEG.
//
            for (int nf = 1; nf <= npt; ++nf)
            {
                double feas = ONE;
                double bigv = ZERO;
                int jsav = 0;
                {
                    int j = 0;

                    PRELIM_80:

                    ++j;
                    if (j <= m && nf >= 2)
                    {
                        double resid = -b[j];
                        for (int i = 1; i <= n; ++i)
                            resid = resid + xpt[nf, i] * amat[i, j];
                        if (resid <= bigv) goto PRELIM_80;
                        bigv = resid;
                        jsav = j;
                        if (resid <= test)
                        {
                            feas = -ONE;
                            goto PRELIM_80;
                        }
                        feas = ZERO;
                    }
                }
                if (feas < ZERO)
                {
                    for (int i = 1; i <= n; ++i)
                        step[i] = xpt[nf, i] + (test - bigv) * amat[i, jsav];
                    for (int k = 1; k <= npt; ++k)
                    {
                        sp[npt + k] = ZERO;
                        for (int j = 1; j <= n; ++j)
                            sp[npt + k] += xpt[k, j] * step[j];
                    }
                    UPDATE(n, npt, xpt, bmat, zmat, idz, sp, step, kbase, ref nf);
                    for (int i = 1; i <= n; ++i)
                        xpt[nf, i] = step[i];
                }
//
//     Calculate the objective function at the current interpolation point,
//       and set KOPT to the index of the first trust region centre.
//
                for (int j = 1; j <= n; ++j)
                    x[j] = xbase[j] + xpt[nf, j];
                double f = feas;
                calfun(n, x, ref f);
                if (iprint == 3)
                {
                    PRINT(logger, PRELIM_140, nf, f, FORMAT("  ", "15:E6", x, 1, n));
                }
                if (nf == 1)
                {
                    kopt = 1;
                }
                else if (f < fval[kopt] && feas > ZERO)
                {
                    kopt = nf;
                }
                fval[nf] = f;
            }
//
//     Set PQ for the first quadratic model.
//
            for (int j = 1; j <= nptm; ++j)
            {
                w[j] = ZERO;
                for (int k = 1; k <= npt; ++k)
                    w[j] += zmat[k, j] * fval[k];
            }
            for (int k = 1; k <= npt; ++k)
            {
                pq[k] = ZERO;
                for (int j = 1; j <= nptm; ++j)
                    pq[k] += zmat[k, j] * w[j];
            }
//
//     Set XOPT, SP, GOPT and HQ for the first quadratic model.
//
            for (int j = 1; j <= n; ++j)
            {
                xopt[j] = xpt[kopt, j];
                xsav[j] = xbase[j] + xopt[j];
                gopt[j] = ZERO;
            }
            for (int k = 1; k <= npt; ++k)
            {
                sp[k] = ZERO;
                for (int j = 1; j <= n; ++j)
                    sp[k] += xpt[k, j] * xopt[j];
                double temp = pq[k] * sp[k];
                for (int j = 1; j <= n; ++j)
                    gopt[j] += fval[k] * bmat[k, j] + temp * xpt[k, j];
            }
            for (int i = 1; i <= (n * n + n) / 2; ++i)
                hq[i] = ZERO;
//
//     Set the initial elements of RESCON.
//
            for (int j = 1; j <= m; ++j)
            {
                double temp = b[j];
                for (int i = 1; i <= n; ++i)
                    temp -= xopt[i] * amat[i, j];
                temp = Math.Max(temp, ZERO);
                if (temp >= rhobeg) temp = -temp;
                rescon[j] = temp;
            }
        }

        private static void QMSTEP(int n, int npt, int m, double[,] amat, double[,] xpt, double[] xopt, int nact,
            int[] iact, double[] rescon, double[,] qfac, int kopt, int knew, double del, double[] step, double[] gl,
            double[] pqw, out int ifeas)
        {
            double[] rstat = new double[1 + m];
            double[] w = new double[1 + n];
//
//     N, NPT, M, AMAT, B, XPT, XOPT, NACT, IACT, RESCON, QFAC, KOPT are the
//       same as the terms with these names in SUBROUTINE LINCOB.
//     KNEW is the index of the interpolation point that is going to be moved.
//     DEL is the current restriction on the length of STEP, which is never
//       greater than the current trust region radius DELTA.
//     STEP will be set to the required step from XOPT to the new point.
//     GL must be set on entry to the gradient of LFUNC at XBASE, where LFUNC
//       is the KNEW-th Lagrange function. It is used also for some other
//       gradients of LFUNC.
//     PQW provides the second derivative parameters of LFUNC.
//     RSTAT and W are used for working space. Their lengths must be at least
//       M and N, respectively. RSTAT(J) is set to -1.0, 0.0, or 1.0 if the
//       J-th constraint is irrelevant, active, or both inactive and relevant,
//       respectively.
//     IFEAS will be set to 0 or 1 if XOPT+STEP is infeasible or feasible.
//
//     STEP is chosen to provide a relatively large value of the modulus of
//       LFUNC(XOPT+STEP), subject to ||STEP|| <= DEL. A projected STEP is
//       calculated too, within the trust region, that does not alter the
//       residuals of the active constraints. The projected step is preferred
//       if its value of | LFUNC(XOPT+STEP) | is at least one fifth of the
//       original one, but the greatest violation of a linear constraint must
//       be at least 0.2*DEL, in order to keep the interpolation points apart.
//       The remedy when the maximum constraint violation is too small is to
//       restore the original step, which is perturbed if necessary so that
//       its maximum constraint violation becomes 0.2*DEL.
//
//     Set some constants.
//
            double test = 0.2 * del;
//
//     Replace GL by the gradient of LFUNC at the trust region centre, and
//       set the elements of RSTAT.
//
            for (int k = 1; k <= npt; ++k)
            {
                double temp = ZERO;
                for (int j = 1; j <= n; ++j)
                    temp += xpt[k, j] * xopt[j];
                temp *= pqw[k];
                for (int i = 1; i <= n; ++i)
                    gl[i] = gl[i] + temp * xpt[k, i];
            }
            if (m > 0)
            {
                for (int j = 1; j <= m; ++j)
                {
                    rstat[j] = ONE;
                    if (Math.Abs(rescon[j]) >= del) rstat[j] = -ONE;
                }
                for (int k = 1; k <= nact; ++k)
                    rstat[iact[k]] = ZERO;
            }
//
//     Find the greatest modulus of LFUNC on a line through XOPT and
//       another interpolation point within the trust region.
//
            int ksav = 0;
            double stpsav = ZERO;
            double vbig = ZERO;
            for (int k = 1; k <= npt; ++k)
            {
                if (k == kopt) continue;
                double ss = ZERO;
                double sp = ZERO;
                for (int i = 1; i <= n; ++i)
                {
                    double temp = xpt[k, i] - xopt[i];
                    ss += temp * temp;
                    sp += gl[i] * temp;
                }
                double stp = -del / Math.Sqrt(ss);
                double vlag;
                if (k == knew)
                {
                    if (sp * (sp - ONE) < ZERO) stp = -stp;
                    vlag = Math.Abs(stp * sp) + stp * stp * Math.Abs(sp - ONE);
                }
                else
                {
                    vlag = Math.Abs(stp * (ONE - stp) * sp);
                }
                if (vlag > vbig)
                {
                    ksav = k;
                    stpsav = stp;
                    vbig = vlag;
                }
            }
//
//     Set STEP to the move that gives the greatest modulus calculated above.
//       This move may be replaced by a steepest ascent step from XOPT.
//
            double gg = ZERO;
            for (int i = 1; i <= n; ++i)
            {
                gg += gl[i] * gl[i];
                step[i] = stpsav * (xpt[ksav, i] - xopt[i]);
            }
            double vgrad = del * Math.Sqrt(gg);
            if (vgrad <= TENTH * vbig) goto QMSTEP_220;
//
//     Make the replacement if it provides a larger value of VBIG.
//
            double ghg = ZERO;
            for (int k = 1; k <= npt; ++k)
            {
                double temp = ZERO;
                for (int j = 1; j <= n; ++j)
                    temp += xpt[k, j] * gl[j];
                ghg += pqw[k] * temp * temp;
            }
            double vnew = vgrad + Math.Abs(HALF * del * del * ghg / gg);
            if (vnew > vbig)
            {
                vbig = vnew;
                double stp = del / Math.Sqrt(gg);
                if (ghg < ZERO) stp = -stp;
                for (int i = 1; i <= n; ++i)
                    step[i] = stp * gl[i];
            }
            if (nact == 0 || nact == n) goto QMSTEP_220;
//
//     Overwrite GL by its projection. Then set VNEW to the greatest
//       value of |LFUNC| on the projected gradient from XOPT subject to
//       the trust region bound. If VNEW is sufficiently large, then STEP
//       may be changed to a move along the projected gradient.
//
            for (int k = nact + 1; k <= n; ++k)
            {
                w[k] = ZERO;
                for (int i = 1; i <= n; ++i)
                    w[k] += gl[i] * qfac[i, k];
            }
            gg = ZERO;
            for (int i = 1; i <= n; ++i)
            {
                gl[i] = ZERO;
                for (int k = nact + 1; k <= n; ++k)
                    gl[i] += qfac[i, k] * w[k];
                gg += gl[i] * gl[i];
            }
            vgrad = del * Math.Sqrt(gg);
            if (vgrad <= TENTH * vbig) goto QMSTEP_220;
            ghg = ZERO;
            for (int k = 1; k <= npt; ++k)
            {
                double temp = ZERO;
                for (int j = 1; j <= n; ++j)
                    temp += xpt[k, j] * gl[j];
                ghg += pqw[k] * temp * temp;
            }
            vnew = vgrad + Math.Abs(HALF * del * del * ghg / gg);
//
//     Set W to the possible move along the projected gradient.
//
            double ww = ZERO;
            {
                double stp = del / Math.Sqrt(gg);
                if (ghg < ZERO) stp = -stp;
                for (int i = 1; i <= n; ++i)
                {
                    w[i] = stp * gl[i];
                    ww += w[i] * w[i];
                }
            }
//
//     Set STEP to W if W gives a sufficiently large value of the modulus
//       of the Lagrange function, and if W either preserves feasibility
//       or gives a constraint violation of at least 0.2*DEL. The purpose
//       of CTOL below is to provide a check on feasibility that includes
//       a tolerance for contributions from computer rounding errors.
//
            if (vnew / vbig >= 0.2)
            {
                ifeas = 1;
                double bigv = ZERO;
                int j = 0;

                QMSTEP_170:

                ++j;
                if (j <= m)
                {
                    if (rstat[j] == ONE)
                    {
                        double temp = -rescon[j];
                        for (int i = 1; i <= n; ++i)
                            temp += w[i] * amat[i, j];
                        bigv = Math.Max(bigv, temp);
                    }
                    if (bigv < test) goto QMSTEP_170;
                    ifeas = 0;
                }
                double ctol = ZERO;
                {
                    double temp = 0.01 * Math.Sqrt(ww);
                    if (bigv > ZERO && bigv < temp)
                    {
                        for (int k = 1; k <= nact; ++k)
                        {
                            j = iact[k];
                            double sum = ZERO;
                            for (int i = 1; i <= n; ++i)
                                sum = sum + w[i] * amat[i, j];
                            ctol = Math.Max(ctol, Math.Abs(sum));
                        }
                    }
                }
                if (bigv <= 10.0 * ctol || bigv >= test)
                {
                    for (int i = 1; i <= n; ++i)
                        step[i] = w[i];
                    return;
                }
            }
//
//     Calculate the greatest constraint violation at XOPT+STEP with STEP at
//       its original value. Modify STEP if this violation is unacceptable.
//
            QMSTEP_220:

            ifeas = 1;
            int jsav = 0;
            double resmax = ZERO;
            {
                int j = 0;
                double bigv = ZERO;

                QMSTEP_230:

                ++j;
                if (j <= m)
                {
                    if (rstat[j] < ZERO) goto QMSTEP_230;
                    double temp = -rescon[j];
                    for (int i = 1; i <= n; ++i)
                        temp += step[i] * amat[i, j];
                    resmax = Math.Max(resmax, temp);
                    if (temp < test)
                    {
                        if (temp <= bigv) goto QMSTEP_230;
                        bigv = temp;
                        jsav = j;
                        ifeas = -1;
                        goto QMSTEP_230;
                    }
                    ifeas = 0;
                }
                if (ifeas == -1)
                {
                    for (int i = 1; i <= n; ++i)
                        step[i] += (test - bigv) * amat[i, jsav];
                    ifeas = 0;
                }
            }
//
//     Return the calculated STEP and the value of IFEAS.
//
        }

        private static void TRSTEP(int n, int npt, int m, double[,] amat, double[,] xpt, double[] hq,
            double[] pq, ref int nact, int[] iact, double[] rescon, double[,] qfac, double[] rfac, ref double snorm,
            double[] step, double[] g)
        {
            double[] resnew = new double[1 + m];
            double[] resact = new double[1 + n];
            double[] d = new double[1 + n];
            double[] dw = new double[1 + n];
            double[] w = new double[1 + Math.Max(m, 2 * n)];
//
//     N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON, QFAC and RFAC
//       are the same as the terms with these names in LINCOB. If RESCON(J)
//       is negative, then |RESCON(J)| must be no less than the trust region
//       radius, so that the J-th constraint can be ignored.
//     SNORM is set to the trust region radius DELTA initially. On the
//       return, however, it is the length of the calculated STEP, which is
//       set to zero if the constraints do not allow a long enough step.
//     STEP is the total calculated step so far from the trust region centre,
//       its final value being given by the sequence of CG iterations, which
//       terminate if the trust region boundary is reached.
//     G must be set on entry to the gradient of the quadratic model at the
//       trust region centre. It is used as working space, however, and is
//       always the gradient of the model at the current STEP, except that
//       on return the value of G(1) is set to ONE instead of to ZERO if
//       and only if GETACT is called more than once.
//     RESNEW, RESACT, D, DW and W are used for working space. A negative
//       value of RESNEW(J) indicates that the J-th constraint does not
//       restrict the CG steps of the current trust region calculation, a
//       zero value of RESNEW(J) indicates that the J-th constraint is active,
//       and otherwise RESNEW(J) is set to the greater of TINY and the actual
//       residual of the J-th constraint for the current STEP. RESACT holds
//       the residuals of the active constraints, which may be positive.
//       D is the search direction of each line search. DW is either another
//       search direction or the change in gradient along D. The length of W
//       must be at least MAX[M,2*N].
//
//     Set some numbers for the conjugate gradient iterations.
//
            double snsq = snorm * snorm;
//
//     Set the initial elements of RESNEW, RESACT and STEP.
//
            if (m > 0)
            {
                for (int j = 1; j <= m; ++j)
                {
                    resnew[j] = rescon[j];
                    if (rescon[j] >= snorm)
                    {
                        resnew[j] = -ONE;
                    }
                    else if (rescon[j] >= ZERO)
                    {
                        resnew[j] = Math.Max(resnew[j], TINY);
                    }
                }
                if (nact > 0)
                {
                    for (int k = 1; k <= nact; ++k)
                    {
                        resact[k] = rescon[iact[k]];
                        resnew[iact[k]] = ZERO;
                    }
                }
            }
            for (int i = 1; i <= n; ++i)
                step[i] = ZERO;
            double ss = ZERO;
            double reduct = ZERO;
            int ncall = 0;
//
//     GETACT picks the active set for the current STEP. It also sets DW to
//       the vector closest to -G that is orthogonal to the normals of the
//       active constraints. DW is scaled to have length 0.2*SNORM, as then
//       a move of DW from STEP is allowed by the linear constraints.
//
            TRSTEP_40:

            ++ncall;
            double dsq = GETACT(n, m, amat, ref nact, iact, qfac, rfac, snorm, resnew, resact, g, dw);
            if (dsq == ZERO) goto TRSTEP_320;
            double scale = 0.2 * snorm / Math.Sqrt(dsq);
            for (int i = 1; i <= n; ++i)
                dw[i] *= scale;
//
//     If the modulus of the residual of an active constraint is substantial,
//       then set D to the shortest move from STEP to the boundaries of the
//       active constraints.
//
            double resmax = ZERO;
            if (nact > 0)
            {
                for (int k = 1; k <= nact; ++k)
                    resmax = Math.Max(resmax, resact[k]);
            }
            double gamma = ZERO;
            if (resmax > 1.0E-4 * snorm)
            {
                int ir = 0;
                for (int k = 1; k <= nact; ++k)
                {
                    double temp = resact[k];
                    if (k >= 2)
                    {
                        for (int i = 1; i <= k - 1; ++i)
                        {
                            ++ir;
                            temp -= rfac[ir] * w[i];
                        }
                    }
                    ++ir;
                    w[k] = temp / rfac[ir];
                }
                for (int i = 1; i <= n; ++i)
                {
                    d[i] = ZERO;
                    for (int k = 1; k <= nact; ++k)
                        d[i] += w[k] * qfac[i, k];
                }
//
//     The vector D that has just been calculated is also the shortest move
//       from STEP+DW to the boundaries of the active constraints. Set GAMMA
//       to the greatest steplength of this move that satisfies the trust
//       region bound.
//
                double rhs = snsq;
                double ds = ZERO;
                double dd = ZERO;
                for (int i = 1; i <= n; ++i)
                {
                    double sum = step[i] + dw[i];
                    rhs -= sum * sum;
                    ds += d[i] * sum;
                    dd += d[i] * d[i];
                }
                if (rhs > ZERO)
                {
                    double temp = Math.Sqrt(ds * ds + dd * rhs);
                    if (ds <= ZERO)
                    {
                        gamma = (temp - ds) / dd;
                    }
                    else
                    {
                        gamma = rhs / (temp + ds);
                    }
                }
//
//     Reduce the steplength GAMMA if necessary so that the move along D
//       also satisfies the linear constraints.
//
                {
                    int j = 0;

                    TRSTEP_110:

                    if (gamma > ZERO)
                    {
                        ++j;
                        if (resnew[j] > ZERO)
                        {
                            double ad = ZERO;
                            double adw = ZERO;
                            for (int i = 1; i <= n; ++i)
                            {
                                ad += amat[i, j] * d[i];
                                adw += amat[i, j] * dw[i];
                            }
                            if (ad > ZERO)
                            {
                                double temp = Math.Max((resnew[j] - adw) / ad, ZERO);
                                gamma = Math.Min(gamma, temp);
                            }
                        }
                        if (j < m) goto TRSTEP_110;
                    }
                }
                gamma = Math.Min(gamma, ONE);
            }
//
//     Set the next direction for seeking a reduction in the model function
//       subject to the trust region bound and the linear constraints.
//
            int icount;
            if (gamma <= ZERO)
            {
                for (int i = 1; i <= n; ++i)
                    d[i] = dw[i];
                icount = nact;
            }
            else
            {
                for (int i = 1; i <= n; ++i)
                    d[i] = dw[i] + gamma * d[i];
                icount = nact - 1;
            }
            double alpbd = ONE;
//
//     Set ALPHA to the steplength from STEP along D to the trust region
//       boundary. Return if the first derivative term of this step is
//       sufficiently small or if no further progress is possible.
//
            TRSTEP_150:

            ++icount;
            double alpha;
            double dg = ZERO;
            {
                double rhs = snsq - ss;
                if (rhs <= ZERO) goto TRSTEP_320;
                double ds = ZERO;
                double dd = ZERO;
                for (int I = 1; I <= n; ++I)
                {
                    dg += d[I] * g[I];
                    ds += d[I] * step[I];
                    dd += d[I] * d[I];
                }
                if (dg >= ZERO) goto TRSTEP_320;
                double temp = Math.Sqrt(rhs * dd + ds * ds);
                if (ds <= ZERO)
                {
                    alpha = (temp - ds) / dd;
                }
                else
                {
                    alpha = rhs / (temp + ds);
                }
                if (-alpha * dg <= CTEST * reduct) goto TRSTEP_320;
            }
//
//     Set DW to the change in gradient along D.
//
            int ih = 0;
            for (int j = 1; j <= n; ++j)
            {
                dw[j] = ZERO;
                for (int i = 1; i <= j; ++i)
                {
                    ++ih;
                    if (i < j) dw[j] += hq[ih] * d[i];
                    dw[i] += hq[ih] * d[j];
                }
            }
            for (int k = 1; k <= npt; ++k)
            {
                double temp = ZERO;
                for (int j = 1; j <= n; ++j)
                    temp += xpt[k, j] * d[j];
                temp *= pq[k];
                for (int i = 1; i <= n; ++i)
                    dw[i] += temp * xpt[k, i];
            }
//
//     Set DGD to the curvature of the model along D. Then reduce ALPHA if
//       necessary to the value that minimizes the model.
//
            double dgd = ZERO;
            for (int i = 1; i <= n; ++i)
                dgd += d[i] * dw[i];
            double alpht = alpha;
            if (dg + alpha * dgd > ZERO)
            {
                alpha = -dg / dgd;
            }
//
//     Make a further reduction in ALPHA if necessary to preserve feasibility,
//       and put some scalar products of D with constraint gradients in W.
//
            double alphm = alpha;
            int jsav = 0;
            if (m > 0)
            {
                for (int j = 1; j <= m; ++j)
                {
                    double ad = ZERO;
                    if (resnew[j] > ZERO)
                    {
                        for (int i = 1; i <= n; ++i)
                            ad += amat[i, j] * d[i];
                        if (alpha * ad > resnew[j])
                        {
                            alpha = resnew[j] / ad;
                            jsav = j;
                        }
                    }
                    w[j] = ad;
                }
            }
            alpha = Math.Max(alpha, alpbd);
            alpha = Math.Min(alpha, alphm);
            if (icount == nact) alpha = Math.Min(alpha, ONE);
//
//     Update STEP, G, RESNEW, RESACT and REDUCT.
//
            ss = ZERO;
            for (int i = 1; i <= n; ++i)
            {
                step[i] += alpha * d[i];
                ss += step[i] * step[i];
                g[i] += alpha * dw[i];
            }
            if (m > 0)
            {
                for (int j = 1; j <= m; ++j)
                {
                    if (resnew[j] > ZERO)
                    {
                        resnew[j] = Math.Max(resnew[j] - alpha * w[j], TINY);
                    }
                }
            }
            if (icount == nact && nact > 0)
            {
                for (int k = 1; k <= nact; ++k)
                    resact[k] *= (ONE - gamma);
            }
            reduct -= alpha * (dg + HALF * alpha * dgd);
//
//     Test for termination. Branch to label 40 if there is a new active
//       constraint and if the distance from STEP to the trust region
//       boundary is at least 0.2*SNORM.
//
            if (alpha == alpht) goto TRSTEP_320;
            {
                double temp = -alphm * (dg + HALF * alphm * dgd);
                if (temp <= CTEST * reduct) goto TRSTEP_320;
            }
            if (jsav > 0)
            {
                if (ss <= 0.64 * snsq) goto TRSTEP_40;
                goto TRSTEP_320;
            }
            if (icount == n) goto TRSTEP_320;
//
//     Calculate the next search direction, which is conjugate to the
//       previous one except in the case ICOUNT=NACT.
//
            if (nact > 0)
            {
                for (int j = nact + 1; j <= n; ++j)
                {
                    w[j] = ZERO;
                    for (int i = 1; i <= n; ++i)
                        w[j] += g[i] * qfac[i, j];
                }
                for (int i = 1; i <= n; ++i)
                {
                    double temp = ZERO;
                    for (int J = nact + 1; J <= n; ++J)
                        temp += qfac[i, J] * w[J];
                    w[n + i] = temp;
                }
            }
            else
            {
                for (int i = 1; i <= n; ++i)
                    w[n + i] = g[i];
            }
            double beta;
            if (icount == nact)
            {
                beta = ZERO;
            }
            else
            {
                double wgd = ZERO;
                for (int i = 1; i <= n; ++i)
                    wgd += w[n + i] * dw[i];
                beta = wgd / dgd;
            }
            for (int i = 1; i <= n; ++i)
                d[i] = -w[n + i] + beta * d[i];
            alpbd = ZERO;
            goto TRSTEP_150;
//
//     Return from the subroutine.
//
            TRSTEP_320:

            snorm = ZERO;
            if (reduct > ZERO) snorm = Math.Sqrt(ss);
            g[1] = ZERO;
            if (ncall > 1) g[1] = ONE;
        }

        private static void UPDATE(int n, int npt, double[,] xpt, double[,] bmat, double[,] zmat, int idz,
            double[] sp, double[] step, int kopt, ref int knew)
        {
            double[] vlag = new double[1 + npt + n];
            double[] w = new double[1 + npt + n];
//
//     The arguments N, NPT, XPT, BMAT, ZMAT, IDZ, NDIM ,SP and STEP are
//       identical to the corresponding arguments in SUBROUTINE LINCOB.
//     KOPT is such that XPT(KOPT,.) is the current trust region centre.
//     KNEW on exit is usually positive, and then it is the index of an
//       interpolation point to be moved to the position XPT(KOPT,.)+STEP(.).
//       It is set on entry either to its final value or to 0. In the latter
//       case, the final value of KNEW is chosen to maximize the denominator
//       of the matrix updating formula times a weighting factor.
//     VLAG and W are used for working space, the first NPT+N elements of
//       both of these vectors being required.
//
//     The arrays BMAT and ZMAT with IDZ are updated, the new matrices being
//       the ones that are suitable after the shift of the KNEW-th point to
//       the new position XPT(KOPT,.)+STEP(.). A return with KNEW set to zero
//       occurs if the calculation fails due to a zero denominator in the
//       updating formula, which should never happen.
//
//     Set some constants.
//
            double nptm = npt - n - 1;
//
//     Calculate VLAG and BETA for the current choice of STEP. The first NPT
//       elements of VLAG are set to the values of the Lagrange functions at
//       XPT(KOPT,.)+STEP(.). The first NPT components of W_check are held
//       in W, where W_check is defined in a paper on the updating method.
//
            for (int k = 1; k <= npt; ++k)
            {
                w[k] = sp[npt + k] * (HALF * sp[npt + k] + sp[k]);
                double sum = ZERO;
                for (int j = 1; j <= n; ++j)
                    sum += bmat[k, j] * step[j];
                vlag[k] = sum;
            }
            double beta = ZERO;
            for (int k = 1; k <= nptm; ++k)
            {
                double sum = ZERO;
                for (int i = 1; i <= npt; ++i)
                    sum = sum + zmat[i, k] * w[i];
                if (k < idz)
                {
                    beta += sum * sum;
                    sum = -sum;
                }
                else
                {
                    beta -= sum * sum;
                }
                for (int i = 1; i <= npt; ++i)
                    vlag[i] = vlag[i] + sum * zmat[i, k];
            }
            double bsum = ZERO;
            double dx = ZERO;
            double ssq = ZERO;
            for (int j = 1; j <= n; ++j)
            {
                double sum = ZERO;
                for (int i = 1; i <= npt; ++i)
                    sum += w[i] * bmat[i, j];
                bsum += sum * step[j];
                int jp = npt + j;
                for (int k = 1; k <= n; ++k)
                    sum += bmat[jp, k] * step[k];
                vlag[jp] = sum;
                bsum += sum * step[j];
                dx += step[j] * xpt[kopt, j];
                ssq += step[j] * step[j];
            }
            beta = dx * dx + ssq * (sp[kopt] + dx + dx + HALF * ssq) + beta - bsum;
            vlag[kopt] += +ONE;
//
//     If KNEW is zero initially, then pick the index of the interpolation
//       point to be deleted, by maximizing the absolute value of the
//       denominator of the updating formula times a weighting factor.
//       
//
            if (knew == 0)
            {
                double denmax = ZERO;
                for (int k = 1; k <= npt; ++k)
                {
                    double hdiag = ZERO;
                    for (int j = 1; j <= nptm; ++j)
                    {
                        double temp = ONE;
                        if (j < idz) temp = -ONE;
                        hdiag += temp * zmat[k, j] * zmat[k, j];
                    }
                    double denabs = Math.Abs(beta * hdiag + vlag[k] * vlag[k]);
                    double distsq = ZERO;
                    for (int j = 1; j <= n; ++j)
                        distsq += Math.Pow(xpt[k, j] - xpt[kopt, j], 2.0);
                    {
                        double temp = denabs * distsq * distsq;
                        if (temp > denmax)
                        {
                            denmax = temp;
                            knew = k;
                        }
                    }
                }
            }
//
//     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
//
            int jl = 1;
            double tempa, tempb = ZERO;
            if (nptm >= 2)
            {
                for (int j = 2; j <= nptm; ++j)
                {
                    if (j == idz)
                    {
                        jl = idz;
                    }
                    else if (zmat[knew, j] != ZERO)
                    {
                        double temp = Math.Sqrt(zmat[knew, jl] * zmat[knew, jl] + zmat[knew, j] * zmat[knew, j]);
                        tempa = zmat[knew, jl] / temp;
                        tempb = zmat[knew, j] / temp;
                        for (int i = 1; i <= npt; ++i)
                        {
                            temp = tempa * zmat[i, jl] + tempb * zmat[i, j];
                            zmat[i, j] = tempa * zmat[i, j] - tempb * zmat[i, jl];
                            zmat[i, jl] = temp;
                        }
                        zmat[knew, j] = ZERO;
                    }
                }
            }
//
//     Put the first NPT components of the KNEW-th column of the Z Z^T matrix
//       into W, and calculate the parameters of the updating formula.
//
            tempa = zmat[knew, 1];
            if (idz >= 2) tempa = -tempa;
            if (jl > 1) tempb = zmat[knew, jl];
            for (int i = 1; i <= npt; ++i)
            {
                w[i] = tempa * zmat[i, 1];
                if (jl > 1) w[i] += tempb * zmat[i, jl];
            }
            double alpha = w[knew];
            double tau = vlag[knew];
            double tausq = tau * tau;
            double denom = alpha * beta + tausq;
            vlag[knew] -= ONE;
            if (denom == ZERO)
            {
                knew = 0;
                return;
            }
            double sqrtdn = Math.Sqrt(Math.Abs(denom));
//
//     Complete the updating of ZMAT when there is only one nonzero element
//       in the KNEW-th row of the new matrix ZMAT. IFLAG is set to one when
//       the value of IDZ is going to be reduced.
//
            int iflag = 0;
            if (jl == 1)
            {
                tempa = tau / sqrtdn;
                tempb = zmat[knew, 1] / sqrtdn;
                for (int i = 1; i <= npt; ++i)
                    zmat[i, 1] = tempa * zmat[i, 1] - tempb * vlag[i];
                if (denom < ZERO)
                {
                    if (idz == 1)
                    {
                        idz = 2;
                    }
                    else
                    {
                        iflag = 1;
                    }
                }
            }
            else
            {
//
//     Complete the updating of ZMAT in the alternative case.
//
                int ja = 1;
                if (beta >= ZERO) ja = jl;
                int jb = jl + 1 - ja;
                double temp = zmat[knew, jb] / denom;
                tempa = temp * beta;
                tempb = temp * tau;
                temp = zmat[knew, ja];
                double scala = ONE / Math.Sqrt(Math.Abs(beta) * temp * temp + tausq);
                double scalb = scala * sqrtdn;
                for (int i = 1; i <= npt; ++i)
                {
                    zmat[i, ja] = scala * (tau * zmat[i, ja] - temp * vlag[i]);
                    zmat[i, jb] = scalb * (zmat[i, jb] - tempa * w[i] - tempb * vlag[i]);
                }
                if (denom <= ZERO)
                {
                    if (beta < ZERO)
                    {
                        idz = idz + 1;
                    }
                    else
                    {
                        iflag = 1;
                    }
                }
            }
//
//     Reduce IDZ when the diagonal part of the ZMAT times Diag(DZ) times
//       ZMAT^T factorization gains another positive element. Then exchange
//       the first and IDZ-th columns of ZMAT.
//
            if (iflag == 1)
            {
                --idz;
                for (int i = 1; i <= npt; ++i)
                {
                    double temp = zmat[i, 1];
                    zmat[i, 1] = zmat[i, idz];
                    zmat[i, idz] = temp;
                }
            }
//
//     Finally, update the matrix BMAT.
//
            for (int j = 1; j <= n; ++j)
            {
                int jp = npt + j;
                w[jp] = bmat[knew, j];
                tempa = (alpha * vlag[jp] - tau * w[jp]) / denom;
                tempb = (-beta * w[jp] - tau * vlag[jp]) / denom;
                for (int i = 1; i <= jp; ++i)
                {
                    bmat[i, j] = bmat[i, j] + tempa * vlag[i] + tempb * w[i];
                    if (i > npt) bmat[jp, i - npt] = bmat[i, j];
                }
            }
        }

        private static void PRINT(TextWriter logger, string format, params object[] args)
        {
            if (logger != null) logger.WriteLine(format, args);
        }

        private static string FORMAT<T>(string separator, string itemFormatter, IEnumerable<T> items, int start, int end)
        {
            return string.Join(
                separator,
                items.Skip(start).Take(end).Select(item => string.Format("{0," + itemFormatter + "}", item)).ToArray());
        }
        // ReSharper restore SuggestUseVarKeywordEvident
        // ReSharper restore InconsistentNaming

        #endregion
    }
}