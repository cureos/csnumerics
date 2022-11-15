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
 *  http://plato.asu.edu/ftp/other_software/bobyqa.zip
 */

using System;
using System.IO;

namespace Cureos.Numerics.Optimizers
{
    // ReSharper disable InconsistentNaming

    #region DELEGATES

    /// <summary>
    /// Delegate for the COBYLA objective function formulation.
    /// </summary>
    /// <param name="n">Number of variables.</param>
    /// <param name="x">Variable array.</param>
    /// <returns>Value of the objective function at <paramref name="x"/>.</returns>
    public delegate double BobyqaObjectiveFunctionDelegate(int n, double[] x);

    #endregion

    /// <summary>
    /// Representation of supported exit statuses from the Bobyqa algorithm.
    /// </summary>
    /// <summary>
    /// C# implementation of Powell�s nonlinear derivative�free bound constrained optimization that uses a quadratic 
    /// approximation approach. The algorithm applies a trust region method that forms quadratic models by interpolation. 
    /// There is usually some freedom in the interpolation conditions, which is taken up by minimizing the Frobenius norm 
    /// of the change to the second derivative of the model, beginning with the zero matrix. 
    /// The values of the variables are constrained by upper and lower bounds.
    /// </summary>
    public class Bobyqa : IQuadraticModelOptimizer
    {
        #region FIELDS

        private readonly int _n;
        private readonly BobyqaObjectiveFunctionDelegate _calfun;

        private readonly double[] _xl;
        private readonly double[] _xu;
        
        private int _npt;
        private int _maxfun;
        private int _iprint;

        private double _rhoend;
        private double _rhobeg;

        private TextWriter _logger;

        private const double INF = 1.0e60;
        private const double INFMIN = -1.0e60;

        private const double ONEMIN = -1.0;
        private const double ZERO = 0.0;
        private const double TENTH = 0.1;
        private const double HALF = 0.5;
        private const double ONE = 1.0;
        private const double TWO = 2.0;
        private const double TEN = 10.0;

        private static readonly string LF = Environment.NewLine;

        private static readonly string InvalidNptText = LF + "Return from BOBYQA because NPT is not in the required interval";
        private static readonly string TooSmallBoundRangeText = LF + "Return from BOBYQA because one of the differences XU[I]-XL[I] is less than 2*RHOBEG.";
        private static readonly string DenominatorCancellationText = LF + "Return from BOBYQA because of much cancellation in a denominator.";
        private static readonly string MaxIterationsText = LF + "Return from BOBYQA because CALFUN has been called MAXFUN times.";
        private static readonly string TrustRegionStepFailureText = LF + "Return from BOBYQA because a trust region step has failed to reduce Q.";

        private static readonly string IterationOutputFormat = LF + "Function number {0,6}    F ={1,18:E10}" + LF + "The corresponding X is: {2}";
        private static readonly string StageCompleteOutputFormat = LF + "Least value of F ={0,15:E9}" + LF + "The corresponding X is: {1}";
        private static readonly string RhoUpdatedFormat = LF + "New RHO ={0,11:E4}" + LF + "Number of function values ={1,6}";
        private static readonly string FinalNumberEvaluationsFormat = LF + "At the return from BOBYQA Number of function values = {0}";

        #endregion

        #region CONSTRUCTORS

        /// <summary>
        /// Constructor for finding a (local) minimum of the objective function <paramref name="calfun"/>, potentially subject to variable
        /// bounds <paramref name="xl"/> and <paramref name="xu"/>.
        /// </summary>
        /// <param name="calfun">Objective function subject to minimization.</param>
        /// <param name="n">Number of optimization variables, must be at least two.</param>
        /// <param name="xl">Lower bounds on the variables. Array is zero-based. If set to null, all variables
        /// are treated as downwards unbounded.</param>
        /// <param name="xu">Upper bounds on the variables. Array is zero-based. If set to null, all variables
        /// are treated as upwards unbounded.</param>
        /// <remarks>The construction of quadratic models requires <paramref name="xl"/> to be strictly less than 
        /// <paramref name="xu"/> for each index I. Further, the contribution to a model from changes to the I-th variable is
        /// damaged severely by rounding errors if difference between upper and lower bound is too small.</remarks>
        public Bobyqa(int n, BobyqaObjectiveFunctionDelegate calfun, double[] xl = null, double[] xu = null)
        {
            _n = n;
            _calfun = calfun;
            _xl = xl;
            _xu = xu;

            _npt = -1;
            _rhobeg = -1.0;
            _rhoend = -1.0;
            _iprint = 1;
            _maxfun = 10000;
            _logger = null;
        }

        #endregion

        #region PROPERTIES

        /// <summary>
        /// Gets or sets the number of maximum function calls.
        /// </summary>
        public int MaximumFunctionCalls
        {
            get
            {
                return this._maxfun;
            }
            set
            {
                this._maxfun = value;
            }
        }

        /// <summary>
        /// Gets or sets the print level to the logger.
        /// </summary>
        public int PrintLevel
        {
            get
            {
                return this._iprint;
            }
            set
            {
                this._iprint = value;
            }
        }

        /// <summary>
        /// Gets or sets the logger to which the optimizer log information should be sent.
        /// </summary>
        public TextWriter Logger
        {
            get
            {
                return this._logger;
            }
            set
            {
                this._logger = value;
            }
        }

        /// <summary>
        /// Gets or sets the number of interpolation conditions.
        /// </summary>
        public int InterpolationConditions
        {
            get
            {
                return this._npt;
            }
            set
            {
                this._npt = value;
            }
        }

        /// <summary>
        /// Gets or sets the final value of the trust region radius.
        /// </summary>
        public double TrustRegionRadiusEnd
        {
            get
            {
                return this._rhoend;
            }
            set
            {
                this._rhoend = value;
            }
        }

        /// <summary>
        /// Gets or sets the start value of the trust region radius.
        /// </summary>
        public double TrustRegionRadiusStart
        {
            get
            {
                return this._rhobeg;
            }
            set
            {
                this._rhobeg = value;
            }
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
            // Verify that the number of variables is greater than 1; BOBYQA does not support 1-D optimization.
            if (_n < 2) return new OptimizationSummary(OptimizationStatus.N_TooSmall, 0, null, Double.NaN);

            // Verify that the number of variables, and bounds if defined, in the respective array is sufficient.
            if (x0.Length < _n || (_xl != null && _xl.Length < _n) || (_xu != null && _xu.Length < _n))
            {
                return new OptimizationSummary(OptimizationStatus.VariableBoundsArrayTooShort, 0, null, Double.NaN);
            }

            // C# arrays are zero-based, whereas BOBYQA methods expect one-based arrays. Therefore define internal matrices
            // to be dispatched to the private BOBYQA methods.
            var xx0 = new double[1 + _n];
            Array.Copy(x0, 0, xx0, 1, _n);

            // If xl and/or xu are null, this is interpreted as that the optimization variables are all unbounded downwards and/or upwards.
            // In that case, assign artificial +/- infinity values to the bounds array(s).
            var xl = new double[1 + _n];
            if (_xl == null)
                for (var i = 1; i <= _n; ++i) xl[i] = INFMIN;
            else
                Array.Copy(_xl, 0, xl, 1, _n);

            var xu = new double[1 + _n];
            if (_xu == null)
                for (var i = 1; i <= _n; ++i) xu[i] = INF;
            else
                Array.Copy(_xu, 0, xu, 1, _n);

            // Verify that all lower bounds are less than upper bounds.
            // If any start value is outside bounds, adjust this value to be within bounds.
            var rng = new double[1 + _n];
            var minrng = Double.MaxValue;
            var maxabsx = 0.0;

            for (var i = 1; i <= _n; ++i)
            {
                if ((rng[i] = xu[i] - xl[i]) <= 0.0)
                    return new OptimizationSummary(OptimizationStatus.InvalidBoundsSpecification, 0, null, Double.NaN);
                minrng = Math.Min(rng[i], minrng);

                if (xx0[i] < xl[i]) xx0[i] = xl[i];
                if (xx0[i] > xu[i]) xx0[i] = xu[i];
                maxabsx = Math.Max(Math.Abs(xx0[i]), maxabsx);
            }

            // If rhobeg is non-positive, set rhobeg based on the absolute values of the variables' start values,
            // using same strategy as R-project BOBYQA wrapper.
            if (_rhobeg <= 0.0) _rhobeg = maxabsx > 0.0 ? Math.Min(0.95, 0.2 * maxabsx) : 0.95;

            // Required that rhobeg is less than half the minimum bounds range; adjust rhobeg if necessary.
            if (_rhobeg > 0.5 * minrng) _rhobeg = 0.2 * minrng;

            // If rhoend is non-negative, set rhoend to one millionth of the rhobeg value (R-project strategy).
            if (_rhoend <= 0.0) _rhoend = 1.0e-6 * _rhobeg;

            // If npt is non-positive, apply default value 2 * n + 1.
            var inpt = _npt > 0 ? _npt : 2 * _n + 1;

            // Define internal calfun to account for that the x vector in the function invocation is one-based.
            var icalfun = new Func<int, double[], double>((n, x) =>
            {
                var xx = new double[_n];
                Array.Copy(x, 1, xx, 0, _n);
                return _calfun(n, xx);
            });

            // Invoke optimization. After completed optimization, transfer the optimized internal variable array to the
            // variable array in the method call.
            return BOBYQA(icalfun, _n, inpt, xx0, xl, xu, _rhobeg, _rhoend, _iprint, _maxfun, _logger);
        }

        #endregion

        #region PRIVATE BOBYQA ALGORITHM METHODS

        private static OptimizationSummary BOBYQA(Func<int, double[], double> calfun, int n, int npt, double[] x,
            double[] xl, double[] xu, double rhobeg, double rhoend, int iprint, int maxfun, TextWriter logger)
        {
            //     This subroutine seeks the least value of a function of many variables,
            //     by applying a trust region method that forms quadratic models by
            //     interpolation. There is usually some freedom in the interpolation
            //     conditions, which is taken up by minimizing the Frobenius norm of
            //     the change to the second derivative of the model, beginning with the
            //     zero matrix. The values of the variables are constrained by upper and
            //     lower bounds. The arguments of the subroutine are as follows.
            //
            //     N must be set to the number of variables and must be at least two.
            //     NPT is the number of interpolation conditions. Its value must be in
            //       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
            //       recommended.
            //     Initial values of the variables must be set in X(1),X(2),...,X(N). They
            //       will be changed to the values that give the least calculated F.
            //     For I=1,2,...,N, XL[I] and XU[I] must provide the lower and upper
            //       bounds, respectively, on X[I]. The construction of quadratic models
            //       requires XL[I] to be strictly less than XU[I] for each I. Further,
            //       the contribution to a model from changes to the I-th variable is
            //       damaged severely by rounding errors if XU[I]-XL[I] is too small.
            //     RHOBEG and RHOEND must be set to the initial and final values of a trust
            //       region radius, so both must be positive with RHOEND no greater than
            //       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
            //       expected change to a variable, while RHOEND should indicate the
            //       accuracy that is required in the final values of the variables. An
            //       error return occurs if any of the differences XU[I]-XL[I], I=1,...,N,
            //       is less than 2*RHOBEG.
            //     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
            //       amount of printing. Specifically, there is no output if IPRINT=0 and
            //       there is output only at the return if IPRINT=1. Otherwise, each new
            //       value of RHO is printed, with the best vector of variables so far and
            //       the corresponding value of the objective function. Further, each new
            //       value of F with its variables are output if IPRINT=3.
            //     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
            //
            //     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
            //     F to the value of the objective function for the current values of the
            //     variables X(1),X(2),...,X(N), which are generated automatically in a
            //     way that satisfies the bounds given in XL and XU.

            //     Return if the value of NPT is unacceptable.

            var np = n + 1;
            if (npt < n + 2 || npt > ((n + 2) * np) / 2)
            {
                if (logger != null) logger.WriteLine(InvalidNptText);
                return new OptimizationSummary(OptimizationStatus.NPT_OutOfRange, 0, null, Double.NaN);
            }

            var ndim = npt + n;

            var sl = new double[1 + n];
            var su = new double[1 + n];

            //     Return if there is insufficient space between the bounds. Modify the
            //     initial X if necessary in order to avoid conflicts between the bounds
            //     and the construction of the first quadratic model. The lower and upper
            //     bounds on moves from the updated X are set now, in the ISL and ISU
            //     partitions of W, in order to provide useful and exact information about
            //     components of X that become within distance RHOBEG from their bounds.

            for (var j = 1; j <= n; ++j)
            {
                double temp = xu[j] - xl[j];
                if (temp < rhobeg + rhobeg)
                {
                    if (logger != null) logger.WriteLine(TooSmallBoundRangeText);
                    return new OptimizationSummary(OptimizationStatus.BoundsRangeTooSmall, 0, null, Double.NaN);
                }
                sl[j] = xl[j] - x[j];
                su[j] = xu[j] - x[j];
                if (sl[j] >= -rhobeg)
                {
                    if (sl[j] >= ZERO)
                    {
                        x[j] = xl[j];
                        sl[j] = ZERO;
                        su[j] = temp;
                    }
                    else
                    {
                        x[j] = xl[j] + rhobeg;
                        sl[j] = -rhobeg;
                        su[j] = Math.Max(xu[j] - x[j], rhobeg);
                    }
                }
                else if (su[j] <= rhobeg)
                {
                    if (su[j] <= ZERO)
                    {
                        x[j] = xu[j];
                        sl[j] = -temp;
                        su[j] = ZERO;
                    }
                    else
                    {
                        x[j] = xu[j] - rhobeg;
                        sl[j] = Math.Min(xl[j] - x[j], -rhobeg);
                        su[j] = rhobeg;
                    }
                }
            }

            //     Make the call of BOBYQB.
            return BOBYQB(calfun, n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, ndim, sl, su, logger);
        }

        private static OptimizationSummary BOBYQB(Func<int, double[], double> calfun, int n, int npt, double[] x, double[] xl,
            double[] xu, double rhobeg, double rhoend, int iprint, int maxfun, int ndim, double[] sl, double[] su,
            TextWriter logger)
        {
            //     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
            //       are identical to the corresponding arguments in SUBROUTINE BOBYQA.
            //     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
            //       All the components of every XOPT are going to satisfy the bounds
            //       SL[I] .LEQ. XOPT[I] .LEQ. SU[I], with appropriate equalities when
            //       XOPT is on a constraint boundary.

            //     Set some constants.

            var np = n + 1;
            var nptm = npt - np;
            var nh = (n * np) / 2;

            //     XBASE holds a shift of origin that should reduce the contributions
            //       from rounding errors to values of the model and Lagrange functions.
            //     XPT is a two-dimensional array that holds the coordinates of the
            //       interpolation points relative to XBASE.
            //     FVAL holds the values of F at the interpolation points.
            //     XOPT is set to the displacement from XBASE of the trust region centre.
            //     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
            //     HQ holds the explicit second derivatives of the quadratic model.
            //     PQ contains the parameters of the implicit second derivatives of the
            //       quadratic model.
            //     BMAT holds the last N columns of H.
            //     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
            //       this factorization being ZMAT times ZMAT^T, which provides both the
            //       correct rank and positive semi-definiteness.
            //     NDIM is the first dimension of BMAT and has the value NPT+N.
            //     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
            //       vector of variables for the next call of CALFUN. XNEW also satisfies
            //       the SL and SU constraints in the way that has just been mentioned.
            //     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
            //       in order to increase the denominator in the updating of UPDATE.
            //     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
            //     VLAG contains the values of the Lagrange functions at a new point X.
            //       They are part of a product that requires VLAG to be of length NDIM.
            //     W is a one-dimensional array that is used for working space. Its length
            //       must be at least 3*NDIM = 3*(NPT+N).

            var xbase = new double[1 + n];
            var xpt = new double[1 + npt,1 + n];
            var fval = new double[1 + npt];
            var xopt = new double[1 + n];
            var gopt = new double[1 + n];
            var hq = new double[1 + n * np / 2];
            var pq = new double[1 + npt];
            var bmat = new double[1 + ndim,1 + n];
            var zmat = new double[1 + npt,1 + npt - np];
            var xnew = new double[1 + n];
            var xalt = new double[1 + n];
            var d = new double[1 + n];
            var vlag = new double[1 + ndim];

            var wn = new double[1 + n];
            var w2npt = new double[1 + 2 * npt];

            var knew = 0;
            var adelt = 0.0;
            var alpha = 0.0;
            var beta = 0.0;
            var cauchy = 0.0;
            var denom = 0.0;
            var diffc = 0.0;
            var ratio = 0.0;
            var f = 0.0;

            double distsq;
            OptimizationStatus status;

            //     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
            //     BMAT and ZMAT for the first iteration, with the corresponding values of
            //     of NF and KOPT, which are the number of calls of CALFUN so far and the
            //     index of the interpolation point at the trust region centre. Then the
            //     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
            //     less than NPT. GOPT will be updated if KOPT is different from KBASE.

            int nf, kopt;
            PRELIM(calfun,n, npt, x, xl, xu, rhobeg, iprint, maxfun, xbase, xpt, fval, gopt, hq, pq, bmat, zmat, ndim, sl, su,
                   out nf, out kopt, logger);

            var xoptsq = ZERO;
            for (var i = 1; i <= n; ++i)
            {
                xopt[i] = xpt[kopt, i];
                xoptsq += xopt[i] * xopt[i];
            }
            var fsave = fval[1];
            if (nf < npt)
            {
                if (iprint > 0 && logger != null) logger.WriteLine(MaxIterationsText);
                status = OptimizationStatus.MAXFUN_Reached;
                goto L_720;
            }
            var kbase = 1;

            //     Complete the settings that are required for the iterative procedure.

            var rho = rhobeg;
            var delta = rho;
            var nresc = nf;
            var ntrits = 0;
            var diffa = ZERO;
            var diffb = ZERO;
            var itest = 0;
            var nfsav = nf;

            //     Update GOPT if necessary before the first iteration and after each
            //     call of RESCUE that makes a call of CALFUN.

            L_20:
            if (kopt != kbase)
            {
                var ih = 0;
                for (var j = 1; j <= n; ++j)
                {
                    for (var i = 1; i <= j; ++i)
                    {
                        ih = ih + 1;
                        if (i < j) gopt[j] += hq[ih] * xopt[i];
                        gopt[i] += hq[ih] * xopt[j];
                    }
                }
                if (nf > npt)
                {
                    for (var k = 1; k <= npt; ++k)
                    {
                        var temp = ZERO;
                        for (var j = 1; j <= n; ++j) temp += xpt[k, j] * xopt[j];
                        temp *= pq[k];
                        for (var i = 1; i <= n; ++i) gopt[i] += temp * xpt[k, i];
                    }
                }
            }

            //     Generate the next point in the trust region that provides a small value
            //     of the quadratic model subject to the constraints on the variables.
            //     The integer NTRITS is set to the number "trust region" iterations that
            //     have occurred since the last "alternative" iteration. If the length
            //     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
            //     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.

            L_60:
            var gnew = new double[1 + n];
            double dsq, crvmin;
            TRSBOX(n, npt, xpt, xopt, gopt, hq, pq, sl, su, delta, xnew, d, gnew, out dsq, out crvmin);

            var dnorm = Math.Min(delta, Math.Sqrt(dsq));
            if (dnorm < HALF * rho)
            {
                ntrits = -1;
                distsq = TEN * TEN * rho * rho;
                if (nf <= nfsav + 2) goto L_650;

                //     The following choice between labels 650 and 680 depends on whether or
                //     not our work with the current RHO seems to be complete. Either RHO is
                //     decreased or termination occurs if the errors in the quadratic model at
                //     the last three interpolation points compare favourably with predictions
                //     of likely improvements to the model within distance HALF*RHO of XOPT.

                var errbig = Math.Max(diffa, Math.Max(diffb, diffc));
                var frhosq = 0.125 * rho * rho;
                if (crvmin > ZERO && errbig > frhosq * crvmin) goto L_650;
                var bdtol = errbig / rho;
                for (var j = 1; j <= n; ++j)
                {
                    var bdtest = bdtol;
                    if (xnew[j] == sl[j]) bdtest = gnew[j];
                    if (xnew[j] == su[j]) bdtest = -gnew[j];
                    if (bdtest < bdtol)
                    {
                        var curv = hq[(j + j * j) / 2];
                        for (var k = 1; k <= npt; ++k)
                        {
                            curv = curv + pq[k] * xpt[k, j] * xpt[k, j];
                        }
                        bdtest = bdtest + HALF * curv * rho;
                        if (bdtest < bdtol) goto L_650;
                    }
                }
                goto L_680;
            }
            ++ntrits;

            //     Severe cancellation is likely to occur if XOPT is too far from XBASE.
            //     If the following test holds, then XBASE is shifted so that XOPT becomes
            //     zero. The appropriate changes are made to BMAT and to the second
            //     derivatives of the current model, beginning with the changes to BMAT
            //     that do not depend on ZMAT. VLAG is used temporarily for working space.

            L_90:
            if (dsq <= 1.0E-3 * xoptsq)
            {
                var fracsq = 0.25 * xoptsq;
                var sumpq = ZERO;
                for (var k = 1; k <= npt; ++k)
                {
                    sumpq += pq[k];
                    var sum = -HALF * xoptsq;
                    for (var i = 1; i <= n; ++i) sum += xpt[k, i] * xopt[i];
                    w2npt[k] = sum;
                    var temp = fracsq - HALF * sum;
                    for (var i = 1; i <= n; ++i)
                    {
                        wn[i] = bmat[k, i];
                        vlag[i] = sum * xpt[k, i] + temp * xopt[i];
                        var ip = npt + i;
                        for (var j = 1; j <= i; ++j) bmat[ip, j] += wn[i] * vlag[j] + vlag[i] * wn[j];
                    }
                }

                //     Then the revisions of BMAT that depend on ZMAT are calculated.

                for (var jj = 1; jj <= nptm; ++jj)
                {
                    var sumz = ZERO;
                    var sumw = ZERO;
                    for (var k = 1; k <= npt; ++k)
                    {
                        sumz += zmat[k, jj];
                        vlag[k] = w2npt[k] * zmat[k, jj];
                        sumw += vlag[k];
                    }
                    for (var j = 1; j <= n; ++j)
                    {
                        var sum = (fracsq * sumz - HALF * sumw) * xopt[j];
                        for (var k = 1; k <= npt; ++k) sum += vlag[k] * xpt[k, j];
                        wn[j] = sum;
                        for (var k = 1; k <= npt; ++k) bmat[k, j] += sum * zmat[k, jj];
                    }
                    for (var i = 1; i <= n; ++i)
                    {
                        var ip = i + npt;
                        var temp = wn[i];
                        for (var j = 1; j <= i; ++j) bmat[ip, j] += temp * wn[j];
                    }
                }

                //     The following instructions complete the shift, including the changes
                //     to the second derivative parameters of the quadratic model.

                var ih = 0;
                for (var j = 1; j <= n; ++j)
                {
                    wn[j] = -HALF * sumpq * xopt[j];
                    for (var k = 1; k <= npt; ++k)
                    {
                        wn[j] += pq[k] * xpt[k, j];
                        xpt[k, j] -= xopt[j];
                    }
                    for (var i = 1; i <= j; ++i)
                    {
                        hq[++ih] += wn[i] * xopt[j] + xopt[i] * wn[j];
                        bmat[npt + i, j] = bmat[npt + j, i];
                    }
                }
                for (var i = 1; i <= n; ++i)
                {
                    xbase[i] += xopt[i];
                    xnew[i] -= xopt[i];
                    sl[i] -= xopt[i];
                    su[i] -= xopt[i];
                    xopt[i] = ZERO;
                }
                xoptsq = ZERO;
            }

            if (ntrits == 0) goto L_210;
            goto L_230;

            //     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
            //     more expensive than the previous shift, because new matrices BMAT and
            //     ZMAT are generated from scratch, which may include the replacement of
            //     interpolation points whose positions seem to be causing near linear
            //     dependence in the interpolation conditions. Therefore RESCUE is called
            //     only if rounding errors have reduced by at least a factor of two the
            //     denominator of the formula for updating the H matrix. It provides a
            //     useful safeguard, but is not invoked in most applications of BOBYQA.

            L_190:
            nfsav = nf;
            kbase = kopt;
            RESCUE(calfun, n, npt, xl, xu, iprint, maxfun, xbase, xpt, fval, xopt, gopt, hq, pq, bmat, zmat, ndim, sl,
                   su, ref nf, delta, ref kopt, vlag, logger);

            //     XOPT is updated now in case the branch below to label 720 is taken.
            //     Any updating of GOPT occurs after the branch below to label 20, which
            //     leads to a trust region iteration as does the branch to label 60.

            xoptsq = ZERO;
            if (kopt != kbase)
            {
                for (var i = 1; i <= n; ++i)
                {
                    xopt[i] = xpt[kopt, i];
                    xoptsq = xoptsq + xopt[i] * xopt[i];
                }
            }
            if (nf < 0)
            {
                nf = maxfun;
                if (iprint > 0 && logger != null) logger.WriteLine(MaxIterationsText);
                status = OptimizationStatus.MAXFUN_Reached;
                goto L_720;
            }
            nresc = nf;
            if (nfsav < nf)
            {
                nfsav = nf;
                goto L_20;
            }
            if (ntrits > 0) goto L_60;

            //     Pick two alternative vectors of variables, relative to XBASE, that
            //     are suitable as new positions of the KNEW-th interpolation point.
            //     Firstly, XNEW is set to the point on a line through XOPT and another
            //     interpolation point that minimizes the predicted value of the next
            //     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
            //     and SU bounds. Secondly, XALT is set to the best feasible point on
            //     a constrained version of the Cauchy step of the KNEW-th Lagrange
            //     function, the corresponding value of the square of this function
            //     being returned in CAUCHY. The choice between these alternatives is
            //     going to be made when the denominator is calculated.

            L_210:
            ALTMOV(n, npt, xpt, xopt, bmat, zmat, sl, su, kopt, knew, adelt, xnew, xalt, out alpha, out cauchy);
            for (var i = 1; i <= n; ++i) d[i] = xnew[i] - xopt[i];

            //     Calculate VLAG and BETA for the current choice of D. The scalar
            //     product of D with XPT(K,.) is going to be held in W(NPT+K) for
            //     use when VQUAD is calculated.

            L_230:
            for (var k = 1; k <= npt; ++k)
            {
                var suma = ZERO;
                var sumb = ZERO;
                var sum = ZERO;
                for (var j = 1; j <= n; ++j)
                {
                    suma += xpt[k, j] * d[j];
                    sumb += xpt[k, j] * xopt[j];
                    sum += bmat[k, j] * d[j];
                }
                w2npt[k] = suma * (HALF * suma + sumb);
                vlag[k] = sum;
                w2npt[npt + k] = suma;
            }

            beta = ZERO;
            for (var jj = 1; jj <= nptm; ++jj)
            {
                var sum = ZERO;
                for (var k = 1; k <= npt; ++k) sum += zmat[k, jj] * w2npt[k];
                beta -= sum * sum;
                for (var k = 1; k <= npt; ++k) vlag[k] += sum * zmat[k, jj];
            }
            dsq = ZERO;
            var bsum = ZERO;
            var dx = ZERO;
            for (var j = 1; j <= n; ++j)
            {
                dsq = dsq + d[j] * d[j];
                var sum = ZERO;
                for (var k = 1; k <= npt; ++k) sum += w2npt[k] * bmat[k, j];
                bsum += sum * d[j];
                var jp = npt + j;
                for (var i = 1; i <= n; ++i) sum += bmat[jp, i] * d[i];
                vlag[jp] = sum;
                bsum += sum * d[j];
                dx += d[j] * xopt[j];
            }
            beta += dx * dx + dsq * (xoptsq + dx + dx + HALF * dsq) - bsum;
            vlag[kopt] += ONE;

            //     If NTRITS is zero, the denominator may be increased by replacing
            //     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
            //     rounding errors have damaged the chosen denominator.

            if (ntrits == 0)
            {
                denom = vlag[knew] * vlag[knew] + alpha * beta;
                if (denom < cauchy && cauchy > ZERO)
                {
                    for (var i = 1; i <= n; ++i)
                    {
                        xnew[i] = xalt[i];
                        d[i] = xnew[i] - xopt[i];
                    }
                    cauchy = ZERO;
                    goto L_230;
                }
                if (denom <= HALF * vlag[knew] * vlag[knew])
                {
                    if (nf > nresc) goto L_190;
                    if (iprint > 0 && logger != null) logger.WriteLine(DenominatorCancellationText);
                    status = OptimizationStatus.DenominatorCancellation;
                    goto L_720;
                }
            }

            //     Alternatively, if NTRITS is positive, then set KNEW to the index of
            //     the next interpolation point to be deleted to make room for a trust
            //     region step. Again RESCUE may be called if rounding errors have damaged
            //     the chosen denominator, which is the reason for attempting to select
            //     KNEW before calculating the next value of the objective function.

            else
            {
                var delsq = delta * delta;
                var scaden = ZERO;
                var biglsq = ZERO;
                knew = 0;
                for (var k = 1; k <= npt; ++k)
                {
                    if (k == kopt) continue;
                    var hdiag = ZERO;
                    for (var jj = 1; jj <= nptm; ++jj) hdiag += zmat[k, jj] * zmat[k, jj];
                    var den = beta * hdiag + vlag[k] * vlag[k];
                    distsq = ZERO;
                    for (var j = 1; j <= n; ++j) distsq += Math.Pow(xpt[k, j] - xopt[j], 2.0);
                    var temp = Math.Max(ONE, Math.Pow(distsq / delsq, 2.0));
                    if (temp * den > scaden)
                    {
                        scaden = temp * den;
                        knew = k;
                        denom = den;
                    }
                    biglsq = Math.Max(biglsq, temp * vlag[k] * vlag[k]);
                }
                if (scaden <= HALF * biglsq)
                {
                    if (nf > nresc) goto L_190;
                    if (iprint > 0 && logger != null) logger.WriteLine(DenominatorCancellationText);
                    status = OptimizationStatus.DenominatorCancellation;
                    goto L_720;
                }
            }

            //     Put the variables for the next calculation of the objective function
            //       in XNEW, with any adjustments for the bounds.

            //     Calculate the value of the objective function at XBASE+XNEW, unless
            //       the limit on the number of calculations of F has been reached.

            L_360:
            for (var i = 1; i <= n; ++i)
            {
                x[i] = Math.Min(Math.Max(xl[i], xbase[i] + xnew[i]), xu[i]);
                if (xnew[i] == sl[i]) x[i] = xl[i];
                if (xnew[i] == su[i]) x[i] = xu[i];
            }
            if (nf >= maxfun)
            {
                if (iprint > 0 && logger != null) logger.WriteLine(MaxIterationsText);
                status = OptimizationStatus.MAXFUN_Reached;
                goto L_720;
            }

            ++nf;
            f = calfun(n, x);

            if (iprint == 3 && logger != null) logger.WriteLine(IterationOutputFormat, nf, f, ToString(x, n));
            if (ntrits == -1)
            {
                fsave = f;
                status = OptimizationStatus.Normal;
                goto L_720;
            }

            //     Use the quadratic model to predict the change in F due to the step D,
            //       and set DIFF to the error of this prediction.

            var fopt = fval[kopt];
            var vquad = ZERO;
            {
                var ih = 0;
                for (var j = 1; j <= n; ++j)
                {
                    vquad += d[j] * gopt[j];
                    for (var i = 1; i <= j; ++i)
                        vquad += hq[++ih] * (i == j ? HALF : ONE) * d[i] * d[j];
                }
            }
            for (var k = 1; k <= npt; ++k) vquad += HALF * pq[k] * w2npt[npt + k] * w2npt[npt + k];

            var diff = f - fopt - vquad;
            diffc = diffb;
            diffb = diffa;
            diffa = Math.Abs(diff);
            if (dnorm > rho) nfsav = nf;

            //     Pick the next value of DELTA after a trust region step.

            if (ntrits > 0)
            {
                if (vquad >= ZERO)
                {
                    if (iprint > 0 && logger != null) logger.WriteLine(TrustRegionStepFailureText);
                    status = OptimizationStatus.TrustRegionStepReductionFailure;
                    goto L_720;
                }
                ratio = (f - fopt) / vquad;

                if (ratio <= TENTH)
                    delta = Math.Min(HALF * delta, dnorm);
                else if (ratio <= 0.7)
                    delta = Math.Max(HALF * delta, dnorm);
                else
                    delta = Math.Max(HALF * delta, dnorm + dnorm);

                if (delta <= 1.5 * rho) delta = rho;

                //     Recalculate KNEW and DENOM if the new F is less than FOPT.

                if (f < fopt)
                {
                    var ksav = knew;
                    var densav = denom;
                    var delsq = delta * delta;
                    var scaden = ZERO;
                    var biglsq = ZERO;
                    knew = 0;
                    for (var k = 1; k <= npt; ++k)
                    {
                        var hdiag = ZERO;
                        for (var jj = 1; jj <= nptm; ++jj) hdiag += zmat[k, jj] * zmat[k, jj];
                        var den = beta * hdiag + vlag[k] * vlag[k];
                        distsq = ZERO;
                        for (var j = 1; j <= n; ++j) distsq += Math.Pow(xpt[k, j] - xnew[j], 2.0);
                        var temp = Math.Max(ONE, Math.Pow(distsq / delsq, 2.0));
                        if (temp * den > scaden)
                        {
                            scaden = temp * den;
                            knew = k;
                            denom = den;
                        }
                        biglsq = Math.Max(biglsq, temp * vlag[k] * vlag[k]);
                    }
                    if (scaden <= HALF * biglsq)
                    {
                        knew = ksav;
                        denom = densav;
                    }
                }
            }

            //     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
            //     moved. Also update the second derivative terms of the model.

            var w = new double[1 + ndim];
            UPDATE(n, npt, bmat, zmat, ndim, vlag, beta, denom, knew, w);

            var pqold = pq[knew];
            pq[knew] = ZERO;
            {
                var ih = 0;
                for (var i = 1; i <= n; ++i)
                {
                    var temp = pqold * xpt[knew, i];
                    for (var j = 1; j <= i; ++j) hq[++ih] += temp * xpt[knew, j];
                }
            }
            for (var jj = 1; jj <= nptm; ++jj)
            {
                var temp = diff * zmat[knew, jj];
                for (var k = 1; k <= npt; ++k) pq[k] += temp * zmat[k, jj];
            }

            //     Include the new interpolation point, and make the changes to GOPT at
            //     the old XOPT that are caused by the updating of the quadratic model.

            fval[knew] = f;
            for (var i = 1; i <= n; ++i)
            {
                xpt[knew, i] = xnew[i];
                wn[i] = bmat[knew, i];
            }
            for (var k = 1; k <= npt; ++k)
            {
                var suma = ZERO;
                for (var jj = 1; jj <= nptm; ++jj) suma += zmat[knew, jj] * zmat[k, jj];
                var sumb = ZERO;
                for (var j = 1; j <= n; ++j) sumb += xpt[k, j] * xopt[j];
                var temp = suma * sumb;
                for (var i = 1; i <= n; ++i) wn[i] += temp * xpt[k, i];
            }
            for (var i = 1; i <= n; ++i) gopt[i] += diff * wn[i];

            //     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.

            if (f < fopt)
            {
                kopt = knew;
                xoptsq = ZERO;
                var ih = 0;
                for (var j = 1; j <= n; ++j)
                {
                    xopt[j] = xnew[j];
                    xoptsq += xopt[j] * xopt[j];
                    for (var i = 1; i <= j; ++i)
                    {
                        ++ih;
                        if (i < j) gopt[j] += +hq[ih] * d[i];
                        gopt[i] += hq[ih] * d[j];
                    }
                }
                for (var k = 1; k <= npt; ++k)
                {
                    var temp = ZERO;
                    for (var j = 1; j <= n; ++j) temp += xpt[k, j] * d[j];
                    temp *= pq[k];
                    for (var i = 1; i <= n; ++i) gopt[i] += temp * xpt[k, i];
                }
            }

            //     Calculate the parameters of the least Frobenius norm interpolant to
            //     the current data, the gradient of this interpolant at XOPT being put
            //     into VLAG(NPT+I), I=1,2,...,N.

            if (ntrits > 0)
            {
                for (var k = 1; k <= npt; ++k)
                {
                    vlag[k] = fval[k] - fval[kopt];
                    w2npt[k] = ZERO;
                }
                for (var j = 1; j <= nptm; ++j)
                {
                    var sum = ZERO;
                    for (var k = 1; k <= npt; ++k) sum += zmat[k, j] * vlag[k];
                    for (var k = 1; k <= npt; ++k) w2npt[k] = w2npt[k] + sum * zmat[k, j];
                }
                for (var k = 1; k <= npt; ++k)
                {
                    var sum = ZERO;
                    for (var j = 1; j <= n; ++j) sum += xpt[k, j] * xopt[j];
                    w2npt[k + npt] = w2npt[k];
                    w2npt[k] *= sum;
                }
                var gqsq = ZERO;
                var gisq = ZERO;
                for (var i = 1; i <= n; ++i)
                {
                    var sum = ZERO;
                    for (var k = 1; k <= npt; ++k) sum += bmat[k, i] * vlag[k] + xpt[k, i] * w2npt[k];
                    if (xopt[i] == sl[i])
                    {
                        gqsq += Math.Pow(Math.Min(ZERO, gopt[i]), 2.0);
                        gisq += Math.Pow(Math.Min(ZERO, sum), 2.0);
                    }
                    else if (xopt[i] == su[i])
                    {
                        gqsq += Math.Pow(Math.Max(ZERO, gopt[i]), 2.0);
                        gisq += Math.Pow(Math.Max(ZERO, sum), 2.0);
                    }
                    else
                    {
                        gqsq += gopt[i] * gopt[i];
                        gisq += sum * sum;
                    }
                    vlag[npt + i] = sum;
                }

                //     Test whether to replace the new quadratic model by the least Frobenius
                //     norm interpolant, making the replacement if the test is satisfied.

                ++itest;
                if (gqsq < TEN * gisq) itest = 0;
                if (itest >= 3)
                {
                    for (var i = 1; i <= Math.Max(npt, nh); ++i)
                    {
                        if (i <= n) gopt[i] = vlag[npt + i];
                        if (i <= npt) pq[i] = w2npt[npt + i];
                        if (i <= nh) hq[i] = ZERO;
                        itest = 0;
                    }
                }
            }

            //     If a trust region step has provided a sufficient decrease in F, then
            //     branch for another trust region calculation. The case NTRITS=0 occurs
            //     when the new interpolation point was reached by an alternative step.

            if (ntrits == 0 || f <= fopt + TENTH * vquad) goto L_60;

            //     Alternatively, find out if the interpolation points are close enough
            //       to the best point so far.

            distsq = Math.Max(TWO * TWO * delta * delta, TEN * TEN * rho * rho);

            L_650:
            knew = 0;
            for (var k = 1; k <= npt; ++k)
            {
                var sum = ZERO;
                for (var j = 1; j <= n; ++j) sum += Math.Pow(xpt[k, j] - xopt[j], 2.0);
                if (sum > distsq)
                {
                    knew = k;
                    distsq = sum;
                }
            }

            //     If KNEW is positive, then ALTMOV finds alternative new positions for
            //     the KNEW-th interpolation point within distance ADELT of XOPT. It is
            //     reached via label 90. Otherwise, there is a branch to label 60 for
            //     another trust region iteration, unless the calculations with the
            //     current RHO are complete.

            if (knew > 0)
            {
                var dist = Math.Sqrt(distsq);
                if (ntrits == -1)
                {
                    delta = Math.Min(TENTH * delta, HALF * dist);
                    if (delta <= 1.5 * rho) delta = rho;
                }
                ntrits = 0;
                adelt = Math.Max(Math.Min(TENTH * dist, delta), rho);
                dsq = adelt * adelt;
                goto L_90;
            }
            if (ntrits == -1) goto L_680;
            if (ratio > ZERO || Math.Max(delta, dnorm) > rho) goto L_60;

            //     The calculations with the current value of RHO are complete. Pick the
            //       next values of RHO and DELTA.

            L_680:
            if (rho > rhoend)
            {
                delta = HALF * rho;
                ratio = rho / rhoend;

                if (ratio <= 16.0)
                    rho = rhoend;
                else if (ratio <= 250.0)
                    rho = Math.Sqrt(ratio) * rhoend;
                else
                    rho = TENTH * rho;

                delta = Math.Max(delta, rho);
                if (iprint >= 2)
                {
                    var bestX = new double[1 + n];
                    for (var i = 1; i <= n; ++i) bestX[i] = xbase[i] + xopt[i];

                    if (logger != null)
                    {
                        logger.WriteLine(RhoUpdatedFormat, rho, nf);
                        logger.WriteLine(StageCompleteOutputFormat, fval[kopt], ToString(bestX, n));
                    }
                }
                ntrits = 0;
                nfsav = nf;
                goto L_60;
            }

            //     Return from the calculation, after another Newton-Raphson step, if
            //       it is too short to have been tried before.

            if (ntrits == -1) goto L_360;
            status = OptimizationStatus.Normal;

            L_720:
            if (fval[kopt] <= fsave)
            {
                for (var i = 1; i <= n; ++i)
                {
                    x[i] = Math.Min(Math.Max(xl[i], xbase[i] + xopt[i]), xu[i]);
                    if (xopt[i] == sl[i]) x[i] = xl[i];
                    if (xopt[i] == su[i]) x[i] = xu[i];
                }
                f = fval[kopt];
            }
            if (iprint >= 1 && logger != null)
            {
                logger.WriteLine(FinalNumberEvaluationsFormat, nf);
                logger.WriteLine(StageCompleteOutputFormat, f, ToString(x, n));
            }

            var xret = new double[n];
            Array.Copy(x, 1, xret, 0, n);

            return new OptimizationSummary(status, nf, xret, f);
        }

        private static void ALTMOV(int n, int npt, double[,] xpt, double[] xopt, double[,] bmat,
            double[,] zmat, double[] sl, double[] su, int kopt, int knew, double adelt,
            double[] xnew, double[] xalt, out double alpha, out double cauchy)
        {
            //     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
            //       the same meanings as the corresponding arguments of BOBYQB.
            //     KOPT is the index of the optimal interpolation point.
            //     KNEW is the index of the interpolation point that is going to be moved.
            //     ADELT is the current trust region bound.
            //     XNEW will be set to a suitable new position for the interpolation point
            //       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
            //       bounds and it should provide a large denominator in the next call of
            //       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
            //       straight lines through XOPT and another interpolation point.
            //     XALT also provides a large value of the modulus of the KNEW-th Lagrange
            //       function subject to the constraints that have been mentioned, its main
            //       difference from XNEW being that XALT-XOPT is a constrained version of
            //       the Cauchy step within the trust region. An exception is that XALT is
            //       not calculated if all components of GLAG (see below) are zero.
            //     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
            //     CAUCHY will be set to the square of the KNEW-th Lagrange function at
            //       the step XALT-XOPT from XOPT for the vector XALT that is returned,
            //       except that CAUCHY is set to zero if XALT is not calculated.

            //     GLAG is a working space vector of length N for the gradient of the
            //       KNEW-th Lagrange function at XOPT.
            //     HCOL is a working space vector of length NPT for the second derivative
            //       coefficients of the KNEW-th Lagrange function.
            //     W is a working space vector of length 2N that is going to hold the
            //       constrained Cauchy step from XOPT of the Lagrange function, followed
            //       by the downhill version of XALT when the uphill step is calculated.

            var glag = new double[1 + n];
            var hcol = new double[1 + npt];
            var w = new double[1 + 2 * n];

            //     Set the first NPT components of W to the leading elements of the
            //     KNEW-th column of the H matrix.

            var @const = ONE + Math.Sqrt(2.0);

            for (var k = 1; k <= npt; ++k) hcol[k] = ZERO;
            for (var j = 1; j <= npt - n - 1; ++j)
            {
                var temp = zmat[knew, j];
                for (var k = 1; k <= npt; ++k) hcol[k] += temp * zmat[k, j];
            }
            alpha = hcol[knew];
            var ha = HALF * alpha;

            //     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
            //
            for (var i = 1; i <= n; ++i) glag[i] = bmat[knew, i];
            for (var k = 1; k <= npt; ++k)
            {
                var temp = ZERO;
                for (var j = 1; j <= n; ++j) temp += xpt[k, j] * xopt[j];
                temp *= hcol[k];
                for (var i = 1; i <= n; ++i) glag[i] += temp * xpt[k, i];
            }

            //     Search for a large denominator along the straight lines through XOPT
            //     and another interpolation point. SLBD and SUBD will be lower and upper
            //     bounds on the step along each of these lines in turn. PREDSQ will be
            //     set to the square of the predicted denominator for each line. PRESAV
            //     will be set to the largest admissible value of PREDSQ that occurs.

            var step = 0.0;
            var ksav = 0;
            var ibdsav = 0;
            var stpsav = 0.0;

            cauchy = ZERO;
            var csave = 0.0;

            var presav = ZERO;
            for (var k = 1; k <= npt; ++k)
            {
                if (k == kopt) continue;
                var dderiv = ZERO;
                var distsq = ZERO;
                for (var i = 1; i <= n; ++i)
                {
                    var temp = xpt[k, i] - xopt[i];
                    dderiv += glag[i] * temp;
                    distsq += temp * temp;
                }
                var subd = adelt / Math.Sqrt(distsq);
                var slbd = -subd;
                var ilbd = 0;
                var iubd = 0;
                var sumin = Math.Min(ONE, subd);

                //     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.

                for (var i = 1; i <= n; ++i)
                {
                    var temp = xpt[k, i] - xopt[i];
                    if (temp > ZERO)
                    {
                        if (slbd * temp < sl[i] - xopt[i])
                        {
                            slbd = (sl[i] - xopt[i]) / temp;
                            ilbd = -i;
                        }
                        if (subd * temp > su[i] - xopt[i])
                        {
                            subd = Math.Max(sumin, (su[i] - xopt[i]) / temp);
                            iubd = i;
                        }
                    }
                    else if (temp < ZERO)
                    {
                        if (slbd * temp > su[i] - xopt[i])
                        {
                            slbd = (su[i] - xopt[i]) / temp;
                            ilbd = i;
                        }
                        if (subd * temp < sl[i] - xopt[i])
                        {
                            subd = Math.Max(sumin, (sl[i] - xopt[i]) / temp);
                            iubd = -i;
                        }
                    }
                }

                //     Seek a large modulus of the KNEW-th Lagrange function when the index
                //     of the other interpolation point on the line through XOPT is KNEW.

                int isbd;
                double vlag;
                if (k == knew)
                {
                    var diff = dderiv - ONE;
                    step = slbd;
                    vlag = slbd * (dderiv - slbd * diff);
                    isbd = ilbd;
                    var temp = subd * (dderiv - subd * diff);
                    if (Math.Abs(temp) > Math.Abs(vlag))
                    {
                        step = subd;
                        vlag = temp;
                        isbd = iubd;
                    }
                    var tempd = HALF * dderiv;
                    var tempa = tempd - diff * slbd;
                    var tempb = tempd - diff * subd;
                    if (tempa * tempb < ZERO)
                    {
                        temp = tempd * tempd / diff;
                        if (Math.Abs(temp) > Math.Abs(vlag))
                        {
                            step = tempd / diff;
                            vlag = temp;
                            isbd = 0;
                        }
                    }

                //     Search along each of the other lines through XOPT and another point.

                }
                else
                {
                    step = slbd;
                    vlag = slbd * (ONE - slbd);
                    isbd = ilbd;
                    var temp = subd * (ONE - subd);
                    if (Math.Abs(temp) > Math.Abs(vlag))
                    {
                        step = subd;
                        vlag = temp;
                        isbd = iubd;
                    }
                    if (subd > HALF)
                    {
                        if (Math.Abs(vlag) < 0.25)
                        {
                            step = HALF;
                            vlag = 0.25;
                            isbd = 0;
                        }
                    }
                    vlag = vlag * dderiv;
                }

                //     Calculate PREDSQ for the current line search and maintain PRESAV.

                {
                    var temp = step * (ONE - step) * distsq;
                    var predsq = vlag * vlag * (vlag * vlag + ha * temp * temp);
                    if (predsq > presav)
                    {
                        presav = predsq;
                        ksav = k;
                        stpsav = step;
                        ibdsav = isbd;
                    }
                }
            }

            //     Construct XNEW in a way that satisfies the bound constraints exactly.
            //
            for (var i = 1; i <= n; ++i)
            {
                var temp = xopt[i] + stpsav * (xpt[ksav, i] - xopt[i]);
                xnew[i] = Math.Max(sl[i], Math.Min(su[i], temp));
            }
            if (ibdsav < 0) xnew[-ibdsav] = sl[-ibdsav];
            if (ibdsav > 0) xnew[ibdsav] = su[ibdsav];
            //
            //     Prepare for the iterative method that assembles the constrained Cauchy
            //     step in W. The sum of squares of the fixed components of W is formed in
            //     WFIXSQ, and the free components of W are set to BIGSTP.
            //
            var bigstp = adelt + adelt;

            for (var iflag = 0; iflag < 2; ++iflag)
            {
                var wfixsq = ZERO;
                var ggfree = ZERO;
                for (var I = 1; I <= n; ++I)
                {
                    w[I] = ZERO;
                    var tempa = Math.Min(xopt[I] - sl[I], glag[I]);
                    var tempb = Math.Max(xopt[I] - su[I], glag[I]);
                    if (tempa > ZERO || tempb < ZERO)
                    {
                        w[I] = bigstp;
                        ggfree += glag[I] * glag[I];
                    }
                }
                if (ggfree == ZERO)
                {
                    cauchy = ZERO;
                    return;
                }

                //     Investigate whether more components of W can be fixed.

                {
                    L_120:
                    var temp = adelt * adelt - wfixsq;

                    if (temp > ZERO)
                    {
                        var wsqsav = wfixsq;
                        step = Math.Sqrt(temp / ggfree);
                        ggfree = ZERO;
                        for (var I = 1; I <= n; ++I)
                        {
                            if (w[I] == bigstp)
                            {
                                temp = xopt[I] - step * glag[I];
                                if (temp <= sl[I])
                                {
                                    w[I] = sl[I] - xopt[I];
                                    wfixsq += w[I] * w[I];
                                }
                                else if (temp >= su[I])
                                {
                                    w[I] = su[I] - xopt[I];
                                    wfixsq += w[I] * w[I];
                                }
                                else
                                {
                                    ggfree += glag[I] * glag[I];
                                }
                            }
                        }
                        if (wfixsq > wsqsav && ggfree > ZERO) goto L_120;
                    }
                }

                //     Set the remaining free components of W and all components of XALT,
                //     except that W may be scaled later.

                var gw = ZERO;
                for (var i = 1; i <= n; ++i)
                {
                    if (w[i] == bigstp)
                    {
                        w[i] = -step * glag[i];
                        xalt[i] = Math.Max(sl[i], Math.Min(su[i], xopt[i] + w[i]));
                    }
                    else if (w[i] == ZERO)
                    {
                        xalt[i] = xopt[i];
                    }
                    else if (glag[i] > ZERO)
                    {
                        xalt[i] = sl[i];
                    }
                    else
                    {
                        xalt[i] = su[i];
                    }
                    gw += glag[i] * w[i];
                }

                //     Set CURV to the curvature of the KNEW-th Lagrange function along W.
                //     Scale W by a factor less than one if that can reduce the modulus of
                //     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
                //     the square of this function.

                var curv = ZERO;
                for (var k = 1; k <= npt; ++k)
                {
                    var temp = ZERO;
                    for (var j = 1; j <= n; ++j) temp += xpt[k, j] * w[j];
                    curv += hcol[k] * temp * temp;
                }
                if (iflag == 1) curv = -curv;
                if (curv > -gw && curv < -@const * gw)
                {
                    var scale = -gw / curv;
                    for (var i = 1; i <= n; ++i)
                    {
                        var temp = xopt[i] + scale * w[i];
                        xalt[i] = Math.Max(sl[i], Math.Min(su[i], temp));
                    }
                    cauchy = Math.Pow(HALF * gw * scale, 2.0);
                }
                else
                {
                    cauchy = Math.Pow(gw + HALF * curv, 2.0);
                }

                //     If IFLAG is zero, then XALT is calculated as before after reversing
                //     the sign of GLAG. Thus two XALT vectors become available. The one that
                //     is chosen is the one that gives the larger value of CAUCHY.

                if (iflag == 0)
                {
                    for (var I = 1; I <= n; ++I)
                    {
                        glag[I] = -glag[I];
                        w[n + I] = xalt[I];
                    }
                    csave = cauchy;
                }
            }

            if (csave > cauchy)
            {
                for (var I = 1; I <= n; ++I) xalt[I] = w[n + I];
                cauchy = csave;
            }
        }

        private static void PRELIM(Func<int, double[], double> calfun, int n, int npt, double[] x, 
            double[] xl, double[] xu, double rhobeg, int iprint, int maxfun, double[] xbase, double[,] xpt, 
            double[] fval, double[] gopt, double[] hq, double[] pq, double[,] bmat, double[,] zmat,
            int ndim, double[] sl, double[] su, out int nf, out int kopt, TextWriter logger)
        {
            //     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
            //       same as the corresponding arguments in SUBROUTINE BOBYQA.
            //     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
            //       are the same as the corresponding arguments in BOBYQB, the elements
            //       of SL and SU being set in BOBYQA.
            //     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
            //       it is set by PRELIM to the gradient of the quadratic model at XBASE.
            //       If XOPT is nonzero, BOBYQB will change it to its usual value later.
            //     NF is maintaned as the number of calls of CALFUN so far.
            //     KOPT will be such that the least calculated value of F so far is at
            //       the point XPT(KOPT,.)+XBASE in the space of the variables.
            //
            //     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
            //     BMAT and ZMAT for the first iteration, and it maintains the values of
            //     NF and KOPT. The vector X is also changed by PRELIM.

            //     Set some constants.

            var rhosq = rhobeg * rhobeg;
            var recip = ONE / rhosq;
            var np = n + 1;

            kopt = 0;

            //     Set XBASE to the initial vector of variables, and set the initial
            //     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.

            for (var j = 1; j <= n; ++j)
            {
                xbase[j] = x[j];
                for (var k = 1; k <= npt; ++k) xpt[k, j] = ZERO;
                for (var i = 1; i <= ndim; ++i) bmat[i, j] = ZERO;
            }
            for (var ih = 1; ih <= n * np / 2; ++ih) hq[ih] = ZERO;
            for (var k = 1; k <= npt; ++k)
            {
                pq[k] = ZERO;
                for (var j = 1; j <= npt - np; ++j) zmat[k, j] = ZERO;
            }

            //     Begin the initialization procedure. NF becomes one more than the number
            //     of function values so far. The coordinates of the displacement of the
            //     next initial interpolation point from XBASE are set in XPT(NF+1,.).

            var ipt = 0;
            var jpt = 0;
            var stepa = 0.0;
            var stepb = 0.0;
            var fbeg = 0.0;

            for (nf = 1; nf < Math.Min(npt, maxfun); ++nf)
            {
                var nfm = nf - 1;
                var nfx = nf - 1 - n;

                if (nfm <= 2 * n)
                {
                    if (nfm >= 1 && nfm <= n)
                    {
                        stepa = rhobeg;
                        if (su[nfm] == ZERO) stepa = -stepa;
                        xpt[nf, nfm] = stepa;
                    }
                    else if (nfm > n)
                    {
                        stepa = xpt[nf - n, nfx];
                        stepb = -rhobeg;
                        if (sl[nfx] == ZERO) stepb = Math.Min(TWO * rhobeg, su[nfx]);
                        if (su[nfx] == ZERO) stepb = Math.Max(-TWO * rhobeg, sl[nfx]);
                        xpt[nf, nfx] = stepb;
                    }
                }
                else
                {
                    var itemp = (nfm - np) / n;
                    jpt = nfm - itemp * n - n;
                    ipt = jpt + itemp;
                    if (ipt > n)
                    {
                        itemp = jpt;
                        jpt = ipt - n;
                        ipt = itemp;
                    }
                    xpt[nf, ipt] = xpt[ipt + 1, ipt];
                    xpt[nf, jpt] = xpt[jpt + 1, jpt];
                }

                //     Calculate the next value of F. The least function value so far and
                //     its index are required.

                for (var j = 1; j <= n; ++j)
                {
                    x[j] = Math.Min(Math.Max(xl[j], xbase[j] + xpt[nf, j]), xu[j]);
                    if (xpt[nf, j] == sl[j]) x[j] = xl[j];
                    if (xpt[nf, j] == su[j]) x[j] = xu[j];
                }

                var f = calfun(n, x);
                if (iprint == 3 && logger != null) logger.WriteLine(IterationOutputFormat, nf, f, ToString(x, n));
                fval[nf] = f;
                if (nf == 1)
                {
                    fbeg = f;
                    kopt = 1;
                }
                else if (f < fval[kopt])
                {
                    kopt = nf;
                }

                //     Set the nonzero initial elements of BMAT and the quadratic model in the
                //     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
                //     of the NF-th and (NF-N)-th interpolation points may be switched, in
                //     order that the function value at the first of them contributes to the
                //     off-diagonal second derivative terms of the initial quadratic model.

                if (nf <= 2 * n + 1)
                {
                    if (nf >= 2 && nf <= n + 1)
                    {
                        gopt[nfm] = (f - fbeg) / stepa;
                        if (npt < nf + n)
                        {
                            bmat[1, nfm] = -ONE / stepa;
                            bmat[nf, nfm] = ONE / stepa;
                            bmat[npt + nfm, nfm] = -HALF * rhosq;
                        }
                    }
                    else if (nf >= n + 2)
                    {
                        var ih = (nfx * (nfx + 1)) / 2;
                        var temp = (f - fbeg) / stepb;
                        var diff = stepb - stepa;
                        hq[ih] = TWO * (temp - gopt[nfx]) / diff;
                        gopt[nfx] = (gopt[nfx] * stepb - temp * stepa) / diff;
                        if (stepa * stepb < ZERO)
                        {
                            if (f < fval[nf - n])
                            {
                                fval[nf] = fval[nf - n];
                                fval[nf - n] = f;
                                if (kopt == nf) kopt = nf - n;
                                xpt[nf - n, nfx] = stepb;
                                xpt[nf, nfx] = stepa;
                            }
                        }
                        bmat[1, nfx] = -(stepa + stepb) / (stepa * stepb);
                        bmat[nf, nfx] = -HALF / xpt[nf - n, nfx];
                        bmat[nf - n, nfx] = -bmat[1, nfx] - bmat[nf, nfx];
                        zmat[1, nfx] = Math.Sqrt(TWO) / (stepa * stepb);
                        zmat[nf, nfx] = Math.Sqrt(HALF) / rhosq;
                        zmat[nf - n, nfx] = -zmat[1, nfx] - zmat[nf, nfx];
                    }

                }
                else
                {
                    //     Set the off-diagonal second derivatives of the Lagrange functions and
                    //     the initial quadratic model.

                    var ih = (ipt * (ipt - 1)) / 2 + jpt;
                    zmat[1, nfx] = recip;
                    zmat[nf, nfx] = recip;
                    zmat[ipt + 1, nfx] = -recip;
                    zmat[jpt + 1, nfx] = -recip;
                    var temp = xpt[nf, ipt] * xpt[nf, jpt];
                    hq[ih] = (fbeg - fval[ipt + 1] - fval[jpt + 1] + f) / temp;
                }
            }
        }

        private static void RESCUE(Func<int, double[], double> calfun,  int n, int npt, double[] xl, double[] xu, int iprint,
            int maxfun, double[] xbase, double[,] xpt, double[] fval, double[] xopt, double[] gopt,
            double[] hq, double[] pq, double[,] bmat, double[,] zmat, int ndim, double[] sl, double[] su,
            ref int nf, double delta, ref int kopt, double[] vlag, TextWriter logger)
        {
            //     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
            //       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
            //       the corresponding arguments of BOBYQB on the entry to RESCUE.
            //     NF is maintained as the number of calls of CALFUN so far, except that
            //       NF is set to -1 if the value of MAXFUN prevents further progress.
            //     KOPT is maintained so that FVAL(KOPT) is the least calculated function
            //       value. Its correct value must be given on entry. It is updated if a
            //       new least function value is found, but the corresponding changes to
            //       XOPT and GOPT have to be made later by the calling program.
            //     DELTA is the current trust region radius.
            //     VLAG is a working space vector that will be used for the values of the
            //       provisional Lagrange functions at each of the interpolation points.
            //       They are part of a product that requires VLAG to be of length NDIM.
            //     The final elements of BMAT and ZMAT are set in a well-conditioned way
            //       to the values that are appropriate for the new interpolation points.
            //     The elements of GOPT, HQ and PQ are also revised to the values that are
            //       appropriate to the final quadratic model.

            //     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and
            //       PTSAUX(2,J) specify the two positions of provisional interpolation
            //       points when a nonzero step is taken along e_J (the J-th coordinate
            //       direction) through XBASE+XOPT, as specified below. Usually these
            //       steps have length DELTA, but other lengths are chosen if necessary
            //       in order to satisfy the given bounds on the variables.
            //     PTSID is also a working space array. It has NPT components that denote
            //       provisional new positions of the original interpolation points, in
            //       case changes are needed to restore the linear independence of the
            //       interpolation conditions. The K-th point is a candidate for change
            //       if and only if PTSID[K] is nonzero. In this case let p and q be the
            //       integer parts of PTSID[K] and (PTSID[K]-p) multiplied by N+1. If p
            //       and q are both positive, the step from XBASE+XOPT to the new K-th
            //       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
            //       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
            //       p=0, respectively.
            //     The first NDIM+NPT elements of the array W are used for working space. 

            var ptsaux = new double[1 + 2, 1 + n];
            var ptsid = new double[1 + npt];
            var w = new double[1 + ndim + npt];

            //     Set some constants.
            //
            var np = n + 1;
            var sfrac = HALF / np;
            var nptm = npt - np;

            //     Shift the interpolation points so that XOPT becomes the origin, and set
            //     the elements of ZMAT to zero. The value of SUMPQ is required in the
            //     updating of HQ below. The squares of the distances from XOPT to the
            //     other interpolation points are set at the end of W. Increments of WINC
            //     may be added later to these squares to balance the consideration of
            //     the choice of point that is going to become current.

            var sumpq = ZERO;
            var winc = ZERO;
            for (var k = 1; k <= npt; ++k)
            {
                var distsq = ZERO;
                for (var j = 1; j <= n; ++j)
                {
                    xpt[k, j] -= xopt[j];
                    distsq += xpt[k, j] * xpt[k, j];
                }
                sumpq += pq[k];
                w[ndim + k] = distsq;
                winc = Math.Max(winc, distsq);
                for (var j = 1; j <= nptm; ++j) zmat[k, j] = ZERO;
            }

            //     Update HQ so that HQ and PQ define the second derivatives of the model
            //     after XBASE has been shifted to the trust region centre.

            {
                var ih = 0;
                for (var j = 1; j <= n; ++j)
                {
                    w[j] = HALF * sumpq * xopt[j];
                    for (var k = 1; k <= npt; ++k) w[j] += pq[k] * xpt[k, j];
                    for (var i = 1; i <= j; ++i) hq[++ih] += w[i] * xopt[j] + w[j] * xopt[i];
                }
            }

            //     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
            //     also set the elements of PTSAUX.

            for (var j = 1; j <= n; ++j)
            {
                xbase[j] += xopt[j];
                sl[j] -= xopt[j];
                su[j] -= xopt[j];
                xopt[j] = ZERO;
                ptsaux[1, j] = Math.Min(delta, su[j]);
                ptsaux[2, j] = Math.Max(-delta, sl[j]);
                if (ptsaux[1, j] + ptsaux[2, j] < ZERO)
                {
                    var temp = ptsaux[1, j];
                    ptsaux[1, j] = ptsaux[2, j];
                    ptsaux[2, j] = temp;
                }
                if (Math.Abs(ptsaux[2, j]) < HALF * Math.Abs(ptsaux[1, j])) ptsaux[2, j] = HALF * ptsaux[1, j];
                for (var i = 1; i <= ndim; ++i) bmat[i, j] = ZERO;
            }
            var fbase = fval[kopt];

            //     Set the identifiers of the artificial interpolation points that are
            //     along a coordinate direction from XOPT, and set the corresponding
            //     nonzero elements of BMAT and ZMAT.

            ptsid[1] = sfrac;
            for (var j = 1; j <= n; ++j)
            {
                var jp = j + 1;
                var jpn = jp + n;
                ptsid[jp] = j + sfrac;
                if (jpn <= npt)
                {
                    ptsid[jpn] = (double)j / np + sfrac;
                    var temp = ONE / (ptsaux[1, j] - ptsaux[2, j]);
                    bmat[jp, j] = -temp + ONE / ptsaux[1, j];
                    bmat[jpn, j] = temp + ONE / ptsaux[2, j];
                    bmat[1, j] = -bmat[jp, j] - bmat[jpn, j];
                    zmat[1, j] = Math.Sqrt(2.0) / Math.Abs(ptsaux[1, j] * ptsaux[2, j]);
                    zmat[jp, j] = zmat[1, j] * ptsaux[2, j] * temp;
                    zmat[jpn, j] = -zmat[1, j] * ptsaux[1, j] * temp;
                }
                else
                {
                    bmat[1, j] = -ONE / ptsaux[1, j];
                    bmat[jp, j] = ONE / ptsaux[1, j];
                    bmat[j + npt, j] = -HALF * ptsaux[1, j] * ptsaux[1, j];
                }
            }

            //     Set any remaining identifiers with their nonzero elements of ZMAT.

            if (npt >= n + np)
            {
                for (var k = 2 * np; k <= npt; ++k)
                {
                    var iw = (int)(k - np - HALF) / n;
                    var ip = k - np - iw * n;
                    var iq = ip + iw;
                    if (iq > n) iq = iq - n;
                    ptsid[k] = ip + (double)iq / np + sfrac;
                    var temp = ONE / (ptsaux[1, ip] * ptsaux[1, iq]);
                    zmat[1, k - np] = temp;
                    zmat[ip + 1, k - np] = -temp;
                    zmat[iq + 1, k - np] = -temp;
                    zmat[k, k - np] = temp;
                }
            }
            var nrem = npt;
            var kold = 1;
            var knew = kopt;

            //     Reorder the provisional points in the way that exchanges PTSID(KOLD)
            //     with PTSID(KNEW).

            var beta = 0.0;
            var denom = 0.0;

            do
            {
                for (var j = 1; j <= n; ++j)
                {
                    var temp = bmat[kold, j];
                    bmat[kold, j] = bmat[knew, j];
                    bmat[knew, j] = temp;
                }
                for (var j = 1; j <= nptm; ++j)
                {
                    var temp = zmat[kold, j];
                    zmat[kold, j] = zmat[knew, j];
                    zmat[knew, j] = temp;
                }
                ptsid[kold] = ptsid[knew];
                ptsid[knew] = ZERO;
                w[ndim + knew] = ZERO;
                --nrem;
                if (knew != kopt)
                {
                    var temp = vlag[kold];
                    vlag[kold] = vlag[knew];
                    vlag[knew] = temp;

                    //     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
                    //     interpolation point can be changed from provisional to original. The
                    //     branch to label 350 occurs if all the original points are reinstated.
                    //     The nonnegative values of W(NDIM+K) are required in the search below.

                    UPDATE(n, npt, bmat, zmat, ndim, vlag, beta, denom, knew, w);
                    if (nrem == 0) return;
                    for (var k = 1; k <= npt; ++k) w[ndim + k] = Math.Abs(w[ndim + k]);
                }

                //     Pick the index KNEW of an original interpolation point that has not
                //     yet replaced one of the provisional interpolation points, giving
                //     attention to the closeness to XOPT and to previous tries with KNEW.

                L_120:
                var dsqmin = ZERO;
                for (var k = 1; k <= npt; ++k)
                {
                    if (w[ndim + k] > ZERO && (dsqmin == ZERO || w[ndim + k] < dsqmin))
                    {
                        knew = k;
                        dsqmin = w[ndim + k];
                    }
                }
                if (dsqmin == ZERO) break;

                //     Form the W-vector of the chosen original interpolation point.

                for (var j = 1; j <= n; ++j) w[npt + j] = xpt[knew, j];
                for (var k = 1; k <= npt; ++k)
                {
                    var sum = ZERO;
                    if (k == kopt) continue;
                    if (ptsid[k] == ZERO)
                        for (var j = 1; j <= n; ++j) sum += w[npt + j] * xpt[k, j];
                    else
                    {
                        var ip = (int)ptsid[k];
                        if (ip > 0) sum = w[npt + ip] * ptsaux[1, ip];
                        var iq = (int)(np * ptsid[k] - ip * np);
                        if (iq > 0)
                        {
                            var iw = ip == 0 ? 2 : 1;
                            sum += w[npt + iq] * ptsaux[iw, iq];
                        }
                    }
                    w[k] = HALF * sum * sum;
                }

                //     Calculate VLAG and BETA for the required updating of the H matrix if
                //     XPT(KNEW,.) is reinstated in the set of interpolation points.
                //
                for (var k = 1; k <= npt; ++k)
                {
                    var sum = ZERO;
                    for (var j = 1; j <= n; ++j) sum += bmat[k, j] * w[npt + j];
                    vlag[k] = sum;
                }
                beta = ZERO;
                for (var j = 1; j <= nptm; ++j)
                {
                    var sum = ZERO;
                    for (var k = 1; k <= npt; ++k) sum += zmat[k, j] * w[k];
                    beta -= sum * sum;
                    for (var k = 1; k <= npt; ++k) vlag[k] = vlag[k] + sum * zmat[k, j];
                }
                {
                    var bsum = ZERO;
                    var distsq = ZERO;
                    for (var j = 1; j <= n; ++j)
                    {
                        var sum = ZERO;
                        for (var k = 1; k <= npt; ++k) sum += bmat[k, j] * w[k];
                        var jp = j + npt;
                        bsum += sum * w[jp];
                        for (var ip = npt + 1; ip <= ndim; ++ip) sum = sum + bmat[ip, j] * w[ip];
                        bsum += sum * w[jp];
                        vlag[jp] = sum;
                        distsq += xpt[knew, j] * xpt[knew, j];
                    }
                    beta += HALF * distsq * distsq - bsum;
                    vlag[kopt] += ONE;
                }

                //     KOLD is set to the index of the provisional interpolation point that is
                //     going to be deleted to make way for the KNEW-th original interpolation
                //     point. The choice of KOLD is governed by the avoidance of a small value
                //     of the denominator in the updating calculation of UPDATE.

                denom = ZERO;
                var vlmxsq = ZERO;
                for (var k = 1; k <= npt; ++k)
                {
                    if (ptsid[k] != ZERO)
                    {
                        var hdiag = ZERO;
                        for (var j = 1; j <= nptm; ++j) hdiag += zmat[k, j] * zmat[k, j];
                        var den = beta * hdiag + vlag[k] * vlag[k];
                        if (den > denom)
                        {
                            kold = k;
                            denom = den;
                        }
                    }
                    vlmxsq = Math.Max(vlmxsq, vlag[k] * vlag[k]);
                }

                if (denom <= 0.01 * vlmxsq)
                {
                    w[ndim + knew] = -w[ndim + knew] - winc;
                    goto L_120;
                }
            } while (true);

            //     When do loop is completed, all the final positions of the interpolation
            //     points have been chosen although any changes have not been included yet
            //     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
            //     from the shift of XBASE, the updating of the quadratic model remains to
            //     be done. The following cycle through the new interpolation points begins
            //     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
            //     except that a RETURN occurs if MAXFUN prohibits another value of F.

            for (var kpt = 1; kpt <= npt; ++kpt)
            {
                if (ptsid[kpt] == ZERO) continue;
                if (nf >= maxfun)
                {
                    nf = -1;
                    return;
                }

                var ih = 0;
                for (var j = 1; j <= n; ++j)
                {
                    w[j] = xpt[kpt, j];
                    xpt[kpt, j] = ZERO;
                    var temp = pq[kpt] * w[j];
                    for (var i = 1; i <= j; ++i) hq[++ih] += temp * w[i];
                }
                pq[kpt] = ZERO;

                var ip = (int)ptsid[kpt];
                var iq = (int)(np * ptsid[kpt] - ip * np);
                var xp = 0.0;
                var xq = 0.0;
                if (ip > 0)
                {
                    xp = ptsaux[1, ip];
                    xpt[kpt, ip] = xp;
                }
                if (iq > 0)
                {
                    xq = ip == 0 ? ptsaux[2, iq] : ptsaux[1, iq];
                    xpt[kpt, iq] = xq;
                }

                //     Set VQUAD to the value of the current model at the new point.

                var vquad = fbase;
                var ihp = 0;
                if (ip > 0)
                {
                    ihp = (ip + ip * ip) / 2;
                    vquad += xp * (gopt[ip] + HALF * xp * hq[ihp]);
                }
                if (iq > 0)
                {
                    var ihq = (iq + iq * iq) / 2;
                    vquad += xq * (gopt[iq] + HALF * xq * hq[ihq]);
                    if (ip > 0)
                    {
                        var iw = Math.Max(ihp, ihq) - Math.Abs(ip - iq);
                        vquad += xp * xq * hq[iw];
                    }
                }
                for (var k = 1; k <= npt; ++k)
                {
                    var temp = ZERO;
                    if (ip > 0) temp += xp * xpt[k, ip];
                    if (iq > 0) temp += xq * xpt[k, iq];
                    vquad += HALF * pq[k] * temp * temp;
                }

                //     Calculate F at the new interpolation point, and set DIFF to the factor
                //     that is going to multiply the KPT-th Lagrange function when the model
                //     is updated to provide interpolation to the new function value.

                for (var i = 1; i <= n; ++i)
                {
                    w[i] = Math.Min(Math.Max(xl[i], xbase[i] + xpt[kpt, i]), xu[i]);
                    if (xpt[kpt, i] == sl[i]) w[i] = xl[i];
                    if (xpt[kpt, i] == su[i]) w[i] = xu[i];
                }

                ++nf;
                var f = calfun(n, w);
                if (iprint == 3 && logger != null) logger.WriteLine(IterationOutputFormat, nf, f, ToString(w, n));

                fval[kpt] = f;
                if (f < fval[kopt]) kopt = kpt;
                var diff = f - vquad;

                //     Update the quadratic model. The RETURN from the subroutine occurs when
                //     all the new interpolation points are included in the model.

                for (var i = 1; i <= n; ++i) gopt[i] += diff * bmat[kpt, i];
                for (var k = 1; k <= npt; ++k)
                {
                    var sum = ZERO;
                    for (var j = 1; j <= nptm; ++j) sum += zmat[k, j] * zmat[kpt, j];
                    var temp = diff * sum;
                    if (ptsid[k] == ZERO)
                        pq[k] += temp;
                    else
                    {
                        ip = (int)ptsid[k];
                        iq = (int)(np * ptsid[k] - ip * np);
                        var ihq = (iq * iq + iq) / 2;
                        if (ip == 0)
                            hq[ihq] += temp * ptsaux[2, iq] * ptsaux[2, iq];
                        else
                        {
                            ihp = (ip * ip + ip) / 2;
                            hq[ihp] += temp * ptsaux[1, ip] * ptsaux[1, ip];
                            if (iq > 0)
                            {
                                hq[ihq] += temp * ptsaux[1, iq] * ptsaux[1, iq];
                                var iw = Math.Max(ihp, ihq) - Math.Abs(iq - ip);
                                hq[iw] += temp * ptsaux[1, ip] * ptsaux[1, iq];
                            }
                        }
                    }
                }
                ptsid[kpt] = ZERO;
            }
        }

        private static void TRSBOX(int n, int npt, double[,] xpt, double[] xopt, double[] gopt,
            double[] hq, double[] pq, double[] sl, double[] su, double delta,
            double[] xnew, double[] d, double[] gnew, out double dsq, out double crvmin)
        {
            //     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
            //       meanings as the corresponding arguments of BOBYQB.
            //     DELTA is the trust region radius for the present calculation, which
            //       seeks a small value of the quadratic model within distance DELTA of
            //       XOPT subject to the bounds on the variables.
            //     XNEW will be set to a new vector of variables that is approximately
            //       the one that minimizes the quadratic model within the trust region
            //       subject to the SL and SU constraints on the variables. It satisfies
            //       as equations the bounds that become active during the calculation.
            //     D is the calculated trial step from XOPT, generated iteratively from an
            //       initial value of zero. Thus XNEW is XOPT+D after the final iteration.
            //     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
            //       when D is updated.
            //     DSQ will be set to the square of the length of XNEW-XOPT.
            //     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
            //       it is set to the least curvature of H that occurs in the conjugate
            //       gradient searches that are not restricted by any constraints. The
            //       value CRVMIN=-1.0D0 is set, however, if all of these searches are
            //       constrained.
            //
            //     A version of the truncated conjugate gradient is applied. If a line
            //     search is restricted by a constraint, then the procedure is restarted,
            //     the values of the variables that are at their bounds being fixed. If
            //     the trust region boundary is reached, then further changes may be made
            //     to D, each one being in the two dimensional space that is spanned
            //     by the current D and the gradient of Q at XOPT+D, staying on the trust
            //     region boundary. Termination occurs when the reduction in Q seems to
            //     be close to the greatest reduction that can be achieved.

            //     The sign of GOPT[I] gives the sign of the change to the I-th variable
            //     that will reduce Q from its value at XOPT. Thus XBDI[I] shows whether
            //     or not to fix the I-th variable at one of its bounds initially, with
            //     NACT being set to the number of fixed variables. D and GNEW are also
            //     set for the first iteration. DELSQ is the upper bound on the sum of
            //     squares of the free variables. QRED is the reduction in Q so far.

            //     XBDI is a working space vector. For I=1,2,...,N, the element XBDI[I] is
            //       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
            //       I-th variable has become fixed at a bound, the bound being SL[I] or
            //       SU[I] in the case XBDI[I]=-1.0 or XBDI[I]=1.0, respectively. This
            //       information is accumulated during the construction of XNEW.
            //     The arrays S, HS and HRED are also used for working space. They hold the
            //       current search direction, and the changes in the gradient of Q along S
            //       and the reduced D, respectively, where the reduced D is the same as D,
            //       except that the components of the fixed variables are zero.

            var xbdi = new double[1 + n];
            var s = new double[1 + n];
            var hs = new double[1 + n];
            var hred = new double[1 + n];

            var iterc = 0;
            var nact = 0;

            for (var i = 1; i <= n; ++i)
            {
                xbdi[i] = ZERO;
                if (xopt[i] <= sl[i])
                    if (gopt[i] >= ZERO) xbdi[i] = ONEMIN;
                    else if (xopt[i] >= su[i])
                        if (gopt[i] <= ZERO) xbdi[i] = ONE;
                if (xbdi[i] != ZERO) ++nact;
                d[i] = ZERO;
                gnew[i] = gopt[i];
            }

            var delsq = delta * delta;
            var qred = ZERO;
            crvmin = ONEMIN;

            //     Set the next search direction of the conjugate gradient method. It is
            //     the steepest descent direction initially and when the iterations are
            //     restarted because a variable has just been fixed by a bound, and of
            //     course the components of the fixed variables are zero. ITERMAX is an
            //     upper bound on the indices of the conjugate gradient iterations.

            double sth;

            var itermax = 0;
            var itcsav = 0;
            var iact = 0;

            var angt = 0.0;
            var angbd = 0.0;
            var dredsq = 0.0;
            var dredg = 0.0;
            var gredsq = 0.0;
            var ggsav = 0.0;
            var rdprev = 0.0;
            var rdnext = 0.0;
            var sredg = 0.0;
            var xsav = 0.0;

            L_20:
            var beta = ZERO;

            L_30:
            var stepsq = ZERO;

            for (var i = 1; i <= n; ++i)
            {
                if (xbdi[i] != ZERO)
                    s[i] = ZERO;
                else if (beta == ZERO)
                    s[i] = -gnew[i];
                else
                    s[i] = beta * s[i] - gnew[i];
                stepsq += s[i] * s[i];
            }
            if (stepsq == ZERO) goto L_190;

            if (beta == ZERO)
            {
                gredsq = stepsq;
                itermax = iterc + n - nact;
            }
            if (gredsq * delsq <= 1.0E-4 * qred * qred) goto L_190;

            //     Multiply the search direction by the second derivative matrix of Q and
            //     calculate some scalars for the choice of steplength. Then set BLEN to
            //     the length of the the step to the trust region boundary and STPLEN to
            //     the steplength, ignoring the simple bounds.

            goto L_210;

            L_50:
            var resid = delsq;
            var ds = ZERO;
            var shs = ZERO;
            for (var i = 1; i <= n; ++i)
            {
                if (xbdi[i] == ZERO)
                {
                    resid -= d[i] * d[i];
                    ds += s[i] * d[i];
                    shs += s[i] * hs[i];
                }
            }
            if (resid <= ZERO) goto L_90;

            var temp = Math.Sqrt(stepsq * resid + ds * ds);
            var blen = ds < ZERO ? (temp - ds) / stepsq : resid / (temp + ds);
            var stplen = shs > ZERO ? Math.Min(blen, gredsq / shs) : blen;

            //     Reduce STPLEN if necessary in order to preserve the simple bounds,
            //     letting IACT be the index of the new constrained variable.

            iact = 0;
            for (var i = 1; i <= n; ++i)
            {
                if (s[i] != ZERO)
                {
                    var xsum = xopt[i] + d[i];
                    if (s[i] > ZERO)
                    {
                        temp = (su[i] - xsum) / s[i];
                    }
                    else
                    {
                        temp = (sl[i] - xsum) / s[i];
                    }
                    if (temp < stplen)
                    {
                        stplen = temp;
                        iact = i;
                    }
                }
            }

            //     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.

            var sdec = ZERO;
            if (stplen > ZERO)
            {
                ++iterc;
                temp = shs / stepsq;
                if (iact == 0 && temp > ZERO)
                {
                    crvmin = Math.Min(crvmin, temp);
                    if (crvmin == ONEMIN) crvmin = temp;
                }
                ggsav = gredsq;
                gredsq = ZERO;
                for (var i = 1; i <= n; ++i)
                {
                    gnew[i] += stplen * hs[i];
                    if (xbdi[i] == ZERO) gredsq += gnew[i] * gnew[i];
                    d[i] += stplen * s[i];
                }
                sdec = Math.Max(stplen * (ggsav - HALF * stplen * shs), ZERO);
                qred += sdec;
            }

            //     Restart the conjugate gradient method if it has hit a new bound.

            if (iact > 0)
            {
                ++nact;
                xbdi[iact] = ONE;
                if (s[iact] < ZERO) xbdi[iact] = ONEMIN;
                delsq -= d[iact] * d[iact];
                if (delsq <= ZERO) goto L_90;
                goto L_20;
            }

            //     If STPLEN is less than BLEN, then either apply another conjugate
            //     gradient iteration or RETURN.

            if (stplen < blen)
            {
                if (iterc == itermax || sdec <= 0.01 * qred) goto L_190;
                beta = gredsq / ggsav;
                goto L_30;
            }

            L_90:
            crvmin = ZERO;

            //     Prepare for the alternative iteration by calculating some scalars
            //     and by multiplying the reduced D by the second derivative matrix of
            //     Q, where S holds the reduced D in the call of GGMULT.

            L_100:
            if (nact >= n - 1) goto L_190;

            dredsq = ZERO;
            dredg = ZERO;
            gredsq = ZERO;
            for (var i = 1; i <= n; ++i)
            {
                if (xbdi[i] == ZERO)
                {
                    dredsq += d[i] * d[i];
                    dredg += d[i] * gnew[i];
                    gredsq += gnew[i] * gnew[i];
                    s[i] = d[i];
                }
                else
                {
                    s[i] = ZERO;
                }
            }
            itcsav = iterc;
            goto L_210;

            //     Let the search direction S be a linear combination of the reduced D
            //     and the reduced G that is orthogonal to the reduced D.

            L_120:
            ++iterc;
            temp = gredsq * dredsq - dredg * dredg;
            if (temp <= 1.0E-4 * qred * qred) goto L_190;

            temp = Math.Sqrt(temp);
            for (var i = 1; i <= n; ++i) s[i] = xbdi[i] == ZERO ? (dredg * d[i] - dredsq * gnew[i]) / temp : ZERO;

            sredg = -temp;

            //     By considering the simple bounds on the variables, calculate an upper
            //     bound on the tangent of half the angle of the alternative iteration,
            //     namely ANGBD, except that, if already a free variable has reached a
            //     bound, there is a branch back to label 100 after fixing that variable.

            angbd = ONE;
            iact = 0;
            for (var i = 1; i <= n; ++i)
            {
                if (xbdi[i] == ZERO)
                {
                    var tempa = xopt[i] + d[i] - sl[i];
                    var tempb = su[i] - xopt[i] - d[i];
                    if (tempa <= ZERO)
                    {
                        ++nact;
                        xbdi[i] = ONEMIN;
                        goto L_100;
                    }
                    if (tempb <= ZERO)
                    {
                        ++nact;
                        xbdi[i] = ONE;
                        goto L_100;
                    }
                    var ssq = d[i] * d[i] + s[i] * s[i];
                    temp = ssq - Math.Pow(xopt[i] - sl[i], 2.0);
                    if (temp > ZERO)
                    {
                        temp = Math.Sqrt(temp) - s[i];
                        if (angbd * temp > tempa)
                        {
                            angbd = tempa / temp;
                            iact = i;
                            xsav = ONEMIN;
                        }
                    }
                    temp = ssq - Math.Pow(su[i] - xopt[i], 2.0);
                    if (temp > ZERO)
                    {
                        temp = Math.Sqrt(temp) + s[i];
                        if (angbd * temp > tempb)
                        {
                            angbd = tempb / temp;
                            iact = i;
                            xsav = ONE;
                        }
                    }
                }
            }

            //     Calculate HHD and some curvatures for the alternative iteration.

            goto L_210;

            L_150:
            shs = ZERO;
            var dhs = ZERO;
            var dhd = ZERO;
            for (var i = 1; i <= n; ++i)
            {
                if (xbdi[i] == ZERO)
                {
                    shs += s[i] * hs[i];
                    dhs += d[i] * hs[i];
                    dhd += d[i] * hred[i];
                }
            }

            //     Seek the greatest reduction in Q for a range of equally spaced values
            //     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
            //     the alternative iteration.

            var redmax = ZERO;
            var isav = 0;
            var redsav = ZERO;
            var iu = (int)(17.0 * angbd + 3.1);
            for (var i = 1; i <= iu; ++i)
            {
                angt = angbd * i / iu;
                sth = (angt + angt) / (ONE + angt * angt);
                temp = shs + angt * (angt * dhd - dhs - dhs);
                var rednew = sth * (angt * dredg - sredg - HALF * sth * temp);
                if (rednew > redmax)
                {
                    redmax = rednew;
                    isav = i;
                    rdprev = redsav;
                }
                else if (i == isav + 1)
                {
                    rdnext = rednew;
                }
                redsav = rednew;
            }

            //     Return if the reduction is zero. Otherwise, set the sine and cosine
            //     of the angle of the alternative iteration, and calculate SDEC.

            if (isav == 0) goto L_190;
            if (isav < iu)
            {
                temp = (rdnext - rdprev) / (redmax + redmax - rdprev - rdnext);
                angt = angbd * (isav + HALF * temp) / iu;
            }
            var cth = (ONE - angt * angt) / (ONE + angt * angt);
            sth = (angt + angt) / (ONE + angt * angt);
            temp = shs + angt * (angt * dhd - dhs - dhs);
            sdec = sth * (angt * dredg - sredg - HALF * sth * temp);
            if (sdec <= ZERO) goto L_190;

            //     Update GNEW, D and HRED. If the angle of the alternative iteration
            //     is restricted by a bound on a free variable, that variable is fixed
            //     at the bound.

            dredg = ZERO;
            gredsq = ZERO;
            for (var i = 1; i <= n; ++i)
            {
                gnew[i] = gnew[i] + (cth - ONE) * hred[i] + sth * hs[i];
                if (xbdi[i] == ZERO)
                {
                    d[i] = cth * d[i] + sth * s[i];
                    dredg += d[i] * gnew[i];
                    gredsq += gnew[i] * gnew[i];
                }
                hred[i] = cth * hred[i] + sth * hs[i];
            }
            qred += sdec;
            if (iact > 0 && isav == iu)
            {
                ++nact;
                xbdi[iact] = xsav;
                goto L_100;
            }

            //     If SDEC is sufficiently small, then RETURN after setting XNEW to
            //     XOPT+D, giving careful attention to the bounds.

            if (sdec > 0.01 * qred) goto L_120;

            L_190:
            dsq = ZERO;
            for (var i = 1; i <= n; ++i)
            {
                xnew[i] = Math.Max(Math.Min(xopt[i] + d[i], su[i]), sl[i]);
                if (xbdi[i] == ONEMIN) xnew[i] = sl[i];
                if (xbdi[i] == ONE) xnew[i] = su[i];
                d[i] = xnew[i] - xopt[i];
                dsq += d[i] * d[i];
            }
            return;

            //     The following instructions multiply the current S-vector by the second
            //     derivative matrix of the quadratic model, putting the product in HS.
            //     They are reached from three different parts of the software above and
            //     they can be regarded as an external subroutine.

            L_210:
            var ih = 0;
            for (var j = 1; j <= n; ++j)
            {
                hs[j] = ZERO;
                for (var i = 1; i <= j; ++i)
                {
                    ++ih;
                    if (i < j) hs[j] += hq[ih] * s[i];
                    hs[i] += hq[ih] * s[j];
                }
            }
            for (var k = 1; k <= npt; ++k)
            {
                if (pq[k] != ZERO)
                {
                    temp = ZERO;
                    for (var j = 1; j <= n; ++j) temp += xpt[k, j] * s[j];
                    temp *= pq[k];
                    for (var i = 1; i <= n; ++i) hs[i] += temp * xpt[k, i];
                }
            }
            if (crvmin != ZERO) goto L_50;
            if (iterc > itcsav) goto L_150;

            for (var i = 1; i <= n; ++i) hred[i] = hs[i];
            goto L_120;
        }

        private static void UPDATE(int n, int npt, double[,] bmat, double[,] zmat, int ndim, double[] vlag, double beta, double denom, int knew, double[] w)
        {
            //     The arrays BMAT and ZMAT are updated, as required by the new position
            //     of the interpolation point that has the index KNEW. The vector VLAG has
            //     N+NPT components, set on entry to the first NPT and last N components
            //     of the product Hw in equation (4.11) of the Powell (2006) paper on
            //     NEWUOA. Further, BETA is set on entry to the value of the parameter
            //     with that name, and DENOM is set to the denominator of the updating
            //     formula. Elements of ZMAT may be treated as zero if their moduli are
            //     at most ZTEST.
            //     The first NDIM elements of W are used for working space.

            //     Set some constants.

            var nptm = npt - n - 1;
            var ztest = ZERO;
            for (var k = 1; k <= npt; ++k)
                for (var j = 1; j <= nptm; ++j)
                    ztest = Math.Max(ztest, Math.Abs(zmat[k, j]));
            ztest = 1.0E-20 * ztest;

            double temp, tempa, tempb;

            //     Apply the rotations that put zeros in the KNEW-th row of ZMAT.

            for (var j = 2; j <= nptm; ++j)
            {
                if (Math.Abs(zmat[knew, j]) > ztest)
                {
                    temp = Math.Sqrt(zmat[knew, 1] * zmat[knew, 1] + zmat[knew, j] * zmat[knew, j]);
                    tempa = zmat[knew, 1] / temp;
                    tempb = zmat[knew, j] / temp;
                    for (var i = 1; i <= npt; ++i)
                    {
                        temp = tempa * zmat[i, 1] + tempb * zmat[i, j];
                        zmat[i, j] = tempa * zmat[i, j] - tempb * zmat[i, 1];
                        zmat[i, 1] = temp;
                    }
                }
                zmat[knew, j] = ZERO;
            }

            //     Put the first NPT components of the KNEW-th column of HLAG into W,
            //     and calculate the parameters of the updating formula.

            for (var i = 1; i <= npt; ++i) w[i] = zmat[knew, 1] * zmat[i, 1];
            var alpha = w[knew];
            var tau = vlag[knew];
            vlag[knew] -= ONE;

            //     Complete the updating of ZMAT.

            temp = Math.Sqrt(denom);
            tempb = zmat[knew, 1] / temp;
            tempa = tau / temp;
            for (var i = 1; i <= npt; ++i) zmat[i, 1] = tempa * zmat[i, 1] - tempb * vlag[i];

            //     Finally, update the matrix BMAT.

            for (var j = 1; j <= n; ++j)
            {
                var jp = npt + j;
                w[jp] = bmat[knew, j];
                tempa = (alpha * vlag[jp] - tau * w[jp]) / denom;
                tempb = (-beta * w[jp] - tau * vlag[jp]) / denom;
                for (var i = 1; i <= jp; ++i)
                {
                    bmat[i, j] = bmat[i, j] + tempa * vlag[i] + tempb * w[i];
                    if (i > npt) bmat[jp, i - npt] = bmat[i, j];
                }
            }
        }

        #endregion

        #region PRIVATE SUPPORT METHODS

        private static string ToString(double[] x, int n)
        {
            var xstr = new string[n];
            for (var i = 0; i < n; ++i) xstr[i] = String.Format("{0,13:F6}", x[1 + i]);
            return String.Concat(xstr);
        }

        #endregion
    }

    // ReSharper restore InconsistentNaming
}
