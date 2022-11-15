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
 *  The original Fortran 77 version of this code was by Michael Powell (M.J.D.Powell @ damtp.cam.ac.uk)
 *  The Fortran 90 version was by Alan Miller (Alan.Miller @ vic.cmis.csiro.au). Latest revision - 30 October 1998
 */

using System;
using System.Collections.Generic;
using System.IO;

namespace Cureos.Numerics.Optimizers
{
    // ReSharper disable InconsistentNaming

    #region DELEGATES

    /// <summary>
    /// Signature for the objective and constraints function evaluation method used in <see cref="Cobyla"/> minimization.
    /// </summary>
    /// <param name="n">Number of variables.</param>
    /// <param name="m">Number of constraints.</param>
    /// <param name="x">Variable values to be employed in function and constraints calculation.</param>
    /// <param name="f">Calculated objective function value.</param>
    /// <param name="con">Calculated function values of the constraints.</param>
    public delegate void CalcfcDelegate(int n, int m, double[] x, out double f, double[] con);

    #endregion

    /// <summary>
    /// Constrained Optimization BY Linear Approximation for .NET
    ///
    /// COBYLA2 is an implementation of Powell�s nonlinear derivative�free constrained optimization that uses 
    /// a linear approximation approach. The algorithm is a sequential trust�region algorithm that employs linear 
    /// approximations to the objective and constraint functions, where the approximations are formed by linear 
    /// interpolation at n + 1 points in the space of the variables and tries to maintain a regular�shaped simplex 
    /// over iterations.
    ///
    /// It solves nonsmooth NLP with a moderate number of variables (about 100). Inequality constraints only.
    ///
    /// The initial point X is taken as one vertex of the initial simplex with zero being another, so, X should 
    /// not be entered as the zero vector.
    /// </summary>
    public class Cobyla : ITrustRegionOptimizer
    {
        #region FIELDS

        private readonly int _n;
        private readonly int _m;
        private readonly CalcfcDelegate _calcfc;

        private int _maxfun;
        private int _iprint;

        private double _rhobeg;
        private double _rhoend;

        private TextWriter _logger;

        private static readonly string LF = Environment.NewLine;
        private static readonly string IterationResultFormatter = LF + "NFVALS = {0,5}   F = {1,13:E6}    MAXCV = {2,13:E6}";

        #endregion

        /// <summary>
        /// Minimizes the objective function F with respect to a set of inequality constraints CON, 
        /// and returns the optimal variable array. F and CON may be non-linear, and should preferably be smooth.
        /// </summary>
        /// <param name="n">Number of variables.</param>
        /// <param name="m">Number of constraints.</param>
        /// <param name="calcfc">Method for calculating objective function and constraints.</param>
        public Cobyla(int n, int m, CalcfcDelegate calcfc)
        {
            _n = n;
            _m = m;
            _calcfc = calcfc;

            _rhobeg = 0.5;
            _rhoend = 1.0e-6;
            _iprint = 1;
            _maxfun = 10000;
            _logger = null;
        }

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
        /// Gets or sets the final value of the trust region radius.
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

        /// <summary>
        /// Gets or sets the start value of the trust region radius.
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

        #endregion

        #region METHODS

        /// <summary>
        /// Find a local minimum of provided objective function satisfying the nonlinear constraints.
        /// </summary>
        /// <param name="x0">Initial variable array.</param>
        /// <returns>Summary of the optimization result.</returns>
        public OptimizationSummary FindMinimum(double[] x0)
        {
            var nf = _maxfun;
            var status = COBYLA(_calcfc, _n, _m, x0, _rhobeg, _rhoend, _iprint, ref nf, _logger);

            double f;
            var con = new double[_m];
            _calcfc(_n, _m, x0, out f, con);

            return new OptimizationSummary(status, nf, x0, f, con);
        }

        private static OptimizationStatus COBYLA(
            CalcfcDelegate calcfc,
            int n,
            int m,
            double[] x,
            double rhobeg,
            double rhoend,
            int iprint,
            ref int iters,
            TextWriter logger)
        {
            //     This subroutine minimizes an objective function F(X) subject to M
            //     inequality constraints on X, where X is a vector of variables that has
            //     N components.  The algorithm employs linear approximations to the
            //     objective and constraint functions, the approximations being formed by
            //     linear interpolation at N+1 points in the space of the variables.
            //     We regard these interpolation points as vertices of a simplex.  The
            //     parameter RHO controls the size of the simplex and it is reduced
            //     automatically from RHOBEG to RHOEND.  For each RHO the subroutine tries
            //     to achieve a good vector of variables for the current size, and then
            //     RHO is reduced until the value RHOEND is reached.  Therefore RHOBEG and
            //     RHOEND should be set to reasonable initial changes to and the required
            //     accuracy in the variables respectively, but this accuracy should be
            //     viewed as a subject for experimentation because it is not guaranteed.
            //     The subroutine has an advantage over many of its competitors, however,
            //     which is that it treats each constraint individually when calculating
            //     a change to the variables, instead of lumping the constraints together
            //     into a single penalty function.  The name of the subroutine is derived
            //     from the phrase Constrained Optimization BY Linear Approximations.

            //     The user must set the values of N, M, RHOBEG and RHOEND, and must
            //     provide an initial vector of variables in X.  Further, the value of
            //     IPRINT should be set to 0, 1, 2 or 3, which controls the amount of
            //     printing during the calculation. Specifically, there is no output if
            //     IPRINT=0 and there is output only at the end of the calculation if
            //     IPRINT=1.  Otherwise each new value of RHO and SIGMA is printed.
            //     Further, the vector of variables and some function information are
            //     given either when RHO is reduced or when each new value of F(X) is
            //     computed in the cases IPRINT=2 or IPRINT=3 respectively. Here SIGMA
            //     is a penalty parameter, it being assumed that a change to X is an
            //     improvement if it reduces the merit function
            //                F(X)+SIGMA*MAX(0.0, - C1(X), - C2(X),..., - CM(X)),
            //     where C1,C2,...,CM denote the constraint functions that should become
            //     nonnegative eventually, at least to the precision of RHOEND. In the
            //     printed output the displayed term that is multiplied by SIGMA is
            //     called MAXCV, which stands for 'MAXimum Constraint Violation'.  The
            //     argument ITERS is an integer variable that must be set by the user to a
            //     limit on the number of calls of CALCFC, the purpose of this routine being
            //     given below.  The value of ITERS will be altered to the number of calls
            //     of CALCFC that are made.

            //     In order to define the objective and constraint functions, we require
            //     a subroutine that has the name and arguments
            //                SUBROUTINE CALCFC (N,M,X,F,CON)
            //                DIMENSION X(:),CON(:)  .
            //     The values of N and M are fixed and have been defined already, while
            //     X is now the current vector of variables. The subroutine should return
            //     the objective and constraint functions at X in F and CON(1),CON(2),
            //     ...,CON(M).  Note that we are trying to adjust X so that F(X) is as
            //     small as possible subject to the constraint functions being nonnegative.

            // Local variables
            var mpp = m + 2;

            // Internal base-1 X array
            var xinout = new double[n + 1];
            Array.Copy(x, 0, xinout, 1, n);

            // Internal representation of the objective and constraints calculation method, 
            // accounting for that X and CON arrays in the cobylb method are base-1 arrays.
            var fcalcfc = new CalcfcDelegate(
                (int nn, int mm, double[] xx, out double f, double[] con) =>
                    {
                        var ixx = new double[nn];
                        Array.Copy(xx, 1, ixx, 0, nn);
                        var ocon = new double[mm];
                        calcfc(nn, mm, ixx, out f, ocon);
                        Array.Copy(ocon, 0, con, 1, mm);
                    });

            var status = COBYLB(fcalcfc, n, m, mpp, xinout, rhobeg, rhoend, iprint, ref iters, logger);

            Array.Copy(xinout, 1, x, 0, n);

            return status;
        }

        private static OptimizationStatus COBYLB(
            CalcfcDelegate calcfc,
            int n,
            int m,
            int mpp,
            double[] x,
            double rhobeg,
            double rhoend,
            int iprint,
            ref int maxfun,
            TextWriter logger)
        {
            // N.B. Arguments CON, SIM, SIMI, DATMAT, A, VSIG, VETA, SIGBAR, DX, W & IACT
            //      have been removed.

            //     Set the initial values of some parameters. The last column of SIM holds
            //     the optimal vertex of the current simplex, and the preceding N columns
            //     hold the displacements from the optimal vertex to the other vertices.
            //     Further, SIMI holds the inverse of the matrix that is contained in the
            //     first N columns of SIM.

            // Local variables

            const double alpha = 0.25;
            const double beta = 2.1;
            const double gamma = 0.5;
            const double delta = 1.1;

            double f, resmax, total;

            var np = n + 1;
            var mp = m + 1;
            var rho = rhobeg;
            var parmu = 0.0;

            var iflag = false;
            var ifull = false;
            var parsig = 0.0;
            var prerec = 0.0;
            var prerem = 0.0;

            var con = new double[1 + mpp];
            var sim = new double[1 + n, 1 + np];
            var simi = new double[1 + n, 1 + n];
            var datmat = new double[1 + mpp, 1 + np];
            var a = new double[1 + n, 1 + mp];
            var vsig = new double[1 + n];
            var veta = new double[1 + n];
            var sigbar = new double[1 + n];
            var dx = new double[1 + n];
            var w = new double[1 + n];

            if (iprint >= 2 && logger != null) logger.WriteLine(LF + "The initial value of RHO is {0,13:F6} and PARMU is set to zero.", rho);

            var nfvals = 0;
            var temp = 1.0 / rho;

            for (var i = 1; i <= n; ++i)
            {
                sim[i, np] = x[i];
                sim[i, i] = rho;
                simi[i, i] = temp;
            }

            var jdrop = np;
            var ibrnch = false;

            OptimizationStatus status;

            //     Make the next call of the user-supplied subroutine CALCFC. These
            //     instructions are also used for calling CALCFC during the iterations of
            //     the algorithm.

            L_40:
            if (nfvals >= maxfun && nfvals > 0)
            {
                if (iprint >= 1 && logger != null) logger.WriteLine(LF + "Return from subroutine COBYLA because the MAXFUN limit has been reached.");
                status = OptimizationStatus.MAXFUN_Reached;
                goto L_600;
            }

            ++nfvals;

            calcfc(n, m, x, out f, con);
            resmax = 0.0;
            for (var k = 1; k <= m; ++k) resmax = Math.Max(resmax, -con[k]);

            if ((nfvals == iprint - 1 || iprint == 3) && logger != null)
            {
                logger.WriteLine(IterationResultFormatter, nfvals, f, resmax);
                logger.WriteLine("X = {0}", x.PART(1, n).FORMAT());
            }

            con[mp] = f;
            con[mpp] = resmax;
            if (ibrnch) goto L_440;

            //     Set the recently calculated function values in a column of DATMAT. This
            //     array has a column for each vertex of the current simplex, the entries of
            //     each column being the values of the constraint functions (if any)
            //     followed by the objective function and the greatest constraint violation
            //     at the vertex.

            for (var i = 1; i <= mpp; ++i) datmat[i, jdrop] = con[i];

            if (nfvals <= np)
            {

                //     Exchange the new vertex of the initial simplex with the optimal vertex if
                //     necessary. Then, if the initial simplex is not complete, pick its next
                //     vertex and calculate the function values there.

                if (jdrop <= n)
                {
                    if (datmat[mp, np] <= f)
                    {
                        x[jdrop] = sim[jdrop, np];
                    }
                    else
                    {
                        sim[jdrop, np] = x[jdrop];
                        for (var k = 1; k <= mpp; ++k)
                        {
                            datmat[k, jdrop] = datmat[k, np];
                            datmat[k, np] = con[k];
                        }
                        for (var k = 1; k <= jdrop; ++k)
                        {
                            sim[jdrop, k] = -rho;
                            temp = 0.0;
                            for (var i = k; i <= jdrop; ++i) temp -= simi[i, k];
                            simi[jdrop, k] = temp;
                        }
                    }
                }
                if (nfvals <= n)
                {
                    jdrop = nfvals;
                    x[jdrop] += rho;
                    goto L_40;
                }
            }

            ibrnch = true;

            //     Identify the optimal vertex of the current simplex.

            L_140:
            var phimin = datmat[mp, np] + parmu * datmat[mpp, np];
            var nbest = np;

            for (var j = 1; j <= n; ++j)
            {
                temp = datmat[mp, j] + parmu * datmat[mpp, j];
                if (temp < phimin)
                {
                    nbest = j;
                    phimin = temp;
                }
                else if (temp == phimin && parmu == 0.0 && datmat[mpp, j] < datmat[mpp, nbest])
                {
                    nbest = j;
                }
            }

            //     Switch the best vertex into pole position if it is not there already,
            //     and also update SIM, SIMI and DATMAT.

            if (nbest <= n)
            {
                for (var i = 1; i <= mpp; ++i)
                {
                    temp = datmat[i, np];
                    datmat[i, np] = datmat[i, nbest];
                    datmat[i, nbest] = temp;
                }
                for (var i = 1; i <= n; ++i)
                {
                    temp = sim[i, nbest];
                    sim[i, nbest] = 0.0;
                    sim[i, np] += temp;

                    var tempa = 0.0;
                    for (var k = 1; k <= n; ++k)
                    {
                        sim[i, k] -= temp;
                        tempa -= simi[k, i];
                    }
                    simi[nbest, i] = tempa;
                }
            }

            //     Make an error return if SIGI is a poor approximation to the inverse of
            //     the leading N by N submatrix of SIG.

            var error = 0.0;
            for (var i = 1; i <= n; ++i)
            {
                for (var j = 1; j <= n; ++j)
                {
                    temp = DOT_PRODUCT(simi.ROW(i).PART(1, n), sim.COL(j).PART(1, n)) - (i == j ? 1.0 : 0.0);
                    error = Math.Max(error, Math.Abs(temp));
                }
            }
            if (error > 0.1)
            {
                if (iprint >= 1 && logger != null)
                    logger.WriteLine(
                        LF + "Return from subroutine COBYLA because rounding errors are becoming damaging.");
                status = OptimizationStatus.X_RoundingErrorsPreventUpdate;
                goto L_600;
            }

            //     Calculate the coefficients of the linear approximations to the objective
            //     and constraint functions, placing minus the objective function gradient
            //     after the constraint gradients in the array A. The vector W is used for
            //     working space.

            for (var k = 1; k <= mp; ++k)
            {
                con[k] = -datmat[k, np];
                for (var j = 1; j <= n; ++j) w[j] = datmat[k, j] + con[k];

                for (var i = 1; i <= n; ++i)
                {
                    a[i, k] = (k == mp ? -1.0 : 1.0) * DOT_PRODUCT(w.PART(1, n), simi.COL(i).PART(1, n));
                }
            }

            //     Calculate the values of sigma and eta, and set IFLAG = 0 if the current
            //     simplex is not acceptable.

            iflag = true;
            parsig = alpha * rho;
            var pareta = beta * rho;

            for (var j = 1; j <= n; ++j)
            {
                var wsig = 0.0;
                for (var k = 1; k <= n; ++k) wsig += simi[j, k] * simi[j, k];
                var weta = 0.0;
                for (var k = 1; k <= n; ++k) weta += sim[k, j] * sim[k, j];
                vsig[j] = 1.0 / Math.Sqrt(wsig);
                veta[j] = Math.Sqrt(weta);
                if (vsig[j] < parsig || veta[j] > pareta) iflag = false;
            }

            //     If a new vertex is needed to improve acceptability, then decide which
            //     vertex to drop from the simplex.

            if (!ibrnch && !iflag)
            {
                jdrop = 0;
                temp = pareta;
                for (var j = 1; j <= n; ++j)
                {
                    if (veta[j] > temp)
                    {
                        jdrop = j;
                        temp = veta[j];
                    }
                }
                if (jdrop == 0)
                {
                    for (var j = 1; j <= n; ++j)
                    {
                        if (vsig[j] < temp)
                        {
                            jdrop = j;
                            temp = vsig[j];
                        }
                    }
                }

                //     Calculate the step to the new vertex and its sign.

                temp = gamma * rho * vsig[jdrop];
                for (var k = 1; k <= n; ++k) dx[k] = temp * simi[jdrop, k];
                var cvmaxp = 0.0;
                var cvmaxm = 0.0;

                total = 0.0;
                for (var k = 1; k <= mp; ++k)
                {
                    total = DOT_PRODUCT(a.COL(k).PART(1, n), dx.PART(1, n));
                    if (k < mp)
                    {
                        temp = datmat[k, np];
                        cvmaxp = Math.Max(cvmaxp, -total - temp);
                        cvmaxm = Math.Max(cvmaxm, total - temp);
                    }
                }
                var dxsign = parmu * (cvmaxp - cvmaxm) > 2.0 * total ? -1.0 : 1.0;

                //     Update the elements of SIM and SIMI, and set the next X.

                temp = 0.0;
                for (var i = 1; i <= n; ++i)
                {
                    dx[i] = dxsign * dx[i];
                    sim[i, jdrop] = dx[i];
                    temp += simi[jdrop, i] * dx[i];
                }
                for (var k = 1; k <= n; ++k) simi[jdrop, k] /= temp;

                for (var j = 1; j <= n; ++j)
                {
                    if (j != jdrop)
                    {
                        temp = DOT_PRODUCT(simi.ROW(j).PART(1, n), dx.PART(1, n));
                        for (var k = 1; k <= n; ++k) simi[j, k] -= temp * simi[jdrop, k];
                    }
                    x[j] = sim[j, np] + dx[j];
                }
                goto L_40;
            }

            //     Calculate DX = x(*)-x(0).
            //     Branch if the length of DX is less than 0.5*RHO.

            TRSTLP(n, m, a, con, rho, dx, out ifull);
            if (!ifull)
            {
                temp = 0.0;
                for (var k = 1; k <= n; ++k) temp += dx[k] * dx[k];
                if (temp < 0.25 * rho * rho)
                {
                    ibrnch = true;
                    goto L_550;
                }
            }

            //     Predict the change to F and the new maximum constraint violation if the
            //     variables are altered from x(0) to x(0) + DX.

            total = 0.0;
            var resnew = 0.0;
            con[mp] = 0.0;
            for (var k = 1; k <= mp; ++k)
            {
                total = con[k] - DOT_PRODUCT(a.COL(k).PART(1, n), dx.PART(1, n));
                if (k < mp) resnew = Math.Max(resnew, total);
            }

            //     Increase PARMU if necessary and branch back if this change alters the
            //     optimal vertex. Otherwise PREREM and PREREC will be set to the predicted
            //     reductions in the merit function and the maximum constraint violation
            //     respectively.

            prerec = datmat[mpp, np] - resnew;
            var barmu = prerec > 0.0 ? total / prerec : 0.0;
            if (parmu < 1.5 * barmu)
            {
                parmu = 2.0 * barmu;
                if (iprint >= 2 && logger != null) logger.WriteLine(LF + "Increase in PARMU to {0,13:F6}", parmu);
                var phi = datmat[mp, np] + parmu * datmat[mpp, np];
                for (var j = 1; j <= n; ++j)
                {
                    temp = datmat[mp, j] + parmu * datmat[mpp, j];
                    if (temp < phi || (temp == phi && parmu == 0.0 && datmat[mpp, j] < datmat[mpp, np])) goto L_140;
                }
            }
            prerem = parmu * prerec - total;

            //     Calculate the constraint and objective functions at x(*).
            //     Then find the actual reduction in the merit function.

            for (var k = 1; k <= n; ++k) x[k] = sim[k, np] + dx[k];
            ibrnch = true;
            goto L_40;

            L_440:
            var vmold = datmat[mp, np] + parmu * datmat[mpp, np];
            var vmnew = f + parmu * resmax;
            var trured = vmold - vmnew;
            if (parmu == 0.0 && f == datmat[mp, np])
            {
                prerem = prerec;
                trured = datmat[mpp, np] - resmax;
            }

            //     Begin the operations that decide whether x(*) should replace one of the
            //     vertices of the current simplex, the change being mandatory if TRURED is
            //     positive. Firstly, JDROP is set to the index of the vertex that is to be
            //     replaced.

            var ratio = trured <= 0.0 ? 1.0 : 0.0;
            jdrop = 0;
            for (var j = 1; j <= n; ++j)
            {
                temp = Math.Abs(DOT_PRODUCT(simi.ROW(j).PART(1, n), dx.PART(1, n)));
                if (temp > ratio)
                {
                    jdrop = j;
                    ratio = temp;
                }
                sigbar[j] = temp * vsig[j];
            }

            //     Calculate the value of ell.

            var edgmax = delta * rho;
            var l = 0;
            for (var j = 1; j <= n; ++j)
            {
                if (sigbar[j] >= parsig || sigbar[j] >= vsig[j])
                {
                    temp = veta[j];
                    if (trured > 0.0)
                    {
                        temp = 0.0;
                        for (var k = 1; k <= n; ++k) temp += Math.Pow(dx[k] - sim[k, j], 2.0);
                        temp = Math.Sqrt(temp);
                    }
                    if (temp > edgmax)
                    {
                        l = j;
                        edgmax = temp;
                    }
                }
            }
            if (l > 0) jdrop = l;
            if (jdrop == 0) goto L_550;

            //     Revise the simplex by updating the elements of SIM, SIMI and DATMAT.

            temp = 0.0;
            for (var i = 1; i <= n; ++i)
            {
                sim[i, jdrop] = dx[i];
                temp += simi[jdrop, i] * dx[i];
            }
            for (var k = 1; k <= n; ++k) simi[jdrop, k] /= temp;
            for (var j = 1; j <= n; ++j)
            {
                if (j != jdrop)
                {
                    temp = DOT_PRODUCT(simi.ROW(j).PART(1, n), dx.PART(1, n));
                    for (var k = 1; k <= n; ++k) simi[j, k] -= temp * simi[jdrop, k];
                }
            }
            for (var k = 1; k <= mpp; ++k) datmat[k, jdrop] = con[k];

            //     Branch back for further iterations with the current RHO.

            if (trured > 0.0 && trured >= 0.1 * prerem) goto L_140;

            L_550:
            if (!iflag)
            {
                ibrnch = false;
                goto L_140;
            }

            //     Otherwise reduce RHO if it is not at its least value and reset PARMU.

            if (rho > rhoend)
            {
                double cmin = 0.0, cmax = 0.0;

                rho *= 0.5;
                if (rho <= 1.5 * rhoend) rho = rhoend;
                if (parmu > 0.0)
                {
                    var denom = 0.0;
                    for (var k = 1; k <= mp; ++k)
                    {
                        cmin = datmat[k, np];
                        cmax = cmin;
                        for (var i = 1; i <= n; ++i)
                        {
                            cmin = Math.Min(cmin, datmat[k, i]);
                            cmax = Math.Max(cmax, datmat[k, i]);
                        }
                        if (k <= m && cmin < 0.5 * cmax)
                        {
                            temp = Math.Max(cmax, 0.0) - cmin;
                            denom = denom <= 0.0 ? temp : Math.Min(denom, temp);
                        }
                    }
                    if (denom == 0.0)
                    {
                        parmu = 0.0;
                    }
                    else if (cmax - cmin < parmu * denom)
                    {
                        parmu = (cmax - cmin) / denom;
                    }
                }
                if (logger != null)
                {
                    if (iprint >= 2) logger.WriteLine(LF + "Reduction in RHO to {0,13:E6}  and PARMU = {1,13:E6}", rho, parmu);
                    if (iprint == 2)
                    {
                        logger.WriteLine(IterationResultFormatter, nfvals, datmat[mp, np], datmat[mpp, np]);
                        logger.WriteLine("X = {0}", sim.COL(np).PART(1, n).FORMAT());
                    }
                }
                goto L_140;
            }

            //     Return the best calculated values of the variables.

            status = OptimizationStatus.Normal;
            if (iprint >= 1 && logger != null) logger.WriteLine(LF + "Normal return from subroutine COBYLA");
            if (ifull) goto L_620;

            L_600:
            for (var k = 1; k <= n; ++k) x[k] = sim[k, np];
            f = datmat[mp, np];
            resmax = datmat[mpp, np];

            L_620:
            if (iprint >= 1 && logger != null)
            {
                logger.WriteLine(IterationResultFormatter, nfvals, f, resmax);
                logger.WriteLine("X = {0}", x.PART(1, n).FORMAT());
            }

            maxfun = nfvals;

            return status;
        }

        private static void TRSTLP(int n, int m, double[,] a, double[] b, double rho, double[] dx, out bool ifull)
        {
            // N.B. Arguments Z, ZDOTA, VMULTC, SDIRN, DXNEW, VMULTD & IACT have been removed.

            //     This subroutine calculates an N-component vector DX by applying the
            //     following two stages. In the first stage, DX is set to the shortest
            //     vector that minimizes the greatest violation of the constraints
            //       A(1,K)*DX(1)+A(2,K)*DX(2)+...+A(N,K)*DX(N) .GE. B(K), K = 2,3,...,M,
            //     subject to the Euclidean length of DX being at most RHO. If its length is
            //     strictly less than RHO, then we use the resultant freedom in DX to
            //     minimize the objective function
            //              -A(1,M+1)*DX(1) - A(2,M+1)*DX(2) - ... - A(N,M+1)*DX(N)
            //     subject to no increase in any greatest constraint violation. This
            //     notation allows the gradient of the objective function to be regarded as
            //     the gradient of a constraint. Therefore the two stages are distinguished
            //     by MCON .EQ. M and MCON .GT. M respectively. It is possible that a
            //     degeneracy may prevent DX from attaining the target length RHO. Then the
            //     value IFULL = 0 would be set, but usually IFULL = 1 on return.

            //     In general NACT is the number of constraints in the active set and
            //     IACT(1),...,IACT(NACT) are their indices, while the remainder of IACT
            //     contains a permutation of the remaining constraint indices.  Further, Z
            //     is an orthogonal matrix whose first NACT columns can be regarded as the
            //     result of Gram-Schmidt applied to the active constraint gradients.  For
            //     J = 1,2,...,NACT, the number ZDOTA(J) is the scalar product of the J-th
            //     column of Z with the gradient of the J-th active constraint.  DX is the
            //     current vector of variables and here the residuals of the active
            //     constraints should be zero. Further, the active constraints have
            //     nonnegative Lagrange multipliers that are held at the beginning of
            //     VMULTC. The remainder of this vector holds the residuals of the inactive
            //     constraints at DX, the ordering of the components of VMULTC being in
            //     agreement with the permutation of the indices of the constraints that is
            //     in IACT. All these residuals are nonnegative, which is achieved by the
            //     shift RESMAX that makes the least residual zero.

            //     Initialize Z and some other variables. The value of RESMAX will be
            //     appropriate to DX = 0, while ICON will be the index of a most violated
            //     constraint if RESMAX is positive. Usually during the first stage the
            //     vector SDIRN gives a search direction that reduces all the active
            //     constraint violations by one simultaneously.

            // Local variables

            ifull = true;

            double temp;

            var nactx = 0;
            var resold = 0.0;

            var z = new double[1 + n, 1 + n];
            var zdota = new double[2 + m];
            var vmultc = new double[2 + m];
            var sdirn = new double[1 + n];
            var dxnew = new double[1 + n];
            var vmultd = new double[2 + m];
            var iact = new int[2 + m];

            var mcon = m;
            var nact = 0;
            for (var i = 1; i <= n; ++i)
            {
                z[i, i] = 1.0;
                dx[i] = 0.0;
            }

            var icon = 0;
            var resmax = 0.0;
            if (m >= 1)
            {
                for (var k = 1; k <= m; ++k)
                {
                    if (b[k] > resmax)
                    {
                        resmax = b[k];
                        icon = k;
                    }
                }
                for (var k = 1; k <= m; ++k)
                {
                    iact[k] = k;
                    vmultc[k] = resmax - b[k];
                }
            }
            if (resmax == 0.0) goto L_480;

            //     End the current stage of the calculation if 3 consecutive iterations
            //     have either failed to reduce the best calculated value of the objective
            //     function or to increase the number of active constraints since the best
            //     value was calculated. This strategy prevents cycling, but there is a
            //     remote possibility that it will cause premature termination.

            L_60:
            var optold = 0.0;
            var icount = 0;

            L_70:
            var optnew = mcon == m ? resmax : -DOT_PRODUCT(dx.PART(1, n), a.COL(mcon).PART(1, n));

            if (icount == 0 || optnew < optold)
            {
                optold = optnew;
                nactx = nact;
                icount = 3;
            }
            else if (nact > nactx)
            {
                nactx = nact;
                icount = 3;
            }
            else
            {
                --icount;
            }
            if (icount == 0) goto L_490;

            //     If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to
            //     the active set. Apply Givens rotations so that the last N-NACT-1 columns
            //     of Z are orthogonal to the gradient of the new constraint, a scalar
            //     product being set to zero if its nonzero value could be due to computer
            //     rounding errors. The array DXNEW is used for working space.

            if (icon <= nact) goto L_260;
            var kk = iact[icon];
            for (var k = 1; k <= n; ++k) dxnew[k] = a[k, kk];
            var tot = 0.0;

            {
                var k = n;
                while (k > nact)
                {
                    var sp = 0.0;
                    var spabs = 0.0;
                    for (var i = 1; i <= n; ++i)
                    {
                        temp = z[i, k] * dxnew[i];
                        sp += temp;
                        spabs += Math.Abs(temp);
                    }
                    var acca = spabs + 0.1 * Math.Abs(sp);
                    var accb = spabs + 0.2 * Math.Abs(sp);
                    if (spabs >= acca || acca >= accb) sp = 0.0;
                    if (tot == 0.0)
                    {
                        tot = sp;
                    }
                    else
                    {
                        var kp = k + 1;
                        temp = Math.Sqrt(sp * sp + tot * tot);
                        var alpha = sp / temp;
                        var beta = tot / temp;
                        tot = temp;
                        for (var i = 1; i <= n; ++i)
                        {
                            temp = alpha * z[i, k] + beta * z[i, kp];
                            z[i, kp] = alpha * z[i, kp] - beta * z[i, k];
                            z[i, k] = temp;
                        }
                    }
                    --k;
                }
            }

            //     Add the new constraint if this can be done without a deletion from the
            //     active set.

            if (tot != 0.0)
            {
                ++nact;
                zdota[nact] = tot;
                vmultc[icon] = vmultc[nact];
                vmultc[nact] = 0.0;
                goto L_210;
            }

            //     The next instruction is reached if a deletion has to be made from the
            //     active set in order to make room for the new active constraint, because
            //     the new constraint gradient is a linear combination of the gradients of
            //     the old active constraints.  Set the elements of VMULTD to the multipliers
            //     of the linear combination.  Further, set IOUT to the index of the
            //     constraint to be deleted, but branch if no suitable index can be found.

            var ratio = -1.0;
            {
                var k = nact;
                do
                {
                    var zdotv = 0.0;
                    var zdvabs = 0.0;

                    for (var i = 1; i <= n; ++i)
                    {
                        temp = z[i, k] * dxnew[i];
                        zdotv = zdotv + temp;
                        zdvabs = zdvabs + Math.Abs(temp);
                    }
                    var acca = zdvabs + 0.1 * Math.Abs(zdotv);
                    var accb = zdvabs + 0.2 * Math.Abs(zdotv);
                    if (zdvabs < acca && acca < accb)
                    {
                        temp = zdotv / zdota[k];
                        if (temp > 0.0 && iact[k] <= m)
                        {
                            var tempa = vmultc[k] / temp;
                            if (ratio < 0.0 || tempa < ratio) ratio = tempa;
                        }

                        if (k >= 2)
                        {
                            var kw = iact[k];
                            for (var i = 1; i <= n; ++i) dxnew[i] -= temp * a[i, kw];
                        }
                        vmultd[k] = temp;
                    }
                    else
                    {
                        vmultd[k] = 0.0;
                    }
                }
                while (--k > 0);
            }
            if (ratio < 0.0) goto L_490;

            //     Revise the Lagrange multipliers and reorder the active constraints so
            //     that the one to be replaced is at the end of the list. Also calculate the
            //     new value of ZDOTA(NACT) and branch if it is not acceptable.

            for (var k = 1; k <= nact; ++k) vmultc[k] = Math.Max(0.0, vmultc[k] - ratio * vmultd[k]);
            if (icon < nact)
            {
                var isave = iact[icon];
                var vsave = vmultc[icon];
                var k = icon;
                do
                {
                    var kp = k + 1;
                    var kw = iact[kp];
                    var sp = DOT_PRODUCT(z.COL(k).PART(1, n), a.COL(kw).PART(1, n));
                    temp = Math.Sqrt(sp * sp + zdota[kp] * zdota[kp]);
                    var alpha = zdota[kp] / temp;
                    var beta = sp / temp;
                    zdota[kp] = alpha * zdota[k];
                    zdota[k] = temp;
                    for (var i = 1; i <= n; ++i)
                    {
                        temp = alpha * z[i, kp] + beta * z[i, k];
                        z[i, kp] = alpha * z[i, k] - beta * z[i, kp];
                        z[i, k] = temp;
                    }
                    iact[k] = kw;
                    vmultc[k] = vmultc[kp];
                    k = kp;
                }
                while (k < nact);
                iact[k] = isave;
                vmultc[k] = vsave;
            }
            temp = DOT_PRODUCT(z.COL(nact).PART(1, n), a.COL(kk).PART(1, n));
            if (temp == 0.0) goto L_490;
            zdota[nact] = temp;
            vmultc[icon] = 0.0;
            vmultc[nact] = ratio;

            //     Update IACT and ensure that the objective function continues to be
            //     treated as the last active constraint when MCON>M.

            L_210:
            iact[icon] = iact[nact];
            iact[nact] = kk;
            if (mcon > m && kk != mcon)
            {
                var k = nact - 1;
                var sp = DOT_PRODUCT(z.COL(k).PART(1, n), a.COL(kk).PART(1, n));
                temp = Math.Sqrt(sp * sp + zdota[nact] * zdota[nact]);
                var alpha = zdota[nact] / temp;
                var beta = sp / temp;
                zdota[nact] = alpha * zdota[k];
                zdota[k] = temp;
                for (var i = 1; i <= n; ++i)
                {
                    temp = alpha * z[i, nact] + beta * z[i, k];
                    z[i, nact] = alpha * z[i, k] - beta * z[i, nact];
                    z[i, k] = temp;
                }
                iact[nact] = iact[k];
                iact[k] = kk;
                temp = vmultc[k];
                vmultc[k] = vmultc[nact];
                vmultc[nact] = temp;
            }

            //     If stage one is in progress, then set SDIRN to the direction of the next
            //     change to the current vector of variables.

            if (mcon > m) goto L_320;
            kk = iact[nact];
            temp = (DOT_PRODUCT(sdirn.PART(1, n), a.COL(kk).PART(1, n)) - 1.0) / zdota[nact];
            for (var k = 1; k <= n; ++k) sdirn[k] -= temp * z[k, nact];
            goto L_340;

            //     Delete the constraint that has the index IACT(ICON) from the active set.

            L_260:
            if (icon < nact)
            {
                var isave = iact[icon];
                var vsave = vmultc[icon];
                var k = icon;
                do
                {
                    var kp = k + 1;
                    kk = iact[kp];
                    var sp = DOT_PRODUCT(z.COL(k).PART(1, n), a.COL(kk).PART(1, n));
                    temp = Math.Sqrt(sp * sp + zdota[kp] * zdota[kp]);
                    var alpha = zdota[kp] / temp;
                    var beta = sp / temp;
                    zdota[kp] = alpha * zdota[k];
                    zdota[k] = temp;
                    for (var i = 1; i <= n; ++i)
                    {
                        temp = alpha * z[i, kp] + beta * z[i, k];
                        z[i, kp] = alpha * z[i, k] - beta * z[i, kp];
                        z[i, k] = temp;
                    }
                    iact[k] = kk;
                    vmultc[k] = vmultc[kp];
                    k = kp;
                }
                while (k < nact);

                iact[k] = isave;
                vmultc[k] = vsave;
            }
            --nact;

            //     If stage one is in progress, then set SDIRN to the direction of the next
            //     change to the current vector of variables.

            if (mcon > m) goto L_320;
            temp = DOT_PRODUCT(sdirn.PART(1, n), z.COL(nact + 1).PART(1, n));
            for (var k = 1; k <= n; ++k) sdirn[k] -= temp * z[k, nact + 1];
            goto L_340;

            //     Pick the next search direction of stage two.

            L_320:
            temp = 1.0 / zdota[nact];
            for (var k = 1; k <= n; ++k) sdirn[k] = temp * z[k, nact];

            //     Calculate the step to the boundary of the trust region or take the step
            //     that reduces RESMAX to zero. The two statements below that include the
            //     factor 1.0E-6 prevent some harmless underflows that occurred in a test
            //     calculation. Further, we skip the step if it could be zero within a
            //     reasonable tolerance for computer rounding errors.

            L_340:
            var dd = rho * rho;
            var sd = 0.0;
            var ss = 0.0;
            for (var i = 1; i <= n; ++i)
            {
                if (Math.Abs(dx[i]) >= 1.0E-6 * rho) dd -= dx[i] * dx[i];
                sd += dx[i] * sdirn[i];
                ss += sdirn[i] * sdirn[i];
            }
            if (dd <= 0.0) goto L_490;
            temp = Math.Sqrt(ss * dd);
            if (Math.Abs(sd) >= 1.0E-6 * temp) temp = Math.Sqrt(ss * dd + sd * sd);
            var stpful = dd / (temp + sd);
            var step = stpful;
            if (mcon == m)
            {
                var acca = step + 0.1 * resmax;
                var accb = step + 0.2 * resmax;
                if (step >= acca || acca >= accb) goto L_480;
                step = Math.Min(step, resmax);
            }

            //     Set DXNEW to the new variables if STEP is the steplength, and reduce
            //     RESMAX to the corresponding maximum residual if stage one is being done.
            //     Because DXNEW will be changed during the calculation of some Lagrange
            //     multipliers, it will be restored to the following value later.

            for (var k = 1; k <= n; ++k) dxnew[k] = dx[k] + step * sdirn[k];
            if (mcon == m)
            {
                resold = resmax;
                resmax = 0.0;
                for (var k = 1; k <= nact; ++k)
                {
                    kk = iact[k];
                    temp = b[kk] - DOT_PRODUCT(a.COL(kk).PART(1, n), dxnew.PART(1, n));
                    resmax = Math.Max(resmax, temp);
                }
            }

            //     Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A
            //     device is included to force VMULTD(K) = 0.0 if deviations from this value
            //     can be attributed to computer rounding errors. First calculate the new
            //     Lagrange multipliers.

            {
                var k = nact;
                L_390:
                var zdotw = 0.0;
                var zdwabs = 0.0;
                for (var i = 1; i <= n; ++i)
                {
                    temp = z[i, k] * dxnew[i];
                    zdotw += temp;
                    zdwabs += Math.Abs(temp);
                }
                var acca = zdwabs + 0.1 * Math.Abs(zdotw);
                var accb = zdwabs + 0.2 * Math.Abs(zdotw);
                if (zdwabs >= acca || acca >= accb) zdotw = 0.0;
                vmultd[k] = zdotw / zdota[k];
                if (k >= 2)
                {
                    kk = iact[k];
                    for (var i = 1; i <= n; ++i) dxnew[i] -= vmultd[k] * a[i, kk];
                    --k;
                    goto L_390;
                }
                if (mcon > m) vmultd[nact] = Math.Max(0.0, vmultd[nact]);
            }

            //     Complete VMULTC by finding the new constraint residuals.

            for (var k = 1; k <= n; ++k) dxnew[k] = dx[k] + step * sdirn[k];
            if (mcon > nact)
            {
                var kl = nact + 1;
                for (var k = kl; k <= mcon; ++k)
                {
                    kk = iact[k];
                    var total = resmax - b[kk];
                    var sumabs = resmax + Math.Abs(b[kk]);
                    for (var i = 1; i <= n; ++i)
                    {
                        temp = a[i, kk] * dxnew[i];
                        total += temp;
                        sumabs += Math.Abs(temp);
                    }
                    var acca = sumabs + 0.1 * Math.Abs(total);
                    var accb = sumabs + 0.2 * Math.Abs(total);
                    if (sumabs >= acca || acca >= accb) total = 0.0;
                    vmultd[k] = total;
                }
            }

            //     Calculate the fraction of the step from DX to DXNEW that will be taken.

            ratio = 1.0;
            icon = 0;
            for (var k = 1; k <= mcon; ++k)
            {
                if (vmultd[k] < 0.0)
                {
                    temp = vmultc[k] / (vmultc[k] - vmultd[k]);
                    if (temp < ratio)
                    {
                        ratio = temp;
                        icon = k;
                    }
                }
            }

            //     Update DX, VMULTC and RESMAX.

            temp = 1.0 - ratio;
            for (var k = 1; k <= n; ++k) dx[k] = temp * dx[k] + ratio * dxnew[k];
            for (var k = 1; k <= mcon; ++k) vmultc[k] = Math.Max(0.0, temp * vmultc[k] + ratio * vmultd[k]);
            if (mcon == m) resmax = resold + ratio * (resmax - resold);

            //     If the full step is not acceptable then begin another iteration.
            //     Otherwise switch to stage two or end the calculation.

            if (icon > 0) goto L_70;
            if (step == stpful) return;

            L_480:
            mcon = m + 1;
            icon = mcon;
            iact[mcon] = mcon;
            vmultc[mcon] = 0.0;
            goto L_60;

            //     We employ any freedom that may be available to reduce the objective
            //     function before returning a DX whose length is less than RHO.

            L_490:
            if (mcon == m) goto L_480;
            ifull = false;
        }

        private static double DOT_PRODUCT(double[] lhs, double[] rhs)
        {
            var sum = 0.0;
            for (var i = 0; i < lhs.Length; ++i) sum += lhs[i] * rhs[i];
            return sum;
        }

        #endregion
    }

    #region EXTENSION METHODS CLASS

    internal static class CobylaExtensionMethods
    {
        internal static T[] ROW<T>(this T[,] src, int rowidx)
        {
            var cols = src.GetLength(1);
            var dest = new T[cols];
            for (var col = 0; col < cols; ++col) dest[col] = src[rowidx, col];
            return dest;
        }

        internal static T[] COL<T>(this T[,] src, int colidx)
        {
            var rows = src.GetLength(0);
            var dest = new T[rows];
            for (var row = 0; row < rows; ++row) dest[row] = src[row, colidx];
            return dest;
        }

        internal static T[] PART<T>(this IList<T> src, int from, int to)
        {
            var dest = new T[to - from + 1];
            var destidx = 0;
            for (var srcidx = from; srcidx <= to; ++srcidx, ++destidx) dest[destidx] = src[srcidx];
            return dest;
        }

        internal static string FORMAT(this double[] x)
        {
            var xStr = new string[x.Length];
            for (var i = 0; i < x.Length; ++i) xStr[i] = String.Format("{0,13:F6}", x[i]);
            return String.Concat(xStr);
        }
    }

    // ReSharper restore InconsistentNaming

    #endregion
}
