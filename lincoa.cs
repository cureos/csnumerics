using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Cureos.Numerics
{
    #region DELEGATES

    internal delegate void CalfunDel(int n, double[] x, ref double f);

    #endregion

    public class Lincoa
    {
        #region INNER TYPES

        public enum Status
        {
            InProgress,
            Success,
            TooFew_N,
            NPT_OutOfRange,
            MAXFUN_NotLargerThan_NPT,
            ConstraintGradientIsZero,
            MAXFUN_Reached,
            X_UpdateHamperedByRoundingErrors,
            UpdatingFormulaDenominatorZero
        }

        #endregion
        
        #region FIELDS

        private const double ZERO=0.0;
        private const double HALF=0.5;
        private const double ONE=1.0;
        private const double TENTH=0.1;
        private const double TINY=1.0E-60;
        private const double CTEST=0.01;

        private static readonly string LF = Environment.NewLine;

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
        private static readonly string LINCOB_590 = "Least value of F = {0,23:E15}         The corresponding X is:" + LF + "{2}";
        private static readonly string LINCOB_620 = "At the return from LINCOA     Number of function values = {0,6:D}";

        #endregion

        #region METHODS

        private static void PRINT(TextWriter logger, string format, params object[] args)
        {
            if (logger != null) logger.WriteLine(format, args);
        }

        private static string FORMAT<T>(string separator, string itemFormatter, IEnumerable<T> items)
        {
            return String.Join(separator, items.Select(item => String.Format("{{0," + itemFormatter + "}}", item)));
        }

        private static Status LINCOA(CalfunDel CALFUN, int N, int NPT, int M, double[,] A, int IA, double[] B, double[] X,
            double RHOBEG, double RHOEND, int IPRINT, int MAXFUN, TextWriter logger)
        {
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
            double SMALLX = 1.0E-6 * RHOEND;
            int NP = N + 1;
            int NPTM = NPT - NP;
            if (N <= 1)
            {
                PRINT(logger, LINCOA_10);
                return Status.TooFew_N;
            }
            if (NPT < N + 2 || NPT > ((N + 2) * NP) / 2)
            {
                PRINT(logger, LINCOA_20);
                return Status.NPT_OutOfRange;
            }
            if (MAXFUN <= NPT)
            {
                PRINT(logger, LINCOA_30);
                return Status.MAXFUN_NotLargerThan_NPT;
            }
//
//     Normalize the constraints, and copy the resultant constraint matrix
//       and right hand sides into working space, after increasing the right
//       hand sides if necessary so that the starting point is feasible.
//
            double[,] AMAT = new double[1+N,1+M];
            double[] BB = new double[1+M];    // B in LINCOB
            int IFLAG = 0;
            if (M > 0)
            {
                for (int J = 1; J <= M; ++J)
                {
                    double SUM = ZERO;
                    double TEMP = ZERO;
                    for (int I = 1; I <= N; ++I)
                    {
                        SUM += A[I, J] * X[I];
                        TEMP += A[I, J] * A[I, J];
                    }
                    if (TEMP == ZERO)
                    {
                        PRINT(logger, LINCOA_50);
                        return Status.ConstraintGradientIsZero;
                    }
                    TEMP = Math.Sqrt(TEMP);
                    if (SUM - B[J] > SMALLX * TEMP) IFLAG = 1;
                    BB[J] = Math.Max(B[J], SUM) / TEMP;
                    for (int I = 1; I <= N; ++I)
                    {
                        AMAT[I,J] = A[I, J] / TEMP;
                    }
                }
            }
            if (IFLAG == 1)
            {
                if (IPRINT > 0) PRINT(logger, LINCOA_70);
            }
//
//     Partition the working space array, so that different parts of it can be
//     treated separately by the subroutine that performs the main calculation.
//
            int IAMAT = Math.Max(M + 3 * N, Math.Max(2 * M + N, 2 * NPT)) + 1;
            int NDIM = NPT + N;

            double[] XBASE = new double[1+N];
            double[,] XPT = new double[1+NPT,1+N];
            double[] FVAL = new double[1+NPT];
            double[] XSAV = new double[1+N];
            double[] XOPT = new double[1+N];
            double[] GOPT = new double[1+N];
            double[] HQ = new double[1+(N * NP) / 2];
            double[] PQ = new double[1+NPT];
            double[,] BMAT = new double[1+NDIM,1+N];
            double[,] ZMAT = new double[1+NPT,1+NPTM];
            double[] STEP = new double[1+N];
            double[] SP = new double[1+ NPT + NPT];
            double[] XNEW = new double[1+N];
            int[] IACT = new int[1+N];
            double[] RESCON = new double[1+M];
            double[,] QFAC = new double[1+N,1+N];
            double[] RFAC = new double[1 + (N * NP) / 2];
            double[] PQW = new double[1+NPT+N];
            double[] W = new double[IAMAT];

//
//     The above settings provide a partition of W for subroutine LINCOB.
//
            return LINCOB(CALFUN, N, NPT, M, AMAT, BB, X, RHOBEG, RHOEND, IPRINT, MAXFUN, XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ,
                PQ, BMAT, ZMAT, NDIM, STEP, SP, XNEW, IACT, RESCON, QFAC, RFAC, PQW, W, logger);
        }

        private static Status LINCOB(CalfunDel CALFUN, int N, int NPT, int M, double[,] AMAT, double[] B, double[] X, double RHOBEG,
            double RHOEND, int IPRINT, int MAXFUN, double[] XBASE, double[,] XPT, double[] FVAL, double[] XSAV,
            double[] XOPT, double[] GOPT, double[] HQ, double[] PQ, double[,] BMAT, double[,] ZMAT, int NDIM,
            double[] STEP, double[] SP, double[] XNEW, int[] IACT, double[] RESCON, double[,] QFAC, double[] RFAC,
            double[] PQW, double[] W, TextWriter logger)
        {
            Status status = Status.InProgress;
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
            int NP = N + 1;
            int NH = (N * NP) / 2;
            int NPTM = NPT - NP;
//
//     Set the elements of XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT,
//       ZMAT and SP for the first iteration. An important feature is that,
//       if the interpolation point XPT(K,.) is not feasible, where K is any
//       integer from [1,NPT], then a change is made to XPT(K,.) if necessary
//       so that the constraint violation is at least 0.2*RHOBEG. Also KOPT
//       is set so that XPT(KOPT,.) is the initial trust region centre.
//
            int KOPT, IDZ;
            PRELIM(N, NPT, M, AMAT, B, X, RHOBEG, IPRINT, XBASE, XPT, FVAL, XSAV, XOPT, GOPT, out KOPT, HQ, PQ, BMAT,
                ZMAT, out IDZ, NDIM, SP, RESCON, STEP, PQW, W);
//
//     Begin the iterative procedure.
//
            int NF = NPT;
            double FOPT = FVAL[1 + KOPT];
            double RHO = RHOBEG;
            double DELTA = RHO;
            int IFEAS = 0;
            int NACT = 0;
            int ITEST = 3;

            int KNEW, NVALA, NVALB;
            double FSAVE, XOPTSQ;

            LINCOB_10:

            KNEW = 0;
            NVALA = 0;
            NVALB = 0;
//
//     Shift XBASE if XOPT may be too far from XBASE. First make the changes
//       to BMAT that do not depend on ZMAT.
//
            LINCOB_20:

            FSAVE = FOPT;
            XOPTSQ = ZERO;
            for (int I = 1; I <= N; ++I)
                XOPTSQ += XOPT[I] * XOPT[I];
            if (XOPTSQ >= 1.0E4 * DELTA * DELTA)
            {
                double QOPTSQ = 0.25 * XOPTSQ;
                for (int K = 1; K <= NPT; ++K)
                {
                    double SUM = ZERO;
                    for (int I = 1; I <= N; ++I)
                        SUM += XPT[K, I] * XOPT[I];
                    SUM -= HALF * XOPTSQ;
                    W[NPT + K] = SUM;
                    SP[K] = ZERO;
                    for (int I = 1; I <= N; ++I)
                    {
                        XPT[K, I] -= HALF * XOPT[I];
                        STEP[I] = BMAT[K, I];
                        W[I] = SUM * XPT[K, I] + QOPTSQ * XOPT[I];
                        int IP = NPT + I;
                        for (int J = 1; J <= I; ++J)
                            BMAT[IP, J] += STEP[I] * W[J] + W[I] * STEP[J];
                    }
                }
//
//     Then the revisions of BMAT that depend on ZMAT are calculated.
//
                for (int K = 1; K <= NPTM; ++K)
                {
                    double SUMZ = ZERO;
                    for (int I = 1; I <= NPT; ++I)
                    {
                        SUMZ += ZMAT[I, K];
                        W[I] = W[NPT + I] * ZMAT[I, K];
                    }
                    for (int J = 1; J <= N; ++J)
                    {
                        double SUM = QOPTSQ * SUMZ * XOPT[J];
                        for (int I = 1; I <= NPT; ++I)
                            SUM += W[I] * XPT[I, J];
                        STEP[J] = SUM;
                        if (K < IDZ) SUM = -SUM;
                        for (int I = 1; I <= NPT; ++I)
                            BMAT[I, J] += SUM * ZMAT[I, K];
                    }
                    for (int I = 1; I <= N; ++I)
                    {
                        int IP = I + NPT;
                        double TEMP = STEP[I];
                        if (K < IDZ) TEMP = -TEMP;
                        for (int J = 1; J <= I; ++J)
                            BMAT[IP, J] += TEMP * STEP[J];
                    }
                }
//
//     Update the right hand sides of the constraints.
//
                if (M > 0)
                {
                    for (int J = 1; J <= M; ++J)
                    {
                        double TEMP = ZERO;
                        for (int I = 1; I <= N; ++I)
                            TEMP += AMAT[I, J] * XOPT[I];
                        B[J] -= TEMP;
                    }
                }
//
//     The following instructions complete the shift of XBASE, including the
//       changes to the parameters of the quadratic model.
//
                for (int IH = 0, J = 1; J <= N; ++J)
                {
                    W[J] = ZERO;
                    for (int K = 1; K <= NPT; ++K)
                    {
                        W[J] += PQ[K] * XPT[K, J];
                        XPT[K, J] -= HALF * XOPT[J];
                    }
                    for (int I = 1; I <= J; ++I)
                    {
                        IH++;
                        HQ[IH] += W[I] * XOPT[J] + XOPT[I] * W[J];
                        BMAT[NPT + I, J] = BMAT[NPT + J, I];
                    }
                }
                for (int J = 1; J <= N; ++J)
                {
                    XBASE[J] += XOPT[J];
                    XOPT[J] = ZERO;
                    XPT[KOPT, J] = ZERO;
                }
            }
//
//     In the case KNEW=0, generate the next trust region step by calling
//       TRSTEP, where SNORM is the current trust region radius initially.
//       The final value of SNORM is the length of the calculated step,
//       except that SNORM is zero on return if the projected gradient is
//       unsuitable for starting the conjugate gradient iterations.
//
            double SNORM;
            double DELSAV = DELTA;
            int KSAVE = KNEW;
            if (KNEW == 0)
            {
                SNORM = DELTA;
                for (int I = 1; I <= N; ++I)
                    XNEW[I] = GOPT[I];
                TRSTEP(N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON, QFAC, RFAC, SNORM, STEP, XNEW, W, W(M + 1),
                    PQW,
                    PQW(NP), W(M + NP));
//
//     A trust region step is applied whenever its length, namely SNORM, is at
//       least HALF*DELTA. It is also applied if its length is at least 0.1999
//       times DELTA and if a line search of TRSTEP has caused a change to the
//       active set. Otherwise there is a branch below to label 530 or 560.
//
                double TEMP = HALF * DELTA;
                if (XNEW[1] >= HALF) TEMP = 0.1999 * DELTA;
                if (SNORM <= TEMP)
                {
                    DELTA *= HALF;
                    if (DELTA <= 1.4 * RHO) DELTA = RHO;
                    ++NVALA;
                    ++NVALB;
                    TEMP = SNORM / RHO;
                    if (DELSAV > RHO) TEMP = ONE;
                    if (TEMP >= HALF) NVALA = 0;
                    if (TEMP >= TENTH) NVALB = 0;
                    if (DELSAV > RHO) goto LINCOB_530;
                    if (NVALA < 5 && NVALB < 3) goto LINCOB_530;
                    if (SNORM > ZERO) KSAVE = -1;
                    goto LINCOB_560;
                }
                NVALA = 0;
                NVALB = 0;
//
//     Alternatively, KNEW is positive. Then the model step is calculated
//       within a trust region of radius DEL, after setting the gradient at
//       XBASE and the second derivative parameters of the KNEW-th Lagrange
//       function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively.
//
            }
            else
            {
                double DEL = Math.Max(TENTH * DELTA, RHO);
                for (int I = 1; I <= N; ++I)
                    W[I] = BMAT[KNEW, I];
                for (int K = 1; K <= NPT; ++K)
                    PQW[K] = ZERO;
                for (int J = 1; J <= NPTM; ++J)
                {
                    double TEMP = ZMAT[KNEW, J];
                    if (J < IDZ) TEMP = -TEMP;
                    for (int K = 1; K <= NPT; K++)
                        PQW[K] = PQW[K] + TEMP * ZMAT[K, J];
                }
                QMSTEP(N, NPT, M, AMAT, B, XPT, XOPT, NACT, IACT, RESCON, QFAC, KOPT, KNEW, DEL, STEP, W, PQW, W(NP),
                    W(NP + M), IFEAS);
            }
//
//     Set VQUAD to the change to the quadratic model when the move STEP is
//       made from XOPT. If STEP is a trust region step, then VQUAD should be
//       negative. If it is nonnegative due to rounding errors in this case,
//       there is a branch to label 530 to try to improve the model.
//
            double VQUAD = ZERO;
            for (int IH = 0, J = 1; J <= N; ++J)
            {
                VQUAD = VQUAD + STEP[J] * GOPT[J];
                for (int I = 1; I <= J; ++I)
                {
                    ++IH;
                    double TEMP = STEP[I] * STEP[J];
                    if (I == J) TEMP *= HALF;
                    VQUAD = VQUAD + TEMP * HQ[IH];
                }
            }
            for (int K = 1; K <= NPT; ++K)
            {
                double TEMP = ZERO;
                for (int J = 1; J <= N; ++J)
                {
                    TEMP = TEMP + XPT[K, J] * STEP[J];
                    SP[NPT + K] = TEMP;
                }
                VQUAD = VQUAD + HALF * PQ[K] * TEMP * TEMP;
            }
            if (KSAVE == 0 && VQUAD >= ZERO) goto LINCOB_530;
//
//     Calculate the next value of the objective function. The difference
//       between the actual new value of F and the value predicted by the
//       model is recorded in DIFF.
//
            LINCOB_220:

            ++NF;
            if (NF > MAXFUN)
            {
                --NF;
                if (IPRINT > 0) PRINT(logger, LINCOB_230);
                status = Status.MAXFUN_Reached;
                goto LINCOB_600;
            }
            double XDIFF = ZERO;
            for (int I = 1; I <= N; ++I)
            {
                XNEW[I] = XOPT[I] + STEP[I];
                X[I] = XBASE[I] + XNEW[I];
                XDIFF += Math.Pow(X[I] - XSAV[I], 2.0);
            }
            XDIFF = Math.Sqrt(XDIFF);
            if (KSAVE == -1) XDIFF = RHO;
            if (XDIFF <= TENTH * RHO || XDIFF >= DELTA + DELTA)
            {
                IFEAS = 0;
                if (IPRINT > 0) PRINT(logger, LINCOB_250);
                status = Status.X_UpdateHamperedByRoundingErrors;
                goto LINCOB_600;
            }
            if (KSAVE <= 0) IFEAS = 1;
            double F = (double)IFEAS;
            CALFUN(N, X, ref F);
            if (IPRINT == 3)
                PRINT(logger, LINCOB_260, NF, F, FORMAT("  ", "15:E6", X));
            if (KSAVE == -1) goto LINCOB_600;
            double DIFF = F - FOPT - VQUAD;
//
//     If X is feasible, then set DFFALT to the difference between the new
//       value of F and the value predicted by the alternative model.
//
            double DFFALT;
            if (IFEAS == 1 && ITEST < 3)
            {
                for (int K = 1; K <= NPT; ++K)
                {
                    PQW[K] = ZERO;
                    W[K] = FVAL[K] - FVAL[KOPT];
                }
                for (int J = 1; J <= NPTM; ++J)
                {
                    double SUM = ZERO;
                    for (int I = 1; I <= NPT; ++I)
                        SUM += W[I] * ZMAT[I, J];
                    if (J < IDZ) SUM = -SUM;
                    for (int K = 1; K <= NPT; ++K)
                        PQW[K] = PQW[K] + SUM * ZMAT[K, J];
                }
                double VQALT = ZERO;
                for (int K = 1; K <= NPT; ++K)
                {
                    double SUM = ZERO;
                    for (int J = 1; J <= N; ++J)
                        SUM += BMAT[K, J] * STEP[J];
                    VQALT = VQALT + SUM * W[K];
                    VQALT += PQW[K] * SP[NPT + K] * (HALF * SP[NPT + K] + SP[K]);
                }
                DFFALT = F - FOPT - VQALT;
            }
            if (ITEST == 3)
            {
                DFFALT = DIFF;
                ITEST = 0;
            }
//
//     Pick the next value of DELTA after a trust region step.
//
            double RATIO;
            if (KSAVE == 0)
            {
                RATIO = (F - FOPT) / VQUAD;
                if (RATIO <= TENTH)
                {
                    DELTA *= HALF;
                }
                else if (RATIO <= 0.7)
                {
                    DELTA = Math.Max(HALF * DELTA, SNORM);
                }
                else
                {
                    double TEMP = Math.Sqrt(2.0) * DELTA;
                    DELTA = Math.Max(HALF * DELTA, SNORM + SNORM);
                    DELTA = Math.Min(DELTA, TEMP);
                }
                if (DELTA <= 1.4 * RHO) DELTA = RHO;
            }
//
//     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
//       can be moved. If STEP is a trust region step, then KNEW is zero at
//       present, but a positive value is picked by subroutine UPDATE.
//
            UPDATE(N, NPT, XPT, BMAT, ZMAT, IDZ, NDIM, SP, STEP, KOPT, KNEW, PQW, W);
            if (KNEW == 0)
            {
                if (IPRINT > 0) PRINT(logger, LINCOB_320);
                status = Status.UpdatingFormulaDenominatorZero;
                goto LINCOB_600;
            }
//
//     If ITEST is increased to 3, then the next quadratic model is the
//       one whose second derivative matrix is least subject to the new
//       interpolation conditions. Otherwise the new model is constructed
//       by the symmetric Broyden method in the usual way.
//
            if (IFEAS == 1)
            {
                ++ITEST;
                if (Math.Abs(DFFALT) >= TENTH * Math.Abs(DIFF)) ITEST = 0;
            }
//
//     Update the second derivatives of the model by the symmetric Broyden
//       method, using PQW for the second derivative parameters of the new
//       KNEW-th Lagrange function. The contribution from the old parameter
//       PQ(KNEW) is included in the second derivative matrix HQ. W is used
//       later for the gradient of the new KNEW-th Lagrange function.       
//
            if (ITEST < 3)
            {
                for (int K = 1; K <= NPT; ++K)
                    PQW[K] = ZERO;
                for (int J = 1; J <= NPTM; ++J)
                {
                    double TEMP = ZMAT[KNEW, J];
                    if (TEMP != ZERO)
                    {
                        if (J < IDZ) TEMP = -TEMP;
                        for ( /*340*/ int K = 1; K <= NPT; ++K)
                            PQW[K] += TEMP * ZMAT[K, J];
                    }
                }
                for (int IH = 0, I = 1; I <= N; ++I)
                {
                    W[I] = BMAT[KNEW, I];
                    double TEMP = PQ[KNEW] * XPT[KNEW, I];
                    for (int J = 1; J <= I; ++J)
                    {
                        ++IH;
                        HQ[IH] += TEMP * XPT[KNEW, J];
                    }
                }
                PQ[KNEW] = ZERO;
                for (int K = 1; K <= NPT; ++K)
                    PQ[K] += DIFF * PQW[K];
            }
//
//     Include the new interpolation point with the corresponding updates of
//       SP. Also make the changes of the symmetric Broyden method to GOPT at
//       the old XOPT if ITEST is less than 3.
//
            FVAL[KNEW] = F;
            SP[KNEW] = SP[KOPT] + SP[NPT + KOPT];
            double SSQ = ZERO;
            for (int I = 1; I <= N; ++I)
            {
                XPT[KNEW, I] = XNEW[I];
                SSQ += STEP[I] * STEP[I];
            }
            SP[NPT + KNEW] = SP[NPT + KOPT] + SSQ;
            if (ITEST < 3)
            {
                for (int K = 1; K <= NPT; ++K)
                {
                    double TEMP = PQW[K] * SP[K];
                    for (int I = 1; I <= N; ++I)
                        W[I] += TEMP * XPT[K, I];
                }
                for (int I = 1; I <= N; ++I)
                    GOPT[I] += +DIFF * W[I];
            }
//
//     Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the
//       least calculated value so far with a feasible vector of variables.
//
            if (F < FOPT && IFEAS == 1)
            {
                FOPT = F;
                for (int J = 1; J <= N; ++J)
                {
                    XSAV[J] = X[J];
                    XOPT[J] = XNEW[J];
                }
                KOPT = KNEW;
                SNORM = Math.Sqrt(SSQ);
                for (int J = 1; J <= M; ++J)
                {
                    if (RESCON[J] >= DELTA + SNORM)
                    {
                        RESCON[J] = SNORM - RESCON[J];
                    }
                    else
                    {
                        RESCON[J] += +SNORM;
                        if (RESCON[J] + DELTA > ZERO)
                        {
                            double TEMP = B[J];
                            for (int I = 1; I <= N; ++I)
                                TEMP -= XOPT[I] * AMAT[I, J];
                            TEMP = Math.Max(TEMP, ZERO);
                            if (TEMP >= DELTA) TEMP = -TEMP;
                            RESCON[J] = TEMP;
                        }
                    }
                }
                for (int K = 1; K <= NPT; ++K)
                    SP[K] = SP[K] + SP[NPT + K];
//
//     Also revise GOPT when symmetric Broyden updating is applied.
//
                if (ITEST < 3)
                {
                    for (int IH = 0, J = 1; J <= N; ++J)
                    {
                        for (int I = 1; I <= J; ++I)
                        {
                            ++IH;
                            if (I < J) GOPT[J] += HQ[IH] * STEP[I];
                            GOPT[I] += HQ[IH] * STEP[J];
                        }
                    }
                    for (int K = 1; K <= NPT; ++K)
                    {
                        double TEMP = PQ[K] * SP[NPT + K];
                        for (int I = 1; I <= N; ++I)
                            GOPT[I] += TEMP * XPT[K, I];
                    }
                }
            }
//
//     Replace the current model by the least Frobenius norm interpolant if
//       this interpolant gives substantial reductions in the predictions
//       of values of F at feasible points.
//
            if (ITEST == 3)
            {
                for (int K = 1; K <= NPT; ++K)
                {
                    PQ[K] = ZERO;
                    W[K] = FVAL[K] - FVAL[KOPT];
                }
                for (int J = 1; J <= NPTM; ++J)
                {
                    double SUM = ZERO;
                    for (int I = 1; I <= NPT; ++I)
                        SUM += W[I] * ZMAT[I, J];
                    if (J < IDZ) SUM = -SUM;
                    for (int K = 1; K <= NPT; ++K)
                        PQ[K] = PQ[K] + SUM * ZMAT[K, J];
                }
                for (int J = 1; J <= N; ++J)
                {
                    GOPT[J] = ZERO;
                    for (int I = 1; I <= NPT; ++I)
                        GOPT[J] += W[I] * BMAT[I, J];
                }
                for (int K = 1; K <= NPT; ++K)
                {
                    double TEMP = PQ[K] * SP[K];
                    for (int I = 1; I <= N; ++I)
                        GOPT[I] += TEMP * XPT[K, I];
                }
                for (int IH = 1; IH <= NH; ++IH)
                    HQ[IH] = ZERO;
            }
//
//     If a trust region step has provided a sufficient decrease in F, then
//       branch for another trust region calculation. Every iteration that
//       takes a model step is followed by an attempt to take a trust region
//       step.
//
            KNEW = 0;
            if (KSAVE > 0) goto LINCOB_20;
            if (RATIO >= TENTH) goto LINCOB_20;
//
//     Alternatively, find out if the interpolation points are close enough
//       to the best point so far.
//
            LINCOB_530:

            double DISTSQ = Math.Max(DELTA * DELTA, 4.0 * RHO * RHO);
            for (int K = 1; K <= NPT; ++K)
            {
                double SUM = ZERO;
                for (int J = 1; J <= N; ++J)
                    SUM += Math.Pow(XPT[K, J] - XOPT[J], 2.0);
                if (SUM > DISTSQ)
                {
                    KNEW = K;
                    DISTSQ = SUM;
                }
            }
//
//     If KNEW is positive, then branch back for the next iteration, which
//       will generate a "model step". Otherwise, if the current iteration
//       has reduced F, or if DELTA was above its lower bound when the last
//       trust region step was calculated, then try a "trust region" step
//       instead.
//
            if (KNEW > 0) goto LINCOB_20;
            KNEW = 0;
            if (FOPT < FSAVE) goto LINCOB_20;
            if (DELSAV > RHO) goto LINCOB_20;
//
//     The calculations with the current value of RHO are complete.
//       Pick the next value of RHO.
//
            LINCOB_560:

            if (RHO > RHOEND)
            {
                DELTA = HALF * RHO;
                if (RHO > 250.0 * RHOEND)
                {
                    RHO *= TENTH;
                }
                else if (RHO <= 16.0 * RHOEND)
                {
                    RHO = RHOEND;
                }
                else
                {
                    RHO = Math.Sqrt(RHO * RHOEND);
                }
                DELTA = Math.Max(DELTA, RHO);
                if (IPRINT >= 2)
                {
                    if (IPRINT >= 3) PRINT(logger, LINCOB_570);
                    PRINT(logger, LINCOB_580, RHO, NF);
                    PRINT(logger, LINCOB_590, FOPT, FORMAT("  ", "15:E6", XBASE.Zip(XOPT, (xb, xo) => xb + xo)));
                }
                goto LINCOB_10;
            }
//
//     Return from the calculation, after branching to label 220 for another
//       Newton-Raphson step if it has not been tried before.
//
            if (KSAVE == -1) goto LINCOB_220;

            status = Status.Success;

            LINCOB_600:

            if (FOPT <= F || IFEAS == 0)
            {
                for (int I = 1; I <= N; ++I)
                    X[I] = XSAV[I];
                F = FOPT;
            }
            if (IPRINT >= 1)
            {
                PRINT(logger, LINCOB_620, NF);
                PRINT(logger, LINCOB_590, F, X);
            }
            W[1] = F;
            W[2] = (double)NF + HALF;

            return status;
        }

        private static void GETACT (int N, int M, double[,] AMAT, double[] B, int NACT, int[] IACT, double[,] QFAC, double[] RFAC, double SNORM, double[] RESNEW, double[] RESACT, double[] G, double[] DW, double[] VLAM, double[] W)
        {
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
            double TDEL = 0.2 * SNORM;
            double DDSAV = ZERO;
            for (int I = 1; I <= N; ++I)
            {
                DDSAV += G[I] * G[I];
                VLAM[I] = ZERO;
            }
            DDSAV *= 2.0;
//
//     Set the initial QFAC to the identity matrix in the case NACT=0.
//
      if (NACT == 0) {
          for (30 int I=1; I<=N; ++I) {
          for (20 int J=1; J<=N; ++J) {
   20     QFAC[I,J]=ZERO
   30     QFAC(I,I)=ONE
          goto 100
      }
//
//     Remove any constraints from the initial active set whose residuals
//       exceed TDEL.
//
      IFLAG=1
      IC=NACT
   40 if (RESACT(IC) > TDEL) goto 800
   50 IC=IC-1
      if (IC > 0) goto 40
//
//     Remove any constraints from the initial active set whose Lagrange
//       multipliers are nonnegative, and set the surviving multipliers.
//
      IFLAG=2
   60 if (NACT == 0) goto 100
      IC=NACT
   70 TEMP=ZERO
      for (80 int I=1; I<=N; ++I) {
   80 TEMP=TEMP+QFAC(I,IC)*G[I]
      IDIAG=(IC*IC+IC)/2
      if (IC < NACT) {
          JW=IDIAG+IC
          for (90 J=IC+1,NACT
          TEMP=TEMP-RFAC(JW)*VLAM[J]
   90     JW=JW+J
      }
      if (TEMP >= ZERO) goto 800
      VLAM(IC)=TEMP/RFAC(IDIAG)
      IC=IC-1
      if (IC > 0) goto 70
//
//     Set the new search direction D. Terminate if the 2-norm of D is zero
//       or does not decrease, or if NACT=N holds. The situation NACT=N
//       occurs for sufficiently large SNORM if the origin is in the convex
//       hull of the constraint gradients.
//
  100 if (NACT == N) goto 290
      for (110 J=NACT+1,N
      W[J]=ZERO
      for (110 int I=1; I<=N; ++I) {
  110 W[J]=W[J]+QFAC[I,J]*G[I]
      DD=ZERO
      for (130 int I=1; I<=N; ++I) {
      DW[I]=ZERO
      for (120 J=NACT+1,N
  120 DW[I]=DW[I]-W[J]*QFAC[I,J]
  130 DD=DD+DW[I]**2
      if (DD >= DDSAV) goto 290
      if (DD == ZERO) goto 300
      DDSAV=DD
      DNORM=Math.Sqrt(DD)
//
//     Pick the next integer L or terminate, a positive value of L being
//       the index of the most violated constraint. The purpose of CTOL
//       below is to estimate whether a positive value of VIOLMX may be
//       due to computer rounding errors.
//
      L=0
      if (M > 0) {
          TEST=DNORM/SNORM
          VIOLMX=ZERO
          for (150 J=1,M
          if (RESNEW[J] > ZERO && RESNEW[J] <= TDEL) {
              SUM=ZERO
              for (140 int I=1; I<=N; ++I) {
  140         SUM=SUM+AMAT[I,J]*DW[I]
              if (SUM > TEST*RESNEW[J]) {
                  if (SUM > VIOLMX) {
                      L=J
                      VIOLMX=SUM
                  }
              }
          }
  150     }
          CTOL=ZERO
          TEMP=0.01*DNORM
          if (VIOLMX > ZERO && VIOLMX < TEMP) {
              if (NACT > 0) {
                  for (170 K=1,NACT
                  J=IACT[K]
                  SUM=ZERO
                  for (160 int I=1; I<=N; ++I) {
  160             SUM=SUM+DW[I]*AMAT[I,J]
  170             CTOL=Math.Max(CTOL,Math.Abs(SUM))
              }
          }
      }
      W[1]=ONE
      if (L == 0) goto 300
      if (VIOLMX <= 10.0*CTOL) goto 300
//
//     Apply Givens rotations to the last (N-NACT) columns of QFAC so that
//       the first (NACT+1) columns of QFAC are the ones required for the
//       addition of the L-th constraint, and add the appropriate column
//       to RFAC.
//
      NACTP=NACT+1
      IDIAG=(NACTP*NACTP-NACTP)/2
      RDIAG=ZERO
      for (200 J=N,1,-1
      SPROD=ZERO
      for (180 int I=1; I<=N; ++I) {
  180 SPROD=SPROD+QFAC[I,J]*AMAT(I,L)
      if (J <= NACT) {
          RFAC(IDIAG+J)=SPROD
      } else {
          if (Math.Abs(RDIAG) <= 1.0E-20*Math.Abs(SPROD)) {
              RDIAG=SPROD
          } else {
              TEMP=Math.Sqrt(SPROD*SPROD+RDIAG*RDIAG)
              COSV=SPROD/TEMP
              SINV=RDIAG/TEMP
              RDIAG=TEMP
              for (190 int I=1; I<=N; ++I) {
              TEMP=COSV*QFAC[I,J]+SINV*QFAC(I,J+1)
              QFAC(I,J+1)=-SINV*QFAC[I,J]+COSV*QFAC(I,J+1)
  190         QFAC[I,J]=TEMP
          }
      }
  200 }
      if (RDIAG < ZERO) {
          for (210 int I=1; I<=N; ++I) {
  210     QFAC(I,NACTP)=-QFAC(I,NACTP)
      }
      RFAC(IDIAG+NACTP)=Math.Abs(RDIAG)
      NACT=NACTP
      IACT(NACT)=L
      RESACT(NACT)=RESNEW(L)
      VLAM(NACT)=ZERO
      RESNEW(L)=ZERO
//
//     Set the components of the vector VMU in W.
//
  220 W(NACT)=ONE/RFAC((NACT*NACT+NACT)/2)**2
      if (NACT > 1) {
          for (240 I=NACT-1,1,-1
          IDIAG=(I*I+I)/2
          JW=IDIAG+I
          SUM=ZERO
          for (230 J=I+1,NACT
          SUM=SUM-RFAC(JW)*W[J]
  230     JW=JW+J
  240     W[I]=SUM/RFAC(IDIAG)
      }
//
//     Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
//
      VMULT=VIOLMX
      IC=0
      J=1
  250 if (J < NACT) {
          if (VLAM[J] >= VMULT*W[J]) {
              IC=J
              VMULT=VLAM[J]/W[J]
          }
          J=J+1
          goto 250
      }
      for (260 J=1,NACT
  260 VLAM[J]=VLAM[J]-VMULT*W[J]
      if (IC > 0) VLAM(IC)=ZERO
      VIOLMX=Math.Max(VIOLMX-VMULT,ZERO)
      if (IC == 0) VIOLMX=ZERO
//
//     Reduce the active set if necessary, so that all components of the
//       new VLAM are negative, with resetting of the residuals of the
//       constraints that become inactive.
//
      IFLAG=3
      IC=NACT
  270 if (VLAM(IC) < ZERO) goto 280
      RESNEW(IACT(IC))=Math.Max(RESACT(IC),TINY)
      goto 800
  280 IC=IC-1
      if (IC > 0) goto 270
//
//     Calculate the next VMU if VIOLMX is positive. Return if NACT=N holds,
//       as then the active constraints imply D=0. Otherwise, go to label
//       100, to calculate the new D and to test for termination.
//
      if (VIOLMX > ZERO) goto 220
      if (NACT < N) goto 100
  290 DD=ZERO
  300 W[1]=DD
      RETURN
//
//     These instructions rearrange the active constraints so that the new
//       value of IACT(NACT) is the old value of IACT(IC). A sequence of
//       Givens rotations is applied to the current QFAC and RFAC. Then NACT
//       is reduced by one.
//
  800 RESNEW(IACT(IC))=Math.Max(RESACT(IC),TINY)
      JC=IC
  810 if (JC < NACT) {
          JCP=JC+1
          IDIAG=JC*JCP/2
          JW=IDIAG+JCP
          TEMP=Math.Sqrt(RFAC(JW-1)**2+RFAC(JW)**2)
          CVAL=RFAC(JW)/TEMP
          SVAL=RFAC(JW-1)/TEMP
          RFAC(JW-1)=SVAL*RFAC(IDIAG)
          RFAC(JW)=CVAL*RFAC(IDIAG)
          RFAC(IDIAG)=TEMP
          if (JCP < NACT) {
              for (820 J=JCP+1,NACT
              TEMP=SVAL*RFAC(JW+JC)+CVAL*RFAC(JW+JCP)
              RFAC(JW+JCP)=CVAL*RFAC(JW+JC)-SVAL*RFAC(JW+JCP)
              RFAC(JW+JC)=TEMP
  820         JW=JW+J
          }
          JDIAG=IDIAG-JC
          for (830 int I=1; I<=N; ++I) {
          if (I < JC) {
              TEMP=RFAC(IDIAG+I)
              RFAC(IDIAG+I)=RFAC(JDIAG+I)
              RFAC(JDIAG+I)=TEMP
          }
          TEMP=SVAL*QFAC(I,JC)+CVAL*QFAC(I,JCP)
          QFAC(I,JCP)=CVAL*QFAC(I,JC)-SVAL*QFAC(I,JCP)
  830     QFAC(I,JC)=TEMP
          IACT(JC)=IACT(JCP)
          RESACT(JC)=RESACT(JCP)
          VLAM(JC)=VLAM(JCP)
          JC=JCP
          goto 810
      }
      NACT=NACT-1
      goto (50,60,280),IFLAG
      END

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% prelim.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      private void PRELIM (N,NPT,M,AMAT,B,X,RHOBEG,IPRINT,XBASE,  XPT,FVAL,XSAV,XOPT,GOPT,KOPT,HQ,PQ,BMAT,ZMAT,IDZ,NDIM,  SP,RESCON,STEP,PQW,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AMAT(N,*),B(*),X(*),XBASE(*),XPT(NPT,*),FVAL(*),  XSAV(*),XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),  SP(*),RESCON(*),STEP(*),PQW(*),W(*)
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
      NPTM=NPT-N-1
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      RECIQ=Math.Sqrt(HALF)/RHOSQ
      TEST=0.2*RHOBEG
      IDZ=1
      KBASE=1
//
//     Set the initial elements of XPT, BMAT, SP and ZMAT to zero. 
//
      for (20 int J=1; J<=N; ++J) {
      XBASE[J]=X[J]
      for (10 int K=1; K<=NPT; ++K) {
   10 XPT[K,J]=ZERO
      for (20 I=1,NDIM
   20 BMAT[I,J]=ZERO
      for (30 int K=1; K<=NPT; ++K) {
      SP[K]=ZERO
      for (30 J=1,NPT-N-1
   30 ZMAT[K,J]=ZERO
//
//     Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT],
//       but they may be altered later to make a constraint violation
//       sufficiently large. The initial nonzero elements of BMAT and of
//       the first min[N,NPT-N-1] columns of ZMAT are set also.
//
      for (40 int J=1; J<=N; ++J) {
      XPT(J+1,J)=RHOBEG
      if (J < NPT-N) {
          JP=N+J+1
          XPT(JP,J)=-RHOBEG
          BMAT(J+1,J)=HALF/RHOBEG
          BMAT(JP,J)=-HALF/RHOBEG
          ZMAT(1,J)=-RECIQ-RECIQ
          ZMAT(J+1,J)=RECIQ
          ZMAT(JP,J)=RECIQ
      } else {
          BMAT(1,J)=-ONE/RHOBEG
          BMAT(J+1,J)=ONE/RHOBEG
          BMAT(NPT+J,J)=-HALF*RHOSQ
      }
   40 }
//
//     Set the remaining initial nonzero elements of XPT and ZMAT when the
//       number of interpolation points exceeds 2*N+1.
//
      if (NPT > 2*N+1) {
          for (50 K=N+1,NPT-N-1
          ITEMP=(K-1)/N
          IPT=K-ITEMP*N
          JPT=IPT+ITEMP
          if (JPT > N) JPT=JPT-N
          XPT(N+K+1,IPT)=RHOBEG
          XPT(N+K+1,JPT)=RHOBEG
          ZMAT(1,K)=RECIP
          ZMAT(IPT+1,K)=-RECIP
          ZMAT(JPT+1,K)=-RECIP
   50     ZMAT(N+K+1,K)=RECIP
      }
//
//     Update the constraint right hand sides to allow for the shift XBASE.
//
      if (M > 0) {
          for (70 J=1,M
          TEMP=ZERO
          for (60 int I=1; I<=N; ++I) {
   60     TEMP=TEMP+AMAT[I,J]*XBASE[I]
   70     B[J]=B[J]-TEMP
      }
//
//     Go through the initial points, shifting every infeasible point if
//       necessary so that its constraint violation is at least 0.2*RHOBEG.
//
      for (150 NF=1,NPT
      FEAS=ONE
      BIGV=ZERO
      J=0
   80 J=J+1
      if (J <= M && NF >= 2) {
          RESID=-B[J]
          for (90 int I=1; I<=N; ++I) {
   90     RESID=RESID+XPT(NF,I)*AMAT[I,J]
          if (RESID <= BIGV) goto 80
          BIGV=RESID
          JSAV=J
          if (RESID <= TEST) {
              FEAS=-ONE
              goto 80
          }
          FEAS=ZERO
      }
      if (FEAS < ZERO) {
          for (100 int I=1; I<=N; ++I) {
  100     STEP[I]=XPT(NF,I)+(TEST-BIGV)*AMAT(I,JSAV)
          for (110 int K=1; K<=NPT; ++K) {
          SP(NPT+K)=ZERO
          for (110 int J=1; J<=N; ++J) {
  110     SP(NPT+K)=SP(NPT+K)+XPT[K,J]*STEP[J]
          CALL UPDATE (N,NPT,XPT,BMAT,ZMAT,IDZ,NDIM,SP,STEP,      KBASE,NF,PQW,W)
          for (120 int I=1; I<=N; ++I) {
  120     XPT(NF,I)=STEP[I]
      }
//
//     Calculate the objective function at the current interpolation point,
//       and set KOPT to the index of the first trust region centre.
//
      for (130 int J=1; J<=N; ++J) {
  130 X[J]=XBASE[J]+XPT(NF,J)
      F=FEAS
      CALL CALFUN (N,X,F)
      if (IPRINT == 3) {
          PRINT 140, NF,F,(X[I],I=1,N)
  140     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,      '    The corresponding X is:'/(2X,5D15.6))
      }
      if (NF == 1) {
          KOPT=1
      } else if (F < FVAL[KOPT] && FEAS > ZERO) {
          KOPT=NF
      }
  150 FVAL(NF)=F
//
//     Set PQ for the first quadratic model.
//
      for (160 int J=1; J<=NPTM; ++J) {
      W[J]=ZERO
      for (160 int K=1; K<=NPT; ++K) {
  160 W[J]=W[J]+ZMAT[K,J]*FVAL[K]
      for (170 int K=1; K<=NPT; ++K) {
      PQ[K]=ZERO
      for (170 int J=1; J<=NPTM; ++J) {
  170 PQ[K]=PQ[K]+ZMAT[K,J]*W[J]
//
//     Set XOPT, SP, GOPT and HQ for the first quadratic model.
//
      for (180 int J=1; J<=N; ++J) {
      XOPT[J]=XPT(KOPT,J)
      XSAV[J]=XBASE[J]+XOPT[J]
  180 GOPT[J]=ZERO
      for (200 int K=1; K<=NPT; ++K) {
      SP[K]=ZERO
      for (190 int J=1; J<=N; ++J) {
  190 SP[K]=SP[K]+XPT[K,J]*XOPT[J]
      TEMP=PQ[K]*SP[K]
      for (200 int J=1; J<=N; ++J) {
  200 GOPT[J]=GOPT[J]+FVAL[K]*BMAT[K,J]+TEMP*XPT[K,J]
      for (210 I=1,(N*N+N)/2
  210 HQ[I]=ZERO
//
//     Set the initial elements of RESCON.
//
      for (230 J=1,M
      TEMP=B[J]
      for (220 int I=1; I<=N; ++I) {
  220 TEMP=TEMP-XOPT[I]*AMAT[I,J]
      TEMP=Math.Max(TEMP,ZERO)
      if (TEMP >= RHOBEG) TEMP=-TEMP
  230 RESCON[J]=TEMP  
      }

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% qmstep.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      private void QMSTEP (N,NPT,M,AMAT,B,XPT,XOPT,NACT,IACT,  RESCON,QFAC,KOPT,KNEW,DEL,STEP,GL,PQW,RSTAT,W,IFEAS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AMAT(N,*),B(*),XPT(NPT,*),XOPT(*),IACT(*),  RESCON(*),QFAC(N,*),STEP(*),GL(*),PQW(*),RSTAT(*),W(*)
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
      TEST=0.2*DEL
//
//     Replace GL by the gradient of LFUNC at the trust region centre, and
//       set the elements of RSTAT.
//
      for (20 int K=1; K<=NPT; ++K) {
      TEMP=ZERO
      for (10 int J=1; J<=N; ++J) {
   10 TEMP=TEMP+XPT[K,J]*XOPT[J]
      TEMP=PQW[K]*TEMP
      for (20 int I=1; I<=N; ++I) {
   20 GL[I]=GL[I]+TEMP*XPT[K,I]
      if (M > 0) {
          for (30 J=1,M
          RSTAT[J]=ONE
   30     if (Math.Abs(RESCON[J]) >= DEL) RSTAT[J]=-ONE
          for (40 K=1,NACT
   40     RSTAT(IACT[K])=ZERO
      }
//
//     Find the greatest modulus of LFUNC on a line through XOPT and
//       another interpolation point within the trust region.
//
      IFLAG=0
      VBIG=ZERO
      for (60 int K=1; K<=NPT; ++K) {
      if (K == KOPT) goto 60
      SS=ZERO
      SP=ZERO
      for (50 int I=1; I<=N; ++I) {
      TEMP=XPT[K,I]-XOPT[I]
      SS=SS+TEMP*TEMP
   50 SP=SP+GL[I]*TEMP
      STP=-DEL/Math.Sqrt(SS)
      if (K == KNEW) {
          if (SP*(SP-ONE) < ZERO) STP=-STP
          VLAG=Math.Abs(STP*SP)+STP*STP*Math.Abs(SP-ONE)
      } else {
          VLAG=Math.Abs(STP*(ONE-STP)*SP)
      }
      if (VLAG > VBIG) {
          KSAV=K
          STPSAV=STP
          VBIG=VLAG
      }
   60 }
//
//     Set STEP to the move that gives the greatest modulus calculated above.
//       This move may be replaced by a steepest ascent step from XOPT.
//
      GG=ZERO
      for (70 int I=1; I<=N; ++I) {
      GG=GG+GL[I]**2
   70 STEP[I]=STPSAV*(XPT(KSAV,I)-XOPT[I])
      VGRAD=DEL*Math.Sqrt(GG)
      if (VGRAD <= TENTH*VBIG) goto 220
//
//     Make the replacement if it provides a larger value of VBIG.
//
      GHG=ZERO
      for (90 int K=1; K<=NPT; ++K) {
      TEMP=ZERO
      for (80 int J=1; J<=N; ++J) {
   80 TEMP=TEMP+XPT[K,J]*GL[J]
   90 GHG=GHG+PQW[K]*TEMP*TEMP
      VNEW=VGRAD+Math.Abs(HALF*DEL*DEL*GHG/GG)
      if (VNEW > VBIG) {
          VBIG=VNEW
          STP=DEL/Math.Sqrt(GG)
          if (GHG < ZERO) STP=-STP
          for (100 int I=1; I<=N; ++I) {
  100     STEP[I]=STP*GL[I]
      }
      if (NACT == 0 || NACT == N) goto 220
//
//     Overwrite GL by its projection. Then set VNEW to the greatest
//       value of |LFUNC| on the projected gradient from XOPT subject to
//       the trust region bound. If VNEW is sufficiently large, then STEP
//       may be changed to a move along the projected gradient.
//
      for (110 K=NACT+1,N
      W[K]=ZERO
      for (110 int I=1; I<=N; ++I) {
  110 W[K]=W[K]+GL[I]*QFAC(I,K)
      GG=ZERO
      for (130 int I=1; I<=N; ++I) {
      GL[I]=ZERO
      for (120 K=NACT+1,N
  120 GL[I]=GL[I]+QFAC(I,K)*W[K]
  130 GG=GG+GL[I]**2
      VGRAD=DEL*Math.Sqrt(GG)
      if (VGRAD <= TENTH*VBIG) goto 220
      GHG=ZERO
      for (150 int K=1; K<=NPT; ++K) {
      TEMP=ZERO
      for (140 int J=1; J<=N; ++J) {
  140 TEMP=TEMP+XPT[K,J]*GL[J]
  150 GHG=GHG+PQW[K]*TEMP*TEMP
      VNEW=VGRAD+Math.Abs(HALF*DEL*DEL*GHG/GG)
//
//     Set W to the possible move along the projected gradient.
//
      STP=DEL/Math.Sqrt(GG)
      if (GHG < ZERO) STP=-STP
      WW=ZERO
      for (160 int I=1; I<=N; ++I) {
      W[I]=STP*GL[I]
  160 WW=WW+W[I]**2
//
//     Set STEP to W if W gives a sufficiently large value of the modulus
//       of the Lagrange function, and if W either preserves feasibility
//       or gives a constraint violation of at least 0.2*DEL. The purpose
//       of CTOL below is to provide a check on feasibility that includes
//       a tolerance for contributions from computer rounding errors.
//
      if (VNEW/VBIG >= 0.2) {
          IFEAS=1
          BIGV=ZERO
          J=0
  170     J=J+1
          if (J <= M) {
              if (RSTAT[J] == ONE) {
                  TEMP=-RESCON[J]
                  for (180 int I=1; I<=N; ++I) {
  180             TEMP=TEMP+W[I]*AMAT[I,J]
                  BIGV=Math.Max(BIGV,TEMP)
              }
              if (BIGV < TEST) goto 170
              IFEAS=0
          }
          CTOL=ZERO
          TEMP=0.01*Math.Sqrt(WW)
          if (BIGV > ZERO && BIGV < TEMP) {
              for (200 K=1,NACT
              J=IACT[K]
              SUM=ZERO
              for (190 int I=1; I<=N; ++I) {
  190         SUM=SUM+W[I]*AMAT[I,J]
  200         CTOL=Math.Max(CTOL,Math.Abs(SUM))
          }
          if (BIGV <= 10.0*CTOL || BIGV >= TEST) {
              for (210 int I=1; I<=N; ++I) {
  210         STEP[I]=W[I]
              goto 260
          }
      }
//
//     Calculate the greatest constraint violation at XOPT+STEP with STEP at
//       its original value. Modify STEP if this violation is unacceptable.
//
  220 IFEAS=1
      BIGV=ZERO
      RESMAX=ZERO
      J=0
  230 J=J+1
      if (J <= M) {
          if (RSTAT[J] < ZERO) goto 230
          TEMP=-RESCON[J]
          for (240 int I=1; I<=N; ++I) {
  240     TEMP=TEMP+STEP[I]*AMAT[I,J]
          RESMAX=Math.Max(RESMAX,TEMP)
          if (TEMP < TEST) {
              if (TEMP <= BIGV) goto 230
              BIGV=TEMP
              JSAV=J
              IFEAS=-1
              goto 230
          }
          IFEAS=0
      }
      if (IFEAS == -1) {
          for (250 int I=1; I<=N; ++I) {
  250     STEP[I]=STEP[I]+(TEST-BIGV)*AMAT(I,JSAV)
          IFEAS=0
      }
//
//     Return the calculated STEP and the value of IFEAS.
//
  260 }

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% trstep.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      private void TRSTEP (N,NPT,M,AMAT,B,XPT,HQ,PQ,NACT,IACT,RESCON,  QFAC,RFAC,SNORM,STEP,G,RESNEW,RESACT,D,DW,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AMAT(N,*),B(*),XPT(NPT,*),HQ(*),PQ(*),IACT(*),  RESCON(*),QFAC(N,*),RFAC(*),STEP(*),G(*),RESNEW(*),RESACT(*),  D(*),DW(*),W(*)
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
      SNSQ=SNORM*SNORM
//
//     Set the initial elements of RESNEW, RESACT and STEP.
//
      if (M > 0) {
          for (10 J=1,M
          RESNEW[J]=RESCON[J]
          if (RESCON[J] >= SNORM) {
              RESNEW[J]=-ONE
          } else if (RESCON[J] >= ZERO) {
              RESNEW[J]=Math.Max(RESNEW[J],TINY)
          }
   10     }
          if (NACT > 0) {
              for (20 K=1,NACT
              RESACT[K]=RESCON(IACT[K])
   20         RESNEW(IACT[K])=ZERO
          }
      }
      for (30 int I=1; I<=N; ++I) {
   30 STEP[I]=ZERO
      SS=ZERO
      REDUCT=ZERO
      NCALL=0
//
//     GETACT picks the active set for the current STEP. It also sets DW to
//       the vector closest to -G that is orthogonal to the normals of the
//       active constraints. DW is scaled to have length 0.2*SNORM, as then
//       a move of DW from STEP is allowed by the linear constraints.
//
   40 NCALL=NCALL+1
      CALL GETACT (N,M,AMAT,B,NACT,IACT,QFAC,RFAC,SNORM,RESNEW,  RESACT,G,DW,W,W(N+1))
      if (W(N+1) == ZERO) goto 320
      SCALE=0.2*SNORM/Math.Sqrt(W(N+1))
      for (50 int I=1; I<=N; ++I) {
   50 DW[I]=SCALE*DW[I]
//
//     If the modulus of the residual of an active constraint is substantial,
//       then set D to the shortest move from STEP to the boundaries of the
//       active constraints.
//
      RESMAX=ZERO
      if (NACT > 0) {
          for (60 K=1,NACT
   60     RESMAX=Math.Max(RESMAX,RESACT[K])
      }
      GAMMA=ZERO
      if (RESMAX > 1.0E-4*SNORM) {
          IR=0
          for (80 K=1,NACT
          TEMP=RESACT[K]
          if (K >= 2) {
              for (70 I=1,K-1
              IR=IR+1
   70         TEMP=TEMP-RFAC(IR)*W[I]
          }
          IR=IR+1
   80     W[K]=TEMP/RFAC(IR)
          for (90 int I=1; I<=N; ++I) {
          D[I]=ZERO
          for (90 K=1,NACT
   90     D[I]=D[I]+W[K]*QFAC(I,K)
//
//     The vector D that has just been calculated is also the shortest move
//       from STEP+DW to the boundaries of the active constraints. Set GAMMA
//       to the greatest steplength of this move that satisfies the trust
//       region bound.
//
          RHS=SNSQ
          DS=ZERO
          DD=ZERO
          for (100 int I=1; I<=N; ++I) {
          SUM=STEP[I]+DW[I]
          RHS=RHS-SUM*SUM
          DS=DS+D[I]*SUM
  100     DD=DD+D[I]**2
          if (RHS > ZERO) {
              TEMP=Math.Sqrt(DS*DS+DD*RHS)
              if (DS <= ZERO) {
                  GAMMA=(TEMP-DS)/DD
              } else {
                  GAMMA=RHS/(TEMP+DS)
              }
          }
//
//     Reduce the steplength GAMMA if necessary so that the move along D
//       also satisfies the linear constraints.
//
          J=0
  110     if (GAMMA > ZERO) {
              J=J+1
              if (RESNEW[J] > ZERO) {
                  AD=ZERO
                  ADW=ZERO
                  for (120 int I=1; I<=N; ++I) {
                  AD=AD+AMAT[I,J]*D[I]
  120             ADW=ADW+AMAT[I,J]*DW[I]
                  if (AD > ZERO) {
                      TEMP=Math.Max((RESNEW[J]-ADW)/AD,ZERO)
                      GAMMA=Math.Min(GAMMA,TEMP)
                  }
              }
              if (J < M) goto 110
          }
          GAMMA=Math.Min(GAMMA,ONE)
      }
//
//     Set the next direction for seeking a reduction in the model function
//       subject to the trust region bound and the linear constraints.
//
      if (GAMMA <= ZERO) {
          for (130 int I=1; I<=N; ++I) {
  130     D[I]=DW[I]
          ICOUNT=NACT
      } else {
          for (140 int I=1; I<=N; ++I) {
  140     D[I]=DW[I]+GAMMA*D[I]
          ICOUNT=NACT-1
      }
      ALPBD=ONE
//
//     Set ALPHA to the steplength from STEP along D to the trust region
//       boundary. Return if the first derivative term of this step is
//       sufficiently small or if no further progress is possible.
//
  150 ICOUNT=ICOUNT+1
      RHS=SNSQ-SS
      if (RHS <= ZERO) goto 320
      DG=ZERO
      DS=ZERO
      DD=ZERO
      for (160 int I=1; I<=N; ++I) {
      DG=DG+D[I]*G[I]
      DS=DS+D[I]*STEP[I]
  160 DD=DD+D[I]**2
      if (DG >= ZERO) goto 320
      TEMP=Math.Sqrt(RHS*DD+DS*DS)
      if (DS <= ZERO) {
          ALPHA=(TEMP-DS)/DD
      } else {
          ALPHA=RHS/(TEMP+DS)
      }
      if (-ALPHA*DG <= CTEST*REDUCT) goto 320
//
//     Set DW to the change in gradient along D.
//
      IH=0
      for (170 int J=1; J<=N; ++J) {
      DW[J]=ZERO
      for (170 int I=1; I<=J; ++I) {
      IH=IH+1
      if (I < J) DW[J]=DW[J]+HQ[IH]*D[I]
  170 DW[I]=DW[I]+HQ[IH]*D[J]
      for (190 int K=1; K<=NPT; ++K) {
      TEMP=ZERO
      for (180 int J=1; J<=N; ++J) {
  180 TEMP=TEMP+XPT[K,J]*D[J]
      TEMP=PQ[K]*TEMP
      for (190 int I=1; I<=N; ++I) {
  190 DW[I]=DW[I]+TEMP*XPT[K,I]
//
//     Set DGD to the curvature of the model along D. Then reduce ALPHA if
//       necessary to the value that minimizes the model.
//
      DGD=ZERO
      for (200 int I=1; I<=N; ++I) {
  200 DGD=DGD+D[I]*DW[I]
      ALPHT=ALPHA
      if (DG+ALPHA*DGD > ZERO) {
          ALPHA=-DG/DGD
      }
//
//     Make a further reduction in ALPHA if necessary to preserve feasibility,
//       and put some scalar products of D with constraint gradients in W.
//
      ALPHM=ALPHA
      JSAV=0
      if (M > 0) {
          for (220 J=1,M
          AD=ZERO
          if (RESNEW[J] > ZERO) {
              for (210 int I=1; I<=N; ++I) {
  210         AD=AD+AMAT[I,J]*D[I]
              if (ALPHA*AD > RESNEW[J]) {
                  ALPHA=RESNEW[J]/AD
                  JSAV=J
              }
          }
  220     W[J]=AD
      }
      ALPHA=Math.Max(ALPHA,ALPBD)
      ALPHA=Math.Min(ALPHA,ALPHM)
      if (ICOUNT == NACT) ALPHA=Math.Min(ALPHA,ONE)
//
//     Update STEP, G, RESNEW, RESACT and REDUCT.
//
      SS=ZERO
      for (230 int I=1; I<=N; ++I) {
      STEP[I]=STEP[I]+ALPHA*D[I]
      SS=SS+STEP[I]**2
  230 G[I]=G[I]+ALPHA*DW[I]
      if (M > 0) {
          for (240 J=1,M
          if (RESNEW[J] > ZERO) {
              RESNEW[J]=Math.Max(RESNEW[J]-ALPHA*W[J],TINY)
          }
  240     }
      }
      if (ICOUNT == NACT && NACT > 0) {
          for (250 K=1,NACT
  250     RESACT[K]=(ONE-GAMMA)*RESACT[K]
      }
      REDUCT=REDUCT-ALPHA*(DG+HALF*ALPHA*DGD)
//
//     Test for termination. Branch to label 40 if there is a new active
//       constraint and if the distance from STEP to the trust region
//       boundary is at least 0.2*SNORM.
//
      if (ALPHA == ALPHT) goto 320
      TEMP=-ALPHM*(DG+HALF*ALPHM*DGD)
      if (TEMP <= CTEST*REDUCT) goto 320
      if (JSAV > 0) {
          if (SS <= 0.64*SNSQ) goto 40
          goto 320
      }
      if (ICOUNT == N) goto 320
//
//     Calculate the next search direction, which is conjugate to the
//       previous one except in the case ICOUNT=NACT.
//
      if (NACT > 0) {
          for (260 J=NACT+1,N
          W[J]=ZERO
          for (260 int I=1; I<=N; ++I) {
  260     W[J]=W[J]+G[I]*QFAC[I,J]
          for (280 int I=1; I<=N; ++I) {
          TEMP=ZERO
          for (270 J=NACT+1,N
  270     TEMP=TEMP+QFAC[I,J]*W[J]
  280     W(N+I)=TEMP
      } else {
          for (290 int I=1; I<=N; ++I) {
  290     W(N+I)=G[I]
      }
      if (ICOUNT == NACT) {
          BETA=ZERO
      } else {
          WGD=ZERO
          for (300 int I=1; I<=N; ++I) {
  300     WGD=WGD+W(N+I)*DW[I]
          BETA=WGD/DGD
      }
      for (310 int I=1; I<=N; ++I) {
  310 D[I]=-W(N+I)+BETA*D[I]
      ALPBD=ZERO
      goto 150
//
//     Return from the subroutine.
//
  320 SNORM=ZERO
      if (REDUCT > ZERO) SNORM=Math.Sqrt(SS)
      G[1]=ZERO
      if (NCALL > 1) G[1]=ONE
      }

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      private void UPDATE (N,NPT,XPT,BMAT,ZMAT,IDZ,NDIM,SP,STEP,  KOPT,KNEW,VLAG,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),SP(*),STEP(*),  VLAG(*),W(*)
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
      NPTM=NPT-N-1
//
//     Calculate VLAG and BETA for the current choice of STEP. The first NPT
//       elements of VLAG are set to the values of the Lagrange functions at
//       XPT(KOPT,.)+STEP(.). The first NPT components of W_check are held
//       in W, where W_check is defined in a paper on the updating method.
//
      for (20 int K=1; K<=NPT; ++K) {
      W[K]=SP(NPT+K)*(HALF*SP(NPT+K)+SP[K])
      SUM=ZERO
      for (10 int J=1; J<=N; ++J) {
   10 SUM=SUM+BMAT[K,J]*STEP[J]
   20 VLAG[K]=SUM
      BETA=ZERO
      for (40 K=1,NPTM
      SUM=ZERO
      for (30 int I=1; I<=NPT; ++I) {
   30 SUM=SUM+ZMAT(I,K)*W[I]
      if (K < IDZ) {
          BETA=BETA+SUM*SUM
          SUM=-SUM
      } else {
          BETA=BETA-SUM*SUM
      }
      for (40 int I=1; I<=NPT; ++I) {
   40 VLAG[I]=VLAG[I]+SUM*ZMAT(I,K)
      BSUM=ZERO
      DX=ZERO
      SSQ=ZERO
      for (70 int J=1; J<=N; ++J) {
      SUM=ZERO
      for (50 int I=1; I<=NPT; ++I) {
   50 SUM=SUM+W[I]*BMAT[I,J]
      BSUM=BSUM+SUM*STEP[J]
      JP=NPT+J
      for (60 K=1,N
   60 SUM=SUM+BMAT(JP,K)*STEP[K]
      VLAG(JP)=SUM
      BSUM=BSUM+SUM*STEP[J]
      DX=DX+STEP[J]*XPT(KOPT,J)
   70 SSQ=SSQ+STEP[J]**2
      BETA=DX*DX+SSQ*(SP[KOPT]+DX+DX+HALF*SSQ)+BETA-BSUM
      VLAG[KOPT]=VLAG[KOPT]+ONE
//
//     If KNEW is zero initially, then pick the index of the interpolation
//       point to be deleted, by maximizing the absolute value of the
//       denominator of the updating formula times a weighting factor.
//       
//
      if (KNEW == 0) {
          DENMAX=ZERO
          for (100 int K=1; K<=NPT; ++K) {
          HDIAG=ZERO
          for (80 int J=1; J<=NPTM; ++J) {
          TEMP=ONE
          if (J < IDZ) TEMP=-ONE
   80     HDIAG=HDIAG+TEMP*ZMAT[K,J]**2
          DENABS=Math.Abs(BETA*HDIAG+VLAG[K]**2)
          DISTSQ=ZERO
          for (90 int J=1; J<=N; ++J) {
   90     DISTSQ=DISTSQ+(XPT[K,J]-XPT(KOPT,J))**2
          TEMP=DENABS*DISTSQ*DISTSQ
          if (TEMP > DENMAX) {
              DENMAX=TEMP
              KNEW=K
          }
  100     }
      }
//
//     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
//
      JL=1
      if (NPTM >= 2) {
          for (120 J=2,NPTM
          if (J == IDZ) {
              JL=IDZ
          } else if (ZMAT(KNEW,J) != ZERO) {
              TEMP=Math.Sqrt(ZMAT(KNEW,JL)**2+ZMAT(KNEW,J)**2)
              TEMPA=ZMAT(KNEW,JL)/TEMP
              TEMPB=ZMAT(KNEW,J)/TEMP
              for (110 int I=1; I<=NPT; ++I) {
              TEMP=TEMPA*ZMAT(I,JL)+TEMPB*ZMAT[I,J]
              ZMAT[I,J]=TEMPA*ZMAT[I,J]-TEMPB*ZMAT(I,JL)
  110         ZMAT(I,JL)=TEMP
              ZMAT(KNEW,J)=ZERO
          }
  120     }
      }
//
//     Put the first NPT components of the KNEW-th column of the Z Z^T matrix
//       into W, and calculate the parameters of the updating formula.
//
      TEMPA=ZMAT(KNEW,1)
      if (IDZ >= 2) TEMPA=-TEMPA
      if (JL > 1) TEMPB=ZMAT(KNEW,JL)
      for (130 int I=1; I<=NPT; ++I) {
      W[I]=TEMPA*ZMAT(I,1)
      if (JL > 1) W[I]=W[I]+TEMPB*ZMAT(I,JL)
  130 }
      ALPHA=W(KNEW)
      TAU=VLAG(KNEW)
      TAUSQ=TAU*TAU
      DENOM=ALPHA*BETA+TAUSQ
      VLAG(KNEW)=VLAG(KNEW)-ONE
      if (DENOM == ZERO) {
          KNEW=0
          goto 180
      }
      SQRTDN=Math.Sqrt(Math.Abs(DENOM))
//
//     Complete the updating of ZMAT when there is only one nonzero element
//       in the KNEW-th row of the new matrix ZMAT. IFLAG is set to one when
//       the value of IDZ is going to be reduced.
//
      IFLAG=0
      if (JL == 1) {
          TEMPA=TAU/SQRTDN
          TEMPB=ZMAT(KNEW,1)/SQRTDN
          for (140 int I=1; I<=NPT; ++I) {
  140     ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG[I]
          if (DENOM < ZERO) {
              if (IDZ == 1) {
                  IDZ=2
              } else {
                  IFLAG=1
              }
          }
      } else {
//
//     Complete the updating of ZMAT in the alternative case.
//
          JA=1
          if (BETA >= ZERO) JA=JL
          JB=JL+1-JA
          TEMP=ZMAT(KNEW,JB)/DENOM
          TEMPA=TEMP*BETA
          TEMPB=TEMP*TAU
          TEMP=ZMAT(KNEW,JA)
          SCALA=ONE/Math.Sqrt(Math.Abs(BETA)*TEMP*TEMP+TAUSQ)
          SCALB=SCALA*SQRTDN
          for (150 int I=1; I<=NPT; ++I) {
          ZMAT(I,JA)=SCALA*(TAU*ZMAT(I,JA)-TEMP*VLAG[I])
  150     ZMAT(I,JB)=SCALB*(ZMAT(I,JB)-TEMPA*W[I]-TEMPB*VLAG[I])
          if (DENOM <= ZERO) {
              if (BETA < ZERO) {
                  IDZ=IDZ+1
              } else {
                  IFLAG=1
              }
          }
      }
//
//     Reduce IDZ when the diagonal part of the ZMAT times Diag(DZ) times
//       ZMAT^T factorization gains another positive element. Then exchange
//       the first and IDZ-th columns of ZMAT.
//
      if (IFLAG == 1) {
          IDZ=IDZ-1
          for (160 int I=1; I<=NPT; ++I) {
          TEMP=ZMAT(I,1)
          ZMAT(I,1)=ZMAT(I,IDZ)
  160     ZMAT(I,IDZ)=TEMP
      }
//
//     Finally, update the matrix BMAT.
//
      for (170 int J=1; J<=N; ++J) {
      JP=NPT+J
      W(JP)=BMAT(KNEW,J)
      TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
      TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
      for (170 int I=1; I<=JP; ++I) {
      BMAT[I,J]=BMAT[I,J]+TEMPA*VLAG[I]+TEMPB*W[I]
      if (I > NPT) BMAT(JP,I-NPT)=BMAT[I,J]
  170 }
  180 }

    #endregion

} // class Lincoa
} // namespace Cureos.Numerics