using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Cureos.Numerics.Optimizers
{

    #region DELEGATES

    public delegate void CalfunDel(int n, double[] x, ref double f);

    #endregion

    public class Lincoa
    {
        #region INNER TYPES

        public enum Status
        {
            Success,
            N_TooSmall,
            NPT_OutOfRange,
            MAXFUN_NotLargerThan_NPT,
            ConstraintGradientIsZero,
            MAXFUN_Reached,
            X_RoundingErrorsPreventUpdate,
            UpdatingFormulaDenominatorZero
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

        #endregion

        #region METHODS

        public static Status LINCOA(CalfunDel CALFUN, int N, int NPT, int M, double[,] A, int IA, double[] B, double[] X,
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
                return Status.N_TooSmall;
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
            double[,] AMAT = new double[1 + N, 1 + M];
            double[] BB = new double[1 + M]; // B in LINCOB
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
                        AMAT[I, J] = A[I, J] / TEMP;
                    }
                }
            }
            if (IFLAG == 1)
            {
                if (IPRINT > 0) PRINT(logger, LINCOA_70);
            }
            return LINCOB(CALFUN, N, NPT, M, AMAT, BB, X, RHOBEG, RHOEND, IPRINT, MAXFUN, logger);
        }

        private static Status LINCOB(CalfunDel CALFUN, int N, int NPT, int M, double[,] AMAT, double[] B, double[] X,
            double RHOBEG, double RHOEND, int IPRINT, int MAXFUN, TextWriter logger)
        {
            Status? status = null;
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
            int IAMAT = Math.Max(M + 3 * N, Math.Max(2 * M + N, 2 * NPT)) + 1;
            int NDIM = NPT + N;
//
//     Partition the working space array, so that different parts of it can be
//     treated separately by the subroutine that performs the main calculation.
//
            double[] XBASE = new double[1 + N];
            double[,] XPT = new double[1 + NPT, 1 + N];
            double[] FVAL = new double[1 + NPT];
            double[] XSAV = new double[1 + N];
            double[] XOPT = new double[1 + N];
            double[] GOPT = new double[1 + N];
            double[] HQ = new double[1 + (N * NP) / 2];
            double[] PQ = new double[1 + NPT];
            double[,] BMAT = new double[1 + NDIM, 1 + N];
            double[,] ZMAT = new double[1 + NPT, 1 + NPTM];
            double[] STEP = new double[1 + N];
            double[] SP = new double[1 + NPT + NPT];
            double[] XNEW = new double[1 + N];
            int[] IACT = new int[1 + N];
            double[] RESCON = new double[1 + M];
            double[,] QFAC = new double[1 + N, 1 + N];
            double[] RFAC = new double[1 + (N * NP) / 2];
            double[] PQW = new double[1 + NPT + N];
            double[] W = new double[IAMAT];
//
//     Set the elements of XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT,
//       ZMAT and SP for the first iteration. An important feature is that,
//       if the interpolation point XPT(K,.) is not feasible, where K is any
//       integer from [1,NPT], then a change is made to XPT(K,.) if necessary
//       so that the constraint violation is at least 0.2*RHOBEG. Also KOPT
//       is set so that XPT(KOPT,.) is the initial trust region centre.
//
            int KOPT, IDZ;
            PRELIM(CALFUN, N, NPT, M, AMAT, B, X, RHOBEG, IPRINT, XBASE, XPT, FVAL, XSAV, XOPT, GOPT, out KOPT, HQ, PQ,
                BMAT, ZMAT, out IDZ, NDIM, SP, RESCON, logger);
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
            double F = ZERO;
            double VQUAD = ZERO;
            double SNORM = ZERO;
            double DELSAV = DELTA;
            int KSAVE = KNEW;
            if (KNEW == 0)
            {
                SNORM = DELTA;
                for (int I = 1; I <= N; ++I)
                    XNEW[I] = GOPT[I];
                TRSTEP(N, NPT, M, AMAT, B, XPT, HQ, PQ, ref NACT, IACT, RESCON, QFAC, RFAC, ref SNORM, STEP, XNEW);
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
                        PQW[K] += TEMP * ZMAT[K, J];
                }
                QMSTEP(N, NPT, M, AMAT, B, XPT, XOPT, NACT, IACT, RESCON, QFAC, KOPT, KNEW, DEL, STEP, W, PQW, ref IFEAS);
            }
//
//     Set VQUAD to the change to the quadratic model when the move STEP is
//       made from XOPT. If STEP is a trust region step, then VQUAD should be
//       negative. If it is nonnegative due to rounding errors in this case,
//       there is a branch to label 530 to try to improve the model.
//
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
                status = Status.X_RoundingErrorsPreventUpdate;
                goto LINCOB_600;
            }
            if (KSAVE <= 0) IFEAS = 1;
            F = (double)IFEAS;
            CALFUN(N, X, ref F);
            if (IPRINT == 3)
                PRINT(logger, LINCOB_260, NF, F, FORMAT("  ", "18:E6", X, 1, N));
            if (KSAVE == -1) goto LINCOB_600;
            double DIFF = F - FOPT - VQUAD;
//
//     If X is feasible, then set DFFALT to the difference between the new
//       value of F and the value predicted by the alternative model.
//
            double DFFALT = ZERO;
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
            double RATIO = ZERO;
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
            UPDATE(N, NPT, XPT, BMAT, ZMAT, IDZ, NDIM, SP, STEP, KOPT, ref KNEW);
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
                    PRINT(logger, LINCOB_590, FOPT, FORMAT("  ", "18:E6", XBASE.Zip(XOPT, (xb, xo) => xb + xo), 1, N));
                }
                goto LINCOB_10;
            }
//
//     Return from the calculation, after branching to label 220 for another
//       Newton-Raphson step if it has not been tried before.
//
            if (KSAVE == -1) goto LINCOB_220;

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
                PRINT(logger, LINCOB_590, F, FORMAT("  ", "18:E6", X, 1, N));
            }
            W[1] = F;
            W[2] = (double)NF + HALF;

            return status.GetValueOrDefault(Status.Success);
        }

        private static double GETACT(int N, int M, double[,] AMAT, double[] B, ref int NACT, int[] IACT, double[,] QFAC,
            double[] RFAC, double SNORM, double[] RESNEW, double[] RESACT, double[] G, double[] DW)
        {
            double[] VLAM = new double[1 + N];
            double[] W = new double[1 + N];
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
            double VIOLMX = ZERO;
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
            if (NACT == 0)
            {
                for (int I = 1; I <= N; ++I)
                {
                    for (int J = 1; J <= N; ++J)
                        QFAC[I, J] = ZERO;
                    QFAC[I, I] = ONE;
                }
                goto GETACT_100;
            }
//
//     Remove any constraints from the initial active set whose residuals
//       exceed TDEL.
//
            int IFLAG = 1;
            int IC = NACT;

            GETACT_40:

            if (RESACT[IC] > TDEL) goto GETACT_800;

            GETACT_50:

            --IC;
            if (IC > 0) goto GETACT_40;
//
//     Remove any constraints from the initial active set whose Lagrange
//       multipliers are nonnegative, and set the surviving multipliers.
//
            IFLAG = 2;

            GETACT_60:

            if (NACT == 0) goto GETACT_100;
            IC = NACT;

            GETACT_70:

            double TEMP = ZERO;
            for (int I = 1; I <= N; ++I)
                TEMP += QFAC[I, IC] * G[I];
            int IDIAG = (IC * IC + IC) / 2;
            if (IC < NACT)
            {
                int JW = IDIAG + IC;
                for (int J = IC + 1; J <= NACT; ++J)
                {
                    TEMP -= RFAC[JW] * VLAM[J];
                    JW += J;
                }
            }
            if (TEMP >= ZERO) goto GETACT_800;
            VLAM[IC] = TEMP / RFAC[IDIAG];
            --IC;
            if (IC > 0) goto GETACT_70;
//
//     Set the new search direction D. Terminate if the 2-norm of D is zero
//       or does not decrease, or if NACT=N holds. The situation NACT=N
//       occurs for sufficiently large SNORM if the origin is in the convex
//       hull of the constraint gradients.
//
            GETACT_100:

            if (NACT == N) goto GETACT_290;
            for (int J = NACT + 1; J <= N; ++J)
            {
                W[J] = ZERO;
                for (int I = 1; I <= N; ++I)
                    W[J] = W[J] + QFAC[I, J] * G[I];
            }
            double DD = ZERO;
            for (int I = 1; I <= N; ++I)
            {
                DW[I] = ZERO;
                for (int J = NACT + 1; J <= N; ++J)
                    DW[I] -= W[J] * QFAC[I, J];
                DD += DW[I] * DW[I];
            }
            if (DD >= DDSAV) goto GETACT_290;
            if (DD == ZERO) goto GETACT_300;
            DDSAV = DD;
            double DNORM = Math.Sqrt(DD);
//
//     Pick the next integer L or terminate, a positive value of L being
//       the index of the most violated constraint. The purpose of CTOL
//       below is to estimate whether a positive value of VIOLMX may be
//       due to computer rounding errors.
//
            int L = 0;
            VIOLMX = ZERO;
            double CTOL = ZERO;
            if (M > 0)
            {
                double TEST = DNORM / SNORM;
                for (int J = 1; J <= M; ++J)
                {
                    if (RESNEW[J] > ZERO && RESNEW[J] <= TDEL)
                    {
                        double SUM = ZERO;
                        for (int I = 1; I <= N; ++I)
                            SUM += AMAT[I, J] * DW[I];
                        if (SUM > TEST * RESNEW[J])
                        {
                            if (SUM > VIOLMX)
                            {
                                L = J;
                                VIOLMX = SUM;
                            }
                        }
                    }
                }
                TEMP = 0.01 * DNORM;
                if (VIOLMX > ZERO && VIOLMX < TEMP)
                {
                    if (NACT > 0)
                    {
                        for (int K = 1; K <= NACT; ++K)
                        {
                            int J = IACT[K];
                            double SUM = ZERO;
                            for (int I = 1; I <= N; ++I)
                                SUM += DW[I] * AMAT[I, J];
                            CTOL = Math.Max(CTOL, Math.Abs(SUM));
                        }
                    }
                }
            }
            W[1] = ONE;
            if (L == 0) goto GETACT_300;
            if (VIOLMX <= 10.0 * CTOL) goto GETACT_300;
//
//     Apply Givens rotations to the last (N-NACT) columns of QFAC so that
//       the first (NACT+1) columns of QFAC are the ones required for the
//       addition of the L-th constraint, and add the appropriate column
//       to RFAC.
//
            int NACTP = NACT + 1;
            IDIAG = (NACTP * NACTP - NACTP) / 2;
            double RDIAG = ZERO;
            for (int J = N; J >= 1; --J)
            {
                double SPROD = ZERO;
                for (int I = 1; I <= N; ++I)
                    SPROD = SPROD + QFAC[I, J] * AMAT[I, L];
                if (J <= NACT)
                {
                    RFAC[IDIAG + J] = SPROD;
                }
                else
                {
                    if (Math.Abs(RDIAG) <= 1.0E-20 * Math.Abs(SPROD))
                    {
                        RDIAG = SPROD;
                    }
                    else
                    {
                        TEMP = Math.Sqrt(SPROD * SPROD + RDIAG * RDIAG);
                        double COSV = SPROD / TEMP;
                        double SINV = RDIAG / TEMP;
                        RDIAG = TEMP;
                        for (int I = 1; I <= N; ++I)
                        {
                            TEMP = COSV * QFAC[I, J] + SINV * QFAC[I, J + 1];
                            QFAC[I, J + 1] = -SINV * QFAC[I, J] + COSV * QFAC[I, J + 1];
                            QFAC[I, J] = TEMP;
                        }
                    }
                }
            }

            if (RDIAG < ZERO)
            {
                for (int I = 1; I <= N; ++I)
                    QFAC[I, NACTP] = -QFAC[I, NACTP];
            }
            RFAC[IDIAG + NACTP] = Math.Abs(RDIAG);
            NACT = NACTP;
            IACT[NACT] = L;
            RESACT[NACT] = RESNEW[L];
            VLAM[NACT] = ZERO;
            RESNEW[L] = ZERO;
//
//     Set the components of the vector VMU in W.
//
            GETACT_220:

            W[NACT] = ONE / Math.Pow(RFAC[(NACT * NACT + NACT) / 2], 2.0);
            if (NACT > 1)
            {
                for (int I = NACT - 1; I >= 1; --I)
                {
                    IDIAG = (I * I + I) / 2;
                    int JW = IDIAG + I;
                    double SUM = ZERO;
                    for (int J = I + 1; J <= NACT; ++J)
                    {
                        SUM -= RFAC[JW] * W[J];
                        JW += +J;
                    }
                    W[I] = SUM / RFAC[IDIAG];
                }
            }
//
//     Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
//
            double VMULT = VIOLMX;
            IC = 0;
            {
                int J = 1;
                while (J < NACT)
                {
                    if (VLAM[J] >= VMULT * W[J])
                    {
                        IC = J;
                        VMULT = VLAM[J] / W[J];
                    }
                    ++J;
                }
            }
            for (int J = 1; J <= NACT; ++J)
                VLAM[J] = VLAM[J] - VMULT * W[J];
            if (IC > 0) VLAM[IC] = ZERO;
            VIOLMX = Math.Max(VIOLMX - VMULT, ZERO);
            if (IC == 0) VIOLMX = ZERO;
//
//     Reduce the active set if necessary, so that all components of the
//       new VLAM are negative, with resetting of the residuals of the
//       constraints that become inactive.
//
            IFLAG = 3;
            IC = NACT;

            GETACT_270:

            if (VLAM[IC] < ZERO) goto GETACT_280;
            RESNEW[IACT[IC]] = Math.Max(RESACT[IC], TINY);
            goto GETACT_800;

            GETACT_280:

            --IC;
            if (IC > 0) goto GETACT_270;
//
//     Calculate the next VMU if VIOLMX is positive. Return if NACT=N holds,
//       as then the active constraints imply D=0. Otherwise, go to label
//       100, to calculate the new D and to test for termination.
//
            if (VIOLMX > ZERO) goto GETACT_220;
            if (NACT < N) goto GETACT_100;

            GETACT_290:

            DD = ZERO;

            GETACT_300:

            return DD;
//
//     These instructions rearrange the active constraints so that the new
//       value of IACT(NACT) is the old value of IACT(IC). A sequence of
//       Givens rotations is applied to the current QFAC and RFAC. Then NACT
//       is reduced by one.
//
            GETACT_800:

            RESNEW[IACT[IC]] = Math.Max(RESACT[IC], TINY);
            int JC = IC;
            while (JC < NACT)
            {
                int JCP = JC + 1;
                IDIAG = JC * JCP / 2;
                int JW = IDIAG + JCP;
                TEMP = Math.Sqrt(RFAC[JW - 1] * RFAC[JW - 1] + RFAC[JW] * RFAC[JW]);
                double CVAL = RFAC[JW] / TEMP;
                double SVAL = RFAC[JW - 1] / TEMP;
                RFAC[JW - 1] = SVAL * RFAC[IDIAG];
                RFAC[JW] = CVAL * RFAC[IDIAG];
                RFAC[IDIAG] = TEMP;
                if (JCP < NACT)
                {
                    for (int J = JCP + 1; J <= NACT; ++J)
                    {
                        TEMP = SVAL * RFAC[JW + JC] + CVAL * RFAC[JW + JCP];
                        RFAC[JW + JCP] = CVAL * RFAC[JW + JC] - SVAL * RFAC[JW + JCP];
                        RFAC[JW + JC] = TEMP;
                        JW += J;
                    }
                }
                int JDIAG = IDIAG - JC;
                for (int I = 1; I <= N; ++I)
                {
                    if (I < JC)
                    {
                        TEMP = RFAC[IDIAG + I];
                        RFAC[IDIAG + I] = RFAC[JDIAG + I];
                        RFAC[JDIAG + I] = TEMP;
                    }
                    TEMP = SVAL * QFAC[I, JC] + CVAL * QFAC[I, JCP];
                    QFAC[I, JCP] = CVAL * QFAC[I, JC] - SVAL * QFAC[I, JCP];
                    QFAC[I, JC] = TEMP;
                }
                IACT[JC] = IACT[JCP];
                RESACT[JC] = RESACT[JCP];
                VLAM[JC] = VLAM[JCP];
                JC = JCP;
            }
            --NACT;
            switch (IFLAG)
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

        private static void PRELIM(CalfunDel CALFUN, int N, int NPT, int M, double[,] AMAT, double[] B, double[] X,
            double RHOBEG, int IPRINT, double[] XBASE, double[,] XPT, double[] FVAL, double[] XSAV, double[] XOPT,
            double[] GOPT, out int KOPT, double[] HQ, double[] PQ, double[,] BMAT, double[,] ZMAT, out int IDZ, int NDIM,
            double[] SP, double[] RESCON, TextWriter logger)
        {
            double[] STEP = new double[1 + N];
            double[] PQW = new double[1 + NPT + N];
            double[] W = new double[1 + NPT + N];
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
            int NPTM = NPT - N - 1;
            double RHOSQ = RHOBEG * RHOBEG;
            double RECIP = ONE / RHOSQ;
            double RECIQ = Math.Sqrt(HALF) / RHOSQ;
            double TEST = 0.2 * RHOBEG;
            KOPT = 0;
            IDZ = 1;
            const int KBASE = 1;
//
//     Set the initial elements of XPT, BMAT, SP and ZMAT to zero. 
//
            for (int J = 1; J <= N; ++J)
            {
                XBASE[J] = X[J];
                for (int K = 1; K <= NPT; ++K)
                    XPT[K, J] = ZERO;
                for (int I = 1; I <= NDIM; ++I)
                    BMAT[I, J] = ZERO;
            }
            for (int K = 1; K <= NPT; ++K)
            {
                SP[K] = ZERO;
                for (int J = 1; J <= NPT - N - 1; ++J)
                    ZMAT[K, J] = ZERO;
            }
//
//     Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT],
//       but they may be altered later to make a constraint violation
//       sufficiently large. The initial nonzero elements of BMAT and of
//       the first min[N,NPT-N-1] columns of ZMAT are set also.
//
            for (int J = 1; J <= N; ++J)
            {
                XPT[J + 1, J] = RHOBEG;
                if (J < NPT - N)
                {
                    int JP = N + J + 1;
                    XPT[JP, J] = -RHOBEG;
                    BMAT[J + 1, J] = HALF / RHOBEG;
                    BMAT[JP, J] = -HALF / RHOBEG;
                    ZMAT[1, J] = -RECIQ - RECIQ;
                    ZMAT[J + 1, J] = RECIQ;
                    ZMAT[JP, J] = RECIQ;
                }
                else
                {
                    BMAT[1, J] = -ONE / RHOBEG;
                    BMAT[J + 1, J] = ONE / RHOBEG;
                    BMAT[NPT + J, J] = -HALF * RHOSQ;
                }
            }
//
//     Set the remaining initial nonzero elements of XPT and ZMAT when the
//       number of interpolation points exceeds 2*N+1.
//
            if (NPT > 2 * N + 1)
            {
                for (int K = N + 1; K <= NPT - N - 1; ++K)
                {
                    int ITEMP = (K - 1) / N;
                    int IPT = K - ITEMP * N;
                    int JPT = IPT + ITEMP;
                    if (JPT > N) JPT -= N;
                    XPT[N + K + 1, IPT] = RHOBEG;
                    XPT[N + K + 1, JPT] = RHOBEG;
                    ZMAT[1, K] = RECIP;
                    ZMAT[IPT + 1, K] = -RECIP;
                    ZMAT[JPT + 1, K] = -RECIP;
                    ZMAT[N + K + 1, K] = RECIP;
                }
            }
//
//     Update the constraint right hand sides to allow for the shift XBASE.
//
            if (M > 0)
            {
                for (int J = 1; J <= M; ++J)
                {
                    double TEMP = ZERO;
                    for (int I = 1; I <= N; ++I)
                        TEMP += AMAT[I, J] * XBASE[I];
                    B[J] -= TEMP;
                }
            }
//
//     Go through the initial points, shifting every infeasible point if
//       necessary so that its constraint violation is at least 0.2*RHOBEG.
//
            for (int NF = 1; NF <= NPT; ++NF)
            {
                double FEAS = ONE;
                double BIGV = ZERO;
                int JSAV = 0;
                {
                    int J = 0;

                    PRELIM_80:

                    ++J;
                    if (J <= M && NF >= 2)
                    {
                        double RESID = -B[J];
                        for (int I = 1; I <= N; ++I)
                            RESID = RESID + XPT[NF, I] * AMAT[I, J];
                        if (RESID <= BIGV) goto PRELIM_80;
                        BIGV = RESID;
                        JSAV = J;
                        if (RESID <= TEST)
                        {
                            FEAS = -ONE;
                            goto PRELIM_80;
                        }
                        FEAS = ZERO;
                    }
                }
                if (FEAS < ZERO)
                {
                    for (int I = 1; I <= N; ++I)
                        STEP[I] = XPT[NF, I] + (TEST - BIGV) * AMAT[I, JSAV];
                    for (int K = 1; K <= NPT; ++K)
                    {
                        SP[NPT + K] = ZERO;
                        for (int J = 1; J <= N; ++J)
                            SP[NPT + K] += XPT[K, J] * STEP[J];
                    }
                    UPDATE(N, NPT, XPT, BMAT, ZMAT, IDZ, NDIM, SP, STEP, KBASE, ref NF);
                    for (int I = 1; I <= N; ++I)
                        XPT[NF, I] = STEP[I];
                }
//
//     Calculate the objective function at the current interpolation point,
//       and set KOPT to the index of the first trust region centre.
//
                for (int J = 1; J <= N; ++J)
                    X[J] = XBASE[J] + XPT[NF, J];
                double F = FEAS;
                CALFUN(N, X, ref F);
                if (IPRINT == 3)
                {
                    PRINT(logger, PRELIM_140, NF, F, FORMAT("  ", "18:E6", X, 1, N));
                }
                if (NF == 1)
                {
                    KOPT = 1;
                }
                else if (F < FVAL[KOPT] && FEAS > ZERO)
                {
                    KOPT = NF;
                }
                FVAL[NF] = F;
            }
//
//     Set PQ for the first quadratic model.
//
            for (int J = 1; J <= NPTM; ++J)
            {
                W[J] = ZERO;
                for (int K = 1; K <= NPT; ++K)
                    W[J] += ZMAT[K, J] * FVAL[K];
            }
            for (int K = 1; K <= NPT; ++K)
            {
                PQ[K] = ZERO;
                for (int J = 1; J <= NPTM; ++J)
                    PQ[K] += ZMAT[K, J] * W[J];
            }
//
//     Set XOPT, SP, GOPT and HQ for the first quadratic model.
//
            for (int J = 1; J <= N; ++J)
            {
                XOPT[J] = XPT[KOPT, J];
                XSAV[J] = XBASE[J] + XOPT[J];
                GOPT[J] = ZERO;
            }
            for (int K = 1; K <= NPT; ++K)
            {
                SP[K] = ZERO;
                for (int J = 1; J <= N; ++J)
                    SP[K] += XPT[K, J] * XOPT[J];
                double TEMP = PQ[K] * SP[K];
                for (int J = 1; J <= N; ++J)
                    GOPT[J] += FVAL[K] * BMAT[K, J] + TEMP * XPT[K, J];
            }
            for (int I = 1; I <= (N * N + N) / 2; ++I)
                HQ[I] = ZERO;
//
//     Set the initial elements of RESCON.
//
            for (int J = 1; J <= M; ++J)
            {
                double TEMP = B[J];
                for (int I = 1; I <= N; ++I)
                    TEMP -= XOPT[I] * AMAT[I, J];
                TEMP = Math.Max(TEMP, ZERO);
                if (TEMP >= RHOBEG) TEMP = -TEMP;
                RESCON[J] = TEMP;
            }
        }

        private static void QMSTEP(int N, int NPT, int M, double[,] AMAT, double[] B, double[,] XPT, double[] XOPT,
            int NACT, int[] IACT, double[] RESCON, double[,] QFAC, int KOPT, int KNEW, double DEL, double[] STEP,
            double[] GL, double[] PQW, ref int IFEAS)
        {
            double[] RSTAT = new double[1 + M]; 
            double[] W = new double[1 + N]; 
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
            double TEST = 0.2 * DEL;
//
//     Replace GL by the gradient of LFUNC at the trust region centre, and
//       set the elements of RSTAT.
//
            for (int K = 1; K <= NPT; ++K)
            {
                double TEMP = ZERO;
                for (int J = 1; J <= N; ++J)
                    TEMP += XPT[K, J] * XOPT[J];
                TEMP *= PQW[K];
                for (int I = 1; I <= N; ++I)
                    GL[I] = GL[I] + TEMP * XPT[K, I];
            }
            if (M > 0)
            {
                for (int J = 1; J <= M; ++J)
                {
                    RSTAT[J] = ONE;
                    if (Math.Abs(RESCON[J]) >= DEL) RSTAT[J] = -ONE;
                }
                for (int K = 1; K <= NACT; ++K)
                    RSTAT[IACT[K]] = ZERO;
            }
//
//     Find the greatest modulus of LFUNC on a line through XOPT and
//       another interpolation point within the trust region.
//
            int KSAV = 0;
            double STPSAV = ZERO;
            double VBIG = ZERO;
            for (int K = 1; K <= NPT; ++K)
            {
                if (K == KOPT) continue;
                double SS = ZERO;
                double SP = ZERO;
                for (int I = 1; I <= N; ++I)
                {
                    double TEMP = XPT[K, I] - XOPT[I];
                    SS += TEMP * TEMP;
                    SP += GL[I] * TEMP;
                }
                double STP = -DEL / Math.Sqrt(SS);
                double VLAG;
                if (K == KNEW)
                {
                    if (SP * (SP - ONE) < ZERO) STP = -STP;
                    VLAG = Math.Abs(STP * SP) + STP * STP * Math.Abs(SP - ONE);
                }
                else
                {
                    VLAG = Math.Abs(STP * (ONE - STP) * SP);
                }
                if (VLAG > VBIG)
                {
                    KSAV = K;
                    STPSAV = STP;
                    VBIG = VLAG;
                }
            }
//
//     Set STEP to the move that gives the greatest modulus calculated above.
//       This move may be replaced by a steepest ascent step from XOPT.
//
            double VGRAD = ZERO;
            double GG = ZERO;
            for (int I = 1; I <= N; ++I)
            {
                GG += GL[I] * GL[I];
                STEP[I] = STPSAV * (XPT[KSAV, I] - XOPT[I]);
            }
            VGRAD = DEL * Math.Sqrt(GG);
            if (VGRAD <= TENTH * VBIG) goto QMSTEP_220;
//
//     Make the replacement if it provides a larger value of VBIG.
//
            double GHG = ZERO;
            for (int K = 1; K <= NPT; ++K)
            {
                double TEMP = ZERO;
                for (int J = 1; J <= N; ++J)
                    TEMP += XPT[K, J] * GL[J];
                GHG += PQW[K] * TEMP * TEMP;
            }
            double VNEW = VGRAD + Math.Abs(HALF * DEL * DEL * GHG / GG);
            if (VNEW > VBIG)
            {
                VBIG = VNEW;
                double STP = DEL / Math.Sqrt(GG);
                if (GHG < ZERO) STP = -STP;
                for (int I = 1; I <= N; ++I)
                    STEP[I] = STP * GL[I];
            }
            if (NACT == 0 || NACT == N) goto QMSTEP_220;
//
//     Overwrite GL by its projection. Then set VNEW to the greatest
//       value of |LFUNC| on the projected gradient from XOPT subject to
//       the trust region bound. If VNEW is sufficiently large, then STEP
//       may be changed to a move along the projected gradient.
//
            for (int K = NACT + 1; K <= N; ++K)
            {
                W[K] = ZERO;
                for (int I = 1; I <= N; ++I)
                    W[K] += GL[I] * QFAC[I, K];
            }
            GG = ZERO;
            for (int I = 1; I <= N; ++I)
            {
                GL[I] = ZERO;
                for (int K = NACT + 1; K <= N; ++K)
                    GL[I] += QFAC[I, K] * W[K];
                GG += GL[I] * GL[I];
            }
            VGRAD = DEL * Math.Sqrt(GG);
            if (VGRAD <= TENTH * VBIG) goto QMSTEP_220;
            GHG = ZERO;
            for (int K = 1; K <= NPT; ++K)
            {
                double TEMP = ZERO;
                for (int J = 1; J <= N; ++J)
                    TEMP += XPT[K, J] * GL[J];
                GHG += PQW[K] * TEMP * TEMP;
            }
            VNEW = VGRAD + Math.Abs(HALF * DEL * DEL * GHG / GG);
//
//     Set W to the possible move along the projected gradient.
//
            double WW = ZERO;
            {
                double STP = DEL / Math.Sqrt(GG);
                if (GHG < ZERO) STP = -STP;
                for (int I = 1; I <= N; ++I)
                {
                    W[I] = STP * GL[I];
                    WW += W[I] * W[I];
                }
            }
//
//     Set STEP to W if W gives a sufficiently large value of the modulus
//       of the Lagrange function, and if W either preserves feasibility
//       or gives a constraint violation of at least 0.2*DEL. The purpose
//       of CTOL below is to provide a check on feasibility that includes
//       a tolerance for contributions from computer rounding errors.
//
            if (VNEW / VBIG >= 0.2)
            {
                IFEAS = 1;
                double BIGV = ZERO;
                int J = 0;

                QMSTEP_170:

                J = J + 1;
                if (J <= M)
                {
                    if (RSTAT[J] == ONE)
                    {
                        double TEMP = -RESCON[J];
                        for (int I = 1; I <= N; ++I)
                            TEMP += W[I] * AMAT[I, J];
                        BIGV = Math.Max(BIGV, TEMP);
                    }
                    if (BIGV < TEST) goto QMSTEP_170;
                    IFEAS = 0;
                }
                double CTOL = ZERO;
                {
                    double TEMP = 0.01 * Math.Sqrt(WW);
                    if (BIGV > ZERO && BIGV < TEMP)
                    {
                        for (int K = 1; K <= NACT; ++K)
                        {
                            J = IACT[K];
                            double SUM = ZERO;
                            for (int I = 1; I <= N; ++I)
                                SUM = SUM + W[I] * AMAT[I, J];
                            CTOL = Math.Max(CTOL, Math.Abs(SUM));
                        }
                    }
                }
                if (BIGV <= 10.0 * CTOL || BIGV >= TEST)
                {
                    for (int I = 1; I <= N; ++I)
                        STEP[I] = W[I];
                    return;
                }
            }
//
//     Calculate the greatest constraint violation at XOPT+STEP with STEP at
//       its original value. Modify STEP if this violation is unacceptable.
//
            QMSTEP_220:

            IFEAS = 1;
            int JSAV = 0;
            double RESMAX = ZERO;
            {
                int J = 0;
                double BIGV = ZERO;

                QMSTEP_230:

                ++J;
                if (J <= M)
                {
                    if (RSTAT[J] < ZERO) goto QMSTEP_230;
                    double TEMP = -RESCON[J];
                    for (int I = 1; I <= N; ++I)
                        TEMP += STEP[I] * AMAT[I, J];
                    RESMAX = Math.Max(RESMAX, TEMP);
                    if (TEMP < TEST)
                    {
                        if (TEMP <= BIGV) goto QMSTEP_230;
                        BIGV = TEMP;
                        JSAV = J;
                        IFEAS = -1;
                        goto QMSTEP_230;
                    }
                    IFEAS = 0;
                }
                if (IFEAS == -1)
                {
                    for (int I = 1; I <= N; ++I)
                        STEP[I] += (TEST - BIGV) * AMAT[I, JSAV];
                    IFEAS = 0;
                }
            }
//
//     Return the calculated STEP and the value of IFEAS.
//
        }

        private static void TRSTEP(int N, int NPT, int M, double[,] AMAT, double[] B, double[,] XPT, double[] HQ,
            double[] PQ, ref int NACT, int[] IACT, double[] RESCON, double[,] QFAC, double[] RFAC, ref double SNORM,
            double[] STEP, double[] G)
        {
            double[] RESNEW = new double[1 + M];
            double[] RESACT = new double[1 + N];
            double[] D = new double[1 + N];
            double[] DW = new double[1 + N];
            double[] W = new double[1 + Math.Max(M, 2 * N)];
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
            double SNSQ = SNORM * SNORM;
//
//     Set the initial elements of RESNEW, RESACT and STEP.
//
            if (M > 0)
            {
                for (int J = 1; J <= M; ++J)
                {
                    RESNEW[J] = RESCON[J];
                    if (RESCON[J] >= SNORM)
                    {
                        RESNEW[J] = -ONE;
                    }
                    else if (RESCON[J] >= ZERO)
                    {
                        RESNEW[J] = Math.Max(RESNEW[J], TINY);
                    }
                }
                if (NACT > 0)
                {
                    for (int K = 1; K <= NACT; ++K)
                    {
                        RESACT[K] = RESCON[IACT[K]];
                        RESNEW[IACT[K]] = ZERO;
                    }
                }
            }
            for (int I = 1; I <= N; ++I)
                STEP[I] = ZERO;
            double SS = ZERO;
            double REDUCT = ZERO;
            int NCALL = 0;
//
//     GETACT picks the active set for the current STEP. It also sets DW to
//       the vector closest to -G that is orthogonal to the normals of the
//       active constraints. DW is scaled to have length 0.2*SNORM, as then
//       a move of DW from STEP is allowed by the linear constraints.
//
            TRSTEP_40:

            ++NCALL;
            double DSQ = GETACT(N, M, AMAT, B, ref NACT, IACT, QFAC, RFAC, SNORM, RESNEW, RESACT, G, DW);
            if (DSQ == ZERO) goto TRSTEP_320;
            double SCALE = 0.2 * SNORM / Math.Sqrt(DSQ);
            for (int I = 1; I <= N; ++I)
                DW[I] *= SCALE;
//
//     If the modulus of the residual of an active constraint is substantial,
//       then set D to the shortest move from STEP to the boundaries of the
//       active constraints.
//
            double RESMAX = ZERO;
            if (NACT > 0)
            {
                for (int K = 1; K <= NACT; ++K)
                    RESMAX = Math.Max(RESMAX, RESACT[K]);
            }
            double GAMMA = ZERO;
            if (RESMAX > 1.0E-4 * SNORM)
            {
                int IR = 0;
                for (int K = 1; K <= NACT; ++K)
                {
                    double TEMP = RESACT[K];
                    if (K >= 2)
                    {
                        for (int I = 1; I <= K - 1; ++I)
                        {
                            ++IR;
                            TEMP -= RFAC[IR] * W[I];
                        }
                    }
                    ++IR;
                    W[K] = TEMP / RFAC[IR];
                }
                for (int I = 1; I <= N; ++I)
                {
                    D[I] = ZERO;
                    for (int K = 1; K <= NACT; ++K)
                        D[I] += W[K] * QFAC[I, K];
                }
//
//     The vector D that has just been calculated is also the shortest move
//       from STEP+DW to the boundaries of the active constraints. Set GAMMA
//       to the greatest steplength of this move that satisfies the trust
//       region bound.
//
                double RHS = SNSQ;
                double DS = ZERO;
                double DD = ZERO;
                for (int I = 1; I <= N; ++I)
                {
                    double SUM = STEP[I] + DW[I];
                    RHS -= SUM * SUM;
                    DS += D[I] * SUM;
                    DD += D[I] * D[I];
                }
                if (RHS > ZERO)
                {
                    double TEMP = Math.Sqrt(DS * DS + DD * RHS);
                    if (DS <= ZERO)
                    {
                        GAMMA = (TEMP - DS) / DD;
                    }
                    else
                    {
                        GAMMA = RHS / (TEMP + DS);
                    }
                }
//
//     Reduce the steplength GAMMA if necessary so that the move along D
//       also satisfies the linear constraints.
//
                {
                    int J = 0;

                    TRSTEP_110:

                    if (GAMMA > ZERO)
                    {
                        ++J;
                        if (RESNEW[J] > ZERO)
                        {
                            double AD = ZERO;
                            double ADW = ZERO;
                            for (int I = 1; I <= N; ++I)
                            {
                                AD += AMAT[I, J] * D[I];
                                ADW += AMAT[I, J] * DW[I];
                            }
                            if (AD > ZERO)
                            {
                                double TEMP = Math.Max((RESNEW[J] - ADW) / AD, ZERO);
                                GAMMA = Math.Min(GAMMA, TEMP);
                            }
                        }
                        if (J < M) goto TRSTEP_110;
                    }
                }
                GAMMA = Math.Min(GAMMA, ONE);
            }
//
//     Set the next direction for seeking a reduction in the model function
//       subject to the trust region bound and the linear constraints.
//
            int ICOUNT;
            if (GAMMA <= ZERO)
            {
                for (int I = 1; I <= N; ++I)
                    D[I] = DW[I];
                ICOUNT = NACT;
            }
            else
            {
                for (int I = 1; I <= N; ++I)
                    D[I] = DW[I] + GAMMA * D[I];
                ICOUNT = NACT - 1;
            }
            double ALPBD = ONE;
//
//     Set ALPHA to the steplength from STEP along D to the trust region
//       boundary. Return if the first derivative term of this step is
//       sufficiently small or if no further progress is possible.
//
            TRSTEP_150:

            ++ICOUNT;
            double ALPHA;
            double DG = ZERO;
            {
                double RHS = SNSQ - SS;
                if (RHS <= ZERO) goto TRSTEP_320;
                double DS = ZERO;
                double DD = ZERO;
                for (int I = 1; I <= N; ++I)
                {
                    DG += D[I] * G[I];
                    DS += D[I] * STEP[I];
                    DD += D[I] * D[I];
                }
                if (DG >= ZERO) goto TRSTEP_320;
                double TEMP = Math.Sqrt(RHS * DD + DS * DS);
                if (DS <= ZERO)
                {
                    ALPHA = (TEMP - DS) / DD;
                }
                else
                {
                    ALPHA = RHS / (TEMP + DS);
                }
                if (-ALPHA * DG <= CTEST * REDUCT) goto TRSTEP_320;
            }
//
//     Set DW to the change in gradient along D.
//
            int IH = 0;
            for (int J = 1; J <= N; ++J)
            {
                DW[J] = ZERO;
                for (int I = 1; I <= J; ++I)
                {
                    ++IH;
                    if (I < J) DW[J] += HQ[IH] * D[I];
                    DW[I] += HQ[IH] * D[J];
                }
            }
            for (int K = 1; K <= NPT; ++K)
            {
                double TEMP = ZERO;
                for (int J = 1; J <= N; ++J)
                    TEMP += XPT[K, J] * D[J];
                TEMP *= PQ[K];
                for (int I = 1; I <= N; ++I)
                    DW[I] += TEMP * XPT[K, I];
            }
//
//     Set DGD to the curvature of the model along D. Then reduce ALPHA if
//       necessary to the value that minimizes the model.
//
            double DGD = ZERO;
            for (int I = 1; I <= N; ++I)
                DGD += D[I] * DW[I];
            double ALPHT = ALPHA;
            if (DG + ALPHA * DGD > ZERO)
            {
                ALPHA = -DG / DGD;
            }
//
//     Make a further reduction in ALPHA if necessary to preserve feasibility,
//       and put some scalar products of D with constraint gradients in W.
//
            double ALPHM = ALPHA;
            int JSAV = 0;
            if (M > 0)
            {
                for (int J = 1; J <= M; ++J)
                {
                    double AD = ZERO;
                    if (RESNEW[J] > ZERO)
                    {
                        for (int I = 1; I <= N; ++I)
                            AD += AMAT[I, J] * D[I];
                        if (ALPHA * AD > RESNEW[J])
                        {
                            ALPHA = RESNEW[J] / AD;
                            JSAV = J;
                        }
                    }
                    W[J] = AD;
                }
            }
            ALPHA = Math.Max(ALPHA, ALPBD);
            ALPHA = Math.Min(ALPHA, ALPHM);
            if (ICOUNT == NACT) ALPHA = Math.Min(ALPHA, ONE);
//
//     Update STEP, G, RESNEW, RESACT and REDUCT.
//
            SS = ZERO;
            for (int I = 1; I <= N; ++I)
            {
                STEP[I] += ALPHA * D[I];
                SS += STEP[I] * STEP[I];
                G[I] += ALPHA * DW[I];
            }
            if (M > 0)
            {
                for (int J = 1; J <= M; ++J)
                {
                    if (RESNEW[J] > ZERO)
                    {
                        RESNEW[J] = Math.Max(RESNEW[J] - ALPHA * W[J], TINY);
                    }
                }
            }
            if (ICOUNT == NACT && NACT > 0)
            {
                for (int K = 1; K <= NACT; ++K)
                    RESACT[K] *= (ONE - GAMMA);
            }
            REDUCT -= ALPHA * (DG + HALF * ALPHA * DGD);
//
//     Test for termination. Branch to label 40 if there is a new active
//       constraint and if the distance from STEP to the trust region
//       boundary is at least 0.2*SNORM.
//
            if (ALPHA == ALPHT) goto TRSTEP_320;
            {
                double TEMP = -ALPHM * (DG + HALF * ALPHM * DGD);
                if (TEMP <= CTEST * REDUCT) goto TRSTEP_320;
            }
            if (JSAV > 0)
            {
                if (SS <= 0.64 * SNSQ) goto TRSTEP_40;
                goto TRSTEP_320;
            }
            if (ICOUNT == N) goto TRSTEP_320;
//
//     Calculate the next search direction, which is conjugate to the
//       previous one except in the case ICOUNT=NACT.
//
            if (NACT > 0)
            {
                for (int J = NACT + 1; J <= N; ++J)
                {
                    W[J] = ZERO;
                    for (int I = 1; I <= N; ++I)
                        W[J] += G[I] * QFAC[I, J];
                }
                for (int I = 1; I <= N; ++I)
                {
                    double TEMP = ZERO;
                    for (int J = NACT + 1; J <= N; ++J)
                        TEMP += QFAC[I, J] * W[J];
                    W[N + I] = TEMP;
                }
            }
            else
            {
                for (int I = 1; I <= N; ++I)
                    W[N + I] = G[I];
            }
            double BETA;
            if (ICOUNT == NACT)
            {
                BETA = ZERO;
            }
            else
            {
                double WGD = ZERO;
                for (int I = 1; I <= N; ++I)
                    WGD += W[N + I] * DW[I];
                BETA = WGD / DGD;
            }
            for (int I = 1; I <= N; ++I)
                D[I] = -W[N + I] + BETA * D[I];
            ALPBD = ZERO;
            goto TRSTEP_150;
//
//     Return from the subroutine.
//
            TRSTEP_320:

            SNORM = ZERO;
            if (REDUCT > ZERO) SNORM = Math.Sqrt(SS);
            G[1] = ZERO;
            if (NCALL > 1) G[1] = ONE;
        }

        private static void UPDATE(int N, int NPT, double[,] XPT, double[,] BMAT, double[,] ZMAT, int IDZ, int NDIM,
            double[] SP, double[] STEP, int KOPT, ref int KNEW)
        {
            double[] VLAG = new double[1 + NPT + N];
            double[] W = new double[1 + NPT + N];
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
            double NPTM = NPT - N - 1;
//
//     Calculate VLAG and BETA for the current choice of STEP. The first NPT
//       elements of VLAG are set to the values of the Lagrange functions at
//       XPT(KOPT,.)+STEP(.). The first NPT components of W_check are held
//       in W, where W_check is defined in a paper on the updating method.
//
            for (int K = 1; K <= NPT; ++K)
            {
                W[K] = SP[NPT + K] * (HALF * SP[NPT + K] + SP[K]);
                double SUM = ZERO;
                for (int J = 1; J <= N; ++J)
                    SUM += BMAT[K, J] * STEP[J];
                VLAG[K] = SUM;
            }
            double BETA = ZERO;
            for (int K = 1; K <= NPTM; ++K)
            {
                double SUM = ZERO;
                for (int I = 1; I <= NPT; ++I)
                    SUM = SUM + ZMAT[I, K] * W[I];
                if (K < IDZ)
                {
                    BETA += SUM * SUM;
                    SUM = -SUM;
                }
                else
                {
                    BETA -= SUM * SUM;
                }
                for (int I = 1; I <= NPT; ++I)
                    VLAG[I] = VLAG[I] + SUM * ZMAT[I, K];
            }
            double BSUM = ZERO;
            double DX = ZERO;
            double SSQ = ZERO;
            for (int J = 1; J <= N; ++J)
            {
                double SUM = ZERO;
                for (int I = 1; I <= NPT; ++I)
                    SUM += W[I] * BMAT[I, J];
                BSUM += SUM * STEP[J];
                int JP = NPT + J;
                for (int K = 1; K <= N; ++K)
                    SUM += BMAT[JP, K] * STEP[K];
                VLAG[JP] = SUM;
                BSUM += SUM * STEP[J];
                DX += STEP[J] * XPT[KOPT, J];
                SSQ += STEP[J] * STEP[J];
            }
            BETA = DX * DX + SSQ * (SP[KOPT] + DX + DX + HALF * SSQ) + BETA - BSUM;
            VLAG[KOPT] += +ONE;
//
//     If KNEW is zero initially, then pick the index of the interpolation
//       point to be deleted, by maximizing the absolute value of the
//       denominator of the updating formula times a weighting factor.
//       
//
            if (KNEW == 0)
            {
                double DENMAX = ZERO;
                for (int K = 1; K <= NPT; ++K)
                {
                    double HDIAG = ZERO;
                    for (int J = 1; J <= NPTM; ++J)
                    {
                        double TEMP = ONE;
                        if (J < IDZ) TEMP = -ONE;
                        HDIAG += TEMP * ZMAT[K, J] * ZMAT[K, J];
                    }
                    double DENABS = Math.Abs(BETA * HDIAG + VLAG[K] * VLAG[K]);
                    double DISTSQ = ZERO;
                    for (int J = 1; J <= N; ++J)
                        DISTSQ += Math.Pow(XPT[K, J] - XPT[KOPT, J], 2.0);
                    {
                        double TEMP = DENABS * DISTSQ * DISTSQ;
                        if (TEMP > DENMAX)
                        {
                            DENMAX = TEMP;
                            KNEW = K;
                        }
                    }
                }
            }
//
//     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
//
            int JL = 1;
            double TEMPA, TEMPB = ZERO;
            if (NPTM >= 2)
            {
                for (int J = 2; J <= NPTM; ++J)
                {
                    if (J == IDZ)
                    {
                        JL = IDZ;
                    }
                    else if (ZMAT[KNEW, J] != ZERO)
                    {
                        double TEMP = Math.Sqrt(ZMAT[KNEW, JL] * ZMAT[KNEW, JL] + ZMAT[KNEW, J] * ZMAT[KNEW, J]);
                        TEMPA = ZMAT[KNEW, JL] / TEMP;
                        TEMPB = ZMAT[KNEW, J] / TEMP;
                        for (int I = 1; I <= NPT; ++I)
                        {
                            TEMP = TEMPA * ZMAT[I, JL] + TEMPB * ZMAT[I, J];
                            ZMAT[I, J] = TEMPA * ZMAT[I, J] - TEMPB * ZMAT[I, JL];
                            ZMAT[I, JL] = TEMP;
                        }
                        ZMAT[KNEW, J] = ZERO;
                    }
                }
            }
//
//     Put the first NPT components of the KNEW-th column of the Z Z^T matrix
//       into W, and calculate the parameters of the updating formula.
//
            TEMPA = ZMAT[KNEW, 1];
            if (IDZ >= 2) TEMPA = -TEMPA;
            if (JL > 1) TEMPB = ZMAT[KNEW, JL];
            for (int I = 1; I <= NPT; ++I)
            {
                W[I] = TEMPA * ZMAT[I, 1];
                if (JL > 1) W[I] += TEMPB * ZMAT[I, JL];
            }
            double ALPHA = W[KNEW];
            double TAU = VLAG[KNEW];
            double TAUSQ = TAU * TAU;
            double DENOM = ALPHA * BETA + TAUSQ;
            VLAG[KNEW] -= ONE;
            if (DENOM == ZERO)
            {
                KNEW = 0;
                return;
            }
            double SQRTDN = Math.Sqrt(Math.Abs(DENOM));
//
//     Complete the updating of ZMAT when there is only one nonzero element
//       in the KNEW-th row of the new matrix ZMAT. IFLAG is set to one when
//       the value of IDZ is going to be reduced.
//
            int IFLAG = 0;
            if (JL == 1)
            {
                TEMPA = TAU / SQRTDN;
                TEMPB = ZMAT[KNEW, 1] / SQRTDN;
                for (int I = 1; I <= NPT; ++I)
                    ZMAT[I, 1] = TEMPA * ZMAT[I, 1] - TEMPB * VLAG[I];
                if (DENOM < ZERO)
                {
                    if (IDZ == 1)
                    {
                        IDZ = 2;
                    }
                    else
                    {
                        IFLAG = 1;
                    }
                }
            }
            else
            {
//
//     Complete the updating of ZMAT in the alternative case.
//
                int JA = 1;
                if (BETA >= ZERO) JA = JL;
                int JB = JL + 1 - JA;
                double TEMP = ZMAT[KNEW, JB] / DENOM;
                TEMPA = TEMP * BETA;
                TEMPB = TEMP * TAU;
                TEMP = ZMAT[KNEW, JA];
                double SCALA = ONE / Math.Sqrt(Math.Abs(BETA) * TEMP * TEMP + TAUSQ);
                double SCALB = SCALA * SQRTDN;
                for (int I = 1; I <= NPT; ++I)
                {
                    ZMAT[I, JA] = SCALA * (TAU * ZMAT[I, JA] - TEMP * VLAG[I]);
                    ZMAT[I, JB] = SCALB * (ZMAT[I, JB] - TEMPA * W[I] - TEMPB * VLAG[I]);
                }
                if (DENOM <= ZERO)
                {
                    if (BETA < ZERO)
                    {
                        IDZ = IDZ + 1;
                    }
                    else
                    {
                        IFLAG = 1;
                    }
                }
            }
//
//     Reduce IDZ when the diagonal part of the ZMAT times Diag(DZ) times
//       ZMAT^T factorization gains another positive element. Then exchange
//       the first and IDZ-th columns of ZMAT.
//
            if (IFLAG == 1)
            {
                --IDZ;
                for (int I = 1; I <= NPT; ++I)
                {
                    double TEMP = ZMAT[I, 1];
                    ZMAT[I, 1] = ZMAT[I, IDZ];
                    ZMAT[I, IDZ] = TEMP;
                }
            }
//
//     Finally, update the matrix BMAT.
//
            for (int J = 1; J <= N; ++J)
            {
                int JP = NPT + J;
                W[JP] = BMAT[KNEW, J];
                TEMPA = (ALPHA * VLAG[JP] - TAU * W[JP]) / DENOM;
                TEMPB = (-BETA * W[JP] - TAU * VLAG[JP]) / DENOM;
                for (int I = 1; I <= JP; ++I)
                {
                    BMAT[I, J] = BMAT[I, J] + TEMPA * VLAG[I] + TEMPB * W[I];
                    if (I > NPT) BMAT[JP, I - NPT] = BMAT[I, J];
                }
            }
        }

        private static void PRINT(TextWriter logger, string format, params object[] args)
        {
            if (logger != null) logger.WriteLine(format, args);
        }

        private static string FORMAT<T>(string separator, string itemFormatter, IEnumerable<T> items, int start, int end)
        {
            return String.Join(separator,
                items.Skip(start).Take(end).Select(item => String.Format("{0," + itemFormatter + "}", item)));
        }

        #endregion
    }
}