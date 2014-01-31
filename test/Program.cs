using System;
using Cureos.Numerics;

namespace test
{
    class Program
    {
        #region FIELDS

        private const double ONE = 1.0;
        private const double TWO = 2.0;
        private const double ZERO = 0.0;
        private static readonly double PI = Math.PI;

        private static double FMAX;

        #endregion

        public static void Main(string[] args)
        {
            //
            //     Calculate the tetrahedron of least volume that encloses the points
            //       (XP(J),YP(J),ZP(J)), J=1,2,...,NP. Our method requires the origin
            //       to be strictly inside the convex hull of these points. There are
            //       twelve variables that define the four faces of each tetrahedron
            //       that is considered. Each face has the form ALPHA*X+BETA*Y+GAMMA*Z=1,
            //       the variables X(3K-2), X(3K-1) and X(3K) being the values of ALPHA,
            //       BETA and GAMMA for the K-th face, K=1,2,3,4. Let the set T contain
            //       all points in three dimensions that can be reached from the origin
            //       without crossing a face. Because the volume of T may be infinite,
            //       the objective function is the smaller of FMAX and the volume of T,
            //       where FMAX is set to an upper bound on the final volume initially.
            //       There are 4*NP linear constraints on the variables, namely that each
            //       of the given points (XP(J),YP(J),ZP(J)) shall be in T. Let XS = min
            //       XP(J), YS = min YP(J), ZS = min ZP(J) and SS = max XP(J)+YP(J)+ZP(J),
            //       where J runs from 1 to NP. The initial values of the variables are
            //       X(1)=1/XS, X(5)=1/YS, X(9)=1/ZS, X(2)=X(3)=X(4)=X(6)=X(7) =X(8)=0
            //       and X(10)=X(11)=X(12)=1/SS, which satisfy the linear constraints,
            //       and which provide the bound FMAX=(SS-XS-YS-ZS)**3/6. Other details
            //       of the test calculation are given below, including the choice of
            //       the data points (XP(J),YP(J),ZP(J)), J=1,2,...,NP. The smaller final
            //       value of the objective function in the case NPT=35 shows that the
            //       problem has local minima.
            //
            //
            //     Set some constants.
            //
            int IA = 12;
            int N = 12;
            //
            //     Set the data points.
            //
            int NP = 50;
            double SUMX = ZERO;
            double SUMY = ZERO;
            double SUMZ = ZERO;

            double[] XP = new double[1 + 50];
            double[] YP = new double[1 + 50];
            double[] ZP = new double[1 + 50];
            for (int J = 1; J <= NP; ++J)
            {
                double THETA = (double)(J - 1) * PI / (double)(NP - 1);
                XP[J] = Math.Cos(THETA) * Math.Cos(TWO * THETA);
                SUMX += XP[J];
                YP[J] = Math.Sin(THETA) * Math.Cos(TWO * THETA);
                SUMY += YP[J];
                ZP[J] = Math.Sin(TWO * THETA);
                SUMZ += ZP[J];
            }
            SUMX /= (double)(NP);
            SUMY /= (double)(NP);
            SUMZ /= (double)(NP);
            for (int J = 1; J <= NP; ++J)
            {
                XP[J] -= SUMX;
                YP[J] -= SUMY;
                ZP[J] -= SUMZ;
            }
            //
            //     Set the linear constraints.
            //
            double[,] A = new double[1 + 12, 1 + 200];
            double[] B = new double[1 + 200];
            int M = 4 * NP;
            for (int K = 1; K <= M; ++K)
            {
                B[K] = ONE;
                for (int I = 1; I <= N; ++I)
                    A[I, K] = ZERO;
            }
            for (int J = 1; J <= NP; ++J)
            {
                for (int I = 1; I <= 4; ++I)
                {
                    int K = 4 * J + I - 4;
                    int IW = 3 * I;
                    A[IW - 2, K] = XP[J];
                    A[IW - 1, K] = YP[J];
                    A[IW, K] = ZP[J];
                }
            }
            //
            //     Set the initial vector of variables. The JCASE=1,6 loop gives six
            //       different choices of NPT when LINCOA is called.
            //
            double[] X = new double[1 + 12];
            double XS = ZERO;
            double YS = ZERO;
            double ZS = ZERO;
            double SS = ZERO;
            for (int J = 1; J <= NP; ++J)
            {
                XS = Math.Min(XS, XP[J]);
                YS = Math.Min(YS, YP[J]);
                ZS = Math.Min(ZS, ZP[J]);
                SS = Math.Max(SS, XP[J] + YP[J] + ZP[J]);
            }
            FMAX = Math.Pow(SS - XS - YS - ZS, 3.0) / 6.0;
            for (int JCASE = 1; JCASE <= 6; ++JCASE)
            {
                for (int I = 2; I <= 8; ++I)
                    X[I] = ZERO;
                X[1] = ONE / XS;
                X[5] = ONE / YS;
                X[9] = ONE / ZS;
                X[10] = ONE / SS;
                X[11] = ONE / SS;
                X[12] = ONE / SS;
                //
                //     Call of LINCOA, which provides the printing given at the end of this
                //       note.
                //
                int NPT = 5 * JCASE + 10;
                double RHOBEG = 1.0;
                double RHOEND = 1.0E-6;
                int IPRINT = 1;
                int MAXFUN = 10000;
                Console.WriteLine("Output from LINCOA with  NPT ={0,4:D}  and  RHOEND ={1,12:E4}", NPT, RHOEND);
                Lincoa.LINCOA(CALFUN, N, NPT, M, A, IA, B, X, RHOBEG, RHOEND, IPRINT, MAXFUN, Console.Out);
                Console.WriteLine();
            }
            Console.ReadLine();
        }

        public static void CALFUN(int N, double[] X, ref double F)
        {
            F = FMAX;
            double V12 = X[1] * X[5] - X[4] * X[2];
            double V13 = X[1] * X[8] - X[7] * X[2];
            double V14 = X[1] * X[11] - X[10] * X[2];
            double V23 = X[4] * X[8] - X[7] * X[5];
            double V24 = X[4] * X[11] - X[10] * X[5];
            double V34 = X[7] * X[11] - X[10] * X[8];
            double DEL1 = V23 * X[12] - V24 * X[9] + V34 * X[6];
            if (DEL1 <= ZERO) return;
            double DEL2 = -V34 * X[3] - V13 * X[12] + V14 * X[9];
            if (DEL2 <= ZERO) return;
            double DEL3 = -V14 * X[6] + V24 * X[3] + V12 * X[12];
            if (DEL3 <= ZERO) return;
            double DEL4 = -V12 * X[9] + V13 * X[6] - V23 * X[3];
            if (DEL4 <= ZERO) return;
            double TEMP = Math.Pow(DEL1 + DEL2 + DEL3 + DEL4, 3.0) / (DEL1 * DEL2 * DEL3 * DEL4);
            F = Math.Min(TEMP / 6.0, FMAX);
        }
    }
}
