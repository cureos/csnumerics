/*
 *  Copyright (c) 2012-2014, Cureos AB.
 *  All rights reserved.
 *  http://www.cureos.com
 *
 *	This file is part of CSNumerics.
 *
 *  CSNumerics is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  CSNumerics is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with CSNumerics.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

using System;
using System.Collections;
using System.IO;
using System.Linq;
using NUnit.Framework;

namespace Cureos.Numerics.Optimizers
{
    [TestFixture]
    public class BobyqaTests
    {
        #region Fields

        private const double TOL = 1.0e-5;

        #endregion

        #region Unit tests

        [Test]
        public void FindMinimum_IntermedRosenbrock_ReturnsValidMinimum()
        {
            var xx = new[] { 1.0, -1.0 };
            Bobyqa.FindMinimum((n, x) => 10.0 * Math.Pow(x[0] * x[0] - x[1], 2.0) + Math.Pow(1.0 + x[0], 2.0), 2, xx);
            Assert.AreEqual(-1.0, xx[0], TOL);
            Assert.AreEqual(1.0, xx[1], TOL);
        }

        [Test]
        public void FindMinimum_HS04_ReturnsValidMinimum()
        {
            var xx = new[] { 1.125, 0.125 };
            Bobyqa.FindMinimum((n, x) => Math.Pow(x[0] + 1.0, 3.0) / 3.0 + x[1], 2, xx, new[] { 1.0, 0.0 }, null, 4);
            Assert.AreEqual(1.0, xx[0], TOL);
            Assert.AreEqual(0.0, xx[1], TOL);
        }

        [Test]
        public void FindMinimum_HS05_ReturnsValidMinimum()
        {
            var xx = new[] { 0.0, 0.0 };
            Bobyqa.FindMinimum(
                (n, x) => Math.Sin(x[0] + x[1]) + Math.Pow(x[0] - x[1], 2.0) - 1.5 * x[0] + 2.5 * x[1] + 1, 2, xx,
                new[] { -1.5, -3.0 }, new[] { 4.0, 3.0 });
            Assert.AreEqual(0.5 - Math.PI / 3.0, xx[0], TOL);
            Assert.AreEqual(-0.5 - Math.PI / 3.0, xx[1], TOL);
        }

        [TestCase(5, 16)]
        [TestCase(5, 21)]
        [TestCase(10, 26)]
        [TestCase(10, 41)]
        public void FindMinimum_BobyqaTestCase_ReturnStatusNormal(int m, int npt)
        {
            //     Test problem for BOBYQA, the objective function being the sum of
            //     the reciprocals of all pairwise distances between the points P_I,
            //     I=1,2,...,M in two dimensions, where M=N/2 and where the components
            //     of P_I are X(2*I-1) and X(2*I). Thus each vector X of N variables
            //     defines the M points P_I. The initial X gives equally spaced points
            //     on a circle. Four different choices of the pairs (N,NPT) are tried,
            //     namely (10,16), (10,21), (20,26) and (20,41). Convergence to a local
            //     minimum that is not global occurs in both the N=10 cases. The details
            //     of the results are highly sensitive to computer rounding errors. The
            //     choice IPRINT=2 provides the current X and optimal F so far whenever
            //     RHO is reduced. The bound constraints of the problem require every
            //     component of X to be in the interval [-1,1].

            var n = 2 * m;
            var x = new double[n];
            var xl = new double[n];
            var xu = new double[n];

            const double bdl = -1.0;
            const double bdu = 1.0;
            for (var i = 0; i < n; ++i)
            {
                xl[i] = bdl;
                xu[i] = bdu;
            }

            Console.WriteLine("{0}2D output with M ={1,4},  N ={2,4}  and  NPT ={3,4}", Environment.NewLine, m, n, npt);

            for (var j = 1; j <= m; ++j)
            {
                var temp = 2.0 * Math.PI * j / m;
                x[2 * j - 2] = Math.Cos(temp);
                x[2 * j - 1] = Math.Sin(temp);
            }

            const int iprint = 2;
            const int maxfun = 500000;
            const double rhobeg = 1.0E-1;
            const double rhoend = 1.0E-6;

            const BobyqaExitStatus expected = BobyqaExitStatus.Normal;
            var actual = Bobyqa.FindMinimum(BobyqaTestCalfun, n, x, xl, xu, npt, rhobeg, rhoend, iprint, maxfun);
            Assert.AreEqual(expected, actual);
        }

        public double BobyqaTestCalfun(int n, double[] x)
        {
            var f = 0.0;
            for (var i = 4; i <= n; i += 2)
            {
                for (var j = 2; j <= i - 2; j += 2)
                {
                    var temp = Math.Max(Math.Pow(x[i - 2] - x[j - 2], 2.0) + Math.Pow(x[i - 1] - x[j - 1], 2.0), 1.0e-6);
                    f += 1.0 / Math.Sqrt(temp);
                }
            }
            return f;
        }

        [TestCase(13, 78)]
        public void FindMinimum_ConstrainedRosenWithAdditionalInterpolationPoints_ReturnsValidMinimum(int n, int maxAdditionalPoints)
        {
            var xl = Enumerable.Repeat(-1.0, n).ToArray();
            var xu = Enumerable.Repeat(2.0, n).ToArray();
            var expected = Enumerable.Repeat(1.0, n).ToArray();

            for (var num = 1; num <= maxAdditionalPoints; ++num)
            {
                Console.WriteLine("\nNumber of additional points = {0}", num);
                var npt = 2 * n + 1 + num;
                var x = Enumerable.Repeat(0.1, n).ToArray();
                Bobyqa.FindMinimum(Rosen, n, x, xl, xu, npt, -1, 1.0e-8, 1, 3000, Console.Out);
                CollectionAssert.AreEqual(expected, x, new DoubleComparer(1.0e-6));
            }
        }
        
        [Test]
        public void FindMinimum_LogOutput_OutputNonEmpty()
        {
            const int n = 9;
            var x = Enumerable.Repeat(0.1, n).ToArray();
            using (var logger = new StringWriter())
            {
                Bobyqa.FindMinimum(Rosen, n, x, null, null, -1, 1.0, 1.0e-8, 1, 2000, logger);
                Assert.Greater(logger.ToString().Length, 0);
            }
        }

        [Test]
        public void FindMinimum_LogOutputToConsole_OutputNonEmpty()
        {
            const int n = 9;
            var x = Enumerable.Repeat(0.1, n).ToArray();
            Bobyqa.FindMinimum(Rosen, n, x, null, null, -1, 1.0, 1.0e-8, 1, 2000, Console.Out);
        }

        public double Rosen(int n, double[] x)
        {
            var f = 0.0;
            for (var i = 0; i < n - 1; ++i)
                f += 1e2 * (x[i] * x[i] - x[i + 1]) * (x[i] * x[i] - x[i + 1]) + (x[i] - 1.0) * (x[i] - 1.0);
            return f;
        }

        #endregion
    }

    internal class DoubleComparer : IComparer
    {
        private readonly double _tol;

        internal DoubleComparer(double tol)
        {
            _tol = tol;
        }

        internal int Compare(double x, double y)
        {
            return Math.Abs(x - y) <= _tol ? 0 : x.CompareTo(y);
        }

        int IComparer.Compare(object x, object y)
        {
            return Compare((double)x, (double)y);
        }
    }
}
