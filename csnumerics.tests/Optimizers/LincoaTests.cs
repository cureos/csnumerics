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
 */

using System;
using NUnit.Framework;

namespace Cureos.Numerics.Optimizers
{
    [TestFixture]
    public class LincoaTests
    {
        #region FIELDS

        private double _fmax;

        #endregion

        #region Unit tests

        [TestCase(15, 0.01)]
        [TestCase(30, 0.001)]
        [TestCase(45, 0.001)]
        public void FindMinimum_PtsinTet_YieldsExpectedValue(int npt, double tol)
        {
            //     Set some constants.
            const int n = 12;

            //     Set the data points.
            const int np = 50;
            var sumx = 0.0;
            var sumy = 0.0;
            var sumz = 0.0;

            var xp = new double[50];
            var yp = new double[50];
            var zp = new double[50];
            for (var j = 0; j < np; ++j)
            {
                var theta = j * Math.PI / (np - 1.0);
                xp[j] = Math.Cos(theta) * Math.Cos(2.0 * theta);
                sumx += xp[j];
                yp[j] = Math.Sin(theta) * Math.Cos(2.0 * theta);
                sumy += yp[j];
                zp[j] = Math.Sin(2.0 * theta);
                sumz += zp[j];
            }
            sumx /= np;
            sumy /= np;
            sumz /= np;
            for (var j = 0; j < np; ++j)
            {
                xp[j] -= sumx;
                yp[j] -= sumy;
                zp[j] -= sumz;
            }

            //     Set the linear constraints.
            const int m = 4 * np;
            var a = new double[m, n];
            var b = new double[m];
            for (var k = 0; k < m; ++k)
            {
                b[k] = 1.0;
                for (var i = 0; i < n; ++i)
                    a[k, i] = 0.0;
            }
            for (var j = 0; j < np; ++j)
            {
                for (var i = 0; i < 4; ++i)
                {
                    var k = 4 * j + i;
                    var iw = 3 * i;
                    a[k, iw] = xp[j];
                    a[k, iw + 1] = yp[j];
                    a[k, iw + 2] = zp[j];
                }
            }

            //     Set the initial vector of variables. 
            var xs = 0.0;
            var ys = 0.0;
            var zs = 0.0;
            var ss = 0.0;
            for (var j = 0; j < np; ++j)
            {
                xs = Math.Min(xs, xp[j]);
                ys = Math.Min(ys, yp[j]);
                zs = Math.Min(zs, zp[j]);
                ss = Math.Max(ss, xp[j] + yp[j] + zp[j]);
            }

            var x = new double[12];
            x[0] = 1.0 / xs;
            x[4] = 1.0 / ys;
            x[8] = 1.0 / zs;
            x[9] = 1.0 / ss;
            x[10] = 1.0 / ss;
            x[11] = 1.0 / ss;

            this._fmax = Math.Pow(ss - xs - ys - zs, 3.0) / 6.0;

            //     Call of LINCOA, which provides the printing given at the end of this note.
            const double rhobeg = 1.0;
            const double rhoend = 1.0E-6;
            const int iprint = 1;
            const int maxfun = 10000;

            var lincoa = new Lincoa(this.PtsinTet, a, b)
                             {
                                 InterpolationConditions = npt,
                                 TrustRegionRadiusStart = rhobeg,
                                 TrustRegionRadiusEnd = rhoend,
                                 MaximumFunctionCalls = maxfun,
                                 PrintLevel = iprint,
                                 Logger = Console.Out
                             };
            var result = lincoa.FindMinimum(x);

            const double expected = 2.761;
            var actual = result.F;
            Assert.AreEqual(expected, actual, tol);
        }

        public double PtsinTet(int n, double[] x, bool isXFeasible)
        {
            var v12 = x[0] * x[4] - x[3] * x[1];
            var v13 = x[0] * x[7] - x[6] * x[1];
            var v14 = x[0] * x[10] - x[9] * x[1];
            var v23 = x[3] * x[7] - x[6] * x[4];
            var v24 = x[3] * x[10] - x[9] * x[4];
            var v34 = x[6] * x[10] - x[9] * x[7];
            var del1 = v23 * x[11] - v24 * x[8] + v34 * x[5];
            if (del1 <= 0.0) return this._fmax;
            var del2 = -v34 * x[2] - v13 * x[11] + v14 * x[8];
            if (del2 <= 0.0) return this._fmax;
            var del3 = -v14 * x[5] + v24 * x[2] + v12 * x[11];
            if (del3 <= 0.0) return this._fmax;
            var del4 = -v12 * x[8] + v13 * x[5] - v23 * x[2];
            if (del4 <= 0.0) return this._fmax;
            var temp = Math.Pow(del1 + del2 + del3 + del4, 3.0) / (del1 * del2 * del3 * del4);
            return Math.Min(temp / 6.0, this._fmax);
        }

        #endregion
        
        [TestCase(2, -100, -300)]
        public void FindMinimum_EasyFun_YieldsExpectedValue(int dim, double minVal, double expect)
        {
            var constraintMatrix = new double[dim,dim];
            var constraintVals = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                constraintMatrix[i, i] = -1;
                constraintVals[i] = -minVal;
            }

            var lincoa = new Lincoa(EasyFun, constraintMatrix, constraintVals)
            {
                TrustRegionRadiusStart = 1,
                TrustRegionRadiusEnd = 0.01,
                MaximumFunctionCalls = 10000,
                PrintLevel = 3,
                Logger = Console.Out
            };
            var result = lincoa.FindMinimum(new double[dim]);

            var actual = result.F;
            Assert.AreEqual(expect, actual, 0.1);

        }

        public double EasyFun(int n, double[] x, bool isXFeasible)
        {
            double total = 0;
            for (var i = 0; i < x.Length; i++)
            {
                total += (i + 1)*x[i];
            }
            return total;
        }
    }
}
