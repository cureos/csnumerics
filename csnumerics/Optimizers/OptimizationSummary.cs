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
 */

namespace Cureos.Numerics.Optimizers
{
    /// <summary>
    /// Summary of output from optimization.
    /// </summary>
    public class OptimizationSummary
    {
        #region CONSTRUCTORS

        /// <summary>
        /// Initializes an instance of the optimization summary.
        /// </summary>
        /// <param name="status">Status of the completed optimization.</param>
        /// <param name="nf">Number of function evaluations.</param>
        /// <param name="x">Optimized variable array.</param>
        /// <param name="f">Optimal value of the objective function.</param>
        /// <param name="g">If defined, values of constraint functions at optimum.</param>
        internal OptimizationSummary(OptimizationStatus status, int nf, double[] x, double f, double[] g = null)
        {
            Status = status;
            Evals = nf;
            X = x;
            F = f;
            G = g;
        }

        #endregion

        #region PROPERTIES

        /// <summary>
        /// Gets the status of the completed optimization.
        /// </summary>
        public OptimizationStatus Status { get; private set; }

        /// <summary>
        /// Gets the number of function evaluations.
        /// </summary>
        public int Evals { get; private set; }

        /// <summary>
        /// Gets the optimized variable array.
        /// </summary>
        public double[] X { get; private set; }

        /// <summary>
        /// Gets the optimal value of the objective function.
        /// </summary>
        public double F { get; private set; }

        /// <summary>
        /// Gets the values of the constraint functions at optimum, if defined.
        /// </summary>
        public double[] G { get; private set; }

        #endregion
    }
}