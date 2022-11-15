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

using System.IO;

namespace Cureos.Numerics.Optimizers
{
    /// <summary>
    /// General interface for optimizers in the CS Numerics class library.
    /// </summary>
    interface IOptimizer
    {
        #region PROPERTIES

        /// <summary>
        /// Gets or sets the number of maximum function calls.
        /// </summary>
        int MaximumFunctionCalls { get; set; }

        /// <summary>
        /// Gets or sets the print level to the logger.
        /// </summary>
        int PrintLevel { get; set; }

        /// <summary>
        /// Gets or sets the logger to which the optimizer log information should be sent.
        /// </summary>
        TextWriter Logger { get; set; }

        #endregion

        #region METHODS

        /// <summary>
        /// Find a local minimum of provided objective function satisfying the provided linear constraints.
        /// </summary>
        /// <param name="x0">Initial variable array.</param>
        /// <returns>Summary of the optimization result.</returns>
        OptimizationSummary FindMinimum(double[] x0);

        #endregion
    }
}