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

namespace Cureos.Numerics.Optimizers
{
    /// <summary>
    /// Interface for derivative-free optimizers following the basic principles of M.J.D. Powell's trust region optimizers.
    /// </summary>
    internal interface ITrustRegionOptimizer : IOptimizer
    {
        /// <summary>
        /// Gets or sets the final value of the trust region radius.
        /// </summary>
        double TrustRegionRadiusStart { get; set; }

        /// <summary>
        /// Gets or sets the start value of the trust region radius.
        /// </summary>
        double TrustRegionRadiusEnd { get; set; }
    }
}