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
    /// Optimization status codes.
    /// </summary>
    public enum OptimizationStatus
    {
// ReSharper disable InconsistentNaming

        /*
         * Common
         */

        /// <summary>
        /// Optimization successfully completed.
        /// </summary>
        Normal,
        N_TooSmall,
        NPT_OutOfRange,
        MAXFUN_NotLargerThan_NPT,

        /// <summary>
        /// Maximum number of iterations (function/constraints evaluations) reached during optimization.
        /// </summary>
        MAXFUN_Reached,

        /// <summary>
        /// Size of rounding error is becoming damaging, terminating prematurely.
        /// </summary>
        X_RoundingErrorsPreventUpdate,

        /*
         * LINCOA specific
         */
        ConstraintGradientIsZero,
        UpdatingFormulaDenominatorZero,

        /*
         * BOBYQA specific
         */
        VariableBoundsArrayTooShort,
        InvalidBoundsSpecification,
        BoundsRangeTooSmall,
        DenominatorCancellation,
        TrustRegionStepReductionFailure

// ReSharper enable InconsistentNaming
    }
}