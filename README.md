<img src="csnumerics.png" alt="CSNumerics logo" height="108" />

# CSNumerics

<b>Portable numerical algorithms in C#</b>

Copyright (c) 2012-2022 Anders Gustafsson, Cureos AB. Licensed under the GNU Lesser General Public License version 3.<br />

[![NuGet](https://img.shields.io/nuget/v/csnumerics.svg)](https://www.nuget.org/packages/csnumerics/)
[![NuGet](https://img.shields.io/nuget/dt/csnumerics.svg)](https://www.nuget.org/packages/csnumerics/)
[![Build status](https://ci.appveyor.com/api/projects/status/2msksh7auurc4hu4?svg=true)](https://ci.appveyor.com/project/anders9ustafsson/csnumerics)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/cureos/csnumerics?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

## Usage

*CSNumerics* is a .NET Standard 1.0 and 2.0 class library of various numerical algorithms written in C#. 

It currently consists of C# implementations of Michael J.D. Powell's optimization algorithms:

* [BOBYQA](https://github.com/cureos/csnumerics/wiki/BOBYQA) for minimizing a nonlinear objective function subject to variable bounds.
* [LINCOA](https://github.com/cureos/csnumerics/wiki/LINCOA) for minimizing a nonlinear objective function subject to linear constraints.
* [COBYLA](https://github.com/cureos/csnumerics/wiki/COBYLA) for minimizing a nonlinear objective function subject to nonlinear constraints.

Use *NuGet* to include [CSNumerics](https://www.nuget.org/packages/csnumerics) in your application.
