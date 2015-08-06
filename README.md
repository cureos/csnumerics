<img src="NuGet/csnumerics.png" alt="CSNumerics logo" height="108" />

# CSNumerics

[![Join the chat at https://gitter.im/cureos/csnumerics](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/cureos/csnumerics?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

<b>Portable numerical algorithms in C#</b>

Copyright (c) 2012-1014 Anders Gustafsson, Cureos AB. Licensed under the General Public License version 3.<br />

[![NuGet](https://img.shields.io/nuget/v/csnumerics.svg)](https://www.nuget.org/packages/csnumerics/)
[![NuGet](https://img.shields.io/nuget/dt/csnumerics.svg)](https://www.nuget.org/packages/csnumerics/)
[![Build status](https://ci.appveyor.com/api/projects/status/2msksh7auurc4hu4?svg=true)](https://ci.appveyor.com/project/anders9ustafsson/csnumerics)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/cureos/csnumerics?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

## Usage

*CSNumerics* is a Portable Class Library (PCL) of various numerical algorithms written in C#. 

It currently consists of C# implementations of Michael J.D. Powell's optimization algorithms:

* [BOBYQA](https://github.com/cureos/csnumerics/wiki/BOBYQA) for minimizing a nonlinear objective function subject to variable bounds.
* [LINCOA](https://github.com/cureos/csnumerics/wiki/LINCOA) for minimizing a nonlinear objective function subject to linear constraints.
* [COBYLA](https://github.com/cureos/csnumerics/wiki/COBYLA) for minimizing a nonlinear objective function subject to nonlinear constraints.

The *CSNumerics* assembly is immediately applicable with the following targets:

* .NET Framework 4 and later
* Windows Store apps (Windows 8 and later)
* Windows Phone Silverlight 8 and 8.1
* Windows Phone 8.1
* Silverlight 5
* Xamarin.iOS
* Xamarin.Android

Use *NuGet* to include [CSNumerics](https://www.nuget.org/packages/csnumerics) in your application.


## Notes on commercial use

The *CSNumerics* class library is published under the General Public License, version 3.
For those interested in using *CSNumerics* without having to adhere to GPL, please contact the copyright holder at

licenses@cureos.com

for commercial licensing alternatives.
