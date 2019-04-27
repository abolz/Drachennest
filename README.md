[![Build Status](https://travis-ci.org/abolz/Grisu.svg?branch=master)](https://travis-ci.org/abolz/Grisu)
[![Build status](https://ci.appveyor.com/api/projects/status/9c62o8sm73tg2k3x?svg=true)](https://ci.appveyor.com/project/abolz/grisu)
[![codecov](https://codecov.io/gh/abolz/Grisu/branch/master/graph/badge.svg)](https://codecov.io/gh/abolz/Grisu)

Converts binary floating-point to decimal floating-point numbers.

---

Contains an implementation of the [Grisu2](https://github.com/abolz/Grisu/blob/master/src/grisu2.h)
and [Grisu3](https://github.com/abolz/Grisu/blob/master/src/grisu3.h) algorithms as described in

* Loitsch, [_Printing Floating-Point Numbers Quickly and Accurately with Integers_](https://dl.acm.org/citation.cfm?id=1806623),

The Grisu3 implementation uses the Dragon4 algorithm as a fallback.

* Steele, White, [_How to Print FloatingPoint Numbers Accurately_](https://dl.acm.org/citation.cfm?id=93559),
* Burger, Dybvig, [_Printing Floating-Point Numbers Quickly and Accurately_](https://dl.acm.org/citation.cfm?id=231397),

Contains an implementation of the [Ryu](https://github.com/abolz/Grisu/blob/master/src/ryu.h)
algorithm as described in

* Adams, [_Ryu: fast float-to-string conversion_](https://dl.acm.org/citation.cfm?id=3192369),

---

Grisu3 and Ryu are optimal, i.e. the output string
1. rounds back to the input number when read in,
2. is as short as possible,
3. is as close to the input number as possible.

Grisu3 and Ryu (currently) assume that the input rounding algorithm uses
_round-to-nearest-even_ to break ties. Grisu2 only is optimal for ~99% of all
floating point numbers, though it guarantees the first property for all of its
inputs, _regardless of how the input rounding mode breaks ties_.
