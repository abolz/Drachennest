[![Build Status](https://travis-ci.org/abolz/Drachennest.svg?branch=master)](https://travis-ci.org/abolz/Drachennest)
[![Build status](https://ci.appveyor.com/api/projects/status/py96h02xct0ycdqs?svg=true)](https://ci.appveyor.com/project/abolz/drachennest)
[![codecov](https://codecov.io/gh/abolz/Drachennest/branch/master/graph/badge.svg)](https://codecov.io/gh/abolz/Drachennest)

Converting binary floating-point to decimal floating-point numbers.

---

Grisu / Dragon
--------------------------------------------------------------------------------

Contains an implementation of the [Grisu2](https://github.com/abolz/Drachennest/blob/master/src/grisu2.h)
and [Grisu3](https://github.com/abolz/Drachennest/blob/master/src/grisu3.h) algorithms as described in

* Loitsch, [_Printing Floating-Point Numbers Quickly and Accurately with Integers_](https://dl.acm.org/citation.cfm?id=1806623),

The Grisu3 implementation uses the [Dragon4](https://github.com/abolz/Drachennest/blob/master/src/dragon4.h)
algorithm as a fallback.

* Steele, White, [_How to Print FloatingPoint Numbers Accurately_](https://dl.acm.org/citation.cfm?id=93559),
* Burger, Dybvig, [_Printing Floating-Point Numbers Quickly and Accurately_](https://dl.acm.org/citation.cfm?id=231397),

Ryu
--------------------------------------------------------------------------------

Contains an implementation of the [Ryu](https://github.com/abolz/Drachennest/blob/master/src/ryu_64.cc)
algorithm as described in

* Adams, [_Ryu: fast float-to-string conversion_](https://dl.acm.org/citation.cfm?id=3192369),

The implemenation also contains a (fast!) `strtod` implementation, which can be
used to convert decimal numbers with at most 17 significant decimal digits back
into binary floating-point numbers. (Note that none of the algorithms here will
ever produce more than 17 significant digits.)

Schubfach
--------------------------------------------------------------------------------

Contains an implementation of the [Schubfach](https://github.com/abolz/Drachennest/blob/master/src/schubfach_64.cc)
algorithm as described in

* Giulietti, [The Schubfach way to render doubles](https://drive.google.com/open?id=1luHhyQF9zKlM8yJ1nebU0OgVYhfC6CBN)

The name of this algorithm "deliberately departs from a long lineage of fabulous drakes".

Dragonbox
--------------------------------------------------------------------------------

Contains a slightly modified version the reference implementation of
Junekey Jeon's [Dragonbox](https://github.com/jk-jeon/dragonbox) algorithm.

---

Grisu3, Ryu, Schubfach, and Dragonbox are optimal, i.e. the output string
1. rounds back to the input number when read in,
2. is as short as possible,
3. is as close to the input number as possible.

These algorithms (currently) assume that the input rounding algorithm uses
_round-to-nearest-even_ to break ties. Grisu2 only is optimal for ~99% of all
floating point numbers, though it guarantees the first property for all of its
inputs, _regardless of how the input rounding mode breaks ties_.

---

Benchmarks
--------------------------------------------------------------------------------

> Benchmarks were run on an Intel Core i7-9750H, using Visual Studio 2019 16.7.7, Clang 10.0, 64-bit.

> Timings are in ns.

---

For this benchmark uniformly distributed random `double`s in the
range `[1,2]` have been generated. These numbers were then rounded to `N`
significant digits and converted to decimal using the given algorithm.

![BenchDigits](/resources/bench_digits.png)

---

Uniformly distributed random numbers in the range `[10^i, 10^(i+1)]` for
`i=-12,...,12`.

![BenchUniform](/resources/bench_uniform.png)

---

Uniformly distributed random numbers in the range `[0, 10^10]`. Each benchmark
is run 10 times (using different numbers each run).

![BenchUniformE10](/resources/bench_uniform_e10.png)

---

Random bit patterns. Each benchmark is run 10 times (using different numbers
each run).

![BenchRandom](/resources/bench_random_bits.png)
