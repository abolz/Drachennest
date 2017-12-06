// Copyright 2017 Alexander Bolz
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <cassert>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

#if 1 && defined(_MSC_VER)
#define GRISU2_INLINE __forceinline
#else
#define GRISU2_INLINE inline
#endif

// If rounding is enabled, Grisu2 produces the optimal output for ~99.8% of all
// input numbers.
// This step is not required for the correctness of the algorithm.
#define GRISU2_ROUND 1

//
// Implements the Grisu2 algorithm for binary to decimal floating-point
// conversion.
//
// This implementation is a slightly modified version of the reference
// implementation by Florian Loitsch which can be obtained from
// http://florian.loitsch.com/publications (bench.tar.gz)
//
// References:
//
// [1]  Loitsch, "Printing Floating-Point Numbers Quickly and Accurately with Integers",
//      Proceedings of the ACM SIGPLAN 2010 Conference on Programming Language Design and Implementation, PLDI 2010
// [2]  Burger, Dybvig, "Printing Floating-Point Numbers Quickly and Accurately",
//      Proceedings of the ACM SIGPLAN 1996 Conference on Programming Language Design and Implementation, PLDI 1996
//
namespace fast_dtoa {

template <typename Float>
struct IEEEFloat
{
    static constexpr bool const is_single
        = std::numeric_limits<Float>::digits == 24 && std::numeric_limits<Float>::max_exponent == 128;
    static constexpr bool const is_double
        = std::numeric_limits<Float>::digits == 53 && std::numeric_limits<Float>::max_exponent == 1024;

    static_assert(std::numeric_limits<Float>::is_iec559 && (is_single || is_double),
        "IEEE-754 single or double precision implementation required");

    using Uint = typename std::conditional<is_single, uint32_t, uint64_t>::type;

    static_assert(sizeof(Float) == sizeof(Uint), "Size mismatch");

    static constexpr int  const kPrecision       = is_single ? 24 : 53; // = p (includes the hidden bit!)
    static constexpr int  const kExponentBias    = is_single ? 0x7F : 0x3FF;
    static constexpr Uint const kHiddenBit       = Uint{1} << (kPrecision - 1);
    static constexpr Uint const kSignMask        = Uint{1} << (is_single ? 31 : 63);
    static constexpr Uint const kExponentMask    = Uint{is_single ? 0xFF : 0x7FF} << (kPrecision - 1);
    static constexpr Uint const kSignificandMask = kHiddenBit - 1;

    union
    {
        Float value;
        Uint bits;
    };

    explicit IEEEFloat(Float value_) : value(value_) {}
    explicit IEEEFloat(Uint bits_) : bits(bits_) {}

    Uint ExponentBits() const {
        return (bits & kExponentMask) >> (kPrecision - 1);
    }

    Uint SignificandBits() const {
        return (bits & kSignificandMask);
    }

    // Returns true if the sign-bit is set.
    bool IsNegative() const {
        return (bits & kSignMask) != 0;
    }

    // Returns true if this value is -0 or +0.
    bool IsZero() const {
        return (bits & ~kSignMask) == 0;
    }

    // Returns true if this value is denormal or 0.
    bool IsDenormal() const {
        return (bits & kExponentMask) == 0;
    }

    // Returns true if this value is NaN
    bool IsNaN() const {
        return (bits & kExponentMask) == kExponentMask && (bits & kSignificandMask) != 0;
    }

    // Returns true if this value is -Inf or +Inf.
    bool IsInf() const {
        return (bits & kExponentMask) == kExponentMask && (bits & kSignificandMask) == 0;
    }

    // Returns this value with the sign-bit cleared.
    Float Abs() const {
        return IEEEFloat(bits & ~kSignMask).value;
    }
};

struct Fp // f * 2^e
{
    static constexpr int const kPrecision = 64; // = q

    uint64_t f;
    int e;

    constexpr Fp() : f(0), e(0) {}
    constexpr Fp(uint64_t f_, int e_) : f(f_), e(e_) {}

    // Returns x - y.
    // Requires: x.e == y.e and x.f >= y.f
    static Fp Sub(Fp x, Fp y);

    // Returns x * y.
    // The result is rounded. (Only the upper q bits are returned.)
    static Fp Mul(Fp x, Fp y);

    // Normalize x such that the significand is >= 2^(q-1).
    // Requires: x.f != 0
    static Fp Normalize(Fp x);

    // Normalize x such that the result has the exponent E.
    // Requires: e >= x.e and the upper e - x.e bits of x.f must be zero.
    static Fp NormalizeTo(Fp x, int e);
};

inline Fp Fp::Sub(Fp x, Fp y)
{
    assert(x.e == y.e);
    assert(x.f >= y.f);

    return Fp(x.f - y.f, x.e);
}

inline Fp Fp::Mul(Fp x, Fp y)
{
    // Computes:
    //  f = round((x.f * y.f) / 2^q)
    //  e = x.e + y.e + q

#if defined(_MSC_VER) && defined(_M_X64)

    uint64_t h = 0;
    uint64_t l = _umul128(x.f, y.f, &h);
    h += l >> 63; // round, ties up: [h, l] += 2^q / 2

    return Fp(h, x.e + y.e + 64);

#elif defined(__GNUC__) && defined(__SIZEOF_INT128__)

    __extension__ using Uint128 = unsigned __int128;

    Uint128 const p = Uint128{x.f} * Uint128{y.f};

    uint64_t h = static_cast<uint64_t>(p >> 64);
    uint64_t l = static_cast<uint64_t>(p);
    h += l >> 63; // round, ties up: [h, l] += 2^q / 2

    return Fp(h, x.e + y.e + 64);

#else

    // Emulate the 64-bit * 64-bit multiplication:
    //
    // p = u * v
    //   = (u_lo + 2^32 u_hi) (v_lo + 2^32 v_hi)
    //   = (u_lo v_lo         ) + 2^32 ((u_lo v_hi         ) + (u_hi v_lo         )) + 2^64 (u_hi v_hi         )
    //   = (p0                ) + 2^32 ((p1                ) + (p2                )) + 2^64 (p3                )
    //   = (p0_lo + 2^32 p0_hi) + 2^32 ((p1_lo + 2^32 p1_hi) + (p2_lo + 2^32 p2_hi)) + 2^64 (p3                )
    //   = (p0_lo             ) + 2^32 (p0_hi + p1_lo + p2_lo                      ) + 2^64 (p1_hi + p2_hi + p3)
    //   = (p0_lo             ) + 2^32 (Q                                          ) + 2^64 (H                 )
    //   = (p0_lo             ) + 2^32 (Q_lo + 2^32 Q_hi                           ) + 2^64 (H                 )
    //
    // (Since Q might be larger than 2^32 - 1)
    //
    //   = (p0_lo + 2^32 Q_lo) + 2^64 (Q_hi + H)
    //
    // (Q_hi + H does not overflow a 64-bit int)
    //
    //   = p_lo + 2^64 p_hi

    uint64_t const u_lo = x.f & 0xFFFFFFFF;
    uint64_t const u_hi = x.f >> 32;
    uint64_t const v_lo = y.f & 0xFFFFFFFF;
    uint64_t const v_hi = y.f >> 32;

    uint64_t const p0 = u_lo * v_lo;
    uint64_t const p1 = u_lo * v_hi;
    uint64_t const p2 = u_hi * v_lo;
    uint64_t const p3 = u_hi * v_hi;

    uint64_t const p0_hi = p0 >> 32;
    uint64_t const p1_lo = p1 & 0xFFFFFFFF;
    uint64_t const p1_hi = p1 >> 32;
    uint64_t const p2_lo = p2 & 0xFFFFFFFF;
    uint64_t const p2_hi = p2 >> 32;

    uint64_t Q = p0_hi + p1_lo + p2_lo;

    // The full product might now be computed as
    //
    // p_hi = p3 + p2_hi + p1_hi + (Q >> 32)
    // p_lo = p0_lo + (Q << 32)
    //
    // But in this particular case here, the full p_lo is not required.
    // Effectively we only need to add the highest bit in p_lo to p_hi (and
    // Q_hi + 1 does not overflow).

    Q += uint64_t{1} << (63 - 32); // round, ties up

    uint64_t const h = p3 + p2_hi + p1_hi + (Q >> 32);

    return Fp(h, x.e + y.e + 64);

#endif
}

// Returns the number of leading zeros of the 64-bit integer n.
// The result is undefined for n = 0.
inline int CountLeadingZeros_64(uint64_t n)
{
    assert(n != 0);

#if defined(_MSC_VER) && defined(_M_X64)

    unsigned long high = 0; // index of most significant 1-bit
    _BitScanReverse64(&high, n);
    return 63 - static_cast<int>(high);

#elif defined(__GNUC__)

    return __builtin_clzll(static_cast<unsigned long long>(n));

#else

    int z = 0;
    while ((n >> 63) == 0)
    {
        z++;
        n <<= 1;
    }

    return z;

#endif
}

inline Fp Fp::Normalize(Fp x)
{
    int const leading_zeros = CountLeadingZeros_64(x.f);
    return Fp(x.f << leading_zeros, x.e - leading_zeros);
}

inline Fp Fp::NormalizeTo(Fp x, int e)
{
    int const delta = x.e - e;

    assert(delta >= 0);
    assert(((x.f << delta) >> delta) == x.f);

    return Fp(x.f << delta, e);
}

#if GRISU2_ROUND
struct FpBoundaries {
    Fp w;
    Fp minus;
    Fp plus;
};
#else
struct FpBoundaries {
    Fp minus;
    Fp plus;
};
#endif

//
// Computes the boundaries m- and m+ of the floating-point value v.
//
// Determine v- and v+, the floating-point predecessor and successor if v,
// respectively.
//
//      v- = v - 2^e        if f != 2^(p-1) or e != e_min                    (A)
//         = v - 2^(e-1)    if f == 2^(p-1) and e > e_min                    (B)
//
//      v+ = v + 2^e
//
// Let m- = (v- + v) / 2 and m+ = (v + v+) / 2. All real numbers _strictly_
// between m- and m+ round to v, regardless of how the input rounding algorithm
// breaks ties.
//
//      ---+-------------+-------------+-------------+-------------+---      (A)
//         v-            m-            v             m+            v+
//
//      -----------------+------+------+-------------+-------------+---      (B)
//                       v-     m-     v             m+            v+
//
// Note that m- and m+ are (by definition) not representable with precision p
// and we therefore need some extra bits of precision.
//
template <typename Float>
GRISU2_INLINE FpBoundaries ComputeBoundaries(Float v_ieee)
{
    using IEEEType = IEEEFloat<Float>;

    //
    // Convert the IEEE representation into a DiyFp.
    //
    // If v is denormal:
    //      value = 0.F * 2^(1 - E_bias) = (          F) * 2^(1 - E_bias - (p-1))
    // If v is normalized:
    //      value = 1.F * 2^(E - E_bias) = (2^(p-1) + F) * 2^(E - E_bias - (p-1))
    //

    IEEEType const v_ieee_bits(v_ieee);

    uint64_t const E = v_ieee_bits.ExponentBits(); // biased exponent
    uint64_t const F = v_ieee_bits.SignificandBits();

    constexpr int const kBias = IEEEType::kExponentBias + (IEEEType::kPrecision - 1);

    Fp const v = (E == 0) // denormal?
        ? Fp(F, 1 - kBias)
        : Fp(IEEEType::kHiddenBit + F, static_cast<int>(E) - kBias);

    //
    // v+ = v + 2^e = (f + 1) * 2^e and therefore
    //
    //      m+ = (v + v+) / 2
    //         = (2*f + 1) * 2^(e-1)
    //
    Fp const m_plus = Fp(2*v.f + 1, v.e - 1);

    //
    // If f != 2^(p-1), then v- = v - 2^e = (f - 1) * 2^e and
    //
    //      m- = (v- + v) / 2
    //         = (2*f - 1) * 2^(e-1)
    //
    // If f = 2^(p-1), then the next smaller _normalized_ floating-point number
    // is actually v- = v - 2^(e-1) = (2^p - 1) * 2^(e-1) and therefore
    //
    //      m- = (4*f - 1) * 2^(e-2)
    //
    // The exception is the smallest normalized floating-point number
    // v = 2^(p-1) * 2^e_min. In this case the predecessor is the largest
    // denormalized floating-point number: v- = (2^(p-1) - 1) * 2^e_min and then
    //
    //      m- = (2*f - 1) * 2^(e-1)
    //
    // If v is denormal, v = f * 2^e_min and v- = v - 2^e = (f - 1) * 2^e and
    // again
    //
    //      m- = (2*f - 1) * 2^(e-1)
    //
    // Note: 0 is not a valid input for Grisu and in case v is denormal:
    // f != 2^(p-1).
    //
    // For IEEE floating-point numbers not equal to 0, the condition f = 2^(p-1)
    // is equivalent to F = 0, and for the smallest normalized number E = 1.
    // For denormals E = 0 (and F != 0).
    //
    Fp const m_minus = (F == 0 && E > 1)
        ? Fp(4*v.f - 1, v.e - 2)
        : Fp(2*v.f - 1, v.e - 1);

    //
    // Determine the normalized w+ = m+.
    //
    Fp const plus = Fp::Normalize(m_plus);

    //
    // Determine w- = m- such that e_(w-) = e_(w+).
    //
    Fp const minus = Fp::NormalizeTo(m_minus, plus.e);

#if GRISU2_ROUND
    return {Fp::Normalize(v), minus, plus};
#else
    return {minus, plus};
#endif
}

//
// Given a (normalized) floating-point number v and its neighbors m- and m+
//
//      ---+---------------------------+---------------------------+---
//         m-                          v                           m+
//
// Grisu first scales the input number w, and its boundaries w- and w+, by an
// approximate power-of-ten c ~= 10^-k (which needs to be precomputed using
// high-precision arithmetic and stored in a table) such that the exponent of
// the products lies within a certain range [alpha, gamma]. It then remains to
// produce the decimal digits of the number M = f * 2^e, where alpha <= e <= gamma.
//
// The choice of alpha and gamma determines the digit generation procedure and
// the size of the look-up table (and/or vice versa...) and depends on the
// extended precision q of the DiyFp's.
//
// In other words, given normalized w, Grisu needs to find a (normalized) cached
// power-of-ten c, such that the exponent of the product c * w = f * 2^e
// satisfies (Definition 3.2 from [1])
//
//      alpha <= e = e_c + e_w + q <= gamma
//
// or
//
//      f_c * f_w * 2^alpha <= f_c 2^(e_c) * f_w 2^(e_w) * 2^q
//                          <= f_c * f_w * 2^gamma
//
// Since c and w are normalized, i.e. 2^(q-1) <= f < 2^q, this implies
//
//      2^(q-1) * 2^(q-1) * 2^alpha <= c * w * 2^q < 2^q * 2^q * 2^gamma
//
// or
//
//      2^(q - 2 + alpha) <= c * w < 2^(q + gamma)
//
// The distance (gamma - alpha) should be as large as possible in order to make
// the table as small as possible, but the digit generation procedure should
// still be efficient.
//
// Assume q = 64 and e < 0. The idea is to cut the number c * w = f * 2^e into
// two parts, which can be processed independently: An integral part p1, and a
// fractional part p2:
//
//      f * 2^e = ( (f div 2^-e) * 2^-e + (f mod 2^-e) ) * 2^e
//              = (f div 2^-e) + (f mod 2^-e) * 2^e
//              = p1 + p2 * 2^e
//
// The conversion of p1 into decimal form requires a series of divisions and
// modulos by (a power of) 10. These operations are faster for 32-bit than for
// 64-bit integers, so p1 should ideally fit into a 32-bit integer. This can be
// achieved by choosing
//
//      -e >= 32   or   e <= -32 := gamma
//
// In order to convert the fractional part
//
//      p2 * 2^e = d[-1] / 10^1 + d[-2] / 10^2 + ... + d[-k] / 10^k + ...
//
// into decimal form, the fraction is repeatedly multiplied by 10 and the digits
// d[-i] are extracted in order:
//
//      (10 * p2) div 2^-e = d[-1]
//      (10 * p2) mod 2^-e = d[-2] / 10^1 + ... + d[-k] / 10^(k-1) + ...
//
// The multiplication by 10 must not overflow. It is sufficient to choose
//
//      10 * p2 < 16 * p2 = 2^4 * p2 <= 2^64.
//
// Since p2 = f mod 2^-e < 2^-e,
//
//      -e <= 60   or   e >= -60 := alpha
//
// On the other hand, if multiplication by 100 does not overflow, two digits of
// p2 might be generated at a time. This requires
//
//      -e <= 64-7   or   e >= -64+7 =: alpha.
//
// This choice requires a slightly larger table of cached powers.
//
constexpr int const kAlpha = -57;
constexpr int const kGamma = -32;

//
// For IEEE double precision floating-point numbers v converted into
// normalized DiyFp's w = f * 2^e, still assuming q = 64,
//
//      e >= -1022      (min IEEE exponent)
//           -52        (IEEE significand size)
//           -52        (possibly normalize denormal IEEE numbers)
//           -11        (normalize the DiyFp)
//         = -1137
//
// and
//
//      e <= +1023      (max IEEE exponent)
//           -52        (IEEE significand size)
//           -11        (normalize the DiyFp)
//         = 960
//
// (For IEEE single precision the exponent range is [-196, 80].)
//
// Now
//
//      alpha <= e_c + e + q <= gamma
//          ==> f_c * 2^alpha <= c * 2^e * 2^q
//
// and since the c's are normalized, 2^(q-1) <= f_c,
//
//          ==> 2^(q - 1 + alpha) <= c * 2^(e + q)
//          ==> 2^(alpha - e - 1) <= c
//
// If c were an exakt power of ten, i.e. c = 10^k, one may determine k as
//
//      k = ceil( log_10( 2^(alpha - e - 1) ) )
//        = ceil( (alpha - e - 1) * log_10(2) )
//
// (From the paper:)
// "In theory the result of the procedure could be wrong since c is rounded,
// and the computation itself is approximated [...]. In practice, however, this
// simple function is sufficient."
//
// The difference gamma - alpha determines the size of the table of precomputed
// powers: The difference of the decimal exponents of adjacent table entries
// must be less than or equal to
//
//      floor( (gamma - alpha) * log_10(2) ) = 7.
//

struct CachedPower { // c = f * 2^e ~= 10^k
    uint64_t f;
    int e;
    int k;
};

inline CachedPower GetCachedPowerForBinaryExponent(int e)
{
    // NB:
    // Actually this function returns c, such that -57 <= e_c + e + 64 <= -34.

    constexpr int const kCachedPowersSize = 91;
    constexpr int const kCachedPowersMinDecExp = -300;
    constexpr int const kCachedPowersDecStep = 7;

    static constexpr CachedPower const kCachedPowers[] = {
        { 0xAB70FE17C79AC6CA, -1060, -300 },
        { 0xCC5FC196FEFD7D0C, -1037, -293 },
        { 0xF3A20279ED56D48A, -1014, -286 },
        { 0x91376C36D99995BE,  -990, -279 },
        { 0xAD1C8EAB5EE43B67,  -967, -272 },
        { 0xCE5D73FF402D98E4,  -944, -265 },
        { 0xF6019DA07F549B2B,  -921, -258 },
        { 0x92A1958A7675175F,  -897, -251 },
        { 0xAECC49914078536D,  -874, -244 },
        { 0xD0601D8EFC57B08C,  -851, -237 },
        { 0xF867241C8CC6D4C1,  -828, -230 },
        { 0x940F4613AE5ED137,  -804, -223 },
        { 0xB080392CC4349DED,  -781, -216 },
        { 0xD267CAA862A12D67,  -758, -209 },
        { 0xFAD2A4B13D1B5D6C,  -735, -202 },
        { 0x9580869F0E7AAC0F,  -711, -195 },
        { 0xB23867FB2A35B28E,  -688, -188 },
        { 0xD47487CC8470652B,  -665, -181 },
        { 0xFD442E4688BD304B,  -642, -174 },
        { 0x96F5600F15A7B7E5,  -618, -167 },
        { 0xB3F4E093DB73A093,  -595, -160 },
        { 0xD686619BA27255A3,  -572, -153 },
        { 0xFFBBCFE994E5C620,  -549, -146 },
        { 0x986DDB5C6B3A76B8,  -525, -139 },
        { 0xB5B5ADA8AAFF80B8,  -502, -132 },
        { 0xD89D64D57A607745,  -479, -125 },
        { 0x811CCC668829B887,  -455, -118 },
        { 0x99EA0196163FA42E,  -432, -111 },
        { 0xB77ADA0617E3BBCB,  -409, -104 },
        { 0xDAB99E59958885C5,  -386,  -97 },
        { 0x825ECC24C8737830,  -362,  -90 },
        { 0x9B69DBE1B548CE7D,  -339,  -83 },
        { 0xB94470938FA89BCF,  -316,  -76 },
        { 0xDCDB1B2798182245,  -293,  -69 },
        { 0x83A3EEEEF9153E89,  -269,  -62 },
        { 0x9CED737BB6C4183D,  -246,  -55 },
        { 0xBB127C53B17EC159,  -223,  -48 },
        { 0xDF01E85F912E37A3,  -200,  -41 }, // ---> single precision
        { 0x84EC3C97DA624AB5,  -176,  -34 }, //
        { 0x9E74D1B791E07E48,  -153,  -27 }, //
        { 0xBCE5086492111AEB,  -130,  -20 }, //
        { 0xE12E13424BB40E13,  -107,  -13 }, //
        { 0x8637BD05AF6C69B6,   -83,   -6 }, //
        { 0xA000000000000000,   -60,    1 }, //
        { 0xBEBC200000000000,   -37,    8 }, //
        { 0xE35FA931A0000000,   -14,   15 }, //
        { 0x878678326EAC9000,    10,   22 }, //
        { 0xA18F07D736B90BE5,    33,   29 }, //
        { 0xC097CE7BC90715B3,    56,   36 }, //
        { 0xE596B7B0C643C719,    79,   43 }, // <--- single precision
        { 0x88D8762BF324CD10,   103,   50 },
        { 0xA321F2D7226895C8,   126,   57 },
        { 0xC2781F49FFCFA6D5,   149,   64 },
        { 0xE7D34C64A9C85D44,   172,   71 },
        { 0x8A2DBF142DFCC7AB,   196,   78 },
        { 0xA4B8CAB1A1563F52,   219,   85 },
        { 0xC45D1DF942711D9A,   242,   92 },
        { 0xEA1575143CF97227,   265,   99 },
        { 0x8B865B215899F46D,   289,  106 },
        { 0xA6539930BF6BFF46,   312,  113 },
        { 0xC646D63501A1511E,   335,  120 },
        { 0xEC5D3FA8CE427B00,   358,  127 },
        { 0x8CE2529E2734BB1D,   382,  134 },
        { 0xA7F26836F282B733,   405,  141 },
        { 0xC83553C5C8965D3D,   428,  148 },
        { 0xEEAABA2E5DBF6785,   451,  155 },
        { 0x8E41ADE9FBEBC27D,   475,  162 },
        { 0xA99541BF57452B28,   498,  169 },
        { 0xCA28A291859BBF93,   521,  176 },
        { 0xF0FDF2D3F3C30B9F,   544,  183 },
        { 0x8FA475791A569D11,   568,  190 },
        { 0xAB3C2FDDEEAAD25B,   591,  197 },
        { 0xCC20CE9BD35C78A5,   614,  204 },
        { 0xF356F7EBF83552FE,   637,  211 },
        { 0x910AB1D4DB9914A0,   661,  218 },
        { 0xACE73CBFDC0BFB7B,   684,  225 },
        { 0xCE1DE40642E3F4B9,   707,  232 },
        { 0xF5B5D7EC8ACB58A3,   730,  239 },
        { 0x92746B9BE2F8552C,   754,  246 },
        { 0xAE9672ABA3D0C321,   777,  253 },
        { 0xD01FEF10A657842C,   800,  260 },
        { 0xF81AA16FDC1B81DB,   823,  267 },
        { 0x93E1AB8252F33B46,   847,  274 },
        { 0xB049DC016ABC5E60,   870,  281 },
        { 0xD226FC195C6A2F8C,   893,  288 },
        { 0xFA856334878FC151,   916,  295 },
        { 0x95527A5202DF0CCB,   940,  302 },
        { 0xB201833B35D63F73,   963,  309 },
        { 0xD433179D9C8CB841,   986,  316 },
        { 0xFCF62C1DEE382C42,  1009,  323 },
        { 0x96C6E0EAB509E64D,  1033,  330 },
    };

    //
    // This computation gives exactly the same results for k as
    //
    //      k = ceil((kAlpha - e - 1) * 0.30102999566398114)
    //
    // for |e| <= 1500, but doesn't require floating-point operations.
    //
    // NB: log_10(2) ~= 78913 / 2^18
    //
    assert(e >= -1500);
    assert(e <=  1500);
    int const f = kAlpha - e - 1;
    int const k = (f * 78913) / (1 << 18) + (f > 0);

    int const index = (-kCachedPowersMinDecExp + k + (kCachedPowersDecStep - 1)) / kCachedPowersDecStep;
    assert(index >= 0);
    assert(index < kCachedPowersSize);
    static_cast<void>(kCachedPowersSize); // Fix warning.

    CachedPower const cached = kCachedPowers[index];
    assert(kAlpha <= cached.e + e + 64);
    assert(kGamma >= cached.e + e + 64);

    return cached;
}

inline char* Itoa100(char* buf, uint32_t digits)
{
    static constexpr char const* const kDigits100 =
        "00010203040506070809"
        "10111213141516171819"
        "20212223242526272829"
        "30313233343536373839"
        "40414243444546474849"
        "50515253545556575859"
        "60616263646566676869"
        "70717273747576777879"
        "80818283848586878889"
        "90919293949596979899";

    assert(digits < 100);
    std::memcpy(buf, kDigits100 + 2*digits, 2);
    return buf + 2;
}

GRISU2_INLINE char* Itoa_32(char* buf, uint32_t n)
{
    uint32_t q;

    if (n >= 1000000000)
    {
//L_10_digits:
        q = n / 100000000;
        n = n % 100000000;
        buf = Itoa100(buf, q);
L_8_digits:
        q = n / 1000000;
        n = n % 1000000;
        buf = Itoa100(buf, q);
L_6_digits:
        q = n / 10000;
        n = n % 10000;
        buf = Itoa100(buf, q);
L_4_digits:
        q = n / 100;
        n = n % 100;
        buf = Itoa100(buf, q);
//L_2_digits:
        return Itoa100(buf, n);
    }

    if (n < 100) {
        if (n >= 10) {
            return Itoa100(buf, n);
        } else {
            buf[0] = static_cast<char>('0' + n);
            return buf + 1;
        }
    }

    if (n < 10000) {
        if (n >= 1000)
            goto L_4_digits;
        else
            goto L_3_digits;
    }

    if (n < 1000000) {
        if (n >= 100000)
            goto L_6_digits;
        else
            goto L_5_digits;
    }

    if (n < 100000000) {
        if (n >= 10000000)
            goto L_8_digits;
        else
            goto L_7_digits;
    }

    //if (n >= 100000000)
    {
//L_9_digits:
        q = n / 10000000;
        n = n % 10000000;
        buf = Itoa100(buf, q);
L_7_digits:
        q = n / 100000;
        n = n % 100000;
        buf = Itoa100(buf, q);
L_5_digits:
        q = n / 1000;
        n = n % 1000;
        buf = Itoa100(buf, q);
L_3_digits:
        q = n / 10;
        n = n % 10;
        buf = Itoa100(buf, q);
//L_1_digit:
        buf[0] = static_cast<char>('0' + n);
        return buf + 1;
    }
}

#if GRISU2_ROUND
GRISU2_INLINE void Grisu2Round(char* buf, int len, uint64_t dist, uint64_t delta, uint64_t rest, uint64_t ten_k)
{
    assert(len >= 1);
    assert(dist <= delta);
    assert(rest <= delta);
    assert(ten_k > 0);

    //
    //               <--------------------------- delta ---->
    //                                  <---- dist --------->
    // --------------[------------------+-------------------]--------------
    //               w-                 w                   w+
    //
    //                                  ten_k
    //                                <------>
    //                                       <---- rest ---->
    // --------------[------------------+----+--------------]--------------
    //                                  w    V
    //                                       = buf * 10^k
    //
    // ulp represents a unit-in-the-last-place in the decimal representation
    // stored in buf.
    // Decrement buf by ulp while this takes buf closer to w.
    //
    // The tests are written in this order to avoid overflow in unsigned
    // integer arithmetic.
    //

    while (rest < dist
        && delta - rest >= ten_k
        && (rest + ten_k < dist || dist - rest > rest + ten_k - dist))
    {
        assert(buf[len - 1] != '0');
        buf[len - 1]--;
        rest += ten_k;
    }
}
#endif

#if GRISU2_ROUND
GRISU2_INLINE void Grisu2DigitGen(char* buffer, int& length, int& decimal_exponent, Fp M_minus, Fp w, Fp M_plus)
#else
GRISU2_INLINE void Grisu2DigitGen(char* buffer, int& length, int& decimal_exponent, Fp M_minus, Fp M_plus)
#endif
{
    //
    // Generates the digits (and the exponent) of a decimal floating-point
    // number V in the range [w-, w+].
    //
    //               <--------------------------- delta ---->
    //                                  <---- dist --------->
    // --------------[------------------+-------------------]--------------
    //               w-                 w                   w+
    //
    // Instead of generating the digits of w, Grisu2 generates the digits
    // of w+ from left to right and stops as soon as V is in [w-,w+].
    //

    static_assert(Fp::kPrecision == 64, "invalid configuration");

    static_assert(kAlpha >= -60, "invalid parameter");
    static_assert(kGamma <= -32, "invalid parameter");

    assert(M_plus.e >= kAlpha);
    assert(M_plus.e <= kGamma);

    uint64_t delta = Fp::Sub(M_plus, M_minus).f; // (significand of (w+ - w-), implicit exponent is e)
#if GRISU2_ROUND
    uint64_t dist  = Fp::Sub(M_plus, w      ).f; // (significand of (w+ - w ), implicit exponent is e)
#endif

    //
    // Split w+ = f * 2^e into two parts p1 and p2 (note: e < 0):
    //
    //      w+ = f * 2^e
    //         = ((f div 2^-e) * 2^-e + (f mod 2^-e)) * 2^e
    //         = ((p1        ) * 2^-e + (p2        )) * 2^e
    //         = p1 + p2 * 2^e
    //

    int      const neg_e = -M_plus.e;
    uint64_t const mod_e = (uint64_t{1} << -M_plus.e) - 1;

    uint32_t p1 = static_cast<uint32_t>(M_plus.f >> neg_e); // p1 = f div 2^-e (Since -e >= 32, p1 fits into a 32-bit int.)
    uint64_t p2 = M_plus.f & mod_e;                         // p2 = f mod 2^-e

    //
    // 1.
    // Generate the digits of the integral part p1 = d[n-1]...d[1]d[0]
    //

    // Since w+ is normalized (f >= 2^(64-1)) and e >= -60, p1 > 0.
    assert(p1 > 0);

    //
    // Now
    //
    //      10^(k-1) <= p1 < 10^k, pow10 = 10^(k-1)
    //
    //      p1 = (p1 div 10^(k-1)) * 10^(k-1) + (p1 mod 10^(k-1))
    //         = (d[k-1]         ) * 10^(k-1) + (p1 mod 10^(k-1))
    //
    //      w+ = p1                                             + p2 * 2^e
    //         = d[k-1] * 10^(k-1) + (p1 mod 10^(k-1))          + p2 * 2^e
    //         = d[k-1] * 10^(k-1) + ((p1 mod 10^(k-1)) * 2^-e + p2) * 2^e
    //         = d[k-1] * 10^(k-1) + (                         rest) * 2^e
    //
    // Now generate the digits d[n] of p1 from left to right (n = k-1,...,0)
    //
    //      p1 = d[k-1]...d[n] * 10^n + d[n-1]...d[0]
    //
    // but stop as soon as
    //
    //      rest * 2^e = (d[n-1]...d[0] * 2^-e + p2) * 2^e <= delta * 2^e
    //
    // The common case here is that all the digits of p1 are required, so
    // optimize for this case and correct the length of the buffer if we have
    // generated too many digits.
    //

    char const* last = Itoa_32(buffer, p1);
    length = static_cast<int>(last - buffer);

    if (p2 <= delta)
    {
        // In this case: Too many digits of p1 might have been generated.
        //
        // Find the largest 0 <= n < k, such that
        //
        //      w+ = (p1 div 10^n) * 10^n + ((p1 mod 10^n) * 2^-e + p2) * 2^e
        //         = (p1 div 10^n) * 10^n + (                     rest) * 2^e
        //
        // and rest <= delta.
        //
        // Compute rest * 2^e = w+ mod 10^n = p1 + p2 * 2^e = (p1 * 2^-e + p2) * 2^e
        // and check if enough digits have been generated:
        //
        //      rest * 2^e <= delta * 2^e
        //
        // This test can be slightly simplified, since
        //
        //      rest = (p1 mod 10^n) * 2^-e + p2 <= delta
        //      <==>    r * 2^-e + p2 <= delta
        //      <==>    r * 2^-e      <= delta - p2 = D = D1 * 2^-e + D2
        //      <==>    r < D1 or (r == D1 and 0 <= D2)
        //      <==>    r <= D1
        //

        uint32_t const D1 = static_cast<uint32_t>((delta - p2) >> neg_e);

        int k = length;
        int n = 0;

        uint32_t r = 0;
        uint32_t pow10 = 1; // 10^n
        for (;;)
        {
            assert(k >= n + 1);
            assert(r <= D1);
            assert(n <= 9);
            assert(static_cast<uint32_t>(buffer[k - (n + 1)] - '0') <= UINT32_MAX / pow10);

            uint32_t r_next = pow10 * static_cast<uint32_t>(buffer[k - (n + 1)] - '0') + r;
            if (r_next > D1)
                break;
            r = r_next;
            n++;
            pow10 *= 10;
        }
        length = k - n;

        //
        // Found V = buffer * 10^n, with w- <= V <= w+.
        //
        decimal_exponent += n;

#if GRISU2_ROUND
        uint64_t const rest = (uint64_t{r} << neg_e) + p2;
        assert(rest <= delta);

        //
        // We may now just stop. But instead look if the buffer could be
        // decremented to bring V closer to w.
        //
        // pow10 = 10^n is now 1 ulp in the decimal representation V.
        // The rounding procedure works with DiyFp's with an implicit
        // exponent of e.
        //
        //      10^n = (10^n * 2^-e) * 2^e = ulp * 2^e
        //
        uint64_t const ten_n = uint64_t{pow10} << neg_e;
        Grisu2Round(buffer, length, dist, delta, rest, ten_n);
#endif
        return;
    }

    //
    // 2.
    // The digits of the integral part have been generated:
    //
    //      w+ = d[k-1]...d[1]d[0] + p2 * 2^e
    //         = buffer            + p2 * 2^e
    //
    // Now generate the digits of the fractional part p2 * 2^e.
    //
    // Note:
    // No decimal point is generated: the exponent is adjusted instead.
    //
    // p2 actually represents the fraction
    //
    //      p2 * 2^e
    //          = p2 / 2^-e
    //          = d[-1] / 10^1 + d[-2] / 10^2 + ...
    //
    // Now generate the digits d[-m] of p1 from left to right (m = 1,2,...)
    //
    //      p2 * 2^e = d[-1]d[-2]...d[-m] * 10^-m
    //                      + 10^-m * (d[-m-1] / 10^1 + d[-m-2] / 10^2 + ...)
    //
    // using
    //
    //      10^m * p2 = ((10^m * p2) div 2^-e) * 2^-e + ((10^m * p2) mod 2^-e)
    //                = (                   d) * 2^-e + (                   r)
    //
    // or
    //      10^m * p2 * 2^e = d + r * 2^e
    //
    // i.e.
    //
    //      w+ = buffer + p2 * 2^e
    //         = buffer + 10^-m * (d + r * 2^e)
    //         = (buffer * 10^m + d) * 10^-m + 10^-m * r * 2^e
    //
    // and stop as soon as 10^-m * r * 2^e <= delta * 2^e
    //

    assert(p2 > delta);
    // (otherwise the loop above would have been exited with rest <= delta)

    int m = 0;
    for (;;)
    {
        //
        // Invariant:
        //      w+ = buffer * 10^-m + 10^-m * (d[-m-1] / 10 + d[-m-2] / 10^2 + ...) * 2^e
        //         = buffer * 10^-m + 10^-m * (p2                                 ) * 2^e
        //         = buffer * 10^-m + 10^-m * (1/100 * (100 * p2)                 ) * 2^e
        //         = buffer * 10^-m + 10^-m * (1/100 * ((100*p2 div 2^-e) * 2^-e + (100*p2 mod 2^-e)) * 2^e
        //

        uint64_t p2_0 = p2;

        assert(p2 <= UINT64_MAX / 100);
        p2 *= 100;

        uint64_t const d = p2 >> neg_e; // d = (100 * p2) div 2^-e
        uint64_t const r = p2 & mod_e;  // r = (100 * p2) mod 2^-e
        //
        //      w+ = buffer * 10^-m + 10^-m * (1/100 * (d * 2^-e + r) * 2^e
        //         = buffer * 10^-m + 10^-m * (1/100 * (d + r * 2^e))
        //         = (buffer * 100 + d) * 10^(-m-2) + 10^(-m-2) * r * 2^e
        //
        Itoa100(buffer + length, static_cast<uint32_t>(d)); // buffer := buffer * 100 + d
        length += 2;
        //
        //      w+ = buffer * 10^(-m-2) + 10^(-m-2) * r * 2^e
        //
        p2 = r;
        m += 2;
        //
        //      w+ = buffer * 10^-m + 10^-m * p2 * 2^e
        //
        // Invariant restored.
        //

        //
        // Check if enough digits have been generated.
        // Compute
        //
        //      10^-m * p2 * 2^e <= delta * 2^e
        //              p2 * 2^e <= 10^m * delta * 2^e
        //                    p2 <= 10^m * delta
        //
        delta *= 100;
#if GRISU2_ROUND
        dist  *= 100;
#endif

        if (p2 <= delta)
        {
            // Almost done.
            // Check if we have generated one digit too much.

            uint64_t const r10 = (10 * p2_0) & mod_e;
            uint64_t const delta10 = delta / 10;

            if (r10 <= delta10) // Only one digit required.
            {
                length--;
                p2 = r10;
                m--;
                delta = delta10;
#if GRISU2_ROUND
                dist /= 10;
#endif
            }

            break;
        }
    }

    decimal_exponent -= m;

#if GRISU2_ROUND
    //
    // 1 ulp in the decimal representation is now 10^-m.
    // Since delta and dist are now scaled by 10^m, we need to do the
    // same with ulp in order to keep the units in sync.
    //
    //      10^m * 10^-m = 1 = 2^-e * 2^e = ten_m * 2^e
    //
    uint64_t const ten_m = uint64_t{1} << neg_e;
    Grisu2Round(buffer, length, dist, delta, p2, ten_m);
#endif

    //
    // By construction this algorithm generates the shortest possible decimal
    // number (Loitsch, Theorem 6.2) which rounds back to w.
    // For an input number of precision p, at least
    //
    //      N = 1 + ceil(p * log_10(2))
    //
    // decimal digits are sufficient to identify all binary floating-point
    // numbers (Matula, "In-and-Out conversions").
    // This implies that the algorithm does not produce more than N decimal
    // digits.
    //
    //      N = 17 for p = 53 (IEEE double precision)
    //      N = 9  for p = 24 (IEEE single precision)
    //
    assert(length <= 17);
}

// v = buf * 10^decimal_exponent
// len is the length of the buffer (number of decimal digits)
#if GRISU2_ROUND
GRISU2_INLINE void Grisu2(char* buf, int& len, int& decimal_exponent, Fp m_minus, Fp v, Fp m_plus)
#else
GRISU2_INLINE void Grisu2(char* buf, int& len, int& decimal_exponent, Fp m_minus, Fp m_plus)
#endif
{
    assert(m_plus.e == m_minus.e);
#if GRISU2_ROUND
    assert(m_plus.e == v.e);
#endif

    //
    //  --------(-----------------------+-----------------------)--------    (A)
    //          m-                      v                       m+
    //
    //  --------------------(-----------+-----------------------)--------    (B)
    //                      m-          v                       m+
    //
    // First scale v (and m- and m+) such that the exponent is in the range
    // [alpha, gamma].
    //

    CachedPower const cached = GetCachedPowerForBinaryExponent(m_plus.e);

    Fp const c_minus_k(cached.f, cached.e); // = c ~= 10^k

    // The exponent of the products is v.e + c_minus_k.e + q
#if GRISU2_ROUND
    Fp const w       = Fp::Mul(v,       c_minus_k);
#endif
    Fp const w_minus = Fp::Mul(m_minus, c_minus_k);
    Fp const w_plus  = Fp::Mul(m_plus,  c_minus_k);

    //
    //  ----(---+---)---------------(---+---)---------------(---+---)----
    //          w-                      w                       w+
    //          = c*m-                  = c*v                   = c*m+
    //
    // Fp::Mul rounds its result and c_minus_k is approximated too. w, w- and
    // w+ are now off by a small amount.
    // In fact:
    //
    //      w - v * 10^k < 1 ulp
    //
    // To account for this inaccuracy, add resp. subtract 1 ulp.
    //
    //  --------+---[---------------(---+---)---------------]---+--------
    //          w-  M-                  w                   M+  w+
    //
    // Now any number in [M-, M+] (bounds included) will round to w when input,
    // regardless of how the input rounding algorithm breaks ties.
    //
    // And DigitGen generates the shortest possible such number in [M-, M+].
    // Note that this does not mean that Grisu2 always generates the shortest
    // possible number in the interval (m-, m+).
    //
    Fp const M_minus = Fp(w_minus.f + 1, w_minus.e);
    Fp const M_plus  = Fp(w_plus.f  - 1, w_plus.e );

    decimal_exponent = -cached.k; // = -(-k) = k

#if GRISU2_ROUND
    Grisu2DigitGen(buf, len, decimal_exponent, M_minus, w, M_plus);
#else
    Grisu2DigitGen(buf, len, decimal_exponent, M_minus, M_plus);
#endif
}

// Returns a pointer to the element following the exponent
GRISU2_INLINE char* AppendExponent(char* buf, int e)
{
    assert(e > -1000);
    assert(e <  1000);

    if (e < 0)
    {
        e = -e;
        *buf++ = '-';
    }
    else
    {
        *buf++ = '+';
    }

    uint32_t k = static_cast<uint32_t>(e);
    if (k < 10)
    {
        buf[0] = static_cast<char>('0' + k);
        return buf + 1;
    }

    if (k < 100)
        return Itoa100(buf, k);

    uint32_t q = k / 100;
    uint32_t r = k % 100;
    buf[0] = static_cast<char>('0' + q);
    return Itoa100(buf + 1, r);
}

GRISU2_INLINE char* FormatBuffer(char* buf, int k, int n)
{
    // v = digits * 10^(n-k)
    // k is the length of the buffer (number of decimal digits)
    // n is the position of the decimal point relative to the start of the buffer.
    //
    // Format the decimal floating-number v in the same way as JavaScript's ToString
    // applied to number type.
    //
    // See:
    // https://tc39.github.io/ecma262/#sec-tostring-applied-to-the-number-type

    if (k <= n && n <= 21)
    {
        // digits[000]

        std::memset(buf + k, '0', static_cast<size_t>(n - k));
        //if (trailing_dot_zero)
        //{
        //    buf[n++] = '.';
        //    buf[n++] = '0';
        //}
        return buf + n;
    }

    if (0 < n && n <= 21)
    {
        // dig.its

        assert(k > n);

        std::memmove(buf + (n + 1), buf + n, static_cast<size_t>(k - n));
        buf[n] = '.';
        return buf + (k + 1);
    }

    if (-6 < n && n <= 0)
    {
        // 0.[000]digits

        std::memmove(buf + (2 + -n), buf, static_cast<size_t>(k));
        buf[0] = '0';
        buf[1] = '.';
        std::memset(buf + 2, '0', static_cast<size_t>(-n));
        return buf + (2 + (-n) + k);
    }

    if (k == 1)
    {
        // dE+123

        buf += 1;
    }
    else
    {
        // d.igitsE+123

        std::memmove(buf + 2, buf + 1, static_cast<size_t>(k - 1));
        buf[1] = '.';
        buf += 1 + k;
    }

    *buf++ = 'e';
    return AppendExponent(buf, n - 1);
}

inline char* StrCopy_unsafe(char* dst, char const* src)
{
    auto const len = strlen(src);
    std::memcpy(dst, src, len);
    return dst + len;
}

//
// Generates a decimal representation of the input floating-point number V in
// BUF.
//
// The result is formatted like JavaScript's ToString applied to a number type.
// Except that:
// An argument representing an infinity is converted to "inf" or "-inf".
// An argument representing a NaN is converted to "nan".
//
// This function never writes more than 25 characters to BUF and returns an
// iterator pointing past-the-end of the decimal representation.
// The result is guaranteed to round-trip (when read back by a correctly
// rounding implementation.)
//
// Note:
// The result is not null-terminated.
//
template <typename Float>
char* ToString(char* next, char* last, Float value)
{
    static constexpr char const* const kNaNString = "NaN";      // assert len <= 25
    static constexpr char const* const kInfString = "Infinity"; // assert len <= 24

    using IEEEType = IEEEFloat<Float>;
    static_assert(Fp::kPrecision >= IEEEType::kPrecision + 3, "insufficient precision");

    assert(last - next >= 25);
    static_cast<void>(last); // unused

    IEEEType const v(value);
    //assert(!v.IsNaN());
    //assert(!v.IsInf());

    if (v.IsNaN())
        return StrCopy_unsafe(next, kNaNString);

    if (v.IsNegative())
        *next++ = '-';

    if (v.IsInf())
        return StrCopy_unsafe(next, kInfString);

    if (v.IsZero())
    {
        *next++ = '0';
        //if (trailing_dot_zero)
        //{
        //    *next++ = '.';
        //    *next++ = '0';
        //}
        return next;
    }

    FpBoundaries const w = ComputeBoundaries(v.Abs());

    // Compute v = buffer * 10^decimal_exponent.
    // The decimal digits are stored in the buffer, which needs to be
    // interpreted as an unsigned decimal integer.
    // len is the length of the buffer, i.e. the number of decimal digits
    int len = 0;
    int decimal_exponent = 0;
#if GRISU2_ROUND
    Grisu2(next, len, decimal_exponent, w.minus, w.w, w.plus);
#else
    Grisu2(next, len, decimal_exponent, w.minus, w.plus);
#endif

    // Compute the position of the decimal point relative to the start of the buffer.
    int const n = decimal_exponent + len;

    return FormatBuffer(next, len, n);
}

} // namespace fast_dtoa

// http://florian.loitsch.com/publications (bench.tar.gz)
//
// Copyright (c) 2009 Florian Loitsch
//
//   Permission is hereby granted, free of charge, to any person
//   obtaining a copy of this software and associated documentation
//   files (the "Software"), to deal in the Software without
//   restriction, including without limitation the rights to use,
//   copy, modify, merge, publish, distribute, sublicense, and/or sell
//   copies of the Software, and to permit persons to whom the
//   Software is furnished to do so, subject to the following
//   conditions:
//
//   The above copyright notice and this permission notice shall be
//   included in all copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
//   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
//   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
//   OTHER DEALINGS IN THE SOFTWARE.
