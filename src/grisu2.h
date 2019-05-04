// Copyright 2019 Alexander Bolz
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include "format_digits.h"
#include "ieee.h"

#include <cassert>
#include <cstdint>
#include <cstring>
#include <limits>
#ifdef _MSC_VER
#include <intrin.h>
#endif

#ifndef GRISU_ASSERT
#define GRISU_ASSERT(X) assert(X)
#endif

#ifndef GRISU_FORCE_INLINE_ATTR
#if __GNUC__
#define GRISU_FORCE_INLINE_ATTR __attribute__((always_inline)) inline
#elif _MSC_VER
#define GRISU_FORCE_INLINE_ATTR __forceinline
#else
#define GRISU_FORCE_INLINE_ATTR inline
#endif
#endif

#ifndef GRISU_INLINE
#define GRISU_INLINE inline
#endif

#ifndef GRISU_FORCE_INLINE
#define GRISU_FORCE_INLINE GRISU_FORCE_INLINE_ATTR
#endif

#ifndef GRISU_SMALL_INT_OPTIMIZATION
#define GRISU_SMALL_INT_OPTIMIZATION 0
#endif

#ifndef GRISU_ROUND
#define GRISU_ROUND 0
#endif

namespace grisu2 {

//==================================================================================================
//
//==================================================================================================

template <typename Float>
struct ToDecimalResult;

template <> struct ToDecimalResult<double> {
    uint64_t digits; // num_digits <= 17
    int exponent;
};

template <> struct ToDecimalResult<float> {
    uint32_t digits; // num_digits <= 9
    int exponent;
};

//==================================================================================================
// DoubleToDigits
//
// Implements the Grisu2 algorithm for (IEEE) binary to decimal floating-point conversion.
//
// References:
//
// [1]  Loitsch, "Printing Floating-Point Numbers Quickly and Accurately with Integers",
//      Proceedings of the ACM SIGPLAN 2010 Conference on Programming Language Design and Implementation, PLDI 2010
// [2]  Burger, Dybvig, "Printing Floating-Point Numbers Quickly and Accurately",
//      Proceedings of the ACM SIGPLAN 1996 Conference on Programming Language Design and Implementation, PLDI 1996
//==================================================================================================
// Constant data: 632 + 200 = 832 bytes

//
// TODO:
// Clean up comments...
//

namespace impl {

struct DiyFp // f * 2^e
{
    static constexpr int SignificandSize = 64; // = q

    uint64_t f = 0;
    int e = 0;

    constexpr DiyFp() = default;
    constexpr DiyFp(uint64_t f_, int e_) : f(f_), e(e_) {}
};

// Returns x - y.
// PRE: x.e == y.e and x.f >= y.f
GRISU_FORCE_INLINE DiyFp Subtract(DiyFp x, DiyFp y)
{
    GRISU_ASSERT(x.e == y.e);
    GRISU_ASSERT(x.f >= y.f);

    return DiyFp(x.f - y.f, x.e);
}

// Returns x * y.
// The result is rounded (ties up). (Only the upper q bits are returned.)
GRISU_FORCE_INLINE DiyFp Multiply(DiyFp x, DiyFp y)
{
    static_assert(DiyFp::SignificandSize == 64, "internal error");

    // Computes:
    //  f = round((x.f * y.f) / 2^q)
    //  e = x.e + y.e + q

#if defined(__SIZEOF_INT128__)
    __extension__ using Uint128 = unsigned __int128;

    const Uint128 p = Uint128{x.f} * Uint128{y.f};

    uint64_t h = static_cast<uint64_t>(p >> 64);
    uint64_t l = static_cast<uint64_t>(p);
    h += l >> 63; // round, ties up: [h, l] += 2^q / 2

    return DiyFp(h, x.e + y.e + 64);
#elif defined(_MSC_VER) && defined(_M_X64)
    uint64_t h = 0;
    uint64_t l = _umul128(x.f, y.f, &h);
    h += l >> 63; // round, ties up: [h, l] += 2^q / 2

    return DiyFp(h, x.e + y.e + 64);
#else
    const uint32_t x_lo = static_cast<uint32_t>(x.f);
    const uint32_t x_hi = static_cast<uint32_t>(x.f >> 32);
    const uint32_t y_lo = static_cast<uint32_t>(y.f);
    const uint32_t y_hi = static_cast<uint32_t>(y.f >> 32);

    const uint64_t b00 = uint64_t{x_lo} * y_lo;
    const uint64_t b01 = uint64_t{x_lo} * y_hi;
    const uint64_t b10 = uint64_t{x_hi} * y_lo;
    const uint64_t b11 = uint64_t{x_hi} * y_hi;

    const uint32_t b00_hi = static_cast<uint32_t>(b00 >> 32);

    const uint64_t mid1 = b10 + b00_hi;
    const uint32_t mid1_lo = static_cast<uint32_t>(mid1);
    const uint32_t mid1_hi = static_cast<uint32_t>(mid1 >> 32);

    const uint64_t mid2 = b01 + mid1_lo;
    const uint32_t mid2_lo = static_cast<uint32_t>(mid2);
    const uint32_t mid2_hi = static_cast<uint32_t>(mid2 >> 32);

    // NB: mid2_lo has the upper 32 bits of the low part of the product.
    const uint32_t r = mid2_lo >> 31;
    const uint64_t h = b11 + mid1_hi + mid2_hi + r;

    return DiyFp(h, x.e + y.e + 64);
#endif
}

// Returns the number of leading 0-bits in x, starting at the most significant bit position.
// If x is 0, the result is undefined.
GRISU_FORCE_INLINE int CountLeadingZeros64(uint64_t x)
{
    GRISU_ASSERT(x != 0);

#if defined(__GNUC__)
    return __builtin_clzll(x);
#elif defined(_MSC_VER) && (defined(_M_ARM) || defined(_M_ARM64))
    return static_cast<int>(_CountLeadingZeros64(x));
#elif defined(_MSC_VER) && defined(_M_X64)
    return static_cast<int>(__lzcnt64(x));
#elif defined(_MSC_VER) && defined(_M_IX86)
    int lz = static_cast<int>( __lzcnt(static_cast<uint32_t>(x >> 32)) );
    if (lz == 32) {
        lz += static_cast<int>( __lzcnt(static_cast<uint32_t>(x)) );
    }
    return lz;
#else
    int lz = 0;
    while ((x >> 63) == 0) {
        x <<= 1;
        ++lz;
    }
    return lz;
#endif
}

// Normalize x such that the significand is >= 2^(q-1).
// PRE: x.f != 0
GRISU_FORCE_INLINE DiyFp Normalize(DiyFp x)
{
    static_assert(DiyFp::SignificandSize == 64, "internal error");

    // For double: lz >= 64-53 = 11
    // For single: lz >= 64-24 = 40
    const int lz = CountLeadingZeros64(x.f);
    return DiyFp(x.f << lz, x.e - lz);
}

// Normalize x such that the result has the exponent E.
// PRE: e >= x.e and the upper e - x.e bits of x.f must be zero.
GRISU_FORCE_INLINE DiyFp NormalizeTo(DiyFp x, int e)
{
    const int delta = x.e - e;

    GRISU_ASSERT(delta >= 0);
    GRISU_ASSERT(((x.f << delta) >> delta) == x.f);

    return DiyFp(x.f << delta, e);
}

// Decomposes `value` into `f * 2^e`.
// The result is not normalized.
// PRE: `value` must be finite and non-negative, i.e. >= +0.0.
template <typename Float>
GRISU_FORCE_INLINE DiyFp DiyFpFromFloat(Float value)
{
    using Fp = dtoa::IEEE<Float>;

    const Fp v(value);
    GRISU_ASSERT(v.IsFinite());
    GRISU_ASSERT(!v.SignBit());

    const auto F = v.PhysicalSignificand();
    const auto E = v.PhysicalExponent();

    // If v is denormal:
    //      value = 0.F * 2^(1 - bias) = (          F) * 2^(1 - bias - (p-1))
    // If v is normalized:
    //      value = 1.F * 2^(E - bias) = (2^(p-1) + F) * 2^(E - bias - (p-1))

    return (E == 0) // denormal?
        ? DiyFp(F, Fp::MinExponent)
        : DiyFp(F | Fp::HiddenBit, static_cast<int>(E) - Fp::ExponentBias);
}

// Compute the boundaries m- and m+ of the floating-point value
// v = f * 2^e.
//
// Determine v- and v+, the floating-point predecessor and successor if v,
// respectively.
//
//      v- = v - 2^e        if f != 2^(p-1) or e == e_min                (A)
//         = v - 2^(e-1)    if f == 2^(p-1) and e > e_min                (B)
//
//      v+ = v + 2^e
//
// Let m- = (v- + v) / 2 and m+ = (v + v+) / 2. All real numbers _strictly_
// between m- and m+ round to v, regardless of how the input rounding
// algorithm breaks ties.
//
//      ---+-------------+-------------+-------------+-------------+---  (A)
//         v-            m-            v             m+            v+
//
//      -----------------+------+------+-------------+-------------+---  (B)
//                       v-     m-     v             m+            v+

struct Boundaries {
#if GRISU_ROUND
    DiyFp v;
#endif
    DiyFp m_minus;
    DiyFp m_plus;
};

// Compute the (normalized) DiyFp representing the input number 'value' and its
// boundaries.
// PRE: 'value' must be finite and positive
template <typename Float>
GRISU_FORCE_INLINE Boundaries ComputeBoundaries(Float value)
{
    using Fp = dtoa::IEEE<Float>;

    GRISU_ASSERT(Fp(value).IsFinite());
    GRISU_ASSERT(value > 0);

    const auto v = DiyFpFromFloat(value);

    // Compute the boundaries of v.
    const bool lower_boundary_is_closer = (v.f == Fp::HiddenBit && v.e > Fp::MinExponent);
    const auto m_minus = DiyFp(4*v.f - 2 + lower_boundary_is_closer, v.e - 2);
    const auto m_plus = DiyFp(4*v.f + 2, v.e - 2);

#if GRISU_ROUND
    // Determine the normalized w = v.
    const auto w = Normalize(v);

    // Determine the normalized w+ = m+.
    // Since e_(w+) == e_(w), one can use NormalizeTo instead of Normalize.
    const auto w_plus = NormalizeTo(m_plus, w.e);

    // Determine w- = m- such that e_(w-) = e_(w+).
    const auto w_minus = NormalizeTo(m_minus, w_plus.e);

    return {w, w_minus, w_plus};
#else
    // Determine the normalized w+ = m+.
    const auto w_plus = Normalize(m_plus);

    // Determine w- = m- such that e_(w-) = e_(w+).
    const auto w_minus = NormalizeTo(m_minus, w_plus.e);

    return {w_minus, w_plus};
#endif
}

// Given normalized DiyFp w, Grisu needs to find a (normalized) cached
// power-of-ten c, such that the exponent of the product c * w = f * 2^e lies
// within a certain range [alpha, gamma] (Definition 3.2 from [1])
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
// The choice of (alpha,gamma) determines the size of the table and the form of
// the digit generation procedure.
//
// Using (alpha,gamma)=(-60,-32) (or any smaller interval contained int [-60,-32]
// works out well in practice:
//
// The idea is to cut the number c * w = f * 2^e into two parts, which can be
// processed independently: An integral part p1, and a fractional part p2:
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
//      p2 * 2^e = p2 / 2^-e = d[-1] / 10^1 + d[-2] / 10^2 + ...
//
// into decimal form, the fraction is repeatedly multiplied by 10 and the digits
// d[-i] are extracted in order:
//
//      (10 * p2) div 2^-e = d[-1]
//      (10 * p2) mod 2^-e = d[-2] / 10^1 + ...
//
// The multiplication by 10 must not overflow. It is sufficient to choose
//
//      10 * p2 < 16 * p2 = 2^4 * p2 <= 2^64.
//
// Since p2 = f mod 2^-e < 2^-e,
//
//      -e <= 60   or   e >= -60 := alpha
//
// A larger alpha would allow to multiply the fractional part by a larger
// power of 10 and therefore allow to generate multiple digits in a single
// iteration.
//
// A smaller gamma would (probably) allow to generate a small amount of
// digits faster.

//
// TODO:
//
// Test different exponent ranges here, e.g., [-57, -43], and adjust DigitGen
// to produce 2 decimal digits per iteration in the (p2 > delta) branch.
//

// Now
//
//      alpha <= e_c + e + q <= gamma                                        (1)
//      ==> f_c * 2^alpha <= c * 2^e * 2^q
//
// and since the c's are normalized, 2^(q-1) <= f_c,
//
//      ==> 2^(q - 1 + alpha) <= c * 2^(e + q)
//      ==> 2^(alpha - e - 1) <= c
//
// If c were an exakt power of ten, i.e. c = 10^k, one may determine k as
//
//      k = ceil( log_10( 2^(alpha - e - 1) ) )
//        = ceil( (alpha - e - 1) * log_10(2) )
//
// From the paper:
// "In theory the result of the procedure could be wrong since c is rounded, and
//  the computation itself is approximated [...]. In practice, however, this
//  simple function is sufficient."
//
// For IEEE double precision floating-point numbers converted into normalized
// DiyFp's w = f * 2^e, with q = 64,
//
//      e >= -1022      (min IEEE exponent)
//           -52        (p - 1)
//           -52        (p - 1, possibly normalize denormal IEEE numbers)
//           -11        (normalize the DiyFp)
//         = -1137
//
// and
//
//      e <= +1023      (max IEEE exponent)
//           -52        (p - 1)
//           -11        (normalize the DiyFp)
//         = 960
//
// For IEEE single-precision the range is [-180, 96].
//
// This binary exponent range [-1137,960] results in a decimal exponent range
// [-307,324]. One does not need to store a cached power for each k in this
// range. For each such k it suffices to find a cached power such that the
// exponent of the product lies in [alpha,gamma].
// This implies that the difference of the decimal exponents of adjacent table
// entries must be less than or equal to
//
//      floor( (gamma - alpha) * log_10(2) ) = 8.
//
// (A smaller distance gamma-alpha would require a larger table.)

// Returns: floor(x / 2^n)
GRISU_FORCE_INLINE int SAR(int x, int n)
{
    // Technically, right-shift of negative integers is implementation defined...
    // Should easily get optimized into SAR (or equivalent) instruction.
#if 1
    return x < 0 ? ~(~x >> n) : (x >> n);
#else
    return x >> n;
#endif
}

// Returns: floor(log_2(10^e))
GRISU_FORCE_INLINE int FloorLog2Pow10(int e)
{
    GRISU_ASSERT(e >= -1233);
    GRISU_ASSERT(e <=  1232);
    return SAR(e * 1741647, 19);
}

// Returns: ceil(log_10(2^e))
GRISU_FORCE_INLINE int CeilLog10Pow2(int e)
{
    GRISU_ASSERT(e >= -2620);
    GRISU_ASSERT(e <=  2620);
    return SAR(e * 315653 + ((1 << 20) - 1), 20);
}

struct CachedPower { // c = f * 2^e ~= 10^k
    uint64_t f;
    int e; // binary exponent
    int k; // decimal exponent
};

constexpr int kAlpha = -50;
constexpr int kGamma = -36;

// For a normalized DiyFp w = f * 2^e, this function returns a (normalized)
// cached power-of-ten c = f_c * 2^e_c, such that the exponent of the product
// w * c satisfies
//
//      kAlpha <= e_c + e + q <= kGamma.
//
GRISU_FORCE_INLINE CachedPower GetCachedPowerForBinaryExponent(int e)
{
    static constexpr int kCachedPowersSize       =  159;
    static constexpr int kCachedPowersMinDecExp  = -304;
    static constexpr int kCachedPowersMaxDecExp  =  328;
    static constexpr int kCachedPowersDecExpStep =    4;

    static constexpr uint64_t kSignificands[] = {
        0x8C71DCD9BA0B4926, // e = -1073, k = -304
        0xAB70FE17C79AC6CA, // e = -1060, k = -300
        0xD1476E2C07286FAA, // e = -1047, k = -296
        0xFF77B1FCBEBCDC4F, // e = -1034, k = -292
        0x9BECCE62836AC577, // e = -1020, k = -288
        0xBE5691EF416BD60C, // e = -1007, k = -284
        0xE858AD248F5C22CA, // e =  -994, k = -280
        0x8DD01FAD907FFC3C, // e =  -980, k = -276
        0xAD1C8EAB5EE43B67, // e =  -967, k = -272
        0xD3515C2831559A83, // e =  -954, k = -268
        0x80FA687F881C7F8E, // e =  -940, k = -264
        0x9D71AC8FADA6C9B5, // e =  -927, k = -260
        0xC0314325637A193A, // e =  -914, k = -256
        0xEA9C227723EE8BCB, // e =  -901, k = -252
        0x8F31CC0937AE58D3, // e =  -887, k = -248
        0xAECC49914078536D, // e =  -874, k = -244
        0xD5605FCDCF32E1D7, // e =  -861, k = -240
        0x823C12795DB6CE57, // e =  -847, k = -236
        0x9EFA548D26E5A6E2, // e =  -834, k = -232
        0xC21094364DFB5637, // e =  -821, k = -228
        0xECE53CEC4A314EBE, // e =  -808, k = -224
        0x9096EA6F3848984F, // e =  -794, k = -220
        0xB080392CC4349DED, // e =  -781, k = -216
        0xD77485CB25823AC7, // e =  -768, k = -212
        0x8380DEA93DA4BC60, // e =  -754, k = -208
        0xA086CFCD97BF97F4, // e =  -741, k = -204
        0xC3F490AA77BD60FD, // e =  -728, k = -200
        0xEF340A98172AACE5, // e =  -715, k = -196
        0x91FF83775423CC06, // e =  -701, k = -192
        0xB23867FB2A35B28E, // e =  -688, k = -188
        0xD98DDAEE19068C76, // e =  -675, k = -184
        0x84C8D4DFD2C63F3B, // e =  -661, k = -180
        0xA21727DB38CB0030, // e =  -648, k = -176
        0xC5DD44271AD3CDBA, // e =  -635, k = -172
        0xF18899B1BC3F8CA2, // e =  -622, k = -168
        0x936B9FCEBB25C996, // e =  -608, k = -164
        0xB3F4E093DB73A093, // e =  -595, k = -160
        0xDBAC6C247D62A584, // e =  -582, k = -156
        0x8613FD0145877586, // e =  -568, k = -152
        0xA3AB66580D5FDAF6, // e =  -555, k = -148
        0xC7CABA6E7C5382C9, // e =  -542, k = -144
        0xF3E2F893DEC3F126, // e =  -529, k = -140
        0x94DB483840B717F0, // e =  -515, k = -136
        0xB5B5ADA8AAFF80B8, // e =  -502, k = -132
        0xDDD0467C64BCE4A1, // e =  -489, k = -128
        0x87625F056C7C4A8B, // e =  -475, k = -124
        0xA54394FE1EEDB8FF, // e =  -462, k = -120
        0xC9BCFF6034C13053, // e =  -449, k = -116
        0xF64335BCF065D37D, // e =  -436, k = -112
        0x964E858C91BA2655, // e =  -422, k = -108
        0xB77ADA0617E3BBCB, // e =  -409, k = -104
        0xDFF9772470297EBD, // e =  -396, k = -100
        0x88B402F7FD75539B, // e =  -382, k =  -96
        0xA6DFBD9FB8E5B88F, // e =  -369, k =  -92
        0xCBB41EF979346BCA, // e =  -356, k =  -88
        0xF8A95FCF88747D94, // e =  -343, k =  -84
        0x97C560BA6B0919A6, // e =  -329, k =  -80
        0xB94470938FA89BCF, // e =  -316, k =  -76
        0xE2280B6C20DD5232, // e =  -303, k =  -72
        0x8A08F0F8BF0F156B, // e =  -289, k =  -68
        0xA87FEA27A539E9A5, // e =  -276, k =  -64
        0xCDB02555653131B6, // e =  -263, k =  -60
        0xFB158592BE068D2F, // e =  -250, k =  -56
        0x993FE2C6D07B7FAC, // e =  -236, k =  -52
        0xBB127C53B17EC159, // e =  -223, k =  -48
        0xE45C10C42A2B3B06, // e =  -210, k =  -44
        0x8B61313BBABCE2C6, // e =  -196, k =  -40
        0xAA242499697392D3, // e =  -183, k =  -36
        0xCFB11EAD453994BA, // e =  -170, k =  -32
        0xFD87B5F28300CA0E, // e =  -157, k =  -28
        0x9ABE14CD44753B53, // e =  -143, k =  -24
        0xBCE5086492111AEB, // e =  -130, k =  -20
        0xE69594BEC44DE15B, // e =  -117, k =  -16
        0x8CBCCC096F5088CC, // e =  -103, k =  -12
        0xABCC77118461CEFD, // e =   -90, k =   -8
        0xD1B71758E219652C, // e =   -77, k =   -4
        0x8000000000000000, // e =   -63, k =    0
        0x9C40000000000000, // e =   -50, k =    4
        0xBEBC200000000000, // e =   -37, k =    8
        0xE8D4A51000000000, // e =   -24, k =   12
        0x8E1BC9BF04000000, // e =   -10, k =   16
        0xAD78EBC5AC620000, // e =     3, k =   20
        0xD3C21BCECCEDA100, // e =    16, k =   24
        0x813F3978F8940984, // e =    30, k =   28
        0x9DC5ADA82B70B59E, // e =    43, k =   32
        0xC097CE7BC90715B3, // e =    56, k =   36
        0xEB194F8E1AE525FD, // e =    69, k =   40
        0x8F7E32CE7BEA5C70, // e =    83, k =   44
        0xAF298D050E4395D7, // e =    96, k =   48
        0xD5D238A4ABE98068, // e =   109, k =   52
        0x82818F1281ED44A0, // e =   123, k =   56
        0x9F4F2726179A2245, // e =   136, k =   60
        0xC2781F49FFCFA6D5, // e =   149, k =   64
        0xED63A231D4C4FB27, // e =   162, k =   68
        0x90E40FBEEA1D3A4B, // e =   176, k =   72
        0xB0DE65388CC8ADA8, // e =   189, k =   76
        0xD7E77A8F87DAF7FC, // e =   202, k =   80
        0x83C7088E1AAB65DB, // e =   216, k =   84
        0xA0DC75F1778E39D6, // e =   229, k =   88
        0xC45D1DF942711D9A, // e =   242, k =   92
        0xEFB3AB16C59B14A3, // e =   255, k =   96
        0x924D692CA61BE758, // e =   269, k =  100
        0xB2977EE300C50FE7, // e =   282, k =  104
        0xDA01EE641A708DEA, // e =   295, k =  108
        0x850FADC09923329E, // e =   309, k =  112
        0xA26DA3999AEF774A, // e =   322, k =  116
        0xC646D63501A1511E, // e =   335, k =  120
        0xF209787BB47D6B85, // e =   348, k =  124
        0x93BA47C980E98CE0, // e =   362, k =  128
        0xB454E4A179DD1877, // e =   375, k =  132
        0xDC21A1171D42645D, // e =   388, k =  136
        0x865B86925B9BC5C2, // e =   402, k =  140
        0xA402B9C5A8D3A6E7, // e =   415, k =  144
        0xC83553C5C8965D3D, // e =   428, k =  148
        0xF46518C2EF5B8CD1, // e =   441, k =  152
        0x952AB45CFA97A0B3, // e =   455, k =  156
        0xB616A12B7FE617AA, // e =   468, k =  160
        0xDE469FBD99A05FE3, // e =   481, k =  164
        0x87AA9AFF79042287, // e =   495, k =  168
        0xA59BC234DB398C25, // e =   508, k =  172
        0xCA28A291859BBF93, // e =   521, k =  176
        0xF6C69A72A3989F5C, // e =   534, k =  180
        0x969EB7C47859E744, // e =   548, k =  184
        0xB7DCBF5354E9BECE, // e =   561, k =  188
        0xE070F78D3927556B, // e =   574, k =  192
        0x88FCF317F22241E2, // e =   588, k =  196
        0xA738C6BEBB12D16D, // e =   601, k =  200
        0xCC20CE9BD35C78A5, // e =   614, k =  204
        0xF92E0C3537826146, // e =   627, k =  208
        0x98165AF37B2153DF, // e =   641, k =  212
        0xB9A74A0637CE2EE1, // e =   654, k =  216
        0xE2A0B5DC971F303A, // e =   667, k =  220
        0x8A5296FFE33CC930, // e =   681, k =  224
        0xA8D9D1535CE3B396, // e =   694, k =  228
        0xCE1DE40642E3F4B9, // e =   707, k =  232
        0xFB9B7CD9A4A7443C, // e =   720, k =  236
        0x9991A6F3D6BF1766, // e =   734, k =  240
        0xBB764C4CA7A44410, // e =   747, k =  244
        0xE4D5E82392A40515, // e =   760, k =  248
        0x8BAB8EEFB6409C1A, // e =   774, k =  252
        0xAA7EEBFB9DF9DE8E, // e =   787, k =  256
        0xD01FEF10A657842C, // e =   800, k =  260
        0xFE0EFB53D30DD4D8, // e =   813, k =  264
        0x9B10A4E5E9913129, // e =   827, k =  268
        0xBD49D14AA79DBC82, // e =   840, k =  272
        0xE7109BFBA19C0C9D, // e =   853, k =  276
        0x8D07E33455637EB3, // e =   867, k =  280
        0xAC2820D9623BF429, // e =   880, k =  284
        0xD226FC195C6A2F8C, // e =   893, k =  288
        0x80444B5E7AA7CF85, // e =   907, k =  292
        0x9C935E00D4B9D8D2, // e =   920, k =  296
        0xBF21E44003ACDD2D, // e =   933, k =  300
        0xE950DF20247C83FD, // e =   946, k =  304
        0x8E679C2F5E44FF8F, // e =   960, k =  308
        0xADD57A27D29339F6, // e =   973, k =  312
        0xD433179D9C8CB841, // e =   986, k =  316
        0x81842F29F2CCE376, // e =  1000, k =  320
        0x9E19DB92B4E31BA9, // e =  1013, k =  324
        0xC0FE908895CF3B44, // e =  1026, k =  328
    };

    GRISU_ASSERT(e >= -1137);
    GRISU_ASSERT(e <=   960);

    const int k = CeilLog10Pow2(kAlpha - e - 1);
    GRISU_ASSERT(k >= kCachedPowersMinDecExp - (kCachedPowersDecExpStep - 1));
    GRISU_ASSERT(k <= kCachedPowersMaxDecExp);

    const unsigned index = static_cast<unsigned>(k - (kCachedPowersMinDecExp - (kCachedPowersDecExpStep - 1))) / kCachedPowersDecExpStep;
    GRISU_ASSERT(index < kCachedPowersSize);
    static_cast<void>(kCachedPowersSize);

    const int k_cached = kCachedPowersMinDecExp + static_cast<int>(index) * kCachedPowersDecExpStep;
    const int e_cached = FloorLog2Pow10(k_cached) + 1 - DiyFp::SignificandSize;

    const CachedPower cached = {kSignificands[index], e_cached, k_cached};
    GRISU_ASSERT(kAlpha <= cached.e + e + DiyFp::SignificandSize);
    GRISU_ASSERT(kGamma >= cached.e + e + DiyFp::SignificandSize);

    return cached;
}

#if GRISU_ROUND
GRISU_INLINE void Grisu2(uint64_t& decimal_digits, int& decimal_exponent, DiyFp m_minus, DiyFp v, DiyFp m_plus)
#else
GRISU_INLINE void Grisu2(uint64_t& decimal_digits, int& decimal_exponent, DiyFp m_minus, DiyFp m_plus)
#endif
{
    //
    // Step 1:
    // Compute rounding interval
    //

    GRISU_ASSERT(m_plus.e == m_minus.e);
#if GRISU_ROUND
    GRISU_ASSERT(m_plus.e == v.e);
#endif
    //GRISU_ASSERT(m_plus.f - m_minus.f >= 1536); // delta >= 2^(q-p-1) + 2^(q-p-2) = 1024 + 512

    //  --------+-----------------------+-----------------------+--------    (A)
    //          m-                      v                       m+
    //
    //  --------------------+-----------+-----------------------+--------    (B)
    //                      m-          v                       m+
    //
    // First scale v (and m- and m+) such that the exponent is in the range
    // [alpha, gamma].

    const auto cached = GetCachedPowerForBinaryExponent(m_plus.e);

    const DiyFp c_minus_k(cached.f, cached.e); // = c ~= 10^-k

#if GRISU_ROUND
    const DiyFp w       = Multiply(v,       c_minus_k);
#endif
    const DiyFp w_minus = Multiply(m_minus, c_minus_k);
    const DiyFp w_plus  = Multiply(m_plus,  c_minus_k);

    // The exponent of the products is = v.e + c_minus_k.e + q and is in the
    // range [alpha, gamma].
    GRISU_ASSERT(w_plus.e >= kAlpha);
    GRISU_ASSERT(w_plus.e <= kGamma);

    // Note:
    // The result of Multiply() is **NOT** neccessarily normalized.
    // But since m+ and c are normalized, w_plus.f >= 2^(q - 2).
    GRISU_ASSERT(w_plus.f >= (uint64_t{1} << (64 - 2)));

    //GRISU_ASSERT(w_plus.f - w_minus.f >= 768); // delta >= 2^(q-p-2) + 2^(q-p-3) = 512 + 256

    //  ----(---+---)---------------(---+---)---------------(---+---)----
    //          w-                      w                       w+
    //          = c*m-                  = c*v                   = c*m+
    //
    // Multiply rounds its result and c_minus_k is approximated too. w, w- and
    // w+ are now off by a small amount.
    // In fact:
    //
    //      w - v * 10^-k < 1 ulp
    //
    // To account for this inaccuracy, add resp. subtract 1 ulp.
    // Note: ulp(w-) = ulp(w) = ulp(w+).
    //
    //  ----(---+---[---------------(---+---)---------------]---+---)----
    //          w-  L                   w                   H   w+
    //
    // Now any number in [L, H] (bounds included) will round to w when input,
    // regardless of how the input rounding algorithm breaks ties.
    //
    // And DigitGen generates the shortest possible such number in [L, H].
    // Note that this does not mean that Grisu2 always generates the shortest
    // possible number in the interval (m-, m+).
    const DiyFp L(w_minus.f + 1, w_minus.e);
    const DiyFp H(w_plus.f  - 1, w_plus.e );

    //
    // Step 2:
    // Generate digits
    //

    static_assert(DiyFp::SignificandSize == 64, "internal error");
    static_assert(kAlpha >= -50, "internal error");
    static_assert(kGamma <= -32, "internal error");

    // Generates the digits (and the exponent) of a decimal floating-point
    // number V = digits * 10^exponent in the range [L, H].
    // The DiyFp's w, L and H share the same exponent e, which satisfies
    // alpha <= e <= gamma.
    //
    //                                  <---- distance ----->
    //              <---------------------------- delta ---->
    //  ----(---+---[---------------(---+---)---------------]---+---)----
    //              L                   w                   H
    //
    // This routine generates the digits of H from left to right and stops as
    // soon as V is in [L, H].

    GRISU_ASSERT(H.e >= kAlpha);
    GRISU_ASSERT(H.e <= kGamma);
    GRISU_ASSERT(H.e == L.e);
#if GRISU_ROUND
    GRISU_ASSERT(H.e == w.e);
#endif // ^^^ GRISU_ROUND ^^^

#if GRISU_ROUND
    uint64_t distance = Subtract(H, w).f; // (significand of (H - w), implicit exponent is e)
#endif // ^^^ GRISU_ROUND ^^^
    uint64_t delta    = Subtract(H, L).f; // (significand of (H - L), implicit exponent is e)
    uint64_t rest;
    uint64_t ten_kappa;

    // Split H = f * 2^e into two parts p1 and p2 (note: e < 0):
    //
    //      H = f * 2^e
    //           = ((f div 2^-e) * 2^-e + (f mod 2^-e)) * 2^e
    //           = ((p1        ) * 2^-e + (p2        )) * 2^e
    //           = p1 + p2 * 2^e

    const DiyFp one(uint64_t{1} << -H.e, H.e); // one = 2^-e * 2^e

    uint32_t p1 = static_cast<uint32_t>(H.f >> -one.e); // p1 = f div 2^-e (Since -e >= 32, p1 fits into a 32-bit int.)
    uint64_t p2 = H.f & (one.f - 1);                    // p2 = f mod 2^-e

    GRISU_ASSERT(p1 >= 4); // (2^(64-2) - 1) >> 60

    // Generate the digits of the integral part p1 = d[n-1]...d[1]d[0]
    //
    //      10^(k-1) <= p1 < 10^k
    //
    //      p1 = (p1 div 10^(k-1)) * 10^(k-1) + (p1 mod 10^(k-1))
    //         = (d[k-1]         ) * 10^(k-1) + (p1 mod 10^(k-1))
    //
    //      H = p1                                             + p2 * 2^e
    //        = d[k-1] * 10^(k-1) + (p1 mod 10^(k-1))          + p2 * 2^e
    //        = d[k-1] * 10^(k-1) + ((p1 mod 10^(k-1)) * 2^-e + p2) * 2^e
    //        = d[k-1] * 10^(k-1) + (                         rest) * 2^e
    //
    // Now generate the digits d[n] of p1 from left to right (n = k-1,...,0)
    //
    //      p1 = d[k-1]...d[n] * 10^n + d[n-1]...d[0]
    //
    // but stop as soon as
    //
    //      rest * 2^e = (d[n-1]...d[0] * 2^-e + p2) * 2^e <= delta * 2^e

    uint64_t digits = p1;
    int exponent = 0;

    if (p2 > delta)
    {
        // We have
        //
        //      H = d[k-1]...d[1]d[0] + p2 * 2^e
        //        = digits            + p2 * 2^e
        //
        // Now generate the digits of the fractional part p2 * 2^e.
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
        //      H = digits + p2 * 2^e
        //        = digits + 10^-m * (d + r * 2^e)
        //        = (digits * 10^m + d) * 10^-m + 10^-m * r * 2^e
        //
        // and stop as soon as 10^-m * r * 2^e <= delta * 2^e

        // unit = 1
        int m = 0;

        auto remove_digits = [&](uint32_t pow10, int e10)
        {
            GRISU_ASSERT(p2 <= 0xFFFFFFFFFFFFFFFFull / pow10);
            const uint64_t s = pow10 * p2;
            const uint32_t d = static_cast<uint32_t>(s >> -one.e); // d = (pow10 * p2) div 2^-e
            const uint64_t r = s & (one.f - 1);                    // r = (pow10 * p2) mod 2^-e
            GRISU_ASSERT(d < pow10);

            // Check if enough digits have been generated.
            if (r <= pow10 * delta)
            {
                digits = pow10 * digits + d;
                p2 = r;
                m += e10;
#if GRISU_ROUND
                delta *= pow10;
                distance *= pow10;
#endif // ^^^ GRISU_ROUND ^^^
                return true;
            }

            return false;
        };

        for (;;)
        {
            // 64 - 42 = 22 bits
            GRISU_ASSERT(digits <= 9999999999999999ull);

            GRISU_ASSERT(p2 <= 0xFFFFFFFFFFFFFFFFull / 10000);
            const uint64_t s = 10000 * p2;
            const uint32_t d = static_cast<uint32_t>(s >> -one.e); // d = (100 * p2) div 2^-e
            const uint64_t r = s & (one.f - 1);                    // r = (100 * p2) mod 2^-e
            GRISU_ASSERT(d <= 9999);

            if (r <= 10000 * delta)
            {
                // We need to remove UP TO 4 digits.
                // But 4 digits might actually be too much.
                if (remove_digits(10, 1) || remove_digits(100, 2) || remove_digits(1000, 3))
                {
                }
                else
                {
                    digits = 10000 * digits + d;
                    p2 = r;
                    m += 4;
#if GRISU_ROUND
                    delta *= 10000;
                    distance *= 10000;
#endif // ^^^ GRISU_ROUND ^^^
                }

                // V = digits * 10^-m, with L <= V <= H.
                exponent = -m;
#if GRISU_ROUND
                rest = p2;
                ten_kappa = one.f; // one.f == 2^-e
#endif // ^^^ GRISU_ROUND ^^^
                break;
            }

            digits = 10000 * digits + d;
            p2 = r;
            m += 4;
            delta *= 10000;
#if GRISU_ROUND
            distance *= 10000;
#endif // ^^^ GRISU_ROUND ^^^
        }
    }
    else // p2 <= delta
    {
        GRISU_ASSERT((uint64_t{p1} << -one.e) + p2 > delta); // Loop terminates.

        // In this case: p1 contains too many digits.
        //
        // Find the largest 0 <= n < k = length, such that
        //
        //      H = (p1 div 10^n) * 10^n + ((p1 mod 10^n) * 2^-e + p2) * 2^e
        //        = (p1 div 10^n) * 10^n + (                     rest) * 2^e
        //
        // and rest <= delta.
        //
        // Compute rest * 2^e = H mod 10^n = p1 + p2 * 2^e = (p1 * 2^-e + p2) * 2^e
        // and check if enough digits have been generated:
        //
        //      rest * 2^e <= delta * 2^e
        //

        // We can work with 32-bit integers here since p1 fits into a 32-bit
        // integer.
        uint32_t d1 = p1;

        rest = p2;

        // 10^n is now 1 ulp in the decimal representation V. The rounding
        // procedure works with DiyFp's with an implicit exponent of e.
        //
        //      10^n = (10^n * 2^-e) * 2^e = ten_kappa * 2^e
        //
        ten_kappa = one.f; // Start with 2^-e

        int n = 0;
#if 1
        for (;;)
        {
            GRISU_ASSERT(rest <= delta);

            // rn = d[n]...d[0] * 2^-e + p2
            const uint32_t q = d1 / 10;
            const uint32_t r = d1 % 10;
            const uint64_t r_next = ten_kappa * r + rest;

            if (r_next > delta)
            {
                digits = d1;
                exponent = n;
                break;
            }

            d1 = q;
            rest = r_next;
            n += 1;
            ten_kappa *= 10;
        }
#else
        auto remove_digits = [&](uint32_t pow10, int e10)
        {
            const uint32_t q = d1 / pow10;
            const uint32_t r = d1 % pow10;
            const uint64_t r_next = ten_kappa * r + rest;

            if (r_next <= delta)
            {
                digits = q;
                exponent = n + e10;
#if GRISU_ROUND
                rest = r_next;
                ten_kappa *= pow10;
#endif // ^^^ GRISU_ROUND ^^^
                return true;
            }

            return false;
        };

        for (;;)
        {
            GRISU_ASSERT(rest <= delta);

            // rn = d[n]...d[0] * 2^-e + p2
            const uint32_t q = d1 / 10000;
            const uint32_t r = d1 % 10000;
            const uint64_t r_next = ten_kappa * r + rest;

            if (r_next > delta)
            {
                if (remove_digits(1000, 3) || remove_digits(100, 2) || remove_digits(10, 1))
                {
                }
                else
                {
                    digits = d1;
                    exponent = n;
                }
                break;
            }

            d1 = q;
            rest = r_next;
            n += 4;
            ten_kappa *= 10000;
        }
#endif
    }

    //
    // Step 3 (optional):
    // Round towards w.
    //

#if GRISU_ROUND
    GRISU_ASSERT(digits >= 1);
    GRISU_ASSERT(distance <= delta);
    GRISU_ASSERT(rest <= delta);
    GRISU_ASSERT(ten_kappa > 0);

    // By generating the digits of H we got the largest (closest to H) value
    // that is still in the interval [L, H]. In the case where w < B <= H we
    // try to decrement this value.
    //
    //                                  <---- distance ----->
    //              <---------------------------- delta ---->
    //                                         <--- rest --->
    //                       <--- ten_kappa --->
    //  ----(---+---[---------------(---+---)--+------------]---+---)----
    //              L                   w      B            H
    //                                         = digits * 10^kappa
    //
    // ten_kappa represents a unit-in-the-last-place in the decimal
    // representation stored in 'digits'.
    //
    // There are three stopping conditions:
    // (The position of the numbers is measured relative to H.)
    //
    //  1)  B is already <= w
    //          rest >= distance
    //
    //  2)  Decrementing B would yield a number B' < L
    //          rest + ten_kappa > delta
    //
    //  3)  Decrementing B would yield a number B' < w and farther away from
    //      w than the current number B: w - B' > B - w
    //          rest + ten_kappa > distance &&
    //          rest + ten_kappa - distance >= distance - rest

    // The tests are written in this order to avoid overflow in unsigned
    // integer arithmetic.

    while (rest < distance
        && delta - rest >= ten_kappa
        && (rest + ten_kappa <= distance || rest + ten_kappa - distance < distance - rest))
    {
        GRISU_ASSERT(digits % 10 != 0);
        digits--;
        rest += ten_kappa;
    }
#endif // ^^^ GRISU_ROUND ^^^

    //
    // Done.
    //

    decimal_digits = digits;
    decimal_exponent = exponent - cached.k;

    // v = decimal_digits * 10^decimal_exponent
}

} // namespace impl

//==================================================================================================
// ToDecimal
//==================================================================================================

GRISU_INLINE ToDecimalResult<double> ToDecimal(double value)
{
    static_assert(grisu2::impl::DiyFp::SignificandSize >= std::numeric_limits<double>::digits + 3,
        "Grisu2 requires at least three extra bits of precision");

    GRISU_ASSERT(dtoa::IEEE<double>(value).IsFinite());
    GRISU_ASSERT(value > 0);

    const auto boundaries = grisu2::impl::ComputeBoundaries(value);

    ToDecimalResult<double> dec;

#if GRISU_ROUND
    grisu2::impl::Grisu2(dec.digits, dec.exponent, boundaries.m_minus, boundaries.v, boundaries.m_plus);
#else
    grisu2::impl::Grisu2(dec.digits, dec.exponent, boundaries.m_minus, boundaries.m_plus);
#endif

    GRISU_ASSERT(dec.digits <= 99999999999999999ull);
    return dec;
}

GRISU_INLINE ToDecimalResult<float> ToDecimal(float value)
{
    //
    // TODO:
    //
    // Test if a specialized implementation for 'float's using a DiyFp with a
    // 32-bit sigificand would be faster...
    //

    static_assert(grisu2::impl::DiyFp::SignificandSize >= std::numeric_limits<float>::digits + 3,
        "Grisu2 requires at least three extra bits of precision");

    GRISU_ASSERT(dtoa::IEEE<float>(value).IsFinite());
    GRISU_ASSERT(value > 0);

    // If the boundaries of 'value' are always computed for double-precision numbers, regardless of
    // type of 'Float', all single-precision numbers can be recovered using strtod (and strtof).
    // However, the resulting decimal representations are not exactly "short".
    //
    // On the other hand, if the boundaries are computed for single-precision numbers, there is a
    // single number (7.0385307e-26f) which can't be recovered using strtod (instead of strtof).
    // The resulting 'double' when cast to 'float' is off by 1 ulp, i.e. for f = 7.0385307e-26f,
    //      f != (float)strtod(ftoa(f))
    // For all other single-precision numbers, equality holds.

    const auto boundaries = grisu2::impl::ComputeBoundaries(value);

    uint64_t decimal_digits;
    int decimal_exponent;

#if GRISU_ROUND
    grisu2::impl::Grisu2(decimal_digits, decimal_exponent, boundaries.m_minus, boundaries.v, boundaries.m_plus);
#else
    grisu2::impl::Grisu2(decimal_digits, decimal_exponent, boundaries.m_minus, boundaries.m_plus);
#endif

    GRISU_ASSERT(decimal_digits <= 999999999);
    return {static_cast<uint32_t>(decimal_digits), decimal_exponent};
}

//==================================================================================================
// ToChars
//==================================================================================================

// Generates a decimal representation of the floating-point number `value` in 'buffer'.
// Note: The result is _not_ null-terminated.
//
// PRE: The buffer must be large enough (32 bytes is sufficient).
template <typename Float>
GRISU_INLINE char* ToChars(char* buffer, Float value, bool force_trailing_dot_zero = false)
{
    using Fp = dtoa::IEEE<Float>;
    const Fp v(value);

    if (!v.IsFinite())
    {
        if (v.IsNaN())
        {
            std::memcpy(buffer, "NaN", 3);
            return buffer + 3;
        }
        if (v.SignBit())
        {
            *buffer++ = '-';
        }
        std::memcpy(buffer, "Infinity", 8);
        return buffer + 8;
    }

    if (v.SignBit())
    {
        value = v.AbsValue();
        *buffer++ = '-';
    }

    if (v.IsZero())
    {
        *buffer++ = '0';
        if (force_trailing_dot_zero)
        {
            *buffer++ = '.';
            *buffer++ = '0';
        }
        return buffer;
    }

    ToDecimalResult<Float> dec;

#if GRISU_SMALL_INT_OPTIMIZATION
    const bool is_small_int = dtoa::impl::ToSmallInt(dec.digits, value);
    if (is_small_int)
    {
        dec.exponent = 0;

        // Move trailing zeros into the exponent.
        //while (dec.digits % 10 == 0)
        //{
        //    dec.digits /= 10;
        //    dec.exponent++;
        //}
    }
    else
#endif
    {
        dec = grisu2::ToDecimal(value);
    }

    return dtoa::FormatDigits(buffer, dec.digits, dec.exponent, force_trailing_dot_zero);
}

} // namespace grisu2

//char* Dtoa(char* buffer, double value)
//{
//    return grisu2::ToChars(buffer, value);
//}

//char* Ftoa(char* buffer, float value)
//{
//    return grisu2::ToChars(buffer, value);
//}
