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

#ifndef GRISU_ROUND
#define GRISU_ROUND 1
#endif

namespace grisu2 {

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

template <typename Dest, typename Source>
GRISU_FORCE_INLINE Dest ReinterpretBits(Source source)
{
    static_assert(sizeof(Dest) == sizeof(Source), "size mismatch");

    Dest dest;
    std::memcpy(&dest, &source, sizeof(Source));
    return dest;
}

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

template <int Precision>
struct BitsType;

template <> struct BitsType<24> { using type = uint32_t; };
template <> struct BitsType<53> { using type = uint64_t; };

template <typename Float>
struct IEEE
{
    // NB:
    // Works for double == long double.
    static_assert(std::numeric_limits<Float>::is_iec559 &&
                  ((std::numeric_limits<Float>::digits == 24 && std::numeric_limits<Float>::max_exponent == 128) ||
                   (std::numeric_limits<Float>::digits == 53 && std::numeric_limits<Float>::max_exponent == 1024)),
        "IEEE-754 single- or double-precision implementation required");

    using value_type = Float;
    using bits_type = typename BitsType<std::numeric_limits<Float>::digits>::type;

    static constexpr int       SignificandSize = std::numeric_limits<value_type>::digits; // = p   (includes the hidden bit)
    static constexpr int       ExponentBias    = std::numeric_limits<value_type>::max_exponent - 1 + (SignificandSize - 1);
    static constexpr int       MaxExponent     = std::numeric_limits<value_type>::max_exponent - 1 - (SignificandSize - 1);
    static constexpr int       MinExponent     = std::numeric_limits<value_type>::min_exponent - 1 - (SignificandSize - 1);
    static constexpr bits_type HiddenBit       = bits_type{1} << (SignificandSize - 1);   // = 2^(p-1)
    static constexpr bits_type SignificandMask = HiddenBit - 1;                           // = 2^(p-1) - 1
    static constexpr bits_type ExponentMask    = bits_type{2 * std::numeric_limits<value_type>::max_exponent - 1} << (SignificandSize - 1);
    static constexpr bits_type SignMask        = ~(~bits_type{0} >> 1);

    bits_type bits;

    explicit IEEE(bits_type bits_) : bits(bits_) {}
    explicit IEEE(value_type value) : bits(ReinterpretBits<bits_type>(value)) {}

    bits_type PhysicalSignificand() const {
        return bits & SignificandMask;
    }

    bits_type PhysicalExponent() const {
        return (bits & ExponentMask) >> (SignificandSize - 1);
    }

    bool IsFinite() const {
        return (bits & ExponentMask) != ExponentMask;
    }

    bool IsInf() const {
        return (bits & ExponentMask) == ExponentMask && (bits & SignificandMask) == 0;
    }

    bool IsNaN() const {
        return (bits & ExponentMask) == ExponentMask && (bits & SignificandMask) != 0;
    }

    bool IsZero() const {
        return (bits & ~SignMask) == 0;
    }

    bool SignBit() const {
        return (bits & SignMask) != 0;
    }

    value_type Value() const {
        return ReinterpretBits<value_type>(bits);
    }

    value_type AbsValue() const {
        return ReinterpretBits<value_type>(bits & ~SignMask);
    }
};

// Decomposes `value` into `f * 2^e`.
// The result is not normalized.
// PRE: `value` must be finite and non-negative, i.e. >= +0.0.
template <typename Float>
GRISU_FORCE_INLINE DiyFp DiyFpFromFloat(Float value)
{
    using Fp = IEEE<Float>;

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
    using Fp = IEEE<Float>;

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

#if 0
constexpr int kAlpha = -60;
constexpr int kGamma = -32;
// k_min = -307
// k_max =  324

constexpr int kCachedPowersSize         =   79;
constexpr int kCachedPowersMinDecExp    = -300;
constexpr int kCachedPowersMaxDecExp    =  324;
constexpr int kCachedPowersDecExpStep   =    8;

// Returns (an approximation) 10^(MinDecExp + index * DecExpStep) in the form f * 2^e.
GRISU_FORCE_INLINE CachedPower GetCachedPower(int index)
{
    static constexpr uint64_t kSignificands[] = {
        0xAB70FE17C79AC6CA, // e = -1060, k = -300
        0xFF77B1FCBEBCDC4F, // e = -1034, k = -292
        0xBE5691EF416BD60C, // e = -1007, k = -284
        0x8DD01FAD907FFC3C, // e =  -980, k = -276
        0xD3515C2831559A83, // e =  -954, k = -268
        0x9D71AC8FADA6C9B5, // e =  -927, k = -260
        0xEA9C227723EE8BCB, // e =  -901, k = -252
        0xAECC49914078536D, // e =  -874, k = -244
        0x823C12795DB6CE57, // e =  -847, k = -236
        0xC21094364DFB5637, // e =  -821, k = -228
        0x9096EA6F3848984F, // e =  -794, k = -220
        0xD77485CB25823AC7, // e =  -768, k = -212
        0xA086CFCD97BF97F4, // e =  -741, k = -204
        0xEF340A98172AACE5, // e =  -715, k = -196
        0xB23867FB2A35B28E, // e =  -688, k = -188
        0x84C8D4DFD2C63F3B, // e =  -661, k = -180
        0xC5DD44271AD3CDBA, // e =  -635, k = -172
        0x936B9FCEBB25C996, // e =  -608, k = -164
        0xDBAC6C247D62A584, // e =  -582, k = -156
        0xA3AB66580D5FDAF6, // e =  -555, k = -148
        0xF3E2F893DEC3F126, // e =  -529, k = -140
        0xB5B5ADA8AAFF80B8, // e =  -502, k = -132
        0x87625F056C7C4A8B, // e =  -475, k = -124
        0xC9BCFF6034C13053, // e =  -449, k = -116
        0x964E858C91BA2655, // e =  -422, k = -108
        0xDFF9772470297EBD, // e =  -396, k = -100
        0xA6DFBD9FB8E5B88F, // e =  -369, k =  -92
        0xF8A95FCF88747D94, // e =  -343, k =  -84
        0xB94470938FA89BCF, // e =  -316, k =  -76
        0x8A08F0F8BF0F156B, // e =  -289, k =  -68
        0xCDB02555653131B6, // e =  -263, k =  -60
        0x993FE2C6D07B7FAC, // e =  -236, k =  -52
        0xE45C10C42A2B3B06, // e =  -210, k =  -44
        0xAA242499697392D3, // e =  -183, k =  -36
        0xFD87B5F28300CA0E, // e =  -157, k =  -28
        0xBCE5086492111AEB, // e =  -130, k =  -20
        0x8CBCCC096F5088CC, // e =  -103, k =  -12
        0xD1B71758E219652C, // e =   -77, k =   -4
        0x9C40000000000000, // e =   -50, k =    4
        0xE8D4A51000000000, // e =   -24, k =   12
        0xAD78EBC5AC620000, // e =     3, k =   20
        0x813F3978F8940984, // e =    30, k =   28
        0xC097CE7BC90715B3, // e =    56, k =   36
        0x8F7E32CE7BEA5C70, // e =    83, k =   44
        0xD5D238A4ABE98068, // e =   109, k =   52
        0x9F4F2726179A2245, // e =   136, k =   60
        0xED63A231D4C4FB27, // e =   162, k =   68
        0xB0DE65388CC8ADA8, // e =   189, k =   76
        0x83C7088E1AAB65DB, // e =   216, k =   84
        0xC45D1DF942711D9A, // e =   242, k =   92
        0x924D692CA61BE758, // e =   269, k =  100
        0xDA01EE641A708DEA, // e =   295, k =  108
        0xA26DA3999AEF774A, // e =   322, k =  116
        0xF209787BB47D6B85, // e =   348, k =  124
        0xB454E4A179DD1877, // e =   375, k =  132
        0x865B86925B9BC5C2, // e =   402, k =  140
        0xC83553C5C8965D3D, // e =   428, k =  148
        0x952AB45CFA97A0B3, // e =   455, k =  156
        0xDE469FBD99A05FE3, // e =   481, k =  164
        0xA59BC234DB398C25, // e =   508, k =  172
        0xF6C69A72A3989F5C, // e =   534, k =  180
        0xB7DCBF5354E9BECE, // e =   561, k =  188
        0x88FCF317F22241E2, // e =   588, k =  196
        0xCC20CE9BD35C78A5, // e =   614, k =  204
        0x98165AF37B2153DF, // e =   641, k =  212
        0xE2A0B5DC971F303A, // e =   667, k =  220
        0xA8D9D1535CE3B396, // e =   694, k =  228
        0xFB9B7CD9A4A7443C, // e =   720, k =  236
        0xBB764C4CA7A44410, // e =   747, k =  244
        0x8BAB8EEFB6409C1A, // e =   774, k =  252
        0xD01FEF10A657842C, // e =   800, k =  260
        0x9B10A4E5E9913129, // e =   827, k =  268
        0xE7109BFBA19C0C9D, // e =   853, k =  276
        0xAC2820D9623BF429, // e =   880, k =  284
        0x80444B5E7AA7CF85, // e =   907, k =  292
        0xBF21E44003ACDD2D, // e =   933, k =  300
        0x8E679C2F5E44FF8F, // e =   960, k =  308
        0xD433179D9C8CB841, // e =   986, k =  316
        0x9E19DB92B4E31BA9, // e =  1013, k =  324
    };

    GRISU_ASSERT(index >= 0);
    GRISU_ASSERT(index < kCachedPowersSize);

    const int k = kCachedPowersMinDecExp + index * kCachedPowersDecExpStep;
    const int e = FloorLog2Pow10(k) + 1 - 64;

    return {kSignificands[index], e, k};
}
#endif
#if 1
constexpr int kAlpha = -57;
constexpr int kGamma = -43;
// k_min = -306
// k_max =  325

constexpr int kCachedPowersSize       =  159;
constexpr int kCachedPowersMinDecExp  = -306;
constexpr int kCachedPowersMaxDecExp  =  326;
constexpr int kCachedPowersDecExpStep =    4;

inline CachedPower GetCachedPower(int index)
{
    static constexpr uint64_t kSignificands[] = {
        0xB3C4F1BA87BC8697, // e = -1080, k = -306
        0xDB71E91432B1A24B, // e = -1067, k = -302
        0x85F0468293F0EB4E, // e = -1053, k = -298
        0xA37FCE126597973D, // e = -1040, k = -294
        0xC795830D75038C1E, // e = -1027, k = -290
        0xF3A20279ED56D48A, // e = -1014, k = -286
        0x94B3A202EB1C3F39, // e = -1000, k = -282
        0xB58547448FFFFB2E, // e =  -987, k = -278
        0xDD95317F31C7FA1D, // e =  -974, k = -274
        0x873E4F75E2224E68, // e =  -960, k = -270
        0xA5178FFF668AE0B6, // e =  -947, k = -266
        0xC987434744AC874F, // e =  -934, k = -262
        0xF6019DA07F549B2B, // e =  -921, k = -258
        0x96267C7535B763B5, // e =  -907, k = -254
        0xB749FAED14125D37, // e =  -894, k = -250
        0xDFBDCECE67006AC9, // e =  -881, k = -246
        0x888F99797A5E012D, // e =  -867, k = -242
        0xA6B34AD8C9DFC070, // e =  -854, k = -238
        0xCB7DDCDDA26DA269, // e =  -841, k = -234
        0xF867241C8CC6D4C1, // e =  -828, k = -230
        0x979CF3CA6CEC5B5B, // e =  -814, k = -226
        0xB913179899F68584, // e =  -801, k = -222
        0xE1EBCE4DC7F16DFC, // e =  -788, k = -218
        0x89E42CAAF9491B61, // e =  -774, k = -214
        0xA8530886B54DBDEC, // e =  -761, k = -210
        0xCD795BE870516656, // e =  -748, k = -206
        0xFAD2A4B13D1B5D6C, // e =  -735, k = -202
        0x991711052D8BF3C5, // e =  -721, k = -198
        0xBAE0A846D2195713, // e =  -708, k = -194
        0xE41F3D6A7377EECA, // e =  -695, k = -190
        0x8B3C113C38F9F37F, // e =  -681, k = -186
        0xA9F6D30A038D1DBC, // e =  -668, k = -182
        0xCF79CC9DB955C2CC, // e =  -655, k = -178
        0xFD442E4688BD304B, // e =  -642, k = -174
        0x9A94DD3E8CF578BA, // e =  -628, k = -170
        0xBCB2B812DB11A5DE, // e =  -615, k = -166
        0xE65829B3046B0AFA, // e =  -602, k = -162
        0x8C974F7383725573, // e =  -588, k = -158
        0xAB9EB47C81F5114F, // e =  -575, k = -154
        0xD17F3B51FCA3A7A1, // e =  -562, k = -150
        0xFFBBCFE994E5C620, // e =  -549, k = -146
        0x9C1661A651213E2D, // e =  -535, k = -142
        0xBE89523386091466, // e =  -522, k = -138
        0xE896A0D7E51E1566, // e =  -509, k = -134
        0x8DF5EFABC5979C90, // e =  -495, k = -130
        0xAD4AB7112EB3929E, // e =  -482, k = -126
        0xD389B47879823479, // e =  -469, k = -122
        0x811CCC668829B887, // e =  -455, k = -118
        0x9D9BA7832936EDC1, // e =  -442, k = -114
        0xC06481FB9BCF8D3A, // e =  -429, k = -110
        0xEADAB0ABA3B2DBE5, // e =  -416, k = -106
        0x8F57FA54C2A9EAB7, // e =  -402, k = -102
        0xAEFAE51477A06B04, // e =  -389, k =  -98
        0xD59944A37C0752A2, // e =  -376, k =  -94
        0x825ECC24C8737830, // e =  -362, k =  -90
        0x9F24B832E6B0F436, // e =  -349, k =  -86
        0xC24452DA229B021C, // e =  -336, k =  -82
        0xED246723473E3813, // e =  -323, k =  -78
        0x90BD77F3483BB9BA, // e =  -309, k =  -74
        0xB0AF48EC79ACE837, // e =  -296, k =  -70
        0xD7ADF884AA879177, // e =  -283, k =  -66
        0x83A3EEEEF9153E89, // e =  -269, k =  -62
        0xA0B19D2AB70E6ED6, // e =  -256, k =  -58
        0xC428D05AA4751E4D, // e =  -243, k =  -54
        0xEF73D256A5C0F77D, // e =  -230, k =  -50
        0x9226712162AB070E, // e =  -216, k =  -46
        0xB267ED1940F1C61C, // e =  -203, k =  -42
        0xD9C7DCED53C72256, // e =  -190, k =  -38
        0x84EC3C97DA624AB5, // e =  -176, k =  -34
        0xA2425FF75E14FC32, // e =  -163, k =  -30
        0xC612062576589DDB, // e =  -150, k =  -26
        0xF1C90080BAF72CB1, // e =  -137, k =  -22
        0x9392EE8E921D5D07, // e =  -123, k =  -18
        0xB424DC35095CD80F, // e =  -110, k =  -14
        0xDBE6FECEBDEDD5BF, // e =   -97, k =  -10
        0x8637BD05AF6C69B6, // e =   -83, k =   -6
        0xA3D70A3D70A3D70A, // e =   -70, k =   -2
        0xC800000000000000, // e =   -57, k =    2
        0xF424000000000000, // e =   -44, k =    6
        0x9502F90000000000, // e =   -30, k =   10
        0xB5E620F480000000, // e =   -17, k =   14
        0xDE0B6B3A76400000, // e =    -4, k =   18
        0x878678326EAC9000, // e =    10, k =   22
        0xA56FA5B99019A5C8, // e =    23, k =   26
        0xC9F2C9CD04674EDF, // e =    36, k =   30
        0xF684DF56C3E01BC7, // e =    49, k =   34
        0x96769950B50D88F4, // e =    63, k =   38
        0xB7ABC627050305AE, // e =    76, k =   42
        0xE0352F62A19E306F, // e =    89, k =   46
        0x88D8762BF324CD10, // e =   103, k =   50
        0xA70C3C40A64E6C52, // e =   116, k =   54
        0xCBEA6F8CEB02BB3A, // e =   129, k =   58
        0xF8EBAD2B84E0D58C, // e =   142, k =   62
        0x97EDD871CFDA3A57, // e =   156, k =   66
        0xB975D6B6EE39E437, // e =   169, k =   70
        0xE264589A4DCDAB15, // e =   182, k =   74
        0x8A2DBF142DFCC7AB, // e =   196, k =   78
        0xA8ACD7C0222311BD, // e =   209, k =   82
        0xCDE6FD5E09ABCF27, // e =   222, k =   86
        0xFB5878494ACE3A5F, // e =   235, k =   90
        0x9968BF6ABBE85F20, // e =   249, k =   94
        0xBB445DA9CA61281F, // e =   262, k =   98
        0xE498F455C38B997A, // e =   275, k =  102
        0x8B865B215899F46D, // e =   289, k =  106
        0xAA51823E34A7EEDF, // e =   302, k =  110
        0xCFE87F7CEF46FF17, // e =   315, k =  114
        0xFDCB4FA002162A63, // e =   328, k =  118
        0x9AE757596946075F, // e =   342, k =  122
        0xBD176620A501FC00, // e =   355, k =  126
        0xE6D3102AD96CEC1E, // e =   368, k =  130
        0x8CE2529E2734BB1D, // e =   382, k =  134
        0xABFA45DA0EDBDE69, // e =   395, k =  138
        0xD1EF0244AF2364FF, // e =   408, k =  142
        0x802221226BE55A65, // e =   422, k =  146
        0x9C69A97284B578D8, // e =   435, k =  150
        0xBEEEFB584AFF8604, // e =   448, k =  154
        0xE912B9D1478CEB17, // e =   461, k =  158
        0x8E41ADE9FBEBC27D, // e =   475, k =  162
        0xADA72CCC20054AEA, // e =   488, k =  166
        0xD3FA922F2D1675F2, // e =   501, k =  170
        0x8161AFB94B44F57D, // e =   515, k =  174
        0x9DEFBF01B061ADAB, // e =   528, k =  178
        0xC0CB28A98FCF3C80, // e =   541, k =  182
        0xEB57FF22FC0C795A, // e =   554, k =  186
        0x8FA475791A569D11, // e =   568, k =  190
        0xAF58416654A6BABB, // e =   581, k =  194
        0xD60B3BD56A5586F2, // e =   594, k =  198
        0x82A45B450226B39D, // e =   608, k =  202
        0x9F79A169BD203E41, // e =   621, k =  206
        0xC2ABF989935DDBFE, // e =   634, k =  210
        0xEDA2EE1C7064130C, // e =   647, k =  214
        0x910AB1D4DB9914A0, // e =   661, k =  218
        0xB10D8E1456105DAD, // e =   674, k =  222
        0xD8210BEFD30EFA5A, // e =   687, k =  226
        0x83EA2B892091E44E, // e =   701, k =  230
        0xA1075A24E4421731, // e =   714, k =  234
        0xC491798A08A2AD4F, // e =   727, k =  238
        0xEFF394DCFF8A948F, // e =   740, k =  242
        0x92746B9BE2F8552C, // e =   754, k =  246
        0xB2C71D5BCA9023F8, // e =   767, k =  250
        0xDA3C0F568CC4F3E9, // e =   780, k =  254
        0x8533285C936B35DF, // e =   794, k =  258
        0xA298F2C501F45F43, // e =   807, k =  262
        0xC67BB4597CE2CE49, // e =   820, k =  266
        0xF24A01A73CF2DCD0, // e =   833, k =  270
        0x93E1AB8252F33B46, // e =   847, k =  274
        0xB484F9DC9641E9DB, // e =   860, k =  278
        0xDC5C5301C56B75F7, // e =   873, k =  282
        0x867F59A9D4BED6C0, // e =   887, k =  286
        0xA42E74F3D032F526, // e =   900, k =  290
        0xC86AB5C39FA63441, // e =   913, k =  294
        0xF4A642E14C6262C9, // e =   926, k =  298
        0x95527A5202DF0CCB, // e =   940, k =  302
        0xB6472E511C81471E, // e =   953, k =  306
        0xDE81E40A034BCF50, // e =   966, k =  310
        0x87CEC76F1C830549, // e =   980, k =  314
        0xA5C7EA73224DEFF3, // e =   993, k =  318
        0xCA5E89B18B602368, // e =  1006, k =  322
        0xF70867153AA2DB39, // e =  1019, k =  326
    };

    GRISU_ASSERT(index >= 0);
    GRISU_ASSERT(index < kCachedPowersSize);

    const int k = kCachedPowersMinDecExp + index * kCachedPowersDecExpStep;
    const int e = FloorLog2Pow10(k) + 1 - 64;

    return {kSignificands[index], e, k};
}
#endif

// For a normalized DiyFp w = f * 2^e, this function returns a (normalized)
// cached power-of-ten c = f_c * 2^e_c, such that the exponent of the product
// w * c satisfies (Definition 3.2 from [1])
//
//      alpha <= e_c + e + q <= gamma.
//
GRISU_FORCE_INLINE CachedPower GetCachedPowerForBinaryExponent(int e)
{
    // For double: -1137 <= e <= 960 ==> -307 <= k <= 324 ==>  0 <= index <= 78
    // For single:  -180 <= e <=  96 ==>  -47 <= k <= 36  ==> 32 <= index <= 42
    GRISU_ASSERT(e >= -1137);
    GRISU_ASSERT(e <=   960);

    const int k = CeilLog10Pow2(kAlpha - e - 1);
    GRISU_ASSERT(k >= kCachedPowersMinDecExp - (kCachedPowersDecExpStep - 1));
    GRISU_ASSERT(k <= kCachedPowersMaxDecExp);

    const int index = static_cast<int>( static_cast<unsigned>(k - (kCachedPowersMinDecExp - (kCachedPowersDecExpStep - 1))) / kCachedPowersDecExpStep );
    GRISU_ASSERT(index >= 0);
    GRISU_ASSERT(index < kCachedPowersSize);
    static_cast<void>(kCachedPowersSize);

    const auto cached = GetCachedPower(index);
    GRISU_ASSERT(kAlpha <= cached.e + e + 64);
    GRISU_ASSERT(kGamma >= cached.e + e + 64);

    return cached;
}

#if GRISU_ROUND
// Modifies the generated digits to approach (round towards) w.
//
// Input:
//  * digits of H/10^kappa
//  * distance    = (H - w) * unit
//  * delta       = (H - L) * unit
//  * rest        = (H - digits * 10^kappa) * unit
//  * ten_kappa   = 10^kappa * unit
GRISU_FORCE_INLINE void Grisu2Round(uint64_t& digits, uint64_t distance, uint64_t delta, uint64_t rest, uint64_t ten_kappa)
{
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
}
#endif

// Generates V = decimal_digits * 10^decimal_exponent, such that L <= V <= H.
// L and H must be normalized and share the same exponent -60 <= e <= -32.
#if GRISU_ROUND
GRISU_FORCE_INLINE void Grisu2DigitGen(uint64_t& decimal_digits, int& decimal_exponent, DiyFp L, DiyFp w, DiyFp H)
#else
GRISU_FORCE_INLINE void Grisu2DigitGen(uint64_t& decimal_digits, int& decimal_exponent, DiyFp L, DiyFp H)
#endif
{
    static_assert(DiyFp::SignificandSize == 64, "internal error");
    static_assert(kAlpha >= -60, "internal error");
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
#endif

#if GRISU_ROUND
    uint64_t distance = Subtract(H, w).f; // (significand of (H - w), implicit exponent is e)
#endif
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
        for (;;)
        {
            GRISU_ASSERT(digits <= 9999999999999999ull);
#if 1
// kAlpha >= -57 --->
            uint64_t p2_prev = p2;

            //
            //      H = digits * 10^-m + 10^-m * (d[-m-1] / 10 + d[-m-2] / 10^2 + ...) * 2^e
            //        = digits * 10^-m + 10^-m * (p2                                 ) * 2^e
            //        = digits * 10^-m + 10^-m * (1/10 * (10 * p2)                   ) * 2^e
            //        = digits * 10^-m + 10^-m * (1/10 * ((10*p2 div 2^-e) * 2^-e + (10*p2 mod 2^-e)) * 2^e
            //

            GRISU_ASSERT(p2 <= 0xFFFFFFFFFFFFFFFFull / 100);
            p2 *= 100;
            const uint32_t d = static_cast<uint32_t>(p2 >> -one.e); // d = (100 * p2) div 2^-e
            const uint64_t r = p2 & (one.f - 1);                    // r = (100 * p2) mod 2^-e
            GRISU_ASSERT(d <= 99);
            //
            //      H = digits * 100^-m + 100^-m * (1/100 * (d * 2^-e + r) * 2^e
            //        = digits * 100^-m + 100^-m * (1/100 * (d + r * 2^e))
            //        = (digits * 100 + d) * 10^(-m-2) + 10^(-m-2) * r * 2^e
            //
            digits = 100 * digits + d;
            //
            //      H = digits * 10^(-m-2) + 10^(-m-2) * r * 2^e
            //
            p2 = r;
            m += 2;
            //
            //      H = digits * 10^-m + 10^-m * p2 * 2^e
            //

            // Keep the units in sync. (unit *= 100)
            delta    *= 100;
#if GRISU_ROUND
            distance *= 100;
#endif

            // Check if enough digits have been generated.
            //
            //      10^-m * p2 * 2^e <= delta * 2^e
            //              p2 * 2^e <= 10^m * delta * 2^e
            //                    p2 <= 10^m * delta
            if (p2 <= delta)
            {
                p2_prev *= 10;
//              const uint32_t d_prev = static_cast<uint32_t>(p2_prev >> -one.e);
                const uint64_t r_prev = p2_prev & (one.f - 1);
//              GRISU_ASSERT(d_prev <= 9);

                const uint64_t delta_prev = delta / 10;
                if (r_prev <= delta_prev)
                {
                    digits /= 10;
                    p2 = r_prev;
                    m--;

                    // Keep the units in sync.
                    delta = delta_prev;
#if GRISU_ROUND
                    distance /= 10;
#endif
                }

                // V = digits * 10^-m, with L <= V <= H.
                exponent = -m;

#if GRISU_ROUND
                rest = p2;

                // 1 ulp in the decimal representation is now 10^-m.
                // Since delta and distance are now scaled by 10^m, we need to do
                // the same with ulp in order to keep the units in sync.
                //
                //      10^m * 10^-m = 1 = 2^-e * 2^e = ten_m * 2^e
                //
                ten_kappa = one.f; // one.f == 2^-e
#endif

                break;
            }
// <--- kAlpha >= -57
#else
            //
            //      H = digits * 10^-m + 10^-m * (d[-m-1] / 10 + d[-m-2] / 10^2 + ...) * 2^e
            //        = digits * 10^-m + 10^-m * (p2                                 ) * 2^e
            //        = digits * 10^-m + 10^-m * (1/10 * (10 * p2)                   ) * 2^e
            //        = digits * 10^-m + 10^-m * (1/10 * ((10*p2 div 2^-e) * 2^-e + (10*p2 mod 2^-e)) * 2^e
            //
            GRISU_ASSERT(p2 <= 0xFFFFFFFFFFFFFFFFull / 10);
            p2 *= 10;
            const uint32_t d = static_cast<uint32_t>(p2 >> -one.e); // d = (10 * p2) div 2^-e
            const uint64_t r = p2 & (one.f - 1);                    // r = (10 * p2) mod 2^-e
            GRISU_ASSERT(d <= 9);
            //
            //      H = digits * 10^-m + 10^-m * (1/10 * (d * 2^-e + r) * 2^e
            //        = digits * 10^-m + 10^-m * (1/10 * (d + r * 2^e))
            //        = (digits * 10 + d) * 10^(-m-1) + 10^(-m-1) * r * 2^e
            //
            digits = 10 * digits + d;
            //
            //      H = digits * 10^(-m-1) + 10^(-m-1) * r * 2^e
            //
            p2 = r;
            m++;
            //
            //      H = digits * 10^-m + 10^-m * p2 * 2^e
            //

            // Keep the units in sync. (unit *= 10)
            delta    *= 10;
#if GRISU_ROUND
            distance *= 10;
#endif

            // Check if enough digits have been generated.
            //
            //      10^-m * p2 * 2^e <= delta * 2^e
            //              p2 * 2^e <= 10^m * delta * 2^e
            //                    p2 <= 10^m * delta
            if (p2 <= delta)
            {
                // V = digits * 10^-m, with L <= V <= H.
                exponent = -m;

#if GRISU_ROUND
                rest = p2;

                // 1 ulp in the decimal representation is now 10^-m.
                // Since delta and distance are now scaled by 10^m, we need to do
                // the same with ulp in order to keep the units in sync.
                //
                //      10^m * 10^-m = 1 = 2^-e * 2^e = ten_m * 2^e
                //
                ten_kappa = one.f; // one.f == 2^-e
#endif

                break;
            }
#endif
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

        rest = p2;

        // 10^n is now 1 ulp in the decimal representation V. The rounding
        // procedure works with DiyFp's with an implicit exponent of e.
        //
        //      10^n = (10^n * 2^-e) * 2^e = ten_kappa * 2^e
        //
        ten_kappa = one.f; // Start with 2^-e

        // We can work with 32-bit integers here since p1 fits into a 32-bit
        // integer.
        uint32_t dn = p1;
        for (int n = 0; /**/; ++n)
        {
            GRISU_ASSERT(rest <= delta);

            // rn = d[n]...d[0] * 2^-e + p2
            const uint32_t qn = dn / 10;
            const uint64_t rn = ten_kappa * (dn - 10 * qn) + rest;

            if (rn > delta)
            {
                digits = dn;
                exponent = n;
                break;
            }

            dn = qn;
            rest = rn;
            ten_kappa *= 10;
        }
    }

    // Now w = digits * 10^exponent.
    // Optionally round the digits towards w.

#if GRISU_ROUND
    Grisu2Round(digits, distance, delta, rest, ten_kappa);
#endif

    decimal_digits = digits;
    decimal_exponent = exponent;
}

#if GRISU_ROUND
/*GRISU_FORCE_INLINE*/GRISU_INLINE void Grisu2(uint64_t& decimal_digits, int& decimal_exponent, DiyFp m_minus, DiyFp v, DiyFp m_plus)
#else
/*GRISU_FORCE_INLINE*/GRISU_INLINE void Grisu2(uint64_t& decimal_digits, int& decimal_exponent, DiyFp m_minus, DiyFp m_plus)
#endif
{
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

#if GRISU_ROUND
    Grisu2DigitGen(decimal_digits, decimal_exponent, L, w, H);
#else
    Grisu2DigitGen(decimal_digits, decimal_exponent, L, H);
#endif
    // w = digits * 10^exponent

    // v = w * 10^k
    decimal_exponent += -cached.k; // cached.k = -k
    // v = digits * 10^exponent
}

} // namespace impl

//==================================================================================================
// ToDecimal
//==================================================================================================

struct F64ToDecimalResult {
    uint64_t digits;
    int exponent;
};

GRISU_INLINE F64ToDecimalResult ToDecimal(double value)
{
    static_assert(grisu2::impl::DiyFp::SignificandSize >= std::numeric_limits<double>::digits + 3,
        "Grisu2 requires at least three extra bits of precision");

    GRISU_ASSERT(grisu2::impl::IEEE<double>(value).IsFinite());
    GRISU_ASSERT(value > 0);

    const auto boundaries = grisu2::impl::ComputeBoundaries(value);

    F64ToDecimalResult dec;

#if GRISU_ROUND
    grisu2::impl::Grisu2(dec.digits, dec.exponent, boundaries.m_minus, boundaries.v, boundaries.m_plus);
#else
    grisu2::impl::Grisu2(dec.digits, dec.exponent, boundaries.m_minus, boundaries.m_plus);
#endif

    GRISU_ASSERT(dec.digits <= 99999999999999999ull);
    return dec;
}

struct F32ToDecimalResult {
    uint32_t digits;
    int exponent;
};

GRISU_INLINE F32ToDecimalResult ToDecimal(float value)
{
    //
    // TODO:
    //
    // Test if a specialized implementation for 'float's using a DiyFp with a
    // 32-bit sigificand would be faster...
    //

    static_assert(grisu2::impl::DiyFp::SignificandSize >= std::numeric_limits<float>::digits + 3,
        "Grisu2 requires at least three extra bits of precision");

    GRISU_ASSERT(grisu2::impl::IEEE<float>(value).IsFinite());
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
    using Fp = grisu2::impl::IEEE<Float>;
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

    const auto dec = grisu2::ToDecimal(value);
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
