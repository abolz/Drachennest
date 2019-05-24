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

#ifndef GRISU_ROUND
#define GRISU_ROUND 0
#endif

namespace grisu2 {

//==================================================================================================
//
//==================================================================================================

template <typename Float>
struct ToDecimalResult;

template <>
struct ToDecimalResult<double>
{
    using SignificandType = uint64_t; // num_digits <= 17
    static constexpr SignificandType MaxDigits = 99999999999999999ull;

    SignificandType digits;
    int exponent;
};

template <>
struct ToDecimalResult<float>
{
    using SignificandType = uint32_t; // num_digits <= 9
    static constexpr SignificandType MaxDigits = 999999999;

    SignificandType digits;
    int exponent;
};

//==================================================================================================
//
//==================================================================================================

namespace impl {

// Returns: floor(x / 2^n)
inline /*__forceinline*/ int SAR(int x, int n)
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
inline /*__forceinline*/ int FloorLog2Pow10(int e)
{
    GRISU_ASSERT(e >= -1233);
    GRISU_ASSERT(e <=  1232);
    return SAR(e * 1741647, 19);
}

// Returns: ceil(log_10(2^e))
inline /*__forceinline*/ int CeilLog10Pow2(int e)
{
    GRISU_ASSERT(e >= -2620);
    GRISU_ASSERT(e <=  2620);
    return SAR(e * 315653 + ((1 << 20) - 1), 20);
}

//==================================================================================================
// Grisu2
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
// Constant data: 159 * 64 bits = 1272 bytes

//
// TODO:
// Clean up comments...
//

struct DiyFp // f * 2^e
{
    static constexpr int SignificandSize = 64; // = q

    uint64_t f = 0;
    int e = 0;

    constexpr DiyFp() = default;
    constexpr DiyFp(uint64_t f_, int e_) : f(f_), e(e_) {}
};

// Returns x * y.
// The result is rounded (ties up). (Only the upper q bits are returned.)
inline /*__forceinline*/ uint64_t MultiplyHighRoundUp(uint64_t x, uint64_t y)
{
    // Computes:
    //  f = round((x.f * y.f) / 2^q)

#if defined(__SIZEOF_INT128__)
    __extension__ using uint128_t = unsigned __int128;

    const uint128_t p = uint128_t{x} * y;

    uint64_t h = static_cast<uint64_t>(p >> 64);
    uint64_t l = static_cast<uint64_t>(p);
    h += l >> 63; // round, ties up: [h, l] += 2^q / 2

    return h;
#elif defined(_MSC_VER) && defined(_M_X64)
    uint64_t h = 0;
    uint64_t l = _umul128(x, y, &h);
    h += l >> 63; // round, ties up: [h, l] += 2^q / 2

    return h;
#else
    const uint32_t x_lo = static_cast<uint32_t>(x);
    const uint32_t x_hi = static_cast<uint32_t>(x >> 32);
    const uint32_t y_lo = static_cast<uint32_t>(y);
    const uint32_t y_hi = static_cast<uint32_t>(y >> 32);

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

    return h;
#endif
}

// Returns the number of leading 0-bits in x, starting at the most significant
// bit position.
// If x is 0, the result is undefined.
inline int CountLeadingZeros64(uint64_t x)
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
// One does not need to store a cached power for each k in this range. For each
// such k it suffices to find a cached power such that the exponent of the
// product lies in [alpha,gamma].
// This implies that the difference of the decimal exponents of adjacent table
// entries must be less than or equal to
//
//      floor( (gamma - alpha) * log_10(2) )
//
// (A smaller distance gamma-alpha would require a larger table.)

struct CachedPower { // c = f * 2^e ~= 10^k
    uint64_t f;
    int e; // binary exponent
    int k; // decimal exponent
};

constexpr int kAlpha = -60;
constexpr int kGamma = -32;
// k_min = -307
// k_max =  324

constexpr int kCachedPowersSize       =   79;
constexpr int kCachedPowersMinDecExp  = -300;
constexpr int kCachedPowersMaxDecExp  =  324;
constexpr int kCachedPowersDecExpStep =    8;

// For a normalized DiyFp w = f * 2^e, this function returns a (normalized)
// cached power-of-ten c = f_c * 2^e_c, such that the exponent of the product
// w * c satisfies
//
//      kAlpha <= e_c + e + q <= kGamma.
//
inline CachedPower GetCachedPowerForBinaryExponent(int e)
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

    GRISU_ASSERT(e >= -1137);
    GRISU_ASSERT(e <=   960);

    const int k = CeilLog10Pow2(kAlpha - e - 1);
    GRISU_ASSERT(k >= kCachedPowersMinDecExp - (kCachedPowersDecExpStep - 1));
    GRISU_ASSERT(k <= kCachedPowersMaxDecExp);

    const unsigned index = static_cast<unsigned>(k - (kCachedPowersMinDecExp - (kCachedPowersDecExpStep - 1))) / kCachedPowersDecExpStep;
    GRISU_ASSERT(index < kCachedPowersSize);

    const int k_cached = kCachedPowersMinDecExp + static_cast<int>(index) * kCachedPowersDecExpStep;
    const int e_cached = FloorLog2Pow10(k_cached) + 1 - 64;

    const CachedPower cached = {kSignificands[index], e_cached, k_cached};
    GRISU_ASSERT(kAlpha <= cached.e + e + 64);
    GRISU_ASSERT(kGamma >= cached.e + e + 64);

    return cached;
}

template <typename Float>
inline void Grisu2(uint64_t& decimal_digits, int& decimal_exponent, Float value)
{
    static_assert(DiyFp::SignificandSize >= std::numeric_limits<double>::digits + 3,
        "Grisu2 requires q >= p + 3");
    static_assert(DiyFp::SignificandSize == 64,
        "This implementation requires q = 64");

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

    using Fp = dtoa::IEEE<Float>;

    GRISU_ASSERT(Fp(value).IsFinite());
    GRISU_ASSERT(value > 0);

    const auto ieee_value = Fp(value);
    const auto ieee_significand = ieee_value.PhysicalSignificand();
    const auto ieee_exponent    = ieee_value.PhysicalExponent();

    int shared_exponent;

#if GRISU_ROUND
    uint64_t m_minus;
    uint64_t v;
    uint64_t m_plus;
    if (ieee_exponent != 0) // normalized floating-point number
    {
        const bool lower_boundary_is_closer = (ieee_significand == 0 && ieee_exponent > 1);

        const auto f2 = ieee_significand | Fp::HiddenBit;
        const auto e2 = static_cast<int>(ieee_exponent) - Fp::ExponentBias;

        const auto fm = 4 * f2 - 2 + (lower_boundary_is_closer ? 1 : 0);
        const auto fv = 4 * f2;
        const auto fp = 4 * f2 + 2;

        const auto shift = DiyFp::SignificandSize - Fp::SignificandSize - 2;

        shared_exponent = e2 - 2 - shift;
        m_minus = uint64_t{fm} << shift;
        v       = uint64_t{fv} << shift;
        m_plus  = uint64_t{fp} << shift;
    }
    else
    {
        const auto f2 = ieee_significand;
        const auto e2 = 1 - Fp::ExponentBias;

        const auto fm = 4 * f2 - 2;
        const auto fv = 4 * f2;
        const auto fp = 4 * f2 + 2;

        const int shift = CountLeadingZeros64(fv);

        shared_exponent = e2 - 2 - shift;
        m_minus = uint64_t{fm} << shift;
        v       = uint64_t{fv} << shift;
        m_plus  = uint64_t{fp} << shift;
    }
#else // ^^^ GRISU_ROUND ^^^
    uint64_t m_minus;
    uint64_t m_plus;
    if (ieee_exponent != 0) // normalized floating-point number
    {
        const bool lower_boundary_is_closer = (ieee_significand == 0 && ieee_exponent > 1);

        const auto f2 = ieee_significand | Fp::HiddenBit;
        const auto e2 = static_cast<int>(ieee_exponent) - Fp::ExponentBias;

        const auto fm = 4 * f2 - 2 + (lower_boundary_is_closer ? 1 : 0);
        const auto fp = 4 * f2 + 2;

        const auto shift = DiyFp::SignificandSize - Fp::SignificandSize - 2;

        shared_exponent = e2 - 2 - shift;
        m_minus = uint64_t{fm} << shift;
        m_plus  = uint64_t{fp} << shift;
    }
    else
    {
        const auto f2 = ieee_significand;
        const auto e2 = 1 - Fp::ExponentBias;

        const auto fm = 4 * f2 - 2;
        const auto fp = 4 * f2 + 2;

        const auto shift = CountLeadingZeros64(fp);

        shared_exponent = e2 - 2 - shift;
        m_minus = uint64_t{fm} << shift;
        m_plus  = uint64_t{fp} << shift;
    }
#endif // ^^^ not GRISU_ROUND ^^^

    //
    // Step 1:
    // Compute rounding interval
    //

    //  --------+-----------------------+-----------------------+--------    (A)
    //          m-                      v                       m+
    //
    //  --------------------+-----------+-----------------------+--------    (B)
    //                      m-          v                       m+
    //
    // First scale v (and m- and m+) such that the exponent is in the range
    // [alpha, gamma].

    const auto cached = GetCachedPowerForBinaryExponent(shared_exponent);

    const uint64_t w_minus = MultiplyHighRoundUp(m_minus, cached.f); // XXX: round down?
#if GRISU_ROUND
    const uint64_t w       = MultiplyHighRoundUp(v,       cached.f); // XXX: compute from w_minus/w_plus?
#endif // ^^^ GRISU_ROUND ^^^
    const uint64_t w_plus  = MultiplyHighRoundUp(m_plus,  cached.f);

    // The exponent of the products is = v.e + cached.e + q and is in the
    // range [alpha, gamma].
    const int e = shared_exponent + cached.e + 64;
    GRISU_ASSERT(e >= kAlpha);
    GRISU_ASSERT(e <= kGamma);

    // Note:
    // The result of Multiply() is **NOT** neccessarily normalized.
    // But since m+ and c are normalized, w+ >= 2^(q - 2).
    GRISU_ASSERT(w_plus >= (uint64_t{1} << (64 - 2)));

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

    const uint64_t L = w_minus + 1;
    const uint64_t H = w_plus  - 1;

    //
    // Step 2:
    // Generate digits
    //

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

#if GRISU_ROUND
    GRISU_ASSERT(H >= w);
    uint64_t distance = H - w; // (significand of (H - w), implicit exponent is H.e)
#endif // ^^^ GRISU_ROUND ^^^
    GRISU_ASSERT(H >= L);
    uint64_t delta    = H - L; // (significand of (H - L), implicit exponent is H.e)
    uint64_t rest;
    uint64_t ten_kappa;

    // Split H = f * 2^e into two parts p1 and p2 (note: e < 0):
    //
    //      H = f * 2^e
    //           = ((f div 2^-e) * 2^-e + (f mod 2^-e)) * 2^e
    //           = ((p1        ) * 2^-e + (p2        )) * 2^e
    //           = p1 + p2 * 2^e

    const DiyFp one(uint64_t{1} << -e, e); // one = 2^-e * 2^e

    uint32_t p1 = static_cast<uint32_t>(H >> -one.e); // p1 = f div 2^-e (Since -e >= 32, p1 fits into a 32-bit int.)
    uint64_t p2 = H & (one.f - 1);                    // p2 = f mod 2^-e

    GRISU_ASSERT(p1 >= 4); // (2^(64-2) - 1) >> 60

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
        // m = 0
        for (;;)
        {
            GRISU_ASSERT(digits <= ToDecimalResult<Float>::MaxDigits);

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
            digits = digits * 10 + d; // digits := digits * 10 + d
            //
            //      H = digits * 10^(-m-1) + 10^(-m-1) * r * 2^e
            //
            p2 = r;
            exponent--; // m += 1
            //
            //      H = digits * 10^-m + 10^-m * p2 * 2^e
            //

            // Keep the units in sync. (unit *= 10)
            delta    *= 10;
#if GRISU_ROUND
            distance *= 10;
#endif // ^^^ GRISU_ROUND ^^^

            // Check if enough digits have been generated.
            //
            //      10^-m * p2 * 2^e <= delta * 2^e
            //              p2 * 2^e <= 10^m * delta * 2^e
            //                    p2 <= 10^m * delta
            if (p2 <= delta)
            {
                // V = digits * 10^-m, with L <= V < H.
                // exponent = -m

#if GRISU_ROUND
                rest = p2;

                // 1 ulp in the decimal representation is now 10^-m.
                // Since delta and distance are now scaled by 10^m, we need to do
                // the same with ulp in order to keep the units in sync.
                //
                //      10^m * 10^-m = 1 = 2^-e * 2^e = ten_m * 2^e
                //
                ten_kappa = one.f; // one.f == 2^-e
#endif // ^^^ GRISU_ROUND ^^^

                break;
            }
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

        // n = 0
        for (;;)
        {
            GRISU_ASSERT(rest <= delta);

            // rn = d[n]...d[0] * 2^-e + p2
            const uint32_t q = p1 / 10;
            const uint32_t r = p1 % 10;
            const uint64_t r_next = ten_kappa * r + rest;

            if (r_next > delta)
            {
                digits = p1;
                break;
            }

            p1 = q;
            exponent += 1; // n += 1
            rest = r_next;
            ten_kappa *= 10;
        }
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

inline ToDecimalResult<double> ToDecimal(double value)
{
    ToDecimalResult<double> dec;

    grisu2::impl::Grisu2(dec.digits, dec.exponent, value);

    GRISU_ASSERT(dec.digits >= 1);
    GRISU_ASSERT(dec.digits <= ToDecimalResult<double>::MaxDigits);

    return dec;
}

inline ToDecimalResult<float> ToDecimal(float value)
{
    //
    // TODO:
    // Test if a specialized implementation for 'float's using a q = 32 would be faster...
    //

    uint64_t decimal_digits;
    int decimal_exponent;

    grisu2::impl::Grisu2(decimal_digits, decimal_exponent, value);

    GRISU_ASSERT(decimal_digits >= 1);
    GRISU_ASSERT(decimal_digits <= ToDecimalResult<float>::MaxDigits);

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
inline char* ToChars(char* buffer, Float value, bool force_trailing_dot_zero = false)
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
