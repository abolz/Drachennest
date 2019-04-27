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

#ifndef GRISU_INLINE
#define GRISU_INLINE inline
#endif

namespace grisu3 {

//==================================================================================================
// Grisu3
//
// Implements the Grisu3 algorithm for (IEEE) binary to decimal floating-point conversion.
//
// References:
//
// [1]  Loitsch, "Printing Floating-Point Numbers Quickly and Accurately with Integers",
//      Proceedings of the ACM SIGPLAN 2010 Conference on Programming Language Design and Implementation, PLDI 2010
// [2]  Burger, Dybvig, "Printing Floating-Point Numbers Quickly and Accurately",
//      Proceedings of the ACM SIGPLAN 1996 Conference on Programming Language Design and Implementation, PLDI 1996
//==================================================================================================
// Constant data: 632 + 200 = 832 bytes

namespace impl {

template <typename Dest, typename Source>
GRISU_INLINE Dest ReinterpretBits(Source source)
{
    static_assert(sizeof(Dest) == sizeof(Source), "size mismatch");

    Dest dest;
    std::memcpy(&dest, &source, sizeof(Source));
    return dest;
}

GRISU_INLINE char* Utoa_2Digits(char* buf, uint32_t digits)
{
    static constexpr char kDigits100[200] = {
        '0','0','0','1','0','2','0','3','0','4','0','5','0','6','0','7','0','8','0','9',
        '1','0','1','1','1','2','1','3','1','4','1','5','1','6','1','7','1','8','1','9',
        '2','0','2','1','2','2','2','3','2','4','2','5','2','6','2','7','2','8','2','9',
        '3','0','3','1','3','2','3','3','3','4','3','5','3','6','3','7','3','8','3','9',
        '4','0','4','1','4','2','4','3','4','4','4','5','4','6','4','7','4','8','4','9',
        '5','0','5','1','5','2','5','3','5','4','5','5','5','6','5','7','5','8','5','9',
        '6','0','6','1','6','2','6','3','6','4','6','5','6','6','6','7','6','8','6','9',
        '7','0','7','1','7','2','7','3','7','4','7','5','7','6','7','7','7','8','7','9',
        '8','0','8','1','8','2','8','3','8','4','8','5','8','6','8','7','8','8','8','9',
        '9','0','9','1','9','2','9','3','9','4','9','5','9','6','9','7','9','8','9','9',
    };

    GRISU_ASSERT(digits < 100);
    std::memcpy(buf, &kDigits100[2 * digits], 2 * sizeof(char));
    return buf + 2;
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
GRISU_INLINE DiyFp Subtract(DiyFp x, DiyFp y)
{
    GRISU_ASSERT(x.e == y.e);
    GRISU_ASSERT(x.f >= y.f);

    return DiyFp(x.f - y.f, x.e);
}

// Returns x * y.
// The result is rounded (ties up). (Only the upper q bits are returned.)
GRISU_INLINE DiyFp Multiply(DiyFp x, DiyFp y)
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
    const uint32_t xLo = static_cast<uint32_t>(x.f);
    const uint32_t xHi = static_cast<uint32_t>(x.f >> 32);
    const uint32_t yLo = static_cast<uint32_t>(y.f);
    const uint32_t yHi = static_cast<uint32_t>(y.f >> 32);

    const uint64_t b00 = uint64_t{xLo} * yLo;
    const uint64_t b01 = uint64_t{xLo} * yHi;
    const uint64_t b10 = uint64_t{xHi} * yLo;
    const uint64_t b11 = uint64_t{xHi} * yHi;

    const uint32_t b00Hi = static_cast<uint32_t>(b00 >> 32);

    const uint64_t mid1 = b10 + b00Hi;
    const uint32_t mid1Lo = static_cast<uint32_t>(mid1);
    const uint32_t mid1Hi = static_cast<uint32_t>(mid1 >> 32);

    const uint64_t mid2 = b01 + mid1Lo;
    const uint32_t mid2Lo = static_cast<uint32_t>(mid2);
    const uint32_t mid2Hi = static_cast<uint32_t>(mid2 >> 32);

    // NB: mid2Lo has the upper 32 bits of the low part of the product.
    const uint32_t r = mid2Lo >> 31;
    const uint64_t h = b11 + mid1Hi + mid2Hi + r;

    return DiyFp(h, x.e + y.e + 64);
#endif
}

// Returns the number of leading 0-bits in x, starting at the most significant bit position.
// If x is 0, the result is undefined.
GRISU_INLINE int CountLeadingZeros64(uint64_t x)
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
GRISU_INLINE DiyFp Normalize(DiyFp x)
{
    static_assert(DiyFp::SignificandSize == 64, "internal error");

    // For double: lz >= 64-53 = 11
    // For single: lz >= 64-24 = 40
    const int lz = CountLeadingZeros64(x.f);
    return DiyFp(x.f << lz, x.e - lz);
}

// Normalize x such that the result has the exponent E.
// PRE: e >= x.e and the upper e - x.e bits of x.f must be zero.
GRISU_INLINE DiyFp NormalizeTo(DiyFp x, int e)
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
GRISU_INLINE DiyFp DiyFpFromFloat(Float value)
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
    DiyFp v;
    DiyFp m_minus;
    DiyFp m_plus;
};

// Compute the (normalized) DiyFp representing the input number 'value' and its
// boundaries.
// PRE: 'value' must be finite and positive
template <typename Float>
GRISU_INLINE Boundaries ComputeBoundaries(Float value)
{
    using Fp = IEEE<Float>;

    GRISU_ASSERT(Fp(value).IsFinite());
    GRISU_ASSERT(value > 0);

    const auto v = DiyFpFromFloat(value);

    // Compute the boundaries of v.
    const bool lower_boundary_is_closer = (v.f == Fp::HiddenBit && v.e > Fp::MinExponent);
    const auto m_minus = DiyFp(4*v.f - 2 + lower_boundary_is_closer, v.e - 2);
    const auto m_plus = DiyFp(4*v.f + 2, v.e - 2);

    // Determine the normalized w = v.
    const auto w = Normalize(v);

    // Determine the normalized w+ = m+.
    // Since e_(w+) == e_(w), one can use NormalizeTo instead of Normalize.
    const auto w_plus = NormalizeTo(m_plus, w.e);

    // Determine w- = m- such that e_(w-) = e_(w+).
    const auto w_minus = NormalizeTo(m_minus, w_plus.e);

    return {w, w_minus, w_plus};
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
// the digit generation procedure. Using (alpha,gamma)=(-60,-32) works out well
// in practice:
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

constexpr int kAlpha = -60;
constexpr int kGamma = -32;

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
GRISU_INLINE int SAR(int x, int n)
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
GRISU_INLINE int FloorLog2Pow10(int e)
{
    GRISU_ASSERT(e >= -1233);
    GRISU_ASSERT(e <=  1232);
    return SAR(e * 1741647, 19);
}

// Returns: ceil(log_10(2^e))
GRISU_INLINE int CeilLog10Pow2(int e)
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

constexpr int kCachedPowersSize         =   79;
constexpr int kCachedPowersMinDecExp    = -300;
constexpr int kCachedPowersMaxDecExp    =  324;
constexpr int kCachedPowersDecExpStep   =    8;

// Returns (an approximation) 10^(MinDecExp + index * DecExpStep) in the form f * 2^e.
GRISU_INLINE CachedPower GetCachedPower(int index)
{
    // Let e = floor(log_2 10^k) + 1 - 64.
    // Negative powers of 10 are stored as: f = round_up(2^-e / 10^-k).
    // Positive powers of 10 are stored as: f = round_up(10^k / 2^e).
    static constexpr uint64_t kSignificands[/*632 bytes*/] = {
        0xAB70FE17C79AC6CA, // e = -1060, k = -300 >>> double-precision
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
        0xE45C10C42A2B3B06, // e =  -210, k =  -44 >>> single-precision
        0xAA242499697392D3, // e =  -183, k =  -36
        0xFD87B5F28300CA0E, // e =  -157, k =  -28
        0xBCE5086492111AEB, // e =  -130, k =  -20
        0x8CBCCC096F5088CC, // e =  -103, k =  -12
        0xD1B71758E219652C, // e =   -77, k =   -4
        0x9C40000000000000, // e =   -50, k =    4
        0xE8D4A51000000000, // e =   -24, k =   12
        0xAD78EBC5AC620000, // e =     3, k =   20
        0x813F3978F8940984, // e =    30, k =   28
        0xC097CE7BC90715B3, // e =    56, k =   36 <<< single-precision
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
        0x9E19DB92B4E31BA9, // e =  1013, k =  324 <<< double-precision
    };

    GRISU_ASSERT(index >= 0);
    GRISU_ASSERT(index < kCachedPowersSize);

    const int k = kCachedPowersMinDecExp + index * kCachedPowersDecExpStep;
    const int e = FloorLog2Pow10(k) + 1 - 64;

    return {kSignificands[index], e, k};
}

// For a normalized DiyFp w = f * 2^e, this function returns a (normalized)
// cached power-of-ten c = f_c * 2^e_c, such that the exponent of the product
// w * c satisfies (Definition 3.2 from [1])
//
//      alpha <= e_c + e + q <= gamma.
//
GRISU_INLINE CachedPower GetCachedPowerForBinaryExponent(int e)
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

    // NB:
    // Actually this function returns c, such that -60 <= e_c + e + 64 <= -34.
    GRISU_ASSERT(-60 <= cached.e + e + 64);
    GRISU_ASSERT(-34 >= cached.e + e + 64);

    return cached;
}

GRISU_INLINE char* GenerateIntegralDigits(char* buf, uint32_t n)
{
//  GRISU_ASSERT(n <= 798336123);
    GRISU_ASSERT(n <= 999999999);

    uint32_t q;
    if (n >= 100000000)
    {
//L_9_digits:
        q = n / 10000000;
        n = n % 10000000;
        buf = Utoa_2Digits(buf, q);
L_7_digits:
        q = n / 100000;
        n = n % 100000;
        buf = Utoa_2Digits(buf, q);
L_5_digits:
        q = n / 1000;
        n = n % 1000;
        buf = Utoa_2Digits(buf, q);
L_3_digits:
        q = n / 10;
        n = n % 10;
        buf = Utoa_2Digits(buf, q);
L_1_digit:
        buf[0] = static_cast<char>('0' + n);
        buf++;
        return buf;
    }

    if (n >= 10000000)
    {
//L_8_digits:
        q = n / 1000000;
        n = n % 1000000;
        buf = Utoa_2Digits(buf, q);
L_6_digits:
        q = n / 10000;
        n = n % 10000;
        buf = Utoa_2Digits(buf, q);
L_4_digits:
        q = n / 100;
        n = n % 100;
        buf = Utoa_2Digits(buf, q);
L_2_digits:
        buf = Utoa_2Digits(buf, n);
        return buf;
    }

#if 1
    if (n >= 100000) { if (n >= 1000000) goto L_7_digits; else goto L_6_digits; }
    if (n >=   1000) { if (n >=   10000) goto L_5_digits; else goto L_4_digits; }
    if (n >=     10) { if (n >=     100) goto L_3_digits; else goto L_2_digits; }
    goto L_1_digit;
#else
    if (n >= 1000000) goto L_7_digits;
    if (n >=  100000) goto L_6_digits;
    if (n >=   10000) goto L_5_digits;
    if (n >=    1000) goto L_4_digits;
    if (n >=     100) goto L_3_digits;
    if (n >=      10) goto L_2_digits;
    goto L_1_digit;
#endif
}

// Modifies the generated digits in the buffer to approach (round towards) w.
//
// Input:
//  * digits of H/10^kappa in [digits, digits + num_digits)
//  * distance    = (H - w) * unit
//  * delta       = (H - L) * unit
//  * rest        = (H - digits * 10^kappa) * unit
//  * ten_kappa   = 10^kappa * unit
GRISU_INLINE bool Grisu3RoundWeed(char* digits, int num_digits, uint64_t distance, uint64_t delta, uint64_t rest, uint64_t ten_kappa, uint64_t unit)
{
    GRISU_ASSERT(num_digits >= 1);
    GRISU_ASSERT(distance <= delta);
    GRISU_ASSERT(rest <= delta);
    GRISU_ASSERT(ten_kappa > 0);
    GRISU_ASSERT(unit > 0);
    GRISU_ASSERT(distance >= unit);
    GRISU_ASSERT(distance <= UINT64_MAX - unit);

    const uint64_t distance_plus  = distance - unit;
    const uint64_t distance_minus = distance + unit;

    // By generating the digits of H we got the largest (closest to H) buffer
    // that is still in the interval [L, H]. In the case where M+ <= B < H we
    // try to decrement the buffer.
    //
    //                                  <------------ distance ----->
    //      <-------------------------------------------- delta ---->
    //                                         <----------- rest --->
    //                       <--- ten_kappa --->
    //  ----[---+---[---------------(---+---)--+------------]---+---)----
    //      L   w-  L+              M-  w   M+ B            H-  w+  H
    //                                         = digits * 10^kappa
    //
    // ten_kappa represents a unit-in-the-last-place in the decimal
    // representation stored in the buffer.
    //
    // There are three stopping conditions:
    // (The position of the numbers is measured relative to H.)
    //
    //  1)  B is already < M+
    //          rest > distance
    //
    //  2)  Decrementing B would yield a number B' <= L
    //          rest + ten_kappa >= delta
    //
    //  3)  Decrementing B would yield a number B' <= M+ and farther away from
    //      M+ than the current number B: M+ - B' > B - M+
    //          rest + ten_kappa > distance &&
    //          rest + ten_kappa - distance >= distance - rest

    // The tests are written in this order to avoid overflow in unsigned
    // integer arithmetic.

    int digit = digits[num_digits - 1] - '0';

    while (rest <= distance_plus
        && delta - rest > ten_kappa
        && (rest + ten_kappa <= distance_plus || rest + ten_kappa - distance_plus < distance_plus - rest))
    {
        GRISU_ASSERT(digit != 0);
        digit--;
        rest += ten_kappa;
    }

    digits[num_digits - 1] = static_cast<char>('0' + digit);

    // Now try to approach M- and check if we might generate a number B' which is closer
    // to M- as B is to M+.
    //
    //  --------(---+-------------------+------+----------------)--------
    //          M-  B'                  w      B                M+
    //
    // If so, there are two representations but Grisu3 is too imprecise to determine which
    // one is actually closer. (Note that w might lie anywhere in the interval (M-, M+).)

    if (rest < distance_minus
        && delta - rest >= ten_kappa
        && (rest + ten_kappa <= distance_minus || rest + ten_kappa - distance_minus < distance_minus - rest))
    {
        return false;
    }

    // Now test if B lies in the safe interval [L+, H-].
    // If it doesn't, Grisu3 is too imprecise and we need to fall back to a more accurate algorithm.
    //
    //      <-------------------------------------------- delta ---->
    //      <--->                              <----------- rest --->
    //       ulp
    //  ----[---+---[--------------------------+------------]---+---)----
    //      L   w-  L+                         B            H-  w+  H

    GRISU_ASSERT(delta >= 4 * unit);

    return 2 * unit <= rest && rest <= delta - 4 * unit;
}

// Generates V = digits * 10^exponent, such that L <= V <= H.
// L and H must be normalized and share the same exponent -60 <= e <= -32.
GRISU_INLINE bool Grisu3DigitGen(char* digits, int& num_digits, int& exponent, DiyFp L, DiyFp w, DiyFp H)
{
    static_assert(DiyFp::SignificandSize == 64, "internal error");
    static_assert(kAlpha >= -60, "internal error");
    static_assert(kGamma <= -32, "internal error");

    // Generates the digits (and the exponent) of a decimal floating-point
    // number V = digits * 10^exponent in the range [L, H).
    // The DiyFp's w, L and H share the same exponent e, which satisfies
    // alpha <= e <= gamma.
    //
    //                                  <------------ distance ----->
    //      <-------------------------------------------- delta ---->
    //  ----(---+---[---------------(---+---)---------------]---+---)----
    //      L   w-  L+              M-  w   M+              H-  w+  H
    //
    // This routine generates the digits of H from left to right and stops as
    // soon as V is in [L, H).

    GRISU_ASSERT(H.e >= kAlpha);
    GRISU_ASSERT(H.e <= kGamma);
    GRISU_ASSERT(H.e == L.e);
    GRISU_ASSERT(H.e == w.e);

    uint64_t distance = Subtract(H, w).f; // (significand of (H - w), implicit exponent is e)
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

    GRISU_ASSERT(p1 >= 4);            // (2^(64-2) - 1) >> 60
    GRISU_ASSERT(p1 <= 798336123);    // Not trivial. Depends on the implementation of GetCachedPowerForBinaryExponent!

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
    //      rest * 2^e = (d[n-1]...d[0] * 2^-e + p2) * 2^e < delta * 2^e

    // The common case is that all the digits of p1 are needed.
    // Optimize for this case and correct later if required.
    num_digits = static_cast<int>(GenerateIntegralDigits(digits, p1) - digits);

    uint64_t unit = 1;
    if (p2 >= delta)
    {
        // The digits of the integral part have been generated (and all of them
        // are significand):
        //
        //      H = d[k-1]...d[1]d[0] + p2 * 2^e
        //        = digits            + p2 * 2^e
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
        //      H = digits + p2 * 2^e
        //        = digits + 10^-m * (d + r * 2^e)
        //        = (digits * 10^m + d) * 10^-m + 10^-m * r * 2^e
        //
        // and stop as soon as 10^-m * r * 2^e < delta * 2^e

        // unit = 1
        int m = 0;
        for (;;)
        {
            // !!! GRISU_ASSERT(num_digits < max_digits10) !!!
            GRISU_ASSERT(num_digits < 17);

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
            digits[num_digits++] = static_cast<char>('0' + d); // digits := digits * 10 + d
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
            distance *= 10;
            unit     *= 10;

            // Check if enough digits have been generated.
            //
            //      10^-m * p2 * 2^e < delta * 2^e
            //              p2 * 2^e < 10^m * delta * 2^e
            //                    p2 < 10^m * delta
            if (p2 < delta)
            {
                // V = digits * 10^-m, with L <= V < H.
                exponent = -m;

                rest = p2;

                // 1 ulp in the decimal representation is now 10^-m.
                // Since delta and distance are now scaled by 10^m, we need to do
                // the same with ulp in order to keep the units in sync.
                //
                //      10^m * 10^-m = 1 = 2^-e * 2^e = ten_m * 2^e
                //
                ten_kappa = one.f; // one.f == 2^-e

                break;
            }
        }
    }
    else // p2 < delta
    {
        GRISU_ASSERT((uint64_t{p1} << -one.e) + p2 >= delta); // Loop terminates.

        // In this case: Too many digits of p1 might have been generated.
        //
        // Find the largest 0 <= n < k = length, such that
        //
        //      H = (p1 div 10^n) * 10^n + ((p1 mod 10^n) * 2^-e + p2) * 2^e
        //        = (p1 div 10^n) * 10^n + (                     rest) * 2^e
        //
        // and rest < delta.
        //
        // Compute rest * 2^e = H mod 10^n = p1 + p2 * 2^e = (p1 * 2^-e + p2) * 2^e
        // and check if enough digits have been generated:
        //
        //      rest * 2^e < delta * 2^e
        //

        const int k = num_digits;
        GRISU_ASSERT(k >= 0);
        GRISU_ASSERT(k <= 9);

        rest = p2;

        // 10^n is now 1 ulp in the decimal representation V. The rounding
        // procedure works with DiyFp's with an implicit exponent of e.
        //
        //      10^n = (10^n * 2^-e) * 2^e = ten_kappa * 2^e
        //
        ten_kappa = one.f; // Start with 2^-e

        for (int n = 0; /**/; ++n)
        {
            GRISU_ASSERT(n <= k - 1);
            GRISU_ASSERT(rest < delta);

            // rn = d[n]...d[0] * 2^-e + p2
            const uint32_t dn = static_cast<uint32_t>(digits[k - 1 - n] - '0');
            const uint64_t rn = dn * ten_kappa + rest;

            if (rn >= delta)
            {
                num_digits = k - n;
                exponent = n;
                break;
            }

            rest = rn;
            ten_kappa *= 10;
        }
    }

    return Grisu3RoundWeed(digits, num_digits, distance, delta, rest, ten_kappa, unit);
}

// v = digits * 10^exponent
// length is the length of the buffer (number of decimal digits)
// The buffer must be large enough, i.e. >= max_digits10.
GRISU_INLINE bool Grisu3(char* digits, int& num_digits, int& exponent, DiyFp m_minus, DiyFp v, DiyFp m_plus)
{
    // For single-precision:  99.172% optimal
    // For double-precision: ~99.45% optimal (uniformly distributed exponent/significands)

    GRISU_ASSERT(m_plus.e == m_minus.e);
    GRISU_ASSERT(m_plus.e == v.e);

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

    const DiyFp w       = Multiply(v,       c_minus_k);
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
    // To account for this inaccuracy, subtract resp. add 1 ulp.
    // Note: ulp(w-) = ulp(w) = ulp(w+).
    //
    //  ----(---+---[---------------(---+---)---------------]---+---)----
    //      L   w-  L+              M-  w   M+              H-  w+  H
    //
    // Now any number in [L+, H-] (bounds included) will round to w when input,
    // regardless of how the input rounding algorithm breaks ties.
    //
    // But since this interval is too narrow, there might be a shorter representation
    // in the interval (w-, w+).
    //
    // Grisu3 now generates the shortest possible number in [L, H), which includes
    // the interval (w-, w+) and also the safe interval [L+, H-].
    //
    // Grisu3 also determines if the result is in the safe interval [L+, H-] and whether
    // the result is unique.
    // If Grisu3 fails, we need to fall back to a more accurate algorithm.
    const DiyFp L(w_minus.f - 1, w_minus.e);
    const DiyFp H(w_plus.f  + 1, w_plus.e );

    const bool result = Grisu3DigitGen(digits, num_digits, exponent, L, w, H);
    // w = digits * 10^exponent

    // v = w * 10^k
    exponent += -cached.k; // cached.k = -k
    // v = digits * 10^exponent

    return result;
}

} // namespace impl

//==================================================================================================
// Dragon4
//
// Implements the Dragon4 algorithm for (IEEE) binary to decimal floating-point conversion.
//
// References:
//
// [1]  Burger, Dybvig, "Printing Floating-Point Numbers Quickly and Accurately",
//      Proceedings of the ACM SIGPLAN 1996 Conference on Programming Language Design and Implementation, PLDI 1996
// [2]  Steele, White, "How to Print FloatingPoint Numbers Accurately",
//      Proceedings of the ACM SIGPLAN 1990 conference on Programming language design and implementation, PLDI 1990
//==================================================================================================
// Constant data: 56 bytes

namespace impl {

GRISU_INLINE int Min(int x, int y) { return y < x ? y : x; }
GRISU_INLINE int Max(int x, int y) { return y < x ? x : y; }

//template <int MaxBits>
struct DiyInt
{
    static constexpr int MaxBits = 1130;
    static constexpr int Capacity = (MaxBits + (32 - 1)) / 32;

    uint32_t bigits[Capacity]; // Significand stored in little-endian form.
    int      size = 0;

    DiyInt() = default;
    DiyInt(DiyInt const&) = delete;             // (not needed here)
    DiyInt& operator=(DiyInt const&) = delete;  // (not needed here)
};

GRISU_INLINE void AssignU32(DiyInt& x, uint32_t value)
{
    x.bigits[0] = value;
    x.size = (value != 0) ? 1 : 0;
}

GRISU_INLINE void AssignU64(DiyInt& x, uint64_t value)
{
    x.bigits[0] = static_cast<uint32_t>(value);
    x.bigits[1] = static_cast<uint32_t>(value >> 32);
    x.size = (x.bigits[1] != 0) ? 2 : ((x.bigits[0] != 0) ? 1 : 0);
}

// x := A * x + B
GRISU_INLINE void MulAddU32(DiyInt& x, uint32_t A, uint32_t B = 0)
{
    GRISU_ASSERT(x.size >= 0);

    if (A == 1 && B == 0)
    {
        return;
    }
    if (A == 0 || x.size <= 0)
    {
        AssignU32(x, B);
        return;
    }

    uint32_t carry = B;
    for (int i = 0; i < x.size; ++i)
    {
        const uint64_t p = uint64_t{x.bigits[i]} * A + carry;
        x.bigits[i]      = static_cast<uint32_t>(p);
        carry            = static_cast<uint32_t>(p >> 32);
    }

    if (carry != 0)
    {
        GRISU_ASSERT(x.size < DiyInt::Capacity);
        x.bigits[x.size++] = carry;
    }
}

// x := x * 2^e2
GRISU_INLINE void MulPow2(DiyInt& x, int e2) // aka left-shift
{
    GRISU_ASSERT(x.size >= 0);
    GRISU_ASSERT(e2 >= 0);

    if (x.size <= 0 || e2 == 0)
        return;

    const int bigit_shift = static_cast<int>(static_cast<unsigned>(e2) / 32);
    const int bit_shift   = static_cast<int>(static_cast<unsigned>(e2) % 32);

    if (bit_shift > 0)
    {
        uint32_t carry = 0;
        for (int i = 0; i < x.size; ++i)
        {
            const uint32_t h = x.bigits[i] >> (32 - bit_shift);
            x.bigits[i]      = x.bigits[i] << bit_shift | carry;
            carry            = h;
        }

        if (carry != 0)
        {
            GRISU_ASSERT(x.size < DiyInt::Capacity);
            x.bigits[x.size++] = carry;
        }
    }

    if (bigit_shift > 0)
    {
        GRISU_ASSERT(bigit_shift <= DiyInt::Capacity);
        GRISU_ASSERT(x.size <= DiyInt::Capacity - bigit_shift);

        std::memmove(x.bigits + bigit_shift, x.bigits, sizeof(uint32_t) * static_cast<uint32_t>(x.size));
        std::memset(x.bigits, 0, sizeof(uint32_t) * static_cast<uint32_t>(bigit_shift));
        x.size += bigit_shift;
    }
}

// x := x * 5^e5
GRISU_INLINE void MulPow5(DiyInt& x, int e5)
{
    // TODO:
    // Optimize for large powers???

    static constexpr uint32_t kPow5_32[] = {
        1, // (unused)
        5,
        25,
        125,
        625,
        3125,
        15625,
        78125,
        390625,
        1953125,
        9765625,
        48828125,
        244140625,
        1220703125, // 5^13
    };

    if (x.size <= 0)
        return;

    GRISU_ASSERT(e5 >= 0);
    while (e5 > 0)
    {
        const int n = Min(e5, 13);
        MulAddU32(x, kPow5_32[n]);
        e5 -= n;
    }
}

GRISU_INLINE void Mul2(DiyInt& x)
{
    uint32_t carry = 0;
    for (int i = 0; i < x.size; ++i)
    {
        const uint32_t h = x.bigits[i] >> 31;
        x.bigits[i]      = x.bigits[i] << 1 | carry;
        carry            = h;
    }
    if (carry != 0)
    {
        GRISU_ASSERT(x.size < DiyInt::Capacity);
        x.bigits[x.size++] = carry;
    }
}

GRISU_INLINE void Mul10(DiyInt& x)
{
#if 0
    MulAddU32(x, 10, 0);
#else
    uint32_t carry = 0;
    for (int i = 0; i < x.size; ++i)
    {
        const uint64_t p = uint64_t{x.bigits[i]} * 10 + carry;
        x.bigits[i]      = static_cast<uint32_t>(p);
        carry            = static_cast<uint32_t>(p >> 32);
    }
    if (carry != 0)
    {
        GRISU_ASSERT(x.size < DiyInt::Capacity);
        x.bigits[x.size++] = carry;
    }
#endif
}

// x := 2^e2
GRISU_INLINE void AssignPow2(DiyInt& x, int e2)
{
    if (e2 == 0)
    {
        AssignU32(x, 1);
        return;
    }

    const int bigit_shift = static_cast<int>(static_cast<uint32_t>(e2)) / 32;
    const int bit_shift   = static_cast<int>(static_cast<uint32_t>(e2)) % 32;

    std::memset(x.bigits, 0, sizeof(uint32_t) * static_cast<uint32_t>(bigit_shift));

    x.bigits[bigit_shift] = 1u << bit_shift;
    x.size = bigit_shift + 1;
}

// x := 5^e5
GRISU_INLINE void AssignPow5(DiyInt& x, int e5)
{
    // TODO:
    // Optimize?!

    AssignU32(x, 1);
    MulPow5(x, e5);
}

// x := 10^e10
GRISU_INLINE void AssignPow10(DiyInt& x, int e10)
{
    // TODO:
    // Optimize?!

    AssignPow5(x, e10);
    MulPow2(x, e10);
}

// x := value * 2^e2
GRISU_INLINE void AssignU64MulPow2(DiyInt& x, uint64_t value, int e2)
{
    GRISU_ASSERT(e2 >= 0);

    if (value == 0 || e2 == 0)
    {
        AssignU64(x, value);
        return;
    }

    const int bigit_shift = static_cast<int>(static_cast<uint32_t>(e2)) / 32;
    const int bit_shift   = static_cast<int>(static_cast<uint32_t>(e2)) % 32;

    std::memset(x.bigits, 0, sizeof(uint32_t) * static_cast<uint32_t>(bigit_shift));

    const uint32_t lo = static_cast<uint32_t>(value);
    const uint32_t hi = static_cast<uint32_t>(value >> 32);
    if (bit_shift == 0)
    {
        // Relax: only write to x.bigits if neccessary?!
        GRISU_ASSERT(DiyInt::Capacity >= bigit_shift + 2);

        x.bigits[bigit_shift + 0] = lo;
        x.bigits[bigit_shift + 1] = hi;
        x.size = bigit_shift + ((hi != 0) ? 2 : 1);
    }
    else
    {
        // Relax: only write to x.bigits if neccessary?!
        GRISU_ASSERT(DiyInt::Capacity >= bigit_shift + 3);

        const uint32_t v0 = lo << bit_shift;
        const uint32_t v1 = hi << bit_shift | lo >> (32 - bit_shift);
        const uint32_t v2 =                   hi >> (32 - bit_shift);
        x.bigits[bigit_shift + 0] = v0;
        x.bigits[bigit_shift + 1] = v1;
        x.bigits[bigit_shift + 2] = v2;
        x.size = bigit_shift + ((v2 != 0) ? 3 : ((v1 != 0) ? 2 : 1));
    }
}

// x := value * 10^e10
GRISU_INLINE void AssignU64MulPow10(DiyInt& x, uint64_t value, int e10)
{
    AssignU64MulPow2(x, value, e10);
    MulPow5(x, e10);
}

// x := 2^e2 * 5^e5
GRISU_INLINE void AssignPow2MulPow5(DiyInt& x, int e2, int e5)
{
    AssignPow2(x, e2);
    MulPow5(x, e5);

    //AssignPow5(x, e5);
    //MulPow2(x, e2);
}

// Returns the number of leading 0-bits in x, starting at the most significant bit position.
// If x is 0, the result is undefined.
GRISU_INLINE int CountLeadingZeros32(uint32_t x)
{
    GRISU_ASSERT(x != 0);

#if defined(__GNUC__)
    return __builtin_clz(x);
#elif defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
    return static_cast<int>(__lzcnt(x));
#else
    int z = 0;
    while ((x >> 31) == 0) {
        x <<= 1;
        ++z;
    }
    return z;
#endif
}

// q, r = divmod(u, v)
// u := r
// return q
// PRE: 0 <= q <= 9
GRISU_INLINE uint32_t DivMod(DiyInt& u, DiyInt const& v)
{
    GRISU_ASSERT(u.size > 0);
    GRISU_ASSERT(v.size > 0);
    GRISU_ASSERT(u.bigits[u.size - 1] != 0);
    GRISU_ASSERT(v.bigits[v.size - 1] != 0);

    const int m = u.size;
    const int n = v.size;
    if (m < n)
    {
        return 0;
    }

    GRISU_ASSERT(m >= n);
    GRISU_ASSERT(n >= 1);

    //--------------------------------------------------------------------------
    // D0.
    //
    // Handle the case of a single digit division first.
    // Note that this step is not only here for performance reasons. The
    // algorithm described below requires at least two digits in the
    // denominator.
    if (n == 1)
    {
        const uint32_t den = v.bigits[0];

        uint32_t q = 0;
        uint32_t r = 0;
        for (int i = m - 1; i >= 0; --i)
        {
            const uint64_t t = (uint64_t{r} << 32) | u.bigits[i];
            q = static_cast<uint32_t>(t / den);
            r = static_cast<uint32_t>(t % den);
        }
        AssignU32(u, r);
        return q;
    }

    GRISU_ASSERT(n >= 2);
    GRISU_ASSERT(DiyInt::Capacity >= m + 1);
    u.bigits[m] = 0;

    //
    // XXX:
    //
    // Is the normalization step only required once in Dragon4 before starting
    // the digit generation procedure?!?!?!
    // It might be more efficient to shift r and s and delta instead of
    // shifting the leading bits in each iteration...
    //

    //--------------------------------------------------------------------------
    // D1. [Normalize.]
    //
    // Set d := b / (v[n-1] + 1). Then set
    //    u' := (u[0] u[1] ... u[m-1] u[m])_b = d * (u[0] u[1] ... u[m-1])_b,
    //    v' := (v[0] v[1] ... v[n-1]     )_b = d * (v[0] v[1] ... v[n-1])_b.
    //
    // Note the introduction of a new digit position u[m] at the right of
    // u[m-1]; if d = 1, all we need to do in this step is set u[m] := 0.
    //
    // On a binary computer it may be preferable to choose d to be a power of 2
    // instead of using the value suggested here; any value of d that results in
    // v[n-1] >= b/2 will suffice.
    //

    // This normalization step is required only to efficiently estimate the
    // quotient q' (see below). Is is not necessary for the other steps of the
    // algorithm.
    // Instead of shifting both u and v into u' and v' resp., the required
    // digits of u' and v' are computed when they are needed.
    //
    // The variables vK here denote v'[n - K], where K = 1, 2, and v' denotes
    // the normalized value d * v.

    uint32_t v1 = v.bigits[n - 1];
    uint32_t v2 = v.bigits[n - 2];

    const int shift = CountLeadingZeros32(v1);
    if (shift > 0)
    {
        const uint32_t v3 = (n >= 3) ? v.bigits[n - 3] : 0;
        v1 = (v1 << shift) | (v2 >> (32 - shift));
        v2 = (v2 << shift) | (v3 >> (32 - shift));
    }
    // v1 and v2 now contain the leading digits of v'.

    //--------------------------------------------------------------------------
    // D2. [Initialize.]
    //
    // Set j := m - n.
    //
    // The loop on j, steps D2 through D7, will be essentially a division of
    // (u[j] u[j+1] ... u[j+n])_b by (v[0] v[1] ... v[n-1])_b to get a single
    // quotient digit.
    //

    //--------------------------------------------------------------------------
    // D3. [Calculate q'.]
    //
    // If u[j+n] = v[n-1], set
    //    q' := b - 1;
    // otherwise set
    //    q' := (u[j+n] * b + u[j+n-1]) / v[n-1].
    // Now test if
    //    q' * v[n-2] > (u[j+n] * b + u[j+n-1] - q' * v[n-1]) * b + u[j+n-2];
    // if so, decrease q' by 1 and repeat this test.
    //
    // The latter test determines at high speed most of the cases in which the
    // trial value q' is one too large, and it eliminates all cases where q' is
    // two too large.
    //

    // The variable uK here denotes u'[j + n - K], where K = 0, 1, 2, and u'
    // denotes the scaled value d * u.

    uint32_t u0 = u.bigits[n];
    uint32_t u1 = u.bigits[n - 1];
    uint32_t u2 = u.bigits[n - 2];

    if (shift > 0)
    {
        GRISU_ASSERT((u0 >> (32 - shift)) == 0);

        const uint32_t u3 = (n >= 3) ? u.bigits[n - 3] : 0;
        u0 = (u0 << shift) | (u1 >> (32 - shift));
        u1 = (u1 << shift) | (u2 >> (32 - shift));
        u2 = (u2 << shift) | (u3 >> (32 - shift));
    }
    // u0, u1 and u2 now contain the leading digits of u'.

    // NB: Use repeated subtraction for division to avoid a 64-bit div.
    uint64_t rp = uint64_t{u0} << 32 | u1;
    uint32_t qp = 0;
    while (rp >= v1)
    {
        rp -= v1;
        qp++;
    }
    GRISU_ASSERT(qp <= 10);

    if (uint64_t{qp} * v2 > (rp << 32 | u2))
    {
        GRISU_ASSERT(qp > 0);
        qp--;
    }
    GRISU_ASSERT(qp <= 9);

    //--------------------------------------------------------------------------
    // D4. [Multiply and subtract.]
    //
    // Replace
    //    (u[j] u[j+1] ... u[j+n])_b
    //      := (u[j] u[j+1] ... u[j+n])_b - q' * (v[0] v[1] ... v[n-1] 0)
    //
    // This step consists of a simple multiplication by a one-place number,
    // combined with subtraction. The digits (u[j] ... u[j+n])_b should be kept
    // positive; if the result of this step is actually negative,
    // (u[j] ... u[j+n])_b should be left as the true value plus b^(n+1), i.e.,
    // as the b's complement of the true value, and a "borrow" to the right
    // should be remembered.
    //

    if (qp == 0)
    {
        // No need to multiply.
        return 0;
    }

    uint32_t borrow = 0;
    for (int i = 0; i < n; ++i)
    {
        const uint32_t ui = u.bigits[i];
        const uint32_t vi = v.bigits[i];
        const uint64_t p  = uint64_t{qp} * vi + borrow;
        const uint32_t si = static_cast<uint32_t>(p);
        borrow            = static_cast<uint32_t>(p >> 32);
        const uint32_t di = ui - si;
        borrow           += di > ui;
        u.bigits[i]       = di;
    }
    // vn = 0:
    const uint32_t un = u.bigits[n];
    const uint32_t dn = un - borrow;
    u.bigits[n] = un;

    //--------------------------------------------------------------------------
    // D5. [Test remainder.]
    //
    // Set q[j] := q'. If the result of step D4 was negative, go to step D6;
    // otherwise go on to step D7.
    //

    const bool was_negative = (dn > un);
    if (was_negative)
    {
        //----------------------------------------------------------------------
        // D6. [Add back.]
        //
        // Decrease q[j] by 1, and add (v[0] v[1] ... v[n-1] 0)_b to
        // (u[j] u[j+1] ... u[j+n])_b. (A carry will occur to the right of
        // u[j+n], and it should be ignored since it cancels with the "borrow"
        // that occurred in D4.)
        //
        // The probability that this step is necessary is very small, on the
        // order of only 2/b; test data that activates this step should
        // therefore be specifically contrived when debugging.
        //

        qp--;

        uint32_t carry = 0;
        for (int i = 0; i < n; ++i)
        {
            const uint32_t ui = u.bigits[i];
            const uint32_t vi = v.bigits[i];
            const uint64_t s  = uint64_t{ui} + vi + carry;
            u.bigits[i]       = static_cast<uint32_t>(s);
            carry             = static_cast<uint32_t>(s >> 32);
        }
        // vn = 0:
        u.bigits[n] += carry;
    }

    //--------------------------------------------------------------------------
    // D7. [Loop on j.]
    //
    // Decrease j by one. Now if j >= 0, go back to D3.
    //

    //--------------------------------------------------------------------------
    // D8. [Unnormalize.]
    //
    // Now (q[0] q[1] ... q[m-n])_b is the desired quotient, and the desired
    // remainder may be obtained by dividing (u[0] u[1] ... u[n-1])_b by d.
    //

    // We didn't multiply in the first place, so we don't need to divide here.

    // Still need to clamp the remainder.
    int k = n;
    for ( ; k > 0 && u.bigits[k - 1] == 0; --k)
    {
    }
    u.size = k;

    return qp;
}

GRISU_INLINE int Compare(DiyInt const& lhs, DiyInt const& rhs)
{
    const int n1 = lhs.size;
    const int n2 = rhs.size;

    if (n1 < n2) return -1;
    if (n1 > n2) return +1;

    for (int i = n1 - 1; i >= 0; --i)
    {
        const uint32_t b1 = lhs.bigits[i];
        const uint32_t b2 = rhs.bigits[i];

        if (b1 < b2) return -1;
        if (b1 > b2) return +1;
    }

    return 0;
}

// Returns Compare(a + b, c)
GRISU_INLINE int CompareAdd(DiyInt const& a, DiyInt const& b, DiyInt const& c)
{
    // NB:
    // This function is only ever called with a >= c, which implies a.size >= c.size.
    GRISU_ASSERT(c.size >= a.size);

    const int na = a.size;
    const int nb = b.size;
    const int nc = c.size;

    const int m = Max(na, nb);
    if (m + 1 < nc)
        return -1; // s = (a + b) cannot be larger or equal to c
    if (m > nc)
        return +1; // max(a, b) > c

    // Perform a (partial) left-to-right subtraction, propagating a borrow digit
    // (base B = 2^32) along to the right, stopping as soon as s > c or s < c.

    uint64_t borrow = 0;
    for (int i = nc - 1; i >= 0; --i)
    {
        // Invariant:
        // The leading digits s[i+1],s[i+2],... of s and the leading digits
        // c[i+1],c[i+2],... (after possibly subtracting a borrow) are equal.

        GRISU_ASSERT(borrow == 0 || borrow == 1);
        const uint64_t ci = borrow << 32 | c.bigits[i];
        const uint32_t ai = i < na ? a.bigits[i] : 0;
        const uint32_t bi = i < nb ? b.bigits[i] : 0;
        const uint64_t si = static_cast<uint64_t>(ai) + bi;
        const uint64_t di = ci - si;
//      if (ci < si)
        if (di > ci)
        {
            // Since all the leading digits are equal, this implies c < s,
            // or a + b > c.
            return +1;
        }
        if (di > 1)
        {
            // In this case, the trailing digits s[i-1],s[i-2],... cannot
            // possibly compensate the difference: therefore c > s, or a + b < c.
            return -1;
        }

        // di == 0 or di == 1.
        // If di == 1, borrow B = 2^32 from ci and add to c[i-1], which restores
        // the invariant.
        //  c:      1   2   9   9  ==>  1   1  19   9
        //  s:      1   1  12   3       1   1  12   3
        //              ^                   ^
        //              i                   i
        borrow = di;
    }

//  return borrow == 0 ? 0 : -1;
    return -static_cast<int>(borrow);
}

GRISU_INLINE int EffectivePrecision(uint64_t f)
{
    GRISU_ASSERT(f != 0);
    return 64 - CountLeadingZeros64(f);
}

GRISU_INLINE int ComputeInitialValuesAndEstimate(DiyInt& r, DiyInt& s, DiyInt& delta, uint64_t f, int e, bool lowerBoundaryIsCloser)
{
    const int boundaryShift = lowerBoundaryIsCloser ? 2 : 1;
    const int p = EffectivePrecision(f);
    GRISU_ASSERT(p >= 1);
    GRISU_ASSERT(p <= 53);
    const int k = CeilLog10Pow2(e + (p - 1));

//  const int cmpf = CompareEstimate((f << boundaryShift) + 1, boundaryShift - e, k);
    if (e >= 0)
    {
        GRISU_ASSERT(e >= 0);
        GRISU_ASSERT(e <= 971);
        GRISU_ASSERT(k >= 0);
        GRISU_ASSERT(k <= 308);

        // r = f * 2^(boundaryShift + e)
        AssignU64MulPow2(r, f << boundaryShift, e);
        // s = 2^boundaryShift * 10^k
        AssignPow2MulPow5(s, boundaryShift + k, k);
        // delta = 2^e
        AssignPow2(delta, e);
    }
    else if (k < 0)
    {
        GRISU_ASSERT(e >= -1074);
        GRISU_ASSERT(e <= -1);
        GRISU_ASSERT(k >= -323);
        GRISU_ASSERT(k <= -1);

        // r = f * 2^boundaryShift * 10^(-k)
        AssignU64MulPow10(r, f << boundaryShift, -k);
        // s = 2^(boundaryShift - e)
        AssignPow2(s, boundaryShift - e);
        // delta = 10^(-k)
        AssignPow10(delta, -k);
    }
    else
    {
        GRISU_ASSERT(e >= -55);
        GRISU_ASSERT(e <= -1);
        GRISU_ASSERT(k >= 0);
        GRISU_ASSERT(k <= 16);

        // r = f * 2^boundaryShift
        AssignU64(r, f << boundaryShift);
        // s = 2^(boundaryShift - e) * 10^k
        AssignPow2MulPow5(s, boundaryShift - e + k, k);
        // delta = 1
        AssignU32(delta, 1);
    }

    return k;
}

GRISU_INLINE char* Dragon4(char* digits, int& num_digits, int& exponent, uint64_t f, int e, bool acceptBounds, bool lowerBoundaryIsCloser)
{
    DiyInt r;
    DiyInt s;
    DiyInt delta;

    //
    // Compute initial values.
    // Estimate k.
    //
    int k = ComputeInitialValuesAndEstimate(r, s, delta, f, e, lowerBoundaryIsCloser);

    //
    // Fixup, in case k is too low.
    //
    const int cmpf = CompareAdd(r, delta, s);
    if (acceptBounds ? (cmpf >= 0) : (cmpf > 0))
    {
        Mul10(s);
        k++;
    }

    //
    // Generate digits from left to right.
    //
    Mul10(r);       // (Move into init step above?)
    Mul10(delta);   // (Move into init step above?)

    int length = 0;
    for (;;)
    {
        GRISU_ASSERT(length < 17);
        GRISU_ASSERT(r.size > 0);

        // q = r / s
        // r = r % s
        uint32_t q = DivMod(r, s);
        GRISU_ASSERT(q <= 9);

        const int cmp1 = Compare(r, delta);
        if (lowerBoundaryIsCloser)
        {
            Mul2(delta);
        }
        const int cmp2 = CompareAdd(r, delta, s);

        const bool tc1 = acceptBounds ? (cmp1 <= 0) : (cmp1 < 0);
        const bool tc2 = acceptBounds ? (cmp2 >= 0) : (cmp2 > 0);
        if (tc1 && tc2)
        {
            // Return the number closer to v.
            // If the two are equidistant from v, use _some_ strategy to break
            // the tie.
            const int cmpr = CompareAdd(r, r, s);
            if (cmpr > 0 || (cmpr == 0 && q % 2 != 0))
            {
                q++;
            }
        }
        else if (!tc1 && tc2)
        {
            q++;
        }

        GRISU_ASSERT(q <= 9);
        digits[length++] = static_cast<char>(q + '0');
        k--;

        if (tc1 || tc2)
            break;

        Mul10(r);
        MulAddU32(delta, lowerBoundaryIsCloser ? 5 : 10);
    }

    num_digits = length;
    exponent = k;

    return digits + length;
}

} // namespace impl

//==================================================================================================
// ToDigits
//==================================================================================================

// v = digits * 10^exponent
// num_digits is the length of the buffer (number of decimal digits)
// PRE: The buffer must be large enough, i.e. >= max_digits10.
// PRE: value must be finite and strictly positive.
template <typename Float>
GRISU_INLINE void ToDigits(char* buffer, int& num_digits, int& exponent, Float value)
{
    using Fp = grisu3::impl::IEEE<Float>;

    static_assert(grisu3::impl::DiyFp::SignificandSize >= std::numeric_limits<Float>::digits + 3,
        "Grisu3 requires at least three extra bits of precision");

    GRISU_ASSERT(grisu3::impl::IEEE<Float>(value).IsFinite());
    GRISU_ASSERT(value > 0);

    // Compute the boundaries of 'value'.
    // These boundaries obviously depend on the type 'Float'.
    //
    // If the boundaries of 'value' are always computed for double-precision numbers, regardless of
    // type of 'Float', all single-precision numbers can be recovered using strtod (and strtof).
    // However, the resulting decimal representations are not exactly "short".
    //
    // On the other hand, if the boundaries are computed for single-precision numbers, there is a
    // single number (7.0385307e-26f) which can't be recovered using strtod (instead of strtof).
    // The resulting 'double' when cast to 'float' is off by 1 ulp, i.e. for f = 7.0385307e-26f,
    //      f != (float)strtod(ftoa(f))
    // For all other single-precision numbers, equality holds.

    const auto boundaries = grisu3::impl::ComputeBoundaries(value);

    const bool ok = grisu3::impl::Grisu3(buffer, num_digits, exponent, boundaries.m_minus, boundaries.v, boundaries.m_plus);
    if (!ok)
    {
        const auto v = grisu3::impl::DiyFpFromFloat(value);

        const bool isEven = (v.f % 2 == 0);
        const bool acceptBounds = isEven;
        const bool lowerBoundaryIsCloser = (v.f == Fp::HiddenBit && v.e > Fp::MinExponent);

        grisu3::impl::Dragon4(buffer, num_digits, exponent, v.f, v.e, acceptBounds, lowerBoundaryIsCloser);
    }

    GRISU_ASSERT(num_digits > 0);
    GRISU_ASSERT(num_digits <= std::numeric_limits<Float>::max_digits10);
}

//==================================================================================================
// ToChars
//==================================================================================================

namespace impl {

// Appends a decimal representation of 'value' to buffer.
// Returns a pointer to the element following the digits.
//
// PRE: -1000 < value < 1000
GRISU_INLINE char* ExponentToString(char* buffer, int value)
{
    GRISU_ASSERT(value > -1000);
    GRISU_ASSERT(value <  1000);

    int n = 0;

    if (value < 0)
    {
        buffer[n++] = '-';
        value = -value;
    }
    else
    {
        buffer[n++] = '+';
    }

    const uint32_t k = static_cast<uint32_t>(value);
    if (k < 10)
    {
        buffer[n++] = static_cast<char>('0' + k);
    }
    else if (k < 100)
    {
        Utoa_2Digits(buffer + n, k);
        n += 2;
    }
    else
    {
        const uint32_t r = k % 10;
        const uint32_t q = k / 10;
        Utoa_2Digits(buffer + n, q);
        n += 2;
        buffer[n++] = static_cast<char>('0' + r);
    }

    return buffer + n;
}

GRISU_INLINE char* FormatFixed(char* buffer, intptr_t num_digits, intptr_t decimal_point, bool force_trailing_dot_zero)
{
    GRISU_ASSERT(buffer != nullptr);
    GRISU_ASSERT(num_digits >= 1);

    if (num_digits <= decimal_point)
    {
        // digits[000]
        // GRISU_ASSERT(buffer_capacity >= decimal_point + (force_trailing_dot_zero ? 2 : 0));

        std::memset(buffer + num_digits, '0', static_cast<size_t>(decimal_point - num_digits));
        buffer += decimal_point;
        if (force_trailing_dot_zero)
        {
            *buffer++ = '.';
            *buffer++ = '0';
        }
        return buffer;
    }
    else if (0 < decimal_point)
    {
        // dig.its
        // GRISU_ASSERT(buffer_capacity >= length + 1);

        std::memmove(buffer + (decimal_point + 1), buffer + decimal_point, static_cast<size_t>(num_digits - decimal_point));
        buffer[decimal_point] = '.';
        return buffer + (num_digits + 1);
    }
    else // decimal_point <= 0
    {
        // 0.[000]digits
        // GRISU_ASSERT(buffer_capacity >= 2 + (-decimal_point) + length);

        std::memmove(buffer + (2 + -decimal_point), buffer, static_cast<size_t>(num_digits));
        buffer[0] = '0';
        buffer[1] = '.';
        std::memset(buffer + 2, '0', static_cast<size_t>(-decimal_point));
        return buffer + (2 + (-decimal_point) + num_digits);
    }
}

GRISU_INLINE char* FormatScientific(char* buffer, intptr_t num_digits, int exponent, bool force_trailing_dot_zero)
{
    GRISU_ASSERT(buffer != nullptr);
    GRISU_ASSERT(num_digits >= 1);

    if (num_digits == 1)
    {
        // dE+123
        // GRISU_ASSERT(buffer_capacity >= num_digits + 5);

        buffer += 1;
        if (force_trailing_dot_zero)
        {
            *buffer++ = '.';
            *buffer++ = '0';
        }
    }
    else
    {
        // d.igitsE+123
        // GRISU_ASSERT(buffer_capacity >= num_digits + 1 + 5);

        std::memmove(buffer + 2, buffer + 1, static_cast<size_t>(num_digits - 1));
        buffer[1] = '.';
        buffer += 1 + num_digits;
    }

    buffer[0] = 'e';
    buffer = ExponentToString(buffer + 1, exponent);

    return buffer;
}

// Format the digits similar to printf's %g style.
GRISU_INLINE char* Format(char* buffer, int num_digits, int exponent, bool force_trailing_dot_zero)
{
    const int decimal_point = num_digits + exponent;

    // NB:
    // These are the values used by JavaScript's ToString applied to Number
    // type. Printf uses the values -4 and max_digits10 resp. (sort of).
    constexpr int kMinExp = -6;
    constexpr int kMaxExp = 21;

    const bool use_fixed = kMinExp < decimal_point && decimal_point <= kMaxExp;

    return use_fixed
        ? grisu3::impl::FormatFixed(buffer, num_digits, decimal_point, force_trailing_dot_zero)
        : grisu3::impl::FormatScientific(buffer, num_digits, decimal_point - 1, force_trailing_dot_zero);
}

} // namespace impl

// Generates a decimal representation of the floating-point number `value` in 'buffer'.
// Note: The result is _not_ null-terminated.
//
// PRE: The buffer must be large enough (32 bytes is sufficient).
template <typename Float>
GRISU_INLINE char* ToChars(char* buffer, Float value, bool force_trailing_dot_zero = false)
{
    using Fp = grisu3::impl::IEEE<Float>;
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

    int num_digits = 0;
    int exponent = 0;
    grisu3::ToDigits(buffer, num_digits, exponent, value);

    return grisu3::impl::Format(buffer, num_digits, exponent, force_trailing_dot_zero);
}

} // namespace grisu3

//char* Dtoa(char* buffer, double value)
//{
//    return grisu3::ToChars(buffer, value);
//}

//char* Ftoa(char* buffer, float value)
//{
//    return grisu3::ToChars(buffer, value);
//}
