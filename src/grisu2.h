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

#ifndef GRISU2_ASSERT
#define GRISU2_ASSERT(X) assert(X)
#endif

#ifndef GRISU2_UNNAMED_NAMESPACE
#define GRISU2_UNNAMED_NAMESPACE 0
#endif

#if GRISU2_UNNAMED_NAMESPACE
namespace {
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

namespace impl {

template <typename Dest, typename Source>
inline Dest ReinterpretBits(Source source)
{
    static_assert(sizeof(Dest) == sizeof(Source), "size mismatch");

    Dest dest;
    std::memcpy(&dest, &source, sizeof(Source));
    return dest;
}

inline char* Utoa_2Digits(char* buf, uint32_t digits)
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

    GRISU2_ASSERT(digits < 100);
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
inline DiyFp Subtract(DiyFp x, DiyFp y)
{
    GRISU2_ASSERT(x.e == y.e);
    GRISU2_ASSERT(x.f >= y.f);

    return DiyFp(x.f - y.f, x.e);
}

// Returns x * y.
// The result is rounded (ties up). (Only the upper q bits are returned.)
inline DiyFp Multiply(DiyFp x, DiyFp y)
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
inline int CountLeadingZeros64(uint64_t x)
{
    GRISU2_ASSERT(x != 0);

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
inline DiyFp Normalize(DiyFp x)
{
    static_assert(DiyFp::SignificandSize == 64, "internal error");

    const int lz = CountLeadingZeros64(x.f);
    return DiyFp(x.f << lz, x.e - lz);
}

// Normalize x such that the result has the exponent E.
// PRE: e >= x.e and the upper e - x.e bits of x.f must be zero.
inline DiyFp NormalizeTo(DiyFp x, int e)
{
    const int delta = x.e - e;

    GRISU2_ASSERT(delta >= 0);
    GRISU2_ASSERT(((x.f << delta) >> delta) == x.f);

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
inline DiyFp DiyFpFromFloat(Float value)
{
    using Fp = IEEE<Float>;

    const Fp v(value);
    GRISU2_ASSERT(v.IsFinite());
    GRISU2_ASSERT(!v.SignBit());

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
inline Boundaries ComputeBoundaries(Float value)
{
    using Fp = IEEE<Float>;

    GRISU2_ASSERT(Fp(value).IsFinite());
    GRISU2_ASSERT(value > 0);

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

// Technically, right-shift of negative integers is implementation defined...
// Portable SAR.
// Should easily be optimized into SAR (or equivalent) instruction.
inline int SAR(int x, int n)
{
//  return x >> n;
    return x < 0 ? ~(~x >> n) : (x >> n);
}

// Returns: floor(log_2(10^e))
inline int FloorLog2Pow10(int e)
{
    GRISU2_ASSERT(e >= -1233);
    GRISU2_ASSERT(e <=  1232);
    return SAR(e * 1741647, 19);
}

// Returns: ceil(log_10(2^e))
inline int CeilLog10Pow2(int e)
{
    GRISU2_ASSERT(e >= -2620);
    GRISU2_ASSERT(e <=  2620);
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
inline CachedPower GetCachedPower(int index)
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

    GRISU2_ASSERT(index >= 0);
    GRISU2_ASSERT(index < kCachedPowersSize);

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
inline CachedPower GetCachedPowerForBinaryExponent(int e)
{
    // For double: -1137 <= e <= 960 ==> -307 <= k <= 324 ==>  0 <= index <= 78
    // For single:  -180 <= e <=  96 ==>  -47 <= k <= 36  ==> 32 <= index <= 42
    GRISU2_ASSERT(e >= -1137);
    GRISU2_ASSERT(e <=   960);

    const int k = CeilLog10Pow2(kAlpha - e - 1);
    GRISU2_ASSERT(k >= kCachedPowersMinDecExp - (kCachedPowersDecExpStep - 1));
    GRISU2_ASSERT(k <= kCachedPowersMaxDecExp);

    const int index = static_cast<int>( static_cast<unsigned>(k - (kCachedPowersMinDecExp - (kCachedPowersDecExpStep - 1))) / kCachedPowersDecExpStep );
    GRISU2_ASSERT(index >= 0);
    GRISU2_ASSERT(index < kCachedPowersSize);
    static_cast<void>(kCachedPowersSize);

    const auto cached = GetCachedPower(index);
    GRISU2_ASSERT(kAlpha <= cached.e + e + 64);
    GRISU2_ASSERT(kGamma >= cached.e + e + 64);

    // NB:
    // Actually this function returns c, such that -60 <= e_c + e + 64 <= -34.
    GRISU2_ASSERT(-60 <= cached.e + e + 64);
    GRISU2_ASSERT(-34 >= cached.e + e + 64);

    return cached;
}

inline char* GenerateIntegralDigits(char* buf, uint32_t n)
{
//  GRISU2_ASSERT(n <= 798336123);
    GRISU2_ASSERT(n <= 999999999);

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

    if (n >= 1000000) goto L_7_digits;
    if (n >=  100000) goto L_6_digits;
    if (n >=   10000) goto L_5_digits;
    if (n >=    1000) goto L_4_digits;
    if (n >=     100) goto L_3_digits;
    if (n >=      10) goto L_2_digits;
    goto L_1_digit;
}

// Modifies the generated digits in the buffer to approach (round towards) w.
//
// Input:
//  * digits of H/10^kappa in [digits, digits + num_digits)
//  * distance    = (H - w) * unit
//  * delta       = (H - L) * unit
//  * rest        = (H - digits * 10^kappa) * unit
//  * ten_kappa   = 10^kappa * unit
inline void Grisu2Round(char* digits, int num_digits, uint64_t distance, uint64_t delta, uint64_t rest, uint64_t ten_kappa)
{
    GRISU2_ASSERT(num_digits >= 1);
    GRISU2_ASSERT(distance <= delta);
    GRISU2_ASSERT(rest <= delta);
    GRISU2_ASSERT(ten_kappa > 0);

    // By generating the digits of H we got the largest (closest to H) buffer
    // that is still in the interval [L, H]. In the case where w < B <= H we
    // try to decrement the buffer.
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
    // representation stored in the buffer.
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

    int digit = digits[num_digits - 1] - '0';

    while (rest < distance
        && delta - rest >= ten_kappa
        && (rest + ten_kappa <= distance || rest + ten_kappa - distance < distance - rest))
    {
        GRISU2_ASSERT(digit != 0);
        digit--;
        rest += ten_kappa;
    }

    digits[num_digits - 1] = static_cast<char>('0' + digit);
}

// Generates V = digits * 10^exponent, such that L <= V <= H.
// L and H must be normalized and share the same exponent -60 <= e <= -32.
inline void Grisu2DigitGen(char* digits, int& num_digits, int& exponent, DiyFp L, DiyFp w, DiyFp H)
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

    GRISU2_ASSERT(H.e >= kAlpha);
    GRISU2_ASSERT(H.e <= kGamma);
    GRISU2_ASSERT(H.e == L.e);
    GRISU2_ASSERT(H.e == w.e);

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

    GRISU2_ASSERT(p1 >= 4);            // (2^(64-2) - 1) >> 60
    GRISU2_ASSERT(p1 <= 798336123);    // Not trivial. Depends on the implementation of GetCachedPowerForBinaryExponent!

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

    // The common case is that all the digits of p1 are needed.
    // Optimize for this case and correct later if required.
    num_digits = static_cast<int>(GenerateIntegralDigits(digits, p1) - digits);

    if (p2 > delta)
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
        // and stop as soon as 10^-m * r * 2^e <= delta * 2^e

        // unit = 1
        int m = 0;
        for (;;)
        {
            // !!! GRISU2_ASSERT(num_digits < max_digits10) !!!
            GRISU2_ASSERT(num_digits < 17);

            //
            //      H = digits * 10^-m + 10^-m * (d[-m-1] / 10 + d[-m-2] / 10^2 + ...) * 2^e
            //        = digits * 10^-m + 10^-m * (p2                                 ) * 2^e
            //        = digits * 10^-m + 10^-m * (1/10 * (10 * p2)                   ) * 2^e
            //        = digits * 10^-m + 10^-m * (1/10 * ((10*p2 div 2^-e) * 2^-e + (10*p2 mod 2^-e)) * 2^e
            //
            GRISU2_ASSERT(p2 <= 0xFFFFFFFFFFFFFFFFull / 10);
            p2 *= 10;
            const uint64_t d = p2 >> -one.e;     // d = (10 * p2) div 2^-e
            const uint64_t r = p2 & (one.f - 1); // r = (10 * p2) mod 2^-e
            GRISU2_ASSERT(d <= 9);
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

            // Check if enough digits have been generated.
            //
            //      10^-m * p2 * 2^e <= delta * 2^e
            //              p2 * 2^e <= 10^m * delta * 2^e
            //                    p2 <= 10^m * delta
            if (p2 <= delta)
            {
                // V = digits * 10^-m, with L <= V <= H.
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
    else // p2 <= delta
    {
        GRISU2_ASSERT((uint64_t{p1} << -one.e) + p2 > delta); // Loop terminates.

        // In this case: Too many digits of p1 might have been generated.
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

        const int k = num_digits;
        GRISU2_ASSERT(k >= 0);
        GRISU2_ASSERT(k <= 9);

        rest = p2;

        // 10^n is now 1 ulp in the decimal representation V. The rounding
        // procedure works with DiyFp's with an implicit exponent of e.
        //
        //      10^n = (10^n * 2^-e) * 2^e = ten_kappa * 2^e
        //
        ten_kappa = one.f; // Start with 2^-e

        for (int n = 0; /**/; ++n)
        {
            GRISU2_ASSERT(n <= k - 1);
            GRISU2_ASSERT(rest <= delta);

            // rn = d[n]...d[0] * 2^-e + p2
            const uint32_t dn = static_cast<uint32_t>(digits[k - 1 - n] - '0');
            const uint64_t rn = dn * ten_kappa + rest;

            if (rn > delta)
            {
                num_digits = k - n;
                exponent = n;
                break;
            }

            rest = rn;
            ten_kappa *= 10;
        }
    }

    // The buffer now contains a correct decimal representation of the input
    // number w = digits * 10^exponent.

    Grisu2Round(digits, num_digits, distance, delta, rest, ten_kappa);
}

// v = digits * 10^exponent
// length is the length of the buffer (number of decimal digits)
// The buffer must be large enough, i.e. >= max_digits10.
inline void Grisu2(char* digits, int& num_digits, int& exponent, DiyFp m_minus, DiyFp v, DiyFp m_plus)
{
    GRISU2_ASSERT(m_plus.e == m_minus.e);
    GRISU2_ASSERT(m_plus.e == v.e);

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
    GRISU2_ASSERT(w_plus.e >= kAlpha);
    GRISU2_ASSERT(w_plus.e <= kGamma);

    // Note:
    // The result of Multiply() is **NOT** neccessarily normalized.
    // But since m+ and c are normalized, w_plus.f >= 2^(q - 2).
    GRISU2_ASSERT(w_plus.f >= (uint64_t{1} << (64 - 2)));

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

    Grisu2DigitGen(digits, num_digits, exponent, L, w, H);
    // w = digits * 10^exponent

    // v = w * 10^k
    exponent += -cached.k; // cached.k = -k
    // v = digits * 10^exponent
}

} // namespace impl

// Maximum number of digits which will be written by DoubleToDigits.
// == numeric_limits<double>::max_digits10
constexpr int kDoubleToDigitsMaxLength = 17;

// v = digits * 10^exponent
// num_digits is the length of the buffer (number of decimal digits)
// PRE: The buffer must be large enough, i.e. >= max_digits10.
// PRE: value must be finite and strictly positive.
template <typename Float>
inline char* DoubleToDigits(char* next, char* last, int& num_digits, int& exponent, Float value)
{
    static_assert(grisu2::impl::DiyFp::SignificandSize >= std::numeric_limits<Float>::digits + 3,
        "Grisu2 requires at least three extra bits of precision");

    GRISU2_ASSERT(last - next >= kDoubleToDigitsMaxLength);
    GRISU2_ASSERT(grisu2::impl::IEEE<Float>(value).IsFinite());
    GRISU2_ASSERT(value > 0);
    static_cast<void>(last); // Fix warning

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
#if 0
    const auto boundaries = grisu2::impl::ComputeBoundaries(static_cast<double>(value));
#else
    const auto boundaries = grisu2::impl::ComputeBoundaries(value);
#endif

    grisu2::impl::Grisu2(next, num_digits, exponent, boundaries.m_minus, boundaries.v, boundaries.m_plus);

    GRISU2_ASSERT(num_digits > 0);
    GRISU2_ASSERT(num_digits <= kDoubleToDigitsMaxLength);

    return next + num_digits;
}

//==================================================================================================
// PositiveDtoa
//==================================================================================================

namespace impl {

// Appends a decimal representation of 'value' to buffer.
// Returns a pointer to the element following the digits.
//
// PRE: -1000 < value < 1000
inline char* ExponentToString(char* buffer, int value)
{
    GRISU2_ASSERT(value > -1000);
    GRISU2_ASSERT(value <  1000);

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

inline char* FormatFixed(char* buffer, intptr_t num_digits, intptr_t decimal_point, bool force_trailing_dot_zero)
{
    GRISU2_ASSERT(buffer != nullptr);
    GRISU2_ASSERT(num_digits >= 1);

    if (num_digits <= decimal_point)
    {
        // digits[000]
        // GRISU2_ASSERT(buffer_capacity >= decimal_point + (force_trailing_dot_zero ? 2 : 0));

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
        // GRISU2_ASSERT(buffer_capacity >= length + 1);

        std::memmove(buffer + (decimal_point + 1), buffer + decimal_point, static_cast<size_t>(num_digits - decimal_point));
        buffer[decimal_point] = '.';
        return buffer + (num_digits + 1);
    }
    else // decimal_point <= 0
    {
        // 0.[000]digits
        // GRISU2_ASSERT(buffer_capacity >= 2 + (-decimal_point) + length);

        std::memmove(buffer + (2 + -decimal_point), buffer, static_cast<size_t>(num_digits));
        buffer[0] = '0';
        buffer[1] = '.';
        std::memset(buffer + 2, '0', static_cast<size_t>(-decimal_point));
        return buffer + (2 + (-decimal_point) + num_digits);
    }
}

inline char* FormatScientific(char* buffer, intptr_t num_digits, int exponent, bool /*force_trailing_dot_zero*/)
{
    GRISU2_ASSERT(buffer != nullptr);
    GRISU2_ASSERT(num_digits >= 1);

    if (num_digits == 1)
    {
        // dE+123
        // GRISU2_ASSERT(buffer_capacity >= num_digits + 5);

        buffer += 1;
#if 0
        if (force_trailing_dot_zero)
        {
            *buffer++ = '.';
            *buffer++ = '0';
        }
#endif
    }
    else
    {
        // d.igitsE+123
        // GRISU2_ASSERT(buffer_capacity >= num_digits + 1 + 5);

        std::memmove(buffer + 2, buffer + 1, static_cast<size_t>(num_digits - 1));
        buffer[1] = '.';
        buffer += 1 + num_digits;
    }

    buffer[0] = 'e';
    buffer = ExponentToString(buffer + 1, exponent);

    return buffer;
}

} // namespace impl

// Maximum number of chars which will be written by PositiveDtoa.
constexpr int kPositiveDtoaMaxLength = 24;

// Generates a decimal representation of the floating-point number `value` in
// the buffer `[next, last)`.
//
// PRE: The input `value` must be strictly positive
// PRE: The buffer must be large enough (>= kPositiveDtoaMaxLength)
//
// Note: The result is _not_ null-terminated
template <typename Float>
inline char* PositiveDtoa(char* next, char* last, Float value, bool force_trailing_dot_zero = false)
{
    GRISU2_ASSERT(last - next >= kPositiveDtoaMaxLength);
    GRISU2_ASSERT(grisu2::impl::IEEE<Float>(value).IsFinite());
    GRISU2_ASSERT(value > 0);

    // Compute v = buffer * 10^exponent.
    // The decimal digits are stored in the buffer, which needs to be
    // interpreted as an unsigned decimal integer.
    // num_digits is the length of the buffer, i.e. the number of decimal digits.
    int num_digits = 0;
    int exponent = 0;
    grisu2::DoubleToDigits(next, last, num_digits, exponent, value);

    // Grisu2 generates at most max_digits10 decimal digits.
    GRISU2_ASSERT(num_digits <= std::numeric_limits<Float>::max_digits10);

    // The position of the decimal point relative to the start of the buffer.
    const int decimal_point = num_digits + exponent;

    // Just appending the exponent would yield a correct decimal representation
    // for the input value.

#if 1
    // Format the digits similar to printf's %g style.
    //
    // NB:
    // These are the values used by JavaScript's ToString applied to Number
    // type. Printf uses the values -4 and max_digits10 resp. (sort of).
    constexpr int kMinExp = -6;
    constexpr int kMaxExp = 21;

    const bool use_fixed = kMinExp < decimal_point && decimal_point <= kMaxExp;
#else
    // NB:
    // Integers <= 2^p = kMaxVal are exactly representable as Float's.
    constexpr auto kMinExp = -6;
    constexpr auto kMaxVal = static_cast<Float>(uint64_t{1} << std::numeric_limits<Float>::digits); // <= 16 digits

    const bool use_fixed = kMinExp < decimal_point && value <= kMaxVal;
#endif

    char* const end = use_fixed
        ? grisu2::impl::FormatFixed(next, num_digits, decimal_point, force_trailing_dot_zero)
        : grisu2::impl::FormatScientific(next, num_digits, decimal_point - 1, force_trailing_dot_zero);

    GRISU2_ASSERT(end - next <= kPositiveDtoaMaxLength);
    return end;
}

//==================================================================================================
// Dtoa
//==================================================================================================

namespace impl {

inline char* StrCopy(char* next, char* last, const char* source)
{
    static_cast<void>(last); // Fix warning

    GRISU2_ASSERT(source != nullptr);

    const auto len = std::strlen(source);
    GRISU2_ASSERT(next <= last);
    GRISU2_ASSERT(static_cast<size_t>(last - next) >= len);

    std::memcpy(next, source, len * sizeof(char));
    return next + len;
}

} // namespace impl

// Maximum number of chars which will be written by Dtoa.
constexpr int kDtoaMaxLength = 1/* minus-sign */ + kPositiveDtoaMaxLength;

// Generates a decimal representation of the floating-point number `value` in
// the buffer `[next, last)`.
//
// PRE: The buffer must be large enough.
// Note:
//   Max(kDtoaMaxLength, strlen(nan_string), 1 + strlen(inf_string))
// is sufficient.
//
// Note: The result is _not_ null-terminated.
template <typename Float>
inline char* Dtoa(
    char*       next,
    char*       last,
    Float       value,
    bool        force_trailing_dot_zero = false,
    const char* nan_string = "NaN",
    const char* inf_string = "Infinity")
{
    using Fp = grisu2::impl::IEEE<Float>;

    GRISU2_ASSERT(last - next >= kDtoaMaxLength);
    GRISU2_ASSERT(nan_string != nullptr);
    GRISU2_ASSERT(inf_string != nullptr);

    const Fp v(value);

    if (!v.IsFinite())
    {
        if (v.IsNaN())
            return grisu2::impl::StrCopy(next, last, nan_string);
        if (v.SignBit())
            *next++ = '-';
        return grisu2::impl::StrCopy(next, last, inf_string);
    }

    if (v.SignBit())
    {
        value = v.AbsValue();
        *next++ = '-';
    }

    if (v.IsZero())
    {
        *next++ = '0';
        if (force_trailing_dot_zero)
        {
            *next++ = '.';
            *next++ = '0';
        }
        return next;
    }

    return grisu2::PositiveDtoa(next, last, value, force_trailing_dot_zero);
}

} // namespace grisu2
#if GRISU2_UNNAMED_NAMESPACE
} // namespace
#endif

//char* Dtoa(char* next, char* last, double value)
//{
//    return grisu2::Dtoa(next, last, value);
//}

//char* Ftoa(char* next, char* last, float value)
//{
//    return grisu2::Dtoa(next, last, value);
//}
