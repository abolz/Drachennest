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
#include <climits>
#include <limits>
#include <type_traits> // conditional

#ifdef _MSC_VER
#include <intrin.h>
#endif

#ifndef STRTOD_ASSERT
#define STRTOD_ASSERT(X) assert(X)
#endif

#ifndef STRTOD_INLINE
#if _MSC_VER
#define STRTOD_INLINE __forceinline
#else
#define STRTOD_INLINE inline
#endif
#endif

#if STRTOD_UNNAMED_NAMESPACE
namespace {
#endif
namespace base_conv {

//==================================================================================================
// DecimalToDouble
//
// Derived from the double-conversion library:
// https://github.com/google/double-conversion
//
// The original license can be found at the end of this file.
//
// [1] Clinger, "How to read floating point numbers accurately",
//     PLDI '90 Proceedings of the ACM SIGPLAN 1990 conference on Programming language design and
//     implementation, Pages 92-101
//==================================================================================================

// Maximum number of significant digits in decimal representation.
//
// The longest possible double in decimal representation is (2^53 - 1) * 5^1074 / 10^1074,
// which has 767 digits.
// If we parse a number whose first digits are equal to a mean of 2 adjacent doubles (that
// could have up to 768 digits) the result must be rounded to the bigger one unless the tail
// consists of zeros, so we don't need to preserve all the digits.
constexpr int kMaxSignificantDigits = 767 + 1;

namespace strtod_impl {

STRTOD_INLINE constexpr int Min(int x, int y) { return y < x ? y : x; }
STRTOD_INLINE constexpr int Max(int x, int y) { return y < x ? x : y; }

STRTOD_INLINE bool IsDigit(char ch)
{
#if 0
    return static_cast<unsigned>(ch - '0') < 10;
#else
    return '0' <= ch && ch <= '9';
#endif
}

STRTOD_INLINE int DigitValue(char ch)
{
    STRTOD_ASSERT(IsDigit(ch));
    return ch - '0';
}

template <typename Dest, typename Source>
STRTOD_INLINE Dest ReinterpretBits(Source source)
{
    static_assert(sizeof(Dest) == sizeof(Source), "size mismatch");

    Dest dest;
    std::memcpy(&dest, &source, sizeof(Source));
    return dest;
}

template <typename Float>
struct IEEE
{
    // NB:
    // Works for double == long double.
    static_assert(std::numeric_limits<Float>::is_iec559 &&
                  ((std::numeric_limits<Float>::digits == 24 && std::numeric_limits<Float>::max_exponent == 128) ||
                   (std::numeric_limits<Float>::digits == 53 && std::numeric_limits<Float>::max_exponent == 1024)),
        "IEEE-754 single- or double-precision implementation required");

    using ieee_type = Float;
    using bits_type = typename std::conditional<std::numeric_limits<Float>::digits == 24, uint32_t, uint64_t>::type;

    static constexpr int       SignificandSize         = std::numeric_limits<ieee_type>::digits;  // = p   (includes the hidden bit)
    static constexpr int       PhysicalSignificandSize = SignificandSize - 1;                     // = p-1 (excludes the hidden bit)
    static constexpr int       UnbiasedMinExponent     = 1;
    static constexpr int       UnbiasedMaxExponent     = 2 * std::numeric_limits<Float>::max_exponent - 1 - 1;
    static constexpr int       ExponentBias            = 2 * std::numeric_limits<Float>::max_exponent / 2 - 1 + (SignificandSize - 1);
    static constexpr int       MinExponent             = UnbiasedMinExponent - ExponentBias;
    static constexpr int       MaxExponent             = UnbiasedMaxExponent - ExponentBias;
    static constexpr bits_type HiddenBit               = bits_type{1} << (SignificandSize - 1);   // = 2^(p-1)
    static constexpr bits_type SignificandMask         = HiddenBit - 1;                           // = 2^(p-1) - 1
    static constexpr bits_type ExponentMask            = bits_type{2 * std::numeric_limits<Float>::max_exponent - 1} << PhysicalSignificandSize;
    static constexpr bits_type SignMask                = ~(~bits_type{0} >> 1);

    bits_type bits;

    explicit IEEE(bits_type bits_) : bits(bits_) {}
    explicit IEEE(ieee_type value) : bits(ReinterpretBits<bits_type>(value)) {}

    bits_type PhysicalSignificand() const {
        return bits & SignificandMask;
    }

    bits_type PhysicalExponent() const {
        return (bits & ExponentMask) >> PhysicalSignificandSize;
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

    ieee_type Value() const {
        return ReinterpretBits<ieee_type>(bits);
    }

    ieee_type AbsValue() const {
        return ReinterpretBits<ieee_type>(bits & ~SignMask);
    }

    ieee_type NextValue() const {
        STRTOD_ASSERT(!SignBit());
        return ReinterpretBits<ieee_type>(IsInf() ? bits : bits + 1);
    }
};

//--------------------------------------------------------------------------------------------------
// StrtodFast
//--------------------------------------------------------------------------------------------------

// Double operations detection based on target architecture.
// Linux uses a 80bit wide floating point stack on x86. This induces double rounding, which in
// turn leads to wrong results.
// An easy way to test if the floating-point operations are correct is to evaluate: 89255.0/1e22.
// If the floating-point stack is 64 bits wide then the result is equal to 89255e-22.
// The best way to test this, is to create a division-function and to compare
// the output of the division with the expected result. (Inlining must be disabled.)
#if defined(_M_X64)              || \
    defined(__x86_64__)          || \
    defined(__ARMEL__)           || \
    defined(__avr32__)           || \
    defined(__hppa__)            || \
    defined(__ia64__)            || \
    defined(__mips__)            || \
    defined(__powerpc__)         || \
    defined(__ppc__)             || \
    defined(__ppc64__)           || \
    defined(_POWER)              || \
    defined(_ARCH_PPC)           || \
    defined(_ARCH_PPC64)         || \
    defined(__sparc__)           || \
    defined(__sparc)             || \
    defined(__s390__)            || \
    defined(__SH4__)             || \
    defined(__alpha__)           || \
    defined(_MIPS_ARCH_MIPS32R2) || \
    defined(__AARCH64EL__)       || \
    defined(__aarch64__)         || \
    defined(__riscv)
#define STRTOD_CORRECT_DOUBLE_OPERATIONS 1
#elif defined(_M_IX86) || defined(__i386__) || defined(__i386)
#ifdef _WIN32
// Windows uses a 64bit wide floating point stack.
#define STRTOD_CORRECT_DOUBLE_OPERATIONS 1
#endif
#endif

// 2^53 = 9007199254740992.
// Any integer with at most 15 decimal digits will hence fit into a double
// (which has a 53bit significand) without loss of precision.
constexpr int kMaxExactDoubleIntegerDecimalDigits = 15;

#if STRTOD_CORRECT_DOUBLE_OPERATIONS

STRTOD_INLINE bool FastPath(double& result, uint64_t digits, int num_digits, int exponent)
{
    static constexpr int kMaxExactPowerOfTen = 22;
    static constexpr double kExactPowersOfTen[] = {
        1.0e+00,
        1.0e+01,
        1.0e+02,
        1.0e+03,
        1.0e+04,
        1.0e+05,
        1.0e+06,
        1.0e+07,
        1.0e+08,
        1.0e+09,
        1.0e+10,
        1.0e+11,
        1.0e+12,
        1.0e+13,
        1.0e+14,
        1.0e+15, // 10^15 < 9007199254740992 = 2^53
        1.0e+16, // 10^16 = 5000000000000000 * 2^1  = (10^15 * 5^1 ) * 2^1
        1.0e+17, // 10^17 = 6250000000000000 * 2^4  = (10^13 * 5^4 ) * 2^4
        1.0e+18, // 10^18 = 7812500000000000 * 2^7  = (10^11 * 5^7 ) * 2^7
        1.0e+19, // 10^19 = 4882812500000000 * 2^11 = (10^8  * 5^11) * 2^11
        1.0e+20, // 10^20 = 6103515625000000 * 2^14 = (10^6  * 5^14) * 2^14
        1.0e+21, // 10^21 = 7629394531250000 * 2^17 = (10^4  * 5^17) * 2^17
        1.0e+22, // 10^22 = 4768371582031250 * 2^21 = (10^1  * 5^21) * 2^21
//      1.0e+23,
    };

    STRTOD_ASSERT(num_digits <= kMaxExactDoubleIntegerDecimalDigits);

    // The significand fits into a double.
    // If 10^exponent (resp. 10^-exponent) fits into a double too then we can
    // compute the result simply by multiplying (resp. dividing) the two
    // numbers.
    // This is possible because IEEE guarantees that floating-point operations
    // return the best possible approximation.

    int const remaining_digits = kMaxExactDoubleIntegerDecimalDigits - num_digits; // 0 <= rd <= 15
    if (-kMaxExactPowerOfTen <= exponent && exponent <= remaining_digits + kMaxExactPowerOfTen)
    {
        double d = static_cast<double>(static_cast<int64_t>(digits));
        if (exponent < 0)
        {
            d /= kExactPowersOfTen[-exponent];
        }
        else if (exponent <= kMaxExactPowerOfTen)
        {
            d *= kExactPowersOfTen[exponent];
        }
        else
        {
            // The buffer is short and we can multiply it with
            // 10^remaining_digits and the remaining exponent fits into a double.
            //
            // Eg. 123 * 10^25 = (123*1000) * 10^22

            d *= kExactPowersOfTen[remaining_digits]; // exact
            d *= kExactPowersOfTen[exponent - remaining_digits];
        }
        result = d;
        return true;
    }

    return false;
}

#else // ^^^ STRTOD_CORRECT_DOUBLE_OPERATIONS

STRTOD_INLINE bool FastPath(double& /*result*/, uint64_t /*digits*/, int /*num_digits*/, int /*exponent*/)
{
    return false;
}

#endif // ^^^ !STRTOD_CORRECT_DOUBLE_OPERATIONS

//--------------------------------------------------------------------------------------------------
// StrtodApprox
//--------------------------------------------------------------------------------------------------

struct DiyFp // f * 2^e
{
    static constexpr int SignificandSize = 64; // = q

    uint64_t f = 0;
    int e = 0;

    constexpr DiyFp() = default;
    constexpr DiyFp(uint64_t f_, int e_) : f(f_), e(e_) {}
};

// Returns whether the given floating point value is normalized.
STRTOD_INLINE bool IsNormalized(DiyFp x)
{
    static_assert(DiyFp::SignificandSize == 64, "internal error");

    return x.f >= (uint64_t{1} << 63);
}

// Returns x - y.
// PRE: x.e == y.e and x.f >= y.f
STRTOD_INLINE DiyFp Subtract(DiyFp x, DiyFp y)
{
    STRTOD_ASSERT(x.e == y.e);
    STRTOD_ASSERT(x.f >= y.f);

    return DiyFp(x.f - y.f, x.e);
}

// Returns x * y.
// The result is rounded (ties up). (Only the upper q bits are returned.)
STRTOD_INLINE DiyFp Multiply(DiyFp x, DiyFp y)
{
    static_assert(DiyFp::SignificandSize == 64, "internal error");

    // Computes:
    //  f = round((x.f * y.f) / 2^q)
    //  e = x.e + y.e + q

#if defined(_MSC_VER) && defined(_M_X64)

    uint64_t h = 0;
    uint64_t l = _umul128(x.f, y.f, &h);
    h += l >> 63; // round, ties up: [h, l] += 2^q / 2

    return DiyFp(h, x.e + y.e + 64);

#elif defined(__GNUC__) && defined(__SIZEOF_INT128__)

    __extension__ using Uint128 = unsigned __int128;

    Uint128 const p = Uint128{x.f} * Uint128{y.f};

    uint64_t h = static_cast<uint64_t>(p >> 64);
    uint64_t l = static_cast<uint64_t>(p);
    h += l >> 63; // round, ties up: [h, l] += 2^q / 2

    return DiyFp(h, x.e + y.e + 64);

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

    // Note:
    // The 32/64-bit casts here help MSVC to avoid calls to the _allmul
    // library function.

    uint32_t const u_lo = static_cast<uint32_t>(x.f /*& 0xFFFFFFFF*/);
    uint32_t const u_hi = static_cast<uint32_t>(x.f >> 32);
    uint32_t const v_lo = static_cast<uint32_t>(y.f /*& 0xFFFFFFFF*/);
    uint32_t const v_hi = static_cast<uint32_t>(y.f >> 32);

    uint64_t const p0 = uint64_t{u_lo} * v_lo;
    uint64_t const p1 = uint64_t{u_lo} * v_hi;
    uint64_t const p2 = uint64_t{u_hi} * v_lo;
    uint64_t const p3 = uint64_t{u_hi} * v_hi;

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

    return DiyFp(h, x.e + y.e + 64);

#endif
}

// Decomposes `value` into `f * 2^e`.
// The result is not normalized.
// PRE: `value` must be finite and non-negative, i.e. >= +0.0.
template <typename Float>
STRTOD_INLINE DiyFp DiyFpFromFloat(Float value)
{
    using Fp = IEEE<Float>;

    auto const v = Fp(value);

    STRTOD_ASSERT(v.IsFinite());
    STRTOD_ASSERT(!v.SignBit());

    auto const F = v.PhysicalSignificand();
    auto const E = v.PhysicalExponent();

    // If v is denormal:
    //      value = 0.F * 2^(1 - bias) = (          F) * 2^(1 - bias - (p-1))
    // If v is normalized:
    //      value = 1.F * 2^(E - bias) = (2^(p-1) + F) * 2^(E - bias - (p-1))

    return (E == 0) // denormal?
        ? DiyFp(F, Fp::MinExponent)
        : DiyFp(F + Fp::HiddenBit, static_cast<int>(E) - Fp::ExponentBias);
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

// Returns the upper boundary of value, i.e. the upper bound of the rounding
// interval for v.
// The result is not normalized.
// PRE: `value` must be finite and non-negative.
template <typename Float>
STRTOD_INLINE DiyFp UpperBoundary(Float value)
{
    auto const v = DiyFpFromFloat(value);
//  return DiyFp(2*v.f + 1, v.e - 1);
    return DiyFp(4*v.f + 2, v.e - 2);
}

// template <typename Float>
// STRTOD_INLINE bool LowerBoundaryIsCloser(Float value)
// {
//     IEEE<Float> const v(value);
//
//     STRTOD_ASSERT(v.IsFinite());
//     STRTOD_ASSERT(!v.SignBit());
//
//     auto const F = v.PhysicalSignificand();
//     auto const E = v.PhysicalExponent();
//     return F == 0 && E > 1;
// }

// Returns the lower boundary of `value`, i.e. the lower bound of the rounding
// interval for `value`.
// The result is not normalized.
// PRE: `value` must be finite and strictly positive.
// template <typename Float>
// STRTOD_INLINE DiyFp LowerBoundary(Float value)
// {
//     STRTOD_ASSERT(IEEE<Float>(value).IsFinite());
//     STRTOD_ASSERT(value > 0);
//
//     auto const v = DiyFpFromFloat(value);
//     return DiyFp(4*v.f - 2 + (LowerBoundaryIsCloser(value) ? 1 : 0), v.e - 2);
// }

struct DiyFpWithError // value = (x.f + delta) * 2^x.e, where |delta| <= error
{
    // We don't want to deal with fractions and therefore work with a common denominator.
    static constexpr int kDenominatorLog = 1;
    static constexpr int kDenominator = 1 << kDenominatorLog;

    DiyFp x;
    uint32_t error = 0;

    constexpr DiyFpWithError() = default;
    constexpr DiyFpWithError(DiyFp x_, uint32_t error_) : x(x_), error(error_) {}
};

// Returns the number of leading 0-bits in x, starting at the most significant bit position.
// If x is 0, the result is undefined.
STRTOD_INLINE int CountLeadingZeros64(uint64_t x)
{
    STRTOD_ASSERT(x != 0);

#if defined(_MSC_VER) && defined(_M_X64)

    return static_cast<int>(__lzcnt64(x));

#elif defined(_MSC_VER) && defined(_M_IX86)

    int lz = static_cast<int>( __lzcnt(static_cast<uint32_t>(x >> 32)) );
    if (lz == 32) {
        lz += static_cast<int>( __lzcnt(static_cast<uint32_t>(x)) );
    }
    return lz;

#elif defined(__GNUC__)

    return __builtin_clzll(x);

#else

    int lz = 0;
    while ((x >> 63) == 0) {
        x <<= 1;
        ++lz;
    }
    return lz;

#endif
}

// Normalize x
// and scale the error, so that the error is in ULP(x)
STRTOD_INLINE void Normalize(DiyFpWithError& num)
{
    int const s = CountLeadingZeros64(num.x.f);

    STRTOD_ASSERT(((num.error << s) >> s) == num.error);

    num.x.f   <<= s;
    num.x.e    -= s;
    num.error <<= s;
}

// 2^64 = 18446744073709551616 > 10^19
// Any integer with at most 19 decimal digits will hence fit into an uint64_t.
constexpr int kMaxUint64DecimalDigits = 19;

template <typename Int>
STRTOD_INLINE Int ReadInt(char const* f, char const* l)
{
    STRTOD_ASSERT(l - f <= std::numeric_limits<Int>::digits10);

    Int value = 0;
#if 1
    for ( ; l - f >= 8; f += 8)
    {
        value = 10 * value + static_cast<unsigned char>(f[0]);
        value = 10 * value + static_cast<unsigned char>(f[1]);
        value = 10 * value + static_cast<unsigned char>(f[2]);
        value = 10 * value + static_cast<unsigned char>(f[3]);
        value = 10 * value + static_cast<unsigned char>(f[4]);
        value = 10 * value + static_cast<unsigned char>(f[5]);
        value = 10 * value + static_cast<unsigned char>(f[6]);
        value = 10 * value + static_cast<unsigned char>(f[7]) - Int{533333328};
    }
#endif
    for ( ; f != l; ++f)
    {
#if 1
        STRTOD_ASSERT(IsDigit(*f));
        value = 10 * value + static_cast<unsigned char>(*f) - '0';
#else
        value = 10 * value + static_cast<uint32_t>(DigitValue(*f));
#endif
    }

    return value;
}

struct CachedPower { // c = f * 2^e ~= 10^k
    uint64_t f;
    int e; // binary exponent
    int k; // decimal exponent
};

constexpr int kCachedPowersSize         =   43;
constexpr int kCachedPowersMinDecExp    = -348;
constexpr int kCachedPowersMaxDecExp    =  324;
constexpr int kCachedPowersDecExpStep   =   16;

// Returns the binary exponent of a cached power for a given decimal exponent.
STRTOD_INLINE int BinaryExponentFromDecimalExponent(int k)
{
    STRTOD_ASSERT(k <=  400);
    STRTOD_ASSERT(k >= -400);

    // log_2(10) ~= [3; 3, 9, 2, 2, 4, 6, 2, 1, 1, 3] = 254370/76573
    // 2^15 * 254370/76573 = 108852.93980907...

//  return (k * 108853) / (1 << 15) - (k < 0) - 63;
    return (k * 108853 - 63 * (1 << 15)) >> 15;
}

STRTOD_INLINE CachedPower GetCachedPower(int index)
{
    static constexpr uint64_t kSignificands[/*340 bytes*/] = {
        0xFA8FD5A0081C0288, // * 2^-1220 <  10^-348
        0x8B16FB203055AC76, // * 2^-1166 <  10^-332
        0x9A6BB0AA55653B2D, // * 2^-1113 <  10^-316
        0xAB70FE17C79AC6CA, // * 2^-1060 <  10^-300
        0xBE5691EF416BD60C, // * 2^-1007 <  10^-284
        0xD3515C2831559A83, // * 2^-954  <  10^-268
        0xEA9C227723EE8BCB, // * 2^-901  <  10^-252
        0x823C12795DB6CE57, // * 2^-847  <  10^-236
        0x9096EA6F3848984F, // * 2^-794  <  10^-220
        0xA086CFCD97BF97F4, // * 2^-741  >  10^-204
        0xB23867FB2A35B28E, // * 2^-688  >  10^-188
        0xC5DD44271AD3CDBA, // * 2^-635  <  10^-172
        0xDBAC6C247D62A584, // * 2^-582  >  10^-156
        0xF3E2F893DEC3F126, // * 2^-529  <  10^-140
        0x87625F056C7C4A8B, // * 2^-475  <  10^-124
        0x964E858C91BA2655, // * 2^-422  <  10^-108
        0xA6DFBD9FB8E5B88F, // * 2^-369  >  10^-92
        0xB94470938FA89BCF, // * 2^-316  >  10^-76
        0xCDB02555653131B6, // * 2^-263  <  10^-60
        0xE45C10C42A2B3B06, // * 2^-210  >  10^-44
        0xFD87B5F28300CA0E, // * 2^-157  >  10^-28
        0x8CBCCC096F5088CC, // * 2^-103  >  10^-12
        0x9C40000000000000, // * 2^-50   == 10^4
        0xAD78EBC5AC620000, // * 2^3     == 10^20
        0xC097CE7BC90715B3, // * 2^56    <  10^36
        0xD5D238A4ABE98068, // * 2^109   <  10^52
        0xED63A231D4C4FB27, // * 2^162   <  10^68
        0x83C7088E1AAB65DB, // * 2^216   <  10^84
        0x924D692CA61BE758, // * 2^269   <  10^100
        0xA26DA3999AEF774A, // * 2^322   >  10^116
        0xB454E4A179DD1877, // * 2^375   <  10^132
        0xC83553C5C8965D3D, // * 2^428   <  10^148
        0xDE469FBD99A05FE3, // * 2^481   <  10^164
        0xF6C69A72A3989F5C, // * 2^534   >  10^180
        0x88FCF317F22241E2, // * 2^588   <  10^196
        0x98165AF37B2153DF, // * 2^641   >  10^212
        0xA8D9D1535CE3B396, // * 2^694   <  10^228
        0xBB764C4CA7A44410, // * 2^747   >  10^244
        0xD01FEF10A657842C, // * 2^800   <  10^260
        0xE7109BFBA19C0C9D, // * 2^853   <  10^276
        0x80444B5E7AA7CF85, // * 2^907   <  10^292
        0x8E679C2F5E44FF8F, // * 2^960   <  10^308
        0x9E19DB92B4E31BA9, // * 2^1013  <  10^324
    };

    STRTOD_ASSERT(index >= 0);
    STRTOD_ASSERT(index < kCachedPowersSize);

    int const k = kCachedPowersMinDecExp + index * kCachedPowersDecExpStep;
    int const e = BinaryExponentFromDecimalExponent(k);

    return {kSignificands[index], e, k};
}

// Returns a cached power of ten x ~= 10^k such that
//  k <= e < k + kCachedPowersDecExpStep.
//
// PRE: e >= kCachedPowersMinDecExp
// PRE: e <  kCachedPowersMaxDecExp + kCachedPowersDecExpStep
STRTOD_INLINE CachedPower GetCachedPowerForDecimalExponent(int e)
{
    STRTOD_ASSERT(e >= kCachedPowersMinDecExp);
    STRTOD_ASSERT(e <  kCachedPowersMaxDecExp + kCachedPowersDecExpStep);

    int const index = static_cast<int>( static_cast<unsigned>(-kCachedPowersMinDecExp + e) / kCachedPowersDecExpStep );
    STRTOD_ASSERT(index >= 0);
    STRTOD_ASSERT(index < kCachedPowersSize);

    auto const cached = GetCachedPower(index);
    STRTOD_ASSERT(e >= cached.k);
    STRTOD_ASSERT(e <  cached.k + kCachedPowersDecExpStep);

    return cached;
}

// Returns 10^k as an exact DiyFp.
// PRE: 1 <= k < kCachedPowersDecExpStep
STRTOD_INLINE DiyFp GetAdjustmentPowerOfTen(int k)
{
    static_assert(kCachedPowersDecExpStep <= 16, "internal error");

    static constexpr uint64_t kSignificands[] = {
        0x8000000000000000, // * 2^-63   == 10^0 (unused)
        0xA000000000000000, // * 2^-60   == 10^1
        0xC800000000000000, // * 2^-57   == 10^2
        0xFA00000000000000, // * 2^-54   == 10^3
        0x9C40000000000000, // * 2^-50   == 10^4
        0xC350000000000000, // * 2^-47   == 10^5
        0xF424000000000000, // * 2^-44   == 10^6
        0x9896800000000000, // * 2^-40   == 10^7
        0xBEBC200000000000, // * 2^-37   == 10^8
        0xEE6B280000000000, // * 2^-34   == 10^9
        0x9502F90000000000, // * 2^-30   == 10^10
        0xBA43B74000000000, // * 2^-27   == 10^11
        0xE8D4A51000000000, // * 2^-24   == 10^12
        0x9184E72A00000000, // * 2^-20   == 10^13
        0xB5E620F480000000, // * 2^-17   == 10^14
        0xE35FA931A0000000, // * 2^-14   == 10^15
    };

    STRTOD_ASSERT(k > 0);
    STRTOD_ASSERT(k < kCachedPowersDecExpStep);

    int const e = BinaryExponentFromDecimalExponent(k);
    return {kSignificands[k], e};
}

// Max double: 1.7976931348623157 * 10^308, which has 309 digits.
// Any x >= 10^309 is interpreted as +infinity.
constexpr int kMaxDecimalPower = 309;

// Min non-zero double: 4.9406564584124654 * 10^-324
// Any x <= 10^-324 is interpreted as 0.
// Note that 2.5e-324 (despite being smaller than the min double) will be read
// as non-zero (equal to the min non-zero double).
constexpr int kMinDecimalPower = -324;

// Returns the significand size for a given order of magnitude.
//
// If v = f * 2^e with 2^(q-1) <= f < 2^q then (q+e) is v's order of magnitude.
// If v = s * 2^e with 1/2 <= s < 1 then e is v's order of magnitude.
//
// This function returns the number of significant binary digits v will have
// once it's encoded into a 'double'. In almost all cases this is equal to
// Double::SignificandSize. The only exceptions are subnormals. They start with
// leading zeroes and their effective significand-size is hence smaller.
STRTOD_INLINE int EffectiveSignificandSize(int order)
{
    using Double = IEEE<double>;

    int const s = order - Double::MinExponent;

    if (s > Double::SignificandSize)
        return Double::SignificandSize;
    if (s < 0)
        return 0;

    return s;
}

// Returns `f * 2^e`.
STRTOD_INLINE double LoadDouble(uint64_t f, int e)
{
    using Double = IEEE<double>;

    STRTOD_ASSERT(f <= Double::HiddenBit + Double::SignificandMask);
    STRTOD_ASSERT(e <= Double::MinExponent || (f & Double::HiddenBit) != 0);

    if (e > Double::MaxExponent)
    {
        return std::numeric_limits<double>::infinity();
    }
    if (e < Double::MinExponent)
    {
        return 0.0;
    }

    uint64_t const exponent = (e == Double::MinExponent && (f & Double::HiddenBit) == 0)
        ? 0 // subnormal
        : static_cast<uint64_t>(e + Double::ExponentBias);

    uint64_t const bits = (exponent << Double::PhysicalSignificandSize) | (f & Double::SignificandMask);

    return ReinterpretBits<double>(bits);
}

// Use DiyFp's to approximate digits * 10^exponent.
//
// If the function returns true then 'result' is the correct double.
// Otherwise 'result' is either the correct double or the double that is just
// below the correct double.
//
// PRE: num_digits + exponent <= kMaxDecimalPower
// PRE: num_digits + exponent >  kMinDecimalPower
STRTOD_INLINE bool StrtodApprox(double& result, char const* digits, int num_digits, int exponent)
{
    using Double = IEEE<double>;

    static_assert(DiyFp::SignificandSize == 64,
        "We use uint64's. This only works if the DiyFp uses uint64's too.");

    STRTOD_ASSERT(num_digits > 0);
    STRTOD_ASSERT(DigitValue(digits[0]) > 0);
//  STRTOD_ASSERT(DigitValue(digits[num_digits - 1]) > 0);
    STRTOD_ASSERT(num_digits + exponent <= kMaxDecimalPower);
    STRTOD_ASSERT(num_digits + exponent >  kMinDecimalPower);

    // Compute an approximation 'input' for B = digits * 10^exponent using DiyFp's.
    // And keep track of the error.
    //
    //                       <-- error -->
    //                               B = digits * 10^exponent
    //  ---------(-----------|-------+---)------------------------------------
    //                       x
    //                       ~= (f * 2^e) * 10^exponent

    constexpr int kULP = DiyFpWithError::kDenominator;

    int const read_digits = Min(num_digits, kMaxUint64DecimalDigits);

    DiyFpWithError input;

    input.x.f = ReadInt<uint64_t>(digits, digits + read_digits);
    input.x.e = 0;
    input.error = 0;

    if (num_digits <= kMaxExactDoubleIntegerDecimalDigits)
    {
        if (FastPath(result, input.x.f, num_digits, exponent))
            return true;
    }

    if (read_digits < num_digits)
    {
        // Round.
        input.x.f += (DigitValue(digits[read_digits]) >= 5);

        // The error is <= 1/2 ULP.
        input.error = kULP / 2;
    }

    // x = f * 2^0

    // Normalize x and scale the error, such that 'error' is in ULP(x).
    Normalize(input);

    // If the input is exact, error == 0.
    // If the input is inexact, we have read 19 digits, i.e., f >= 10^(19-1) > 2^59.
    // The scaling factor in the normalization step above therefore is <= 2^(63-59) = 2^4.
    STRTOD_ASSERT(input.error <= 16 * (kULP / 2));

    // Move the remaining decimals into the (decimal) exponent.
    exponent += num_digits - read_digits;

    // Let x and y be normalized floating-point numbers
    //
    //      x = f_x * 2^e_x,    2^(q-1) <= f_x < 2^q
    //      y = f_y * 2^e_y,    2^(q-1) <= f_y < 2^q
    //
    // Then
    //
    //      z = Multiply(x,y) = f_z * 2^e_z
    //
    // returns the floating-point number closest to the product x*y. The result z is
    // not neccessarily normalized, but the error is bounded by 1/2 ulp, i.e.,
    //
    //      |x*y - z| <= 1/2 ulp
    //
    // or
    //
    //      x*y = (f_z + eps_z) * 2^e_z,    |eps_z| <= 1/2, e_z = e_x + e_y + q.
    //
    // If x and y are approximations to real numbers X and Y, i.e.,
    //
    //      X = (f_x + eps_x) * 2^e_x,      |eps_x| <= err_x,
    //      Y = (f_y + eps_y) * 2^e_y,      |eps_y| <= err_y,
    //
    // then the error introduced by a multiplication Multiply(x,y) is (see [1])
    //
    //      |X*Y - z| <= 1/2 + err_x + err_y + (err_x * err_y - err_x - err_y) / 2^q
    //
    // And if err_x < 1 (or err_y < 1), then
    //
    //      |X*Y - z| <= 1/2 + (err_x + err_y)

    auto const cached = GetCachedPowerForDecimalExponent(exponent);
    auto const cached_power = DiyFp(cached.f, cached.e);

    // Not all powers-of-ten are cached.
    // If cached.k != exponent we need to multiply 'x' by the difference first.
    // This may introduce an additional error.

    if (cached.k != exponent)
    {
        auto const adjustment_exponent = exponent - cached.k;
        auto const adjustment_power = GetAdjustmentPowerOfTen(adjustment_exponent);

        STRTOD_ASSERT(IsNormalized(input.x));
        STRTOD_ASSERT(IsNormalized(adjustment_power));

        input.x = Multiply(input.x, adjustment_power);
        // x ~= digits * 10^adjustment_exponent

        // Adjust error.
        // The adjustment_power is exact (err_y = 0).
        // There is hence only an additional error of (at most) 1/2.

        if (num_digits + adjustment_exponent <= kMaxUint64DecimalDigits)
        {
            // x and adjustment_power are exact.
            // The product (digits * 10^adjustment_exponent) fits into an uint64_t.
            // x * adjustment_power is therefore exact, too, and there is no additional error.
        }
        else
        {
            input.error += kULP / 2;

            STRTOD_ASSERT(input.error <= 17 * (kULP / 2));
        }

        // The result of the multiplication might not be normalized.
        // Normalize 'x' again and scale the error.
        Normalize(input);

        // Since both factors are normalized, input.f >= 2^(q-2), and the scaling
        // factor in the normalization step above is bounded by 2^1.
        STRTOD_ASSERT(input.error <= 34 * (kULP / 2));
    }

    STRTOD_ASSERT(IsNormalized(input.x));
    STRTOD_ASSERT(IsNormalized(cached_power));

    input.x = Multiply(input.x, cached_power);
    // x ~= digits * 10^exponent

    // Adjust the error.
    // Since all cached powers have an error of less than 1/2 ulp, err_y = 1/2,
    // and the error is therefore less than 1/2 + (err_x + err_y).

#if 1
    input.error += static_cast<unsigned>(kULP / 2 + (0 <= exponent && exponent <= 27 ? 0 : kULP / 2));
#else
    input.error += kULP / 2 + kULP / 2;
#endif

    STRTOD_ASSERT(input.error <= 36 * (kULP / 2));

    // The result of the multiplication might not be normalized.
    // Normalize 'x' again and scale the error.
    Normalize(input);

    // Since both factors were normalized, the scaling factor in the
    // normalization step above is again bounded by 2^1.
    STRTOD_ASSERT(input.error <= 72 * (kULP / 2));

    // We now have an approximation x = f * 2^e ~= digits * 10^exponent.
    //
    //                       <-- error -->
    //                               B = digits * 10^exponent
    //  ---------(-----------|-------+---)------------------------------------
    //                       x
    //                       ~= digits * 10^exponent
    //
    // B = (x.f + delta) * 2^x.e, where |delta| <= error / kULP
    //
    // When converting f * 2^e, which has a q-bit significand, into an IEEE
    // double-precision number, we need to drop some 'excess_bits' bits of
    // precision.

    int const prec = EffectiveSignificandSize(DiyFp::SignificandSize + input.x.e);
    STRTOD_ASSERT(prec >= 0);
    STRTOD_ASSERT(prec <= 53);

    int const excess_bits = DiyFp::SignificandSize - prec;

    // n = excess_bits
    //
    // f = (f div 2^n) * 2^n + (f mod 2^n)
    //   = (p1       ) * 2^n + (p2       )
    //
    //                             f = p1 * 2^n + p2
    //   <--- p2 ------------------>
    //                 <-- error --+-- error -->
    // --|-------------(-----------+------|----)---------------------------|--
    //   p1 * 2^n                                                 (p1 + 1) * 2^n
    //   <------------- half ------------->
    //                  = 2^n / 2
    //
    // The correct double now is either p1 * 2^(e + n) or (p1 + 1) * 2^(e + n).
    // See [1], Theorem 11.
    //
    // In case p2 + error < half, we can safely round down. If p2 - error > half
    // we can safely round up. Otherwise, we are too inaccurate. In this case
    // we round down, so the returned double is either the correct double or the
    // double just below the correct double. In this case we return false, so
    // that the we can fall back to a more precise algorithm.

    STRTOD_ASSERT(excess_bits >= 11);
    STRTOD_ASSERT(excess_bits <= 64);

    uint64_t const p2 = (excess_bits < 64) ? (input.x.f & ((uint64_t{1} << excess_bits) - 1)) : input.x.f;
    uint64_t const half = uint64_t{1} << (excess_bits - 1);

    // Truncate the significand to p = q - n bits and move the discarded bits into the (binary) exponent.
    // (Right shift of >= bit-width is undefined.)
    input.x.f = (excess_bits < 64) ? (input.x.f >> excess_bits) : 0;
    input.x.e += excess_bits;

    // Split up error into high (integral) and low (fractional) parts,
    // since half * kULP might overflow.
    uint64_t const error_hi = input.error / kULP;
    uint64_t const error_lo = input.error % kULP;

    STRTOD_ASSERT(input.error > 0);
    STRTOD_ASSERT(half >= error_hi && half - error_hi <= UINT64_MAX / kULP && (half - error_hi) * kULP >= error_lo);
    STRTOD_ASSERT(half <= UINT64_MAX - error_hi);
    static_cast<void>(error_lo);

    // Note:
    // Since error is non-zero, we can safely use '<=' and '>=' in the comparisons below.

    bool success;
#if 1
    // p2 * U >= half * U + error
    // <=> p2 * U >= half * U + (error_hi * U + error_lo)
    // <=> p2 * U >= (half + error_hi) * U + error_lo
    // <=> p2 >= (half + error_hi) + error_lo / U
    if (p2 > half + error_hi)
#else
    if (p2 * kULP >= half * kULP + input.error)
#endif
    {
        // Round up.
        success = true;

        ++input.x.f;

        // Rounding up may overflow the p-bit significand.
        // But in this case the significand is 2^53 and we don't loose any
        // bits by normalizing 'input' (we just move a factor of 2 into the
        // binary exponent).
        if (input.x.f > Double::HiddenBit + Double::SignificandMask)
        {
            STRTOD_ASSERT(input.x.f == (Double::HiddenBit << 1));

            input.x.f >>= 1;
            input.x.e  += 1;
        }
    }
#if 1
    // p2 * U <= half * U - error
    // <=> half * U >= p2 * U + error
    // <=> half * U >= p2 * U + (error_hi * U + error_lo)
    // <=> half * U >= (p2 + error_hi) * U + error_lo
    // <=> half >= (p2 + error_hi) + error_lo / U
    else if (half > p2 + error_hi)
#else
    else if (p2 * kULP <= half * kULP - input.error)
#endif
    {
        // Round down.
        success = true;
    }
    else
    {
        // Too imprecise.
        // Round down and return false, so that we can fall back to a more
        // precise algorithm.
        success = false;
    }

    result = LoadDouble(input.x.f, input.x.e);
    return success;
}

STRTOD_INLINE bool ComputeGuess(double& result, char const* digits, int num_digits, int exponent)
{
    STRTOD_ASSERT(num_digits > 0);
    STRTOD_ASSERT(num_digits <= kMaxSignificantDigits);
    STRTOD_ASSERT(DigitValue(digits[0]) > 0);
//  STRTOD_ASSERT(DigitValue(digits[num_digits - 1]) > 0);

    // Any v >= 10^309 is interpreted as +Infinity.
    if (num_digits + exponent > kMaxDecimalPower)
    {
        // Overflow.
        result = std::numeric_limits<double>::infinity();
        return true;
    }

    // Any v <= 10^-324 is interpreted as 0.
    if (num_digits + exponent <= kMinDecimalPower)
    {
        // Underflow.
        result = 0;
        return true;
    }

    return StrtodApprox(result, digits, num_digits, exponent);
}

//--------------------------------------------------------------------------------------------------
// StrtodBignum
//--------------------------------------------------------------------------------------------------

struct DiyInt // bigits * 2^exponent
{
    static constexpr int MaxBits = 64 + 2536 /*log_2(5^(324 - 1 + 769))*/ + 32;
    static constexpr int BigitSize = 32;
    static constexpr int Capacity = (MaxBits + (BigitSize - 1)) / BigitSize;

    uint32_t bigits[Capacity]; // Significand stored in little-endian form.
    int      size = 0;
    int      exponent = 0;

    DiyInt() = default;
    DiyInt(DiyInt const&) = delete;             // (not needed here)
    DiyInt& operator=(DiyInt const&) = delete;  // (not needed here)
};

STRTOD_INLINE void AssignZero(DiyInt& x)
{
    x.size = 0;
    x.exponent = 0;
}

STRTOD_INLINE void AssignU32(DiyInt& x, uint32_t value)
{
    AssignZero(x);

    if (value == 0)
        return;

    x.bigits[0] = value;
    x.size = 1;
}

STRTOD_INLINE void AssignU64(DiyInt& x, uint64_t value)
{
    AssignZero(x);

    if (value == 0)
        return;

    x.bigits[0] = static_cast<uint32_t>(value);
    x.bigits[1] = static_cast<uint32_t>(value >> DiyInt::BigitSize);
    x.size = (x.bigits[1] == 0) ? 1 : 2;
}

// x := A * x + B
STRTOD_INLINE void MulAddU32(DiyInt& x, uint32_t A, uint32_t B = 0)
{
    STRTOD_ASSERT(B == 0 || x.exponent == 0);

    if (A == 1 && B == 0)
    {
        return;
    }
    if (A == 0 || x.size == 0)
    {
        AssignU32(x, B);
        return;
    }

    uint32_t carry = B;
    for (int i = 0; i < x.size; ++i)
    {
        uint64_t const p = uint64_t{x.bigits[i]} * A + carry;
        x.bigits[i]      = static_cast<uint32_t>(p);
        carry            = static_cast<uint32_t>(p >> DiyInt::BigitSize);
    }

    if (carry != 0)
    {
        STRTOD_ASSERT(x.size < DiyInt::Capacity);
        x.bigits[x.size++] = carry;
    }
}

STRTOD_INLINE void AssignDecimalDigits(DiyInt& x, char const* digits, int num_digits)
{
    static constexpr uint32_t kPow10[] = {
        1, // (unused)
        10,
        100,
        1000,
        10000,
        100000,
        1000000,
        10000000,
        100000000,
        1000000000, // 10^9
    };

    AssignZero(x);

    while (num_digits > 0)
    {
        int const n = Min(num_digits, 9);
        MulAddU32(x, kPow10[n], ReadInt<uint32_t>(digits, digits + n));
        digits     += n;
        num_digits -= n;
    }
}

STRTOD_INLINE void MulPow2(DiyInt& x, int exp) // aka left-shift
{
    STRTOD_ASSERT(exp >= 0);

    if (x.size == 0)
        return;
    if (exp == 0)
        return;

    int const bigit_shift = exp / DiyInt::BigitSize;
    int const bit_shift   = exp % DiyInt::BigitSize;

    if (bit_shift > 0)
    {
        uint32_t carry = 0;
        for (int i = 0; i < x.size; ++i)
        {
            uint32_t const h = x.bigits[i] >> (DiyInt::BigitSize - bit_shift);
            x.bigits[i]      = x.bigits[i] << bit_shift | carry;
            carry            = h;
        }

        if (carry != 0)
        {
            STRTOD_ASSERT(x.size < DiyInt::Capacity);
            x.bigits[x.size++] = carry;
        }
    }

    x.exponent += bigit_shift;
}

STRTOD_INLINE void MulPow5(DiyInt& x, int exp)
{
    static constexpr uint32_t kPow5[] = {
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

    if (x.size == 0)
        return;

    STRTOD_ASSERT(exp >= 0);
    if (exp == 0)
        return;

    while (exp > 0)
    {
        int const n = Min(exp, 13);
        MulAddU32(x, kPow5[n]);
        exp -= n;
    }
}

STRTOD_INLINE int Compare(DiyInt const& lhs, DiyInt const& rhs)
{
    int const e1 = lhs.exponent;
    int const e2 = rhs.exponent;
    int const n1 = lhs.size + e1;
    int const n2 = rhs.size + e2;

    if (n1 < n2) return -1;
    if (n1 > n2) return +1;

    for (int i = n1 - 1; i >= Min(e1, e2); --i)
    {
        uint32_t const b1 = (i - e1) >= 0 ? lhs.bigits[i - e1] : 0;
        uint32_t const b2 = (i - e2) >= 0 ? rhs.bigits[i - e2] : 0;

        if (b1 < b2) return -1;
        if (b1 > b2) return +1;
    }

    return 0;
}

// Compare digits * 10^exponent with v = f * 2^e.
//
// PRE: num_digits + exponent <= kMaxDecimalPower
// PRE: num_digits + exponent >  kMinDecimalPower
// PRE: num_digits            <= kMaxSignificantDigits
STRTOD_INLINE int CompareBufferWithDiyFp(char const* digits, int num_digits, int exponent, bool nonzero_tail, DiyFp v)
{
    STRTOD_ASSERT(num_digits > 0);
    STRTOD_ASSERT(num_digits + exponent <= kMaxDecimalPower);
    STRTOD_ASSERT(num_digits + exponent >  kMinDecimalPower);
    STRTOD_ASSERT(num_digits            <= kMaxSignificantDigits);

    DiyInt lhs;
    DiyInt rhs;

    AssignDecimalDigits(lhs, digits, num_digits);
    if (nonzero_tail)
    {
        MulAddU32(lhs, 10, 1);
        exponent--;
    }
    AssignU64(rhs, v.f);

    STRTOD_ASSERT(lhs.size <= (2555 + 31) / 32); // bits <= log_2(10^769) = 2555
    STRTOD_ASSERT(rhs.size <= (  64 + 31) / 32); // bits <= 64

    int lhs_exp5 = 0;
    int rhs_exp5 = 0;
    int lhs_exp2 = 0;
    int rhs_exp2 = 0;

    if (exponent >= 0)
    {
        lhs_exp5 += exponent;
        lhs_exp2 += exponent;
    }
    else
    {
        rhs_exp5 -= exponent;
        rhs_exp2 -= exponent;
    }

    if (v.e >= 0)
    {
        rhs_exp2 += v.e;
    }
    else
    {
        lhs_exp2 -= v.e;
    }

#if 1
    if (lhs_exp5 > 0 || rhs_exp5 > 0)
    {
        MulPow5((lhs_exp5 > 0) ? lhs : rhs, (lhs_exp5 > 0) ? lhs_exp5 : rhs_exp5);
    }
#else
    if (lhs_exp5 > 0) // rhs >= digits
    {
        MulPow5(lhs, lhs_exp5);

        // num_digits + exponent <= kMaxDecimalPower + 1
        STRTOD_ASSERT(lhs.size <= (1030 + 31) / 32);  // 1030 = log_2(10^(309 + 1)))
        STRTOD_ASSERT(rhs.size <= (  64 + 31) / 32);
    }
    else if (rhs_exp5 > 0)
    {
        MulPow5(rhs, rhs_exp5);

        // kMinDecimalPower + 1 <= num_digits + exponent <= kMaxDecimalPower + 1
        // rhs_exp5 = -exponent <= -kMinDecimalPower - 1 + num_digits = 324 - 1 + num_digits <= 324 - 1 + 769
        // rhs_exp5 = -exponent >= -kMaxDecimalPower - 1 + num_digits
        STRTOD_ASSERT(lhs.size <= (2555        + 31) / 32);
        STRTOD_ASSERT(rhs.size <= (  64 + 2536 + 31) / 32); // 2536 = log_2(5^(324 - 1 + 769)) ---- XXX: 2504
    }
#endif

#if 0
    // Cancel common factors of 2.
    int const min_exp2 = Min(lhs_exp2, rhs_exp2);
    lhs_exp2 -= min_exp2;
    rhs_exp2 -= min_exp2;

    MulPow2(lhs, lhs_exp2);
    MulPow2(rhs, rhs_exp2);
#else
#if 1
    int const diff_exp2 = lhs_exp2 - rhs_exp2;
    if (diff_exp2 != 0)
    {
        MulPow2((diff_exp2 > 0) ? lhs : rhs, (diff_exp2 > 0) ? diff_exp2 : -diff_exp2);
    }
#else
    if (lhs_exp2 > rhs_exp2)
    {
        MulPow2(lhs, lhs_exp2 - rhs_exp2);
    }
    else if (rhs_exp2 > lhs_exp2)
    {
        MulPow2(rhs, rhs_exp2 - lhs_exp2);
    }
#endif
#endif

    STRTOD_ASSERT(lhs.size <= (2555        + 32 + 31) / 32);
    STRTOD_ASSERT(rhs.size <= (  64 + 2536 + 32 + 31) / 32);

    return Compare(lhs, rhs);
}

//--------------------------------------------------------------------------------------------------
// DecimalToDouble
//--------------------------------------------------------------------------------------------------

// Returns whether the significand f of v = f * 2^e is even.
STRTOD_INLINE bool SignificandIsEven(double v)
{
    return (IEEE<double>(v).PhysicalSignificand() & 1) == 0;
}

// Returns the next larger double-precision value.
// If v is +Infinity returns v.
STRTOD_INLINE double NextFloat(double v)
{
    return IEEE<double>(v).NextValue();
}

// Convert the decimal representation 'digits * 10^exponent' into an IEEE
// double-precision number.
//
// PRE: digits must contain only ASCII characters in the range '0'...'9'.
// PRE: num_digits >= 0
// PRE: num_digits + exponent must not overflow.
inline double DecimalToDouble(char const* digits, int num_digits, int exponent, bool nonzero_tail)
{
    STRTOD_ASSERT(num_digits >= 0);
    STRTOD_ASSERT(exponent <= INT_MAX - num_digits);

    // Ignore leading zeros
    while (num_digits > 0 && digits[0] == '0')
    {
        digits++;
        num_digits--;
    }

    // Move trailing zeros into the exponent
    while (num_digits > 0 && digits[num_digits - 1] == '0')
    {
        num_digits--;
        exponent++;
    }

    if (num_digits > kMaxSignificantDigits)
    {
        STRTOD_ASSERT(DigitValue(digits[num_digits - 1]) > 0); // since trailing zeros have been trimmed above.

        nonzero_tail = true;

        // Discard insignificant digits.
        exponent += num_digits - kMaxSignificantDigits;
        num_digits = kMaxSignificantDigits;

#if 0
        // Move trailing zeros into the exponent
        while (num_digits > 0 && digits[num_digits - 1] == '0')
        {
            num_digits--;
            exponent++;
        }
#endif
    }

    if (num_digits == 0)
    {
        return 0;
    }

    double v;
    if (ComputeGuess(v, digits, num_digits, exponent))
    {
        return v;
    }

    // Now v is either the correct or the next-lower double (i.e. the correct double is v+).
    // Compare B = buffer * 10^exponent with v's upper boundary m+.
    //
    //     v             m+            v+
    //  ---+--------+----+-------------+---
    //              B

    int const cmp = CompareBufferWithDiyFp(digits, num_digits, exponent, nonzero_tail, UpperBoundary(v));
    if (cmp < 0 || (cmp == 0 && SignificandIsEven(v)))
    {
        return v;
    }
    return NextFloat(v);
}

} // namespace strtod_impl

// Convert the decimal representation 'digits * 10^exponent' into an IEEE
// double-precision number.
//
// PRE: digits must contain only ASCII characters in the range '0'...'9'.
// PRE: num_digits >= 0
// PRE: num_digits + exponent must not overflow.
inline double DecimalToDouble(char const* digits, int num_digits, int exponent, bool nonzero_tail = false)
{
    return base_conv::strtod_impl::DecimalToDouble(digits, num_digits, exponent, nonzero_tail);
}

//==================================================================================================
// Strtod
//==================================================================================================

enum class StrtodStatus {
    success,
    input_too_large,
    no_digits,
    syntax_error,
    // XXX: exponent_too_large,
    // XXX: overflow,
    // XXX: underflow,
};

inline StrtodStatus Strtod(double& result, char const*& next, char const* last)
{
    using base_conv::strtod_impl::IsDigit;
    using base_conv::strtod_impl::DigitValue;

    // Inputs larger than kMaxInt (currently) can not be handled.
    // To avoid overflow in integer arithmetic.
    constexpr int const kMaxInt = INT_MAX / 4;

    StrtodStatus status = StrtodStatus::success;
    double       value  = 0;
    char const*  curr   = next;

    char digits[kMaxSignificantDigits];
    int  num_digits   = 0;
    int  exponent     = 0;
    bool nonzero_tail = false;
    bool is_neg       = false;

    if (last - curr >= kMaxInt)
    {
        status = StrtodStatus::input_too_large;
        goto L_done;
    }

    if (curr == last)
    {
        status = StrtodStatus::no_digits;
        goto L_done;
    }

    is_neg = (*curr == '-');
    if (is_neg)
    {
        ++curr;
    }
    else if (/*allow_leading_plus &&*/ *curr == '+')
    {
        ++curr;
    }

    if (curr == last)
    {
        status = StrtodStatus::syntax_error;
        goto L_done;
    }

    if (*curr == '0')
    {
        ++curr;
        if (curr == last)
        {
            goto L_done;
        }
    }
    else if (IsDigit(*curr))
    {
        for (;;)
        {
            if (num_digits < kMaxSignificantDigits)
            {
                digits[num_digits++] = *curr;
            }
            else
            {
                ++exponent;
                nonzero_tail = nonzero_tail || *curr != '0';
            }
            ++curr;
            if (curr == last)
            {
                goto L_convert;
            }
            if (!IsDigit(*curr))
            {
                break;
            }
        }
    }
    else if (/*allow_leading_dot &&*/ *curr == '.')
    {
        // Do nothing.
        // Will be parsed again below.
    }
    else
    {
        //
        // TODO:
        // Parse NaN and Infinity here.
        //
        goto L_done;
    }

    if (*curr == '.')
    {
        ++curr;
        if (curr == last)
        {
            status = (num_digits > 0) ? StrtodStatus::success : StrtodStatus::syntax_error;
            goto L_convert;
        }

        if (num_digits == 0)
        {
            // Integer part consists of 0 (or is absent).
            // Significant digits start after leading zeros (if any).
            while (*curr == '0')
            {
                ++curr;
                if (curr == last)
                {
                    goto L_done;
                }

                // Move this 0 into the exponent.
                --exponent;
            }
        }

        // There is a fractional part.
        // We don't emit a '.', but adjust the exponent instead.
        while (IsDigit(*curr))
        {
            if (num_digits < kMaxSignificantDigits)
            {
                digits[num_digits++] = *curr;
                --exponent;
            }
            else
            {
                nonzero_tail = nonzero_tail || *curr != '0';
            }
            ++curr;
            if (curr == last)
            {
                goto L_convert;
            }
        }
    }

    // Parse exponential part.
    if (*curr == 'e' || *curr == 'E')
    {
        ++curr;
        if (curr == last)
        {
            status = StrtodStatus::syntax_error;
            goto L_done;
        }

        bool const exp_is_neg = (*curr == '-');

        if (exp_is_neg || *curr == '+')
        {
            ++curr;
            if (curr == last)
            {
                status = StrtodStatus::syntax_error;
                goto L_done;
            }
        }

        if (!IsDigit(*curr))
        {
            status = StrtodStatus::syntax_error;
            goto L_done;
        }

        int num = 0;
        for (;;)
        {
            int const digit = DigitValue(*curr);

//          if (num > kMaxInt / 10 || digit > kMaxInt - 10 * num)
            if (num > kMaxInt / 10 - 9)
            {
                //status = StrtodStatus::exponent_too_large;
                num = kMaxInt;
                break;
            }

            num = num * 10 + digit;
            ++curr;
            if (curr == last)
            {
                break;
            }
            if (!IsDigit(*curr))
            {
                break;
            }
        }

        // Skip the rest of the exponent (ignored).
        for ( ; curr != last && IsDigit(*curr); ++curr)
        {
        }

        exponent += exp_is_neg ? -num : num;
    }

L_convert:
    value = base_conv::strtod_impl::DecimalToDouble(digits, num_digits, exponent, nonzero_tail);

L_done:
    result = is_neg ? -value : value;
    next = curr;

    return status;
}

STRTOD_INLINE double Strtod(char const* first, char const* last)
{
    double d = 0.0;
    base_conv::Strtod(d, first, last);
    return d;
}

} // namespace base_conv
#if STRTOD_UNNAMED_NAMESPACE
} // namespace
#endif

/*
Copyright 2006-2011, the V8 project authors. All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.
    * Neither the name of Google Inc. nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
