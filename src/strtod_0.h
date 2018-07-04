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

#ifndef STRTOD_OPTIMIZE_SIZE
#define STRTOD_OPTIMIZE_SIZE 1
#endif

namespace base_conv {

namespace strtod_impl {

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

} // namespace strtod_impl

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

    Uint128 const p = Uint128{x.f} * y.f + (uint64_t{1} << 63);
    return DiyFp(static_cast<uint64_t>(p >> 64), x.e + y.e + 64);

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

#if STRTOD_OPTIMIZE_SIZE
// sizeof(tables) = 340 + 128 = 468 bytes

constexpr int kCachedPowersSize         =   43;
constexpr int kCachedPowersMinDecExp    = -348;
constexpr int kCachedPowersMaxDecExp    =  324;
constexpr int kCachedPowersDecExpStep   =   16; // XXX: 27?!

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

#else // ^^^ STRTOD_OPTIMIZE_SIZE
// sizeof(tables) = 5384 bytes

constexpr int kCachedPowersSize         =  673;
constexpr int kCachedPowersMinDecExp    = -348;
constexpr int kCachedPowersMaxDecExp    =  324;
constexpr int kCachedPowersDecExpStep   =    1;

STRTOD_INLINE CachedPower GetCachedPower(int index)
{
    static constexpr uint64_t kSignificands[/*5384 bytes*/] = {
        0xFA8FD5A0081C0288, // * 2^-1220 <  10^-348
        0x9C99E58405118195, // * 2^-1216 <  10^-347
        0xC3C05EE50655E1FA, // * 2^-1213 <  10^-346
        0xF4B0769E47EB5A79, // * 2^-1210 >  10^-345
        0x98EE4A22ECF3188C, // * 2^-1206 >  10^-344
        0xBF29DCABA82FDEAE, // * 2^-1203 <  10^-343
        0xEEF453D6923BD65A, // * 2^-1200 <  10^-342
        0x9558B4661B6565F8, // * 2^-1196 <  10^-341
        0xBAAEE17FA23EBF76, // * 2^-1193 <  10^-340
        0xE95A99DF8ACE6F54, // * 2^-1190 >  10^-339
        0x91D8A02BB6C10594, // * 2^-1186 <  10^-338
        0xB64EC836A47146FA, // * 2^-1183 >  10^-337
        0xE3E27A444D8D98B8, // * 2^-1180 >  10^-336
        0x8E6D8C6AB0787F73, // * 2^-1176 >  10^-335
        0xB208EF855C969F50, // * 2^-1173 >  10^-334
        0xDE8B2B66B3BC4724, // * 2^-1170 >  10^-333
        0x8B16FB203055AC76, // * 2^-1166 <  10^-332
        0xADDCB9E83C6B1794, // * 2^-1163 >  10^-331
        0xD953E8624B85DD79, // * 2^-1160 >  10^-330
        0x87D4713D6F33AA6C, // * 2^-1156 >  10^-329
        0xA9C98D8CCB009506, // * 2^-1153 <  10^-328
        0xD43BF0EFFDC0BA48, // * 2^-1150 <  10^-327
        0x84A57695FE98746D, // * 2^-1146 <  10^-326
        0xA5CED43B7E3E9188, // * 2^-1143 <  10^-325
        0xCF42894A5DCE35EA, // * 2^-1140 <  10^-324
        0x818995CE7AA0E1B2, // * 2^-1136 <  10^-323
        0xA1EBFB4219491A1F, // * 2^-1133 <  10^-322
        0xCA66FA129F9B60A7, // * 2^-1130 >  10^-321
        0xFD00B897478238D1, // * 2^-1127 >  10^-320
        0x9E20735E8CB16382, // * 2^-1123 <  10^-319
        0xC5A890362FDDBC63, // * 2^-1120 >  10^-318
        0xF712B443BBD52B7C, // * 2^-1117 >  10^-317
        0x9A6BB0AA55653B2D, // * 2^-1113 <  10^-316
        0xC1069CD4EABE89F9, // * 2^-1110 >  10^-315
        0xF148440A256E2C77, // * 2^-1107 >  10^-314
        0x96CD2A865764DBCA, // * 2^-1103 <  10^-313
        0xBC807527ED3E12BD, // * 2^-1100 >  10^-312
        0xEBA09271E88D976C, // * 2^-1097 >  10^-311
        0x93445B8731587EA3, // * 2^-1093 <  10^-310
        0xB8157268FDAE9E4C, // * 2^-1090 <  10^-309
        0xE61ACF033D1A45DF, // * 2^-1087 <  10^-308
        0x8FD0C16206306BAC, // * 2^-1083 >  10^-307
        0xB3C4F1BA87BC8697, // * 2^-1080 >  10^-306
        0xE0B62E2929ABA83C, // * 2^-1077 <  10^-305
        0x8C71DCD9BA0B4926, // * 2^-1073 >  10^-304
        0xAF8E5410288E1B6F, // * 2^-1070 <  10^-303
        0xDB71E91432B1A24B, // * 2^-1067 >  10^-302
        0x892731AC9FAF056F, // * 2^-1063 >  10^-301
        0xAB70FE17C79AC6CA, // * 2^-1060 <  10^-300
        0xD64D3D9DB981787D, // * 2^-1057 <  10^-299
        0x85F0468293F0EB4E, // * 2^-1053 <  10^-298
        0xA76C582338ED2622, // * 2^-1050 >  10^-297
        0xD1476E2C07286FAA, // * 2^-1047 <  10^-296
        0x82CCA4DB847945CA, // * 2^-1043 <  10^-295
        0xA37FCE126597973D, // * 2^-1040 >  10^-294
        0xCC5FC196FEFD7D0C, // * 2^-1037 <  10^-293
        0xFF77B1FCBEBCDC4F, // * 2^-1034 <  10^-292
        0x9FAACF3DF73609B1, // * 2^-1030 <  10^-291
        0xC795830D75038C1E, // * 2^-1027 >  10^-290
        0xF97AE3D0D2446F25, // * 2^-1024 <  10^-289
        0x9BECCE62836AC577, // * 2^-1020 <  10^-288
        0xC2E801FB244576D5, // * 2^-1017 <  10^-287
        0xF3A20279ED56D48A, // * 2^-1014 <  10^-286
        0x9845418C345644D7, // * 2^-1010 >  10^-285
        0xBE5691EF416BD60C, // * 2^-1007 <  10^-284
        0xEDEC366B11C6CB8F, // * 2^-1004 <  10^-283
        0x94B3A202EB1C3F39, // * 2^-1000 <  10^-282
        0xB9E08A83A5E34F08, // * 2^-997  >  10^-281
        0xE858AD248F5C22CA, // * 2^-994  >  10^-280
        0x91376C36D99995BE, // * 2^-990  <  10^-279
        0xB58547448FFFFB2E, // * 2^-987  >  10^-278
        0xE2E69915B3FFF9F9, // * 2^-984  <  10^-277
        0x8DD01FAD907FFC3C, // * 2^-980  >  10^-276
        0xB1442798F49FFB4B, // * 2^-977  >  10^-275
        0xDD95317F31C7FA1D, // * 2^-974  <  10^-274
        0x8A7D3EEF7F1CFC52, // * 2^-970  <  10^-273
        0xAD1C8EAB5EE43B67, // * 2^-967  >  10^-272
        0xD863B256369D4A41, // * 2^-964  >  10^-271
        0x873E4F75E2224E68, // * 2^-960  <  10^-270
        0xA90DE3535AAAE202, // * 2^-957  <  10^-269
        0xD3515C2831559A83, // * 2^-954  <  10^-268
        0x8412D9991ED58092, // * 2^-950  >  10^-267
        0xA5178FFF668AE0B6, // * 2^-947  <  10^-266
        0xCE5D73FF402D98E4, // * 2^-944  >  10^-265
        0x80FA687F881C7F8E, // * 2^-940  <  10^-264
        0xA139029F6A239F72, // * 2^-937  <  10^-263
        0xC987434744AC874F, // * 2^-934  >  10^-262
        0xFBE9141915D7A922, // * 2^-931  <  10^-261
        0x9D71AC8FADA6C9B5, // * 2^-927  <  10^-260
        0xC4CE17B399107C23, // * 2^-924  >  10^-259
        0xF6019DA07F549B2B, // * 2^-921  <  10^-258
        0x99C102844F94E0FB, // * 2^-917  <  10^-257
        0xC0314325637A193A, // * 2^-914  >  10^-256
        0xF03D93EEBC589F88, // * 2^-911  <  10^-255
        0x96267C7535B763B5, // * 2^-907  <  10^-254
        0xBBB01B9283253CA3, // * 2^-904  >  10^-253
        0xEA9C227723EE8BCB, // * 2^-901  <  10^-252
        0x92A1958A7675175F, // * 2^-897  <  10^-251
        0xB749FAED14125D37, // * 2^-894  >  10^-250
        0xE51C79A85916F485, // * 2^-891  >  10^-249
        0x8F31CC0937AE58D3, // * 2^-887  >  10^-248
        0xB2FE3F0B8599EF08, // * 2^-884  >  10^-247
        0xDFBDCECE67006AC9, // * 2^-881  <  10^-246
        0x8BD6A141006042BE, // * 2^-877  >  10^-245
        0xAECC49914078536D, // * 2^-874  <  10^-244
        0xDA7F5BF590966849, // * 2^-871  >  10^-243
        0x888F99797A5E012D, // * 2^-867  <  10^-242
        0xAAB37FD7D8F58179, // * 2^-864  >  10^-241
        0xD5605FCDCF32E1D7, // * 2^-861  >  10^-240
        0x855C3BE0A17FCD26, // * 2^-857  <  10^-239
        0xA6B34AD8C9DFC070, // * 2^-854  >  10^-238
        0xD0601D8EFC57B08C, // * 2^-851  >  10^-237
        0x823C12795DB6CE57, // * 2^-847  <  10^-236
        0xA2CB1717B52481ED, // * 2^-844  <  10^-235
        0xCB7DDCDDA26DA269, // * 2^-841  >  10^-234
        0xFE5D54150B090B03, // * 2^-838  >  10^-233
        0x9EFA548D26E5A6E2, // * 2^-834  >  10^-232
        0xC6B8E9B0709F109A, // * 2^-831  <  10^-231
        0xF867241C8CC6D4C1, // * 2^-828  >  10^-230
        0x9B407691D7FC44F8, // * 2^-824  <  10^-229
        0xC21094364DFB5637, // * 2^-821  >  10^-228
        0xF294B943E17A2BC4, // * 2^-818  <  10^-227
        0x979CF3CA6CEC5B5B, // * 2^-814  >  10^-226
        0xBD8430BD08277231, // * 2^-811  <  10^-225
        0xECE53CEC4A314EBE, // * 2^-808  >  10^-224
        0x940F4613AE5ED137, // * 2^-804  >  10^-223
        0xB913179899F68584, // * 2^-801  <  10^-222
        0xE757DD7EC07426E5, // * 2^-798  <  10^-221
        0x9096EA6F3848984F, // * 2^-794  <  10^-220
        0xB4BCA50B065ABE63, // * 2^-791  <  10^-219
        0xE1EBCE4DC7F16DFC, // * 2^-788  >  10^-218
        0x8D3360F09CF6E4BD, // * 2^-784  <  10^-217
        0xB080392CC4349DED, // * 2^-781  >  10^-216
        0xDCA04777F541C568, // * 2^-778  >  10^-215
        0x89E42CAAF9491B61, // * 2^-774  >  10^-214
        0xAC5D37D5B79B6239, // * 2^-771  <  10^-213
        0xD77485CB25823AC7, // * 2^-768  <  10^-212
        0x86A8D39EF77164BD, // * 2^-764  >  10^-211
        0xA8530886B54DBDEC, // * 2^-761  >  10^-210
        0xD267CAA862A12D67, // * 2^-758  >  10^-209
        0x8380DEA93DA4BC60, // * 2^-754  <  10^-208
        0xA46116538D0DEB78, // * 2^-751  <  10^-207
        0xCD795BE870516656, // * 2^-748  <  10^-206
        0x806BD9714632DFF6, // * 2^-744  <  10^-205
        0xA086CFCD97BF97F4, // * 2^-741  >  10^-204
        0xC8A883C0FDAF7DF0, // * 2^-738  <  10^-203
        0xFAD2A4B13D1B5D6C, // * 2^-735  <  10^-202
        0x9CC3A6EEC6311A64, // * 2^-731  >  10^-201
        0xC3F490AA77BD60FD, // * 2^-728  >  10^-200
        0xF4F1B4D515ACB93C, // * 2^-725  >  10^-199
        0x991711052D8BF3C5, // * 2^-721  <  10^-198
        0xBF5CD54678EEF0B7, // * 2^-718  >  10^-197
        0xEF340A98172AACE5, // * 2^-715  >  10^-196
        0x9580869F0E7AAC0F, // * 2^-711  >  10^-195
        0xBAE0A846D2195713, // * 2^-708  >  10^-194
        0xE998D258869FACD7, // * 2^-705  <  10^-193
        0x91FF83775423CC06, // * 2^-701  <  10^-192
        0xB67F6455292CBF08, // * 2^-698  <  10^-191
        0xE41F3D6A7377EECA, // * 2^-695  <  10^-190
        0x8E938662882AF53E, // * 2^-691  <  10^-189
        0xB23867FB2A35B28E, // * 2^-688  >  10^-188
        0xDEC681F9F4C31F31, // * 2^-685  <  10^-187
        0x8B3C113C38F9F37F, // * 2^-681  >  10^-186
        0xAE0B158B4738705F, // * 2^-678  >  10^-185
        0xD98DDAEE19068C76, // * 2^-675  <  10^-184
        0x87F8A8D4CFA417CA, // * 2^-671  >  10^-183
        0xA9F6D30A038D1DBC, // * 2^-668  <  10^-182
        0xD47487CC8470652B, // * 2^-665  <  10^-181
        0x84C8D4DFD2C63F3B, // * 2^-661  <  10^-180
        0xA5FB0A17C777CF0A, // * 2^-658  >  10^-179
        0xCF79CC9DB955C2CC, // * 2^-655  <  10^-178
        0x81AC1FE293D599C0, // * 2^-651  >  10^-177
        0xA21727DB38CB0030, // * 2^-648  >  10^-176
        0xCA9CF1D206FDC03C, // * 2^-645  >  10^-175
        0xFD442E4688BD304B, // * 2^-642  >  10^-174
        0x9E4A9CEC15763E2F, // * 2^-638  >  10^-173
        0xC5DD44271AD3CDBA, // * 2^-635  <  10^-172
        0xF7549530E188C129, // * 2^-632  >  10^-171
        0x9A94DD3E8CF578BA, // * 2^-628  >  10^-170
        0xC13A148E3032D6E8, // * 2^-625  >  10^-169
        0xF18899B1BC3F8CA2, // * 2^-622  >  10^-168
        0x96F5600F15A7B7E5, // * 2^-618  <  10^-167
        0xBCB2B812DB11A5DE, // * 2^-615  <  10^-166
        0xEBDF661791D60F56, // * 2^-612  <  10^-165
        0x936B9FCEBB25C996, // * 2^-608  >  10^-164
        0xB84687C269EF3BFB, // * 2^-605  <  10^-163
        0xE65829B3046B0AFA, // * 2^-602  <  10^-162
        0x8FF71A0FE2C2E6DC, // * 2^-598  <  10^-161
        0xB3F4E093DB73A093, // * 2^-595  <  10^-160
        0xE0F218B8D25088B8, // * 2^-592  <  10^-159
        0x8C974F7383725573, // * 2^-588  <  10^-158
        0xAFBD2350644EEAD0, // * 2^-585  >  10^-157
        0xDBAC6C247D62A584, // * 2^-582  >  10^-156
        0x894BC396CE5DA772, // * 2^-578  <  10^-155
        0xAB9EB47C81F5114F, // * 2^-575  <  10^-154
        0xD686619BA27255A3, // * 2^-572  >  10^-153
        0x8613FD0145877586, // * 2^-568  >  10^-152
        0xA798FC4196E952E7, // * 2^-565  <  10^-151
        0xD17F3B51FCA3A7A1, // * 2^-562  >  10^-150
        0x82EF85133DE648C5, // * 2^-558  >  10^-149
        0xA3AB66580D5FDAF6, // * 2^-555  >  10^-148
        0xCC963FEE10B7D1B3, // * 2^-552  <  10^-147
        0xFFBBCFE994E5C620, // * 2^-549  >  10^-146
        0x9FD561F1FD0F9BD4, // * 2^-545  >  10^-145
        0xC7CABA6E7C5382C9, // * 2^-542  >  10^-144
        0xF9BD690A1B68637B, // * 2^-539  <  10^-143
        0x9C1661A651213E2D, // * 2^-535  <  10^-142
        0xC31BFA0FE5698DB8, // * 2^-532  <  10^-141
        0xF3E2F893DEC3F126, // * 2^-529  <  10^-140
        0x986DDB5C6B3A76B8, // * 2^-525  >  10^-139
        0xBE89523386091466, // * 2^-522  >  10^-138
        0xEE2BA6C0678B597F, // * 2^-519  <  10^-137
        0x94DB483840B717F0, // * 2^-515  >  10^-136
        0xBA121A4650E4DDEC, // * 2^-512  >  10^-135
        0xE896A0D7E51E1566, // * 2^-509  <  10^-134
        0x915E2486EF32CD60, // * 2^-505  <  10^-133
        0xB5B5ADA8AAFF80B8, // * 2^-502  <  10^-132
        0xE3231912D5BF60E6, // * 2^-499  <  10^-131
        0x8DF5EFABC5979C90, // * 2^-495  >  10^-130
        0xB1736B96B6FD83B4, // * 2^-492  >  10^-129
        0xDDD0467C64BCE4A1, // * 2^-489  >  10^-128
        0x8AA22C0DBEF60EE4, // * 2^-485  <  10^-127
        0xAD4AB7112EB3929E, // * 2^-482  >  10^-126
        0xD89D64D57A607745, // * 2^-479  >  10^-125
        0x87625F056C7C4A8B, // * 2^-475  <  10^-124
        0xA93AF6C6C79B5D2E, // * 2^-472  >  10^-123
        0xD389B47879823479, // * 2^-469  <  10^-122
        0x843610CB4BF160CC, // * 2^-465  >  10^-121
        0xA54394FE1EEDB8FF, // * 2^-462  >  10^-120
        0xCE947A3DA6A9273E, // * 2^-459  <  10^-119
        0x811CCC668829B887, // * 2^-455  <  10^-118
        0xA163FF802A3426A9, // * 2^-452  >  10^-117
        0xC9BCFF6034C13053, // * 2^-449  >  10^-116
        0xFC2C3F3841F17C68, // * 2^-446  >  10^-115
        0x9D9BA7832936EDC1, // * 2^-442  >  10^-114
        0xC5029163F384A931, // * 2^-439  <  10^-113
        0xF64335BCF065D37D, // * 2^-436  <  10^-112
        0x99EA0196163FA42E, // * 2^-432  <  10^-111
        0xC06481FB9BCF8D3A, // * 2^-429  >  10^-110
        0xF07DA27A82C37088, // * 2^-426  <  10^-109
        0x964E858C91BA2655, // * 2^-422  <  10^-108
        0xBBE226EFB628AFEB, // * 2^-419  >  10^-107
        0xEADAB0ABA3B2DBE5, // * 2^-416  <  10^-106
        0x92C8AE6B464FC96F, // * 2^-412  <  10^-105
        0xB77ADA0617E3BBCB, // * 2^-409  <  10^-104
        0xE55990879DDCAABE, // * 2^-406  >  10^-103
        0x8F57FA54C2A9EAB7, // * 2^-402  >  10^-102
        0xB32DF8E9F3546564, // * 2^-399  <  10^-101
        0xDFF9772470297EBD, // * 2^-396  <  10^-100
        0x8BFBEA76C619EF36, // * 2^-392  <  10^-99
        0xAEFAE51477A06B04, // * 2^-389  >  10^-98
        0xDAB99E59958885C5, // * 2^-386  >  10^-97
        0x88B402F7FD75539B, // * 2^-382  <  10^-96
        0xAAE103B5FCD2A882, // * 2^-379  >  10^-95
        0xD59944A37C0752A2, // * 2^-376  <  10^-94
        0x857FCAE62D8493A5, // * 2^-372  <  10^-93
        0xA6DFBD9FB8E5B88F, // * 2^-369  >  10^-92
        0xD097AD07A71F26B2, // * 2^-366  <  10^-91
        0x825ECC24C8737830, // * 2^-362  >  10^-90
        0xA2F67F2DFA90563B, // * 2^-359  <  10^-89
        0xCBB41EF979346BCA, // * 2^-356  <  10^-88
        0xFEA126B7D78186BD, // * 2^-353  >  10^-87
        0x9F24B832E6B0F436, // * 2^-349  <  10^-86
        0xC6EDE63FA05D3144, // * 2^-346  >  10^-85
        0xF8A95FCF88747D94, // * 2^-343  <  10^-84
        0x9B69DBE1B548CE7D, // * 2^-339  >  10^-83
        0xC24452DA229B021C, // * 2^-336  >  10^-82
        0xF2D56790AB41C2A3, // * 2^-333  >  10^-81
        0x97C560BA6B0919A6, // * 2^-329  >  10^-80
        0xBDB6B8E905CB600F, // * 2^-326  <  10^-79
        0xED246723473E3813, // * 2^-323  <  10^-78
        0x9436C0760C86E30C, // * 2^-319  >  10^-77
        0xB94470938FA89BCF, // * 2^-316  >  10^-76
        0xE7958CB87392C2C3, // * 2^-313  >  10^-75
        0x90BD77F3483BB9BA, // * 2^-309  >  10^-74
        0xB4ECD5F01A4AA828, // * 2^-306  <  10^-73
        0xE2280B6C20DD5232, // * 2^-303  <  10^-72
        0x8D590723948A535F, // * 2^-299  <  10^-71
        0xB0AF48EC79ACE837, // * 2^-296  <  10^-70
        0xDCDB1B2798182245, // * 2^-293  >  10^-69
        0x8A08F0F8BF0F156B, // * 2^-289  <  10^-68
        0xAC8B2D36EED2DAC6, // * 2^-286  >  10^-67
        0xD7ADF884AA879177, // * 2^-283  <  10^-66
        0x86CCBB52EA94BAEB, // * 2^-279  >  10^-65
        0xA87FEA27A539E9A5, // * 2^-276  <  10^-64
        0xD29FE4B18E88640F, // * 2^-273  >  10^-63
        0x83A3EEEEF9153E89, // * 2^-269  <  10^-62
        0xA48CEAAAB75A8E2B, // * 2^-266  <  10^-61
        0xCDB02555653131B6, // * 2^-263  <  10^-60
        0x808E17555F3EBF12, // * 2^-259  >  10^-59
        0xA0B19D2AB70E6ED6, // * 2^-256  <  10^-58
        0xC8DE047564D20A8C, // * 2^-253  >  10^-57
        0xFB158592BE068D2F, // * 2^-250  >  10^-56
        0x9CED737BB6C4183D, // * 2^-246  <  10^-55
        0xC428D05AA4751E4D, // * 2^-243  >  10^-54
        0xF53304714D9265E0, // * 2^-240  >  10^-53
        0x993FE2C6D07B7FAC, // * 2^-236  >  10^-52
        0xBF8FDB78849A5F97, // * 2^-233  >  10^-51
        0xEF73D256A5C0F77D, // * 2^-230  >  10^-50
        0x95A8637627989AAE, // * 2^-226  >  10^-49
        0xBB127C53B17EC159, // * 2^-223  <  10^-48
        0xE9D71B689DDE71B0, // * 2^-220  >  10^-47
        0x9226712162AB070E, // * 2^-216  >  10^-46
        0xB6B00D69BB55C8D1, // * 2^-213  <  10^-45
        0xE45C10C42A2B3B06, // * 2^-210  >  10^-44
        0x8EB98A7A9A5B04E3, // * 2^-206  <  10^-43
        0xB267ED1940F1C61C, // * 2^-203  <  10^-42
        0xDF01E85F912E37A3, // * 2^-200  <  10^-41
        0x8B61313BBABCE2C6, // * 2^-196  <  10^-40
        0xAE397D8AA96C1B78, // * 2^-193  >  10^-39
        0xD9C7DCED53C72256, // * 2^-190  >  10^-38
        0x881CEA14545C7575, // * 2^-186  <  10^-37
        0xAA242499697392D3, // * 2^-183  >  10^-36
        0xD4AD2DBFC3D07788, // * 2^-180  >  10^-35
        0x84EC3C97DA624AB5, // * 2^-176  >  10^-34
        0xA6274BBDD0FADD62, // * 2^-173  >  10^-33
        0xCFB11EAD453994BA, // * 2^-170  <  10^-32
        0x81CEB32C4B43FCF5, // * 2^-166  >  10^-31
        0xA2425FF75E14FC32, // * 2^-163  >  10^-30
        0xCAD2F7F5359A3B3E, // * 2^-160  <  10^-29
        0xFD87B5F28300CA0E, // * 2^-157  >  10^-28
        0x9E74D1B791E07E48, // * 2^-153  <  10^-27
        0xC612062576589DDB, // * 2^-150  >  10^-26
        0xF79687AED3EEC551, // * 2^-147  <  10^-25
        0x9ABE14CD44753B53, // * 2^-143  >  10^-24
        0xC16D9A0095928A27, // * 2^-140  <  10^-23
        0xF1C90080BAF72CB1, // * 2^-137  <  10^-22
        0x971DA05074DA7BEF, // * 2^-133  >  10^-21
        0xBCE5086492111AEB, // * 2^-130  >  10^-20
        0xEC1E4A7DB69561A5, // * 2^-127  <  10^-19
        0x9392EE8E921D5D07, // * 2^-123  <  10^-18
        0xB877AA3236A4B449, // * 2^-120  <  10^-17
        0xE69594BEC44DE15B, // * 2^-117  <  10^-16
        0x901D7CF73AB0ACD9, // * 2^-113  <  10^-15
        0xB424DC35095CD80F, // * 2^-110  <  10^-14
        0xE12E13424BB40E13, // * 2^-107  <  10^-13
        0x8CBCCC096F5088CC, // * 2^-103  >  10^-12
        0xAFEBFF0BCB24AAFF, // * 2^-100  >  10^-11
        0xDBE6FECEBDEDD5BF, // * 2^-97   >  10^-10
        0x89705F4136B4A597, // * 2^-93   <  10^-9
        0xABCC77118461CEFD, // * 2^-90   >  10^-8
        0xD6BF94D5E57A42BC, // * 2^-87   <  10^-7
        0x8637BD05AF6C69B6, // * 2^-83   >  10^-6
        0xA7C5AC471B478423, // * 2^-80   <  10^-5
        0xD1B71758E219652C, // * 2^-77   >  10^-4
        0x83126E978D4FDF3B, // * 2^-73   <  10^-3
        0xA3D70A3D70A3D70A, // * 2^-70   <  10^-2
        0xCCCCCCCCCCCCCCCD, // * 2^-67   >  10^-1
        0x8000000000000000, // * 2^-63   == 10^0
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
        0x8E1BC9BF04000000, // * 2^-10   == 10^16
        0xB1A2BC2EC5000000, // * 2^-7    == 10^17
        0xDE0B6B3A76400000, // * 2^-4    == 10^18
        0x8AC7230489E80000, // * 2^0     == 10^19
        0xAD78EBC5AC620000, // * 2^3     == 10^20
        0xD8D726B7177A8000, // * 2^6     == 10^21
        0x878678326EAC9000, // * 2^10    == 10^22
        0xA968163F0A57B400, // * 2^13    == 10^23
        0xD3C21BCECCEDA100, // * 2^16    == 10^24
        0x84595161401484A0, // * 2^20    == 10^25
        0xA56FA5B99019A5C8, // * 2^23    == 10^26
        0xCECB8F27F4200F3A, // * 2^26    == 10^27
        0x813F3978F8940984, // * 2^30    <  10^28
        0xA18F07D736B90BE5, // * 2^33    <  10^29
        0xC9F2C9CD04674EDF, // * 2^36    >  10^30
        0xFC6F7C4045812296, // * 2^39    <  10^31
        0x9DC5ADA82B70B59E, // * 2^43    >  10^32
        0xC5371912364CE305, // * 2^46    <  10^33
        0xF684DF56C3E01BC7, // * 2^49    >  10^34
        0x9A130B963A6C115C, // * 2^53    <  10^35
        0xC097CE7BC90715B3, // * 2^56    <  10^36
        0xF0BDC21ABB48DB20, // * 2^59    <  10^37
        0x96769950B50D88F4, // * 2^63    <  10^38
        0xBC143FA4E250EB31, // * 2^66    <  10^39
        0xEB194F8E1AE525FD, // * 2^69    <  10^40
        0x92EFD1B8D0CF37BE, // * 2^73    <  10^41
        0xB7ABC627050305AE, // * 2^76    >  10^42
        0xE596B7B0C643C719, // * 2^79    <  10^43
        0x8F7E32CE7BEA5C70, // * 2^83    >  10^44
        0xB35DBF821AE4F38C, // * 2^86    >  10^45
        0xE0352F62A19E306F, // * 2^89    >  10^46
        0x8C213D9DA502DE45, // * 2^93    <  10^47
        0xAF298D050E4395D7, // * 2^96    >  10^48
        0xDAF3F04651D47B4C, // * 2^99    <  10^49
        0x88D8762BF324CD10, // * 2^103   >  10^50
        0xAB0E93B6EFEE0054, // * 2^106   >  10^51
        0xD5D238A4ABE98068, // * 2^109   <  10^52
        0x85A36366EB71F041, // * 2^113   <  10^53
        0xA70C3C40A64E6C52, // * 2^116   >  10^54
        0xD0CF4B50CFE20766, // * 2^119   >  10^55
        0x82818F1281ED44A0, // * 2^123   >  10^56
        0xA321F2D7226895C8, // * 2^126   >  10^57
        0xCBEA6F8CEB02BB3A, // * 2^129   >  10^58
        0xFEE50B7025C36A08, // * 2^132   <  10^59
        0x9F4F2726179A2245, // * 2^136   <  10^60
        0xC722F0EF9D80AAD6, // * 2^139   <  10^61
        0xF8EBAD2B84E0D58C, // * 2^142   >  10^62
        0x9B934C3B330C8577, // * 2^146   <  10^63
        0xC2781F49FFCFA6D5, // * 2^149   <  10^64
        0xF316271C7FC3908B, // * 2^152   >  10^65
        0x97EDD871CFDA3A57, // * 2^156   >  10^66
        0xBDE94E8E43D0C8EC, // * 2^159   <  10^67
        0xED63A231D4C4FB27, // * 2^162   <  10^68
        0x945E455F24FB1CF9, // * 2^166   >  10^69
        0xB975D6B6EE39E437, // * 2^169   >  10^70
        0xE7D34C64A9C85D44, // * 2^172   <  10^71
        0x90E40FBEEA1D3A4B, // * 2^176   >  10^72
        0xB51D13AEA4A488DD, // * 2^179   <  10^73
        0xE264589A4DCDAB15, // * 2^182   >  10^74
        0x8D7EB76070A08AED, // * 2^186   >  10^75
        0xB0DE65388CC8ADA8, // * 2^189   <  10^76
        0xDD15FE86AFFAD912, // * 2^192   <  10^77
        0x8A2DBF142DFCC7AB, // * 2^196   <  10^78
        0xACB92ED9397BF996, // * 2^199   <  10^79
        0xD7E77A8F87DAF7FC, // * 2^202   >  10^80
        0x86F0AC99B4E8DAFD, // * 2^206   <  10^81
        0xA8ACD7C0222311BD, // * 2^209   >  10^82
        0xD2D80DB02AABD62C, // * 2^212   >  10^83
        0x83C7088E1AAB65DB, // * 2^216   <  10^84
        0xA4B8CAB1A1563F52, // * 2^219   <  10^85
        0xCDE6FD5E09ABCF27, // * 2^222   >  10^86
        0x80B05E5AC60B6178, // * 2^226   <  10^87
        0xA0DC75F1778E39D6, // * 2^229   <  10^88
        0xC913936DD571C84C, // * 2^232   <  10^89
        0xFB5878494ACE3A5F, // * 2^235   <  10^90
        0x9D174B2DCEC0E47B, // * 2^239   <  10^91
        0xC45D1DF942711D9A, // * 2^242   <  10^92
        0xF5746577930D6501, // * 2^245   >  10^93
        0x9968BF6ABBE85F20, // * 2^249   <  10^94
        0xBFC2EF456AE276E9, // * 2^252   >  10^95
        0xEFB3AB16C59B14A3, // * 2^255   >  10^96
        0x95D04AEE3B80ECE6, // * 2^259   >  10^97
        0xBB445DA9CA61281F, // * 2^262   <  10^98
        0xEA1575143CF97227, // * 2^265   >  10^99
        0x924D692CA61BE758, // * 2^269   <  10^100
        0xB6E0C377CFA2E12E, // * 2^272   <  10^101
        0xE498F455C38B997A, // * 2^275   <  10^102
        0x8EDF98B59A373FEC, // * 2^279   <  10^103
        0xB2977EE300C50FE7, // * 2^282   <  10^104
        0xDF3D5E9BC0F653E1, // * 2^285   <  10^105
        0x8B865B215899F46D, // * 2^289   >  10^106
        0xAE67F1E9AEC07188, // * 2^292   >  10^107
        0xDA01EE641A708DEA, // * 2^295   >  10^108
        0x884134FE908658B2, // * 2^299   <  10^109
        0xAA51823E34A7EEDF, // * 2^302   >  10^110
        0xD4E5E2CDC1D1EA96, // * 2^305   <  10^111
        0x850FADC09923329E, // * 2^309   <  10^112
        0xA6539930BF6BFF46, // * 2^312   >  10^113
        0xCFE87F7CEF46FF17, // * 2^315   >  10^114
        0x81F14FAE158C5F6E, // * 2^319   <  10^115
        0xA26DA3999AEF774A, // * 2^322   >  10^116
        0xCB090C8001AB551C, // * 2^325   <  10^117
        0xFDCB4FA002162A63, // * 2^328   <  10^118
        0x9E9F11C4014DDA7E, // * 2^332   <  10^119
        0xC646D63501A1511E, // * 2^335   >  10^120
        0xF7D88BC24209A565, // * 2^338   <  10^121
        0x9AE757596946075F, // * 2^342   <  10^122
        0xC1A12D2FC3978937, // * 2^345   <  10^123
        0xF209787BB47D6B85, // * 2^348   >  10^124
        0x9745EB4D50CE6333, // * 2^352   >  10^125
        0xBD176620A501FC00, // * 2^355   >  10^126
        0xEC5D3FA8CE427B00, // * 2^358   >  10^127
        0x93BA47C980E98CE0, // * 2^362   >  10^128
        0xB8A8D9BBE123F018, // * 2^365   >  10^129
        0xE6D3102AD96CEC1E, // * 2^368   >  10^130
        0x9043EA1AC7E41393, // * 2^372   >  10^131
        0xB454E4A179DD1877, // * 2^375   <  10^132
        0xE16A1DC9D8545E95, // * 2^378   >  10^133
        0x8CE2529E2734BB1D, // * 2^382   <  10^134
        0xB01AE745B101E9E4, // * 2^385   <  10^135
        0xDC21A1171D42645D, // * 2^388   <  10^136
        0x899504AE72497EBA, // * 2^392   <  10^137
        0xABFA45DA0EDBDE69, // * 2^395   <  10^138
        0xD6F8D7509292D603, // * 2^398   <  10^139
        0x865B86925B9BC5C2, // * 2^402   <  10^140
        0xA7F26836F282B733, // * 2^405   >  10^141
        0xD1EF0244AF2364FF, // * 2^408   <  10^142
        0x8335616AED761F1F, // * 2^412   <  10^143
        0xA402B9C5A8D3A6E7, // * 2^415   <  10^144
        0xCD036837130890A1, // * 2^418   <  10^145
        0x802221226BE55A65, // * 2^422   >  10^146
        0xA02AA96B06DEB0FE, // * 2^425   >  10^147
        0xC83553C5C8965D3D, // * 2^428   <  10^148
        0xFA42A8B73ABBF48D, // * 2^431   >  10^149
        0x9C69A97284B578D8, // * 2^435   >  10^150
        0xC38413CF25E2D70E, // * 2^438   >  10^151
        0xF46518C2EF5B8CD1, // * 2^441   <  10^152
        0x98BF2F79D5993803, // * 2^445   >  10^153
        0xBEEEFB584AFF8604, // * 2^448   >  10^154
        0xEEAABA2E5DBF6785, // * 2^451   >  10^155
        0x952AB45CFA97A0B3, // * 2^455   >  10^156
        0xBA756174393D88E0, // * 2^458   >  10^157
        0xE912B9D1478CEB17, // * 2^461   <  10^158
        0x91ABB422CCB812EF, // * 2^465   >  10^159
        0xB616A12B7FE617AA, // * 2^468   <  10^160
        0xE39C49765FDF9D95, // * 2^471   >  10^161
        0x8E41ADE9FBEBC27D, // * 2^475   <  10^162
        0xB1D219647AE6B31C, // * 2^478   <  10^163
        0xDE469FBD99A05FE3, // * 2^481   <  10^164
        0x8AEC23D680043BEE, // * 2^485   <  10^165
        0xADA72CCC20054AEA, // * 2^488   >  10^166
        0xD910F7FF28069DA4, // * 2^491   <  10^167
        0x87AA9AFF79042287, // * 2^495   >  10^168
        0xA99541BF57452B28, // * 2^498   <  10^169
        0xD3FA922F2D1675F2, // * 2^501   <  10^170
        0x847C9B5D7C2E09B7, // * 2^505   <  10^171
        0xA59BC234DB398C25, // * 2^508   <  10^172
        0xCF02B2C21207EF2F, // * 2^511   >  10^173
        0x8161AFB94B44F57D, // * 2^515   <  10^174
        0xA1BA1BA79E1632DC, // * 2^518   <  10^175
        0xCA28A291859BBF93, // * 2^521   <  10^176
        0xFCB2CB35E702AF78, // * 2^524   <  10^177
        0x9DEFBF01B061ADAB, // * 2^528   <  10^178
        0xC56BAEC21C7A1916, // * 2^531   <  10^179
        0xF6C69A72A3989F5C, // * 2^534   >  10^180
        0x9A3C2087A63F6399, // * 2^538   <  10^181
        0xC0CB28A98FCF3C80, // * 2^541   >  10^182
        0xF0FDF2D3F3C30B9F, // * 2^544   <  10^183
        0x969EB7C47859E744, // * 2^548   >  10^184
        0xBC4665B596706115, // * 2^551   >  10^185
        0xEB57FF22FC0C795A, // * 2^554   >  10^186
        0x9316FF75DD87CBD8, // * 2^558   <  10^187
        0xB7DCBF5354E9BECE, // * 2^561   <  10^188
        0xE5D3EF282A242E82, // * 2^564   >  10^189
        0x8FA475791A569D11, // * 2^568   >  10^190
        0xB38D92D760EC4455, // * 2^571   <  10^191
        0xE070F78D3927556B, // * 2^574   >  10^192
        0x8C469AB843B89563, // * 2^578   >  10^193
        0xAF58416654A6BABB, // * 2^581   <  10^194
        0xDB2E51BFE9D0696A, // * 2^584   <  10^195
        0x88FCF317F22241E2, // * 2^588   <  10^196
        0xAB3C2FDDEEAAD25B, // * 2^591   >  10^197
        0xD60B3BD56A5586F2, // * 2^594   >  10^198
        0x85C7056562757457, // * 2^598   >  10^199
        0xA738C6BEBB12D16D, // * 2^601   >  10^200
        0xD106F86E69D785C8, // * 2^604   >  10^201
        0x82A45B450226B39D, // * 2^608   >  10^202
        0xA34D721642B06084, // * 2^611   <  10^203
        0xCC20CE9BD35C78A5, // * 2^614   <  10^204
        0xFF290242C83396CE, // * 2^617   <  10^205
        0x9F79A169BD203E41, // * 2^621   <  10^206
        0xC75809C42C684DD1, // * 2^624   <  10^207
        0xF92E0C3537826146, // * 2^627   >  10^208
        0x9BBCC7A142B17CCC, // * 2^631   >  10^209
        0xC2ABF989935DDBFE, // * 2^634   <  10^210
        0xF356F7EBF83552FE, // * 2^637   <  10^211
        0x98165AF37B2153DF, // * 2^641   >  10^212
        0xBE1BF1B059E9A8D6, // * 2^644   <  10^213
        0xEDA2EE1C7064130C, // * 2^647   <  10^214
        0x9485D4D1C63E8BE8, // * 2^651   >  10^215
        0xB9A74A0637CE2EE1, // * 2^654   <  10^216
        0xE8111C87C5C1BA9A, // * 2^657   >  10^217
        0x910AB1D4DB9914A0, // * 2^661   <  10^218
        0xB54D5E4A127F59C8, // * 2^664   <  10^219
        0xE2A0B5DC971F303A, // * 2^667   <  10^220
        0x8DA471A9DE737E24, // * 2^671   <  10^221
        0xB10D8E1456105DAD, // * 2^674   <  10^222
        0xDD50F1996B947519, // * 2^677   >  10^223
        0x8A5296FFE33CC930, // * 2^681   >  10^224
        0xACE73CBFDC0BFB7B, // * 2^684   <  10^225
        0xD8210BEFD30EFA5A, // * 2^687   <  10^226
        0x8714A775E3E95C78, // * 2^691   <  10^227
        0xA8D9D1535CE3B396, // * 2^694   <  10^228
        0xD31045A8341CA07C, // * 2^697   <  10^229
        0x83EA2B892091E44E, // * 2^701   >  10^230
        0xA4E4B66B68B65D61, // * 2^704   >  10^231
        0xCE1DE40642E3F4B9, // * 2^707   <  10^232
        0x80D2AE83E9CE78F4, // * 2^711   >  10^233
        0xA1075A24E4421731, // * 2^714   >  10^234
        0xC94930AE1D529CFD, // * 2^717   >  10^235
        0xFB9B7CD9A4A7443C, // * 2^720   <  10^236
        0x9D412E0806E88AA6, // * 2^724   >  10^237
        0xC491798A08A2AD4F, // * 2^727   >  10^238
        0xF5B5D7EC8ACB58A3, // * 2^730   >  10^239
        0x9991A6F3D6BF1766, // * 2^734   >  10^240
        0xBFF610B0CC6EDD3F, // * 2^737   <  10^241
        0xEFF394DCFF8A948F, // * 2^740   >  10^242
        0x95F83D0A1FB69CD9, // * 2^744   <  10^243
        0xBB764C4CA7A44410, // * 2^747   >  10^244
        0xEA53DF5FD18D5514, // * 2^750   >  10^245
        0x92746B9BE2F8552C, // * 2^754   <  10^246
        0xB7118682DBB66A77, // * 2^757   <  10^247
        0xE4D5E82392A40515, // * 2^760   <  10^248
        0x8F05B1163BA6832D, // * 2^764   <  10^249
        0xB2C71D5BCA9023F8, // * 2^767   <  10^250
        0xDF78E4B2BD342CF7, // * 2^770   >  10^251
        0x8BAB8EEFB6409C1A, // * 2^774   <  10^252
        0xAE9672ABA3D0C321, // * 2^777   >  10^253
        0xDA3C0F568CC4F3E9, // * 2^780   >  10^254
        0x8865899617FB1871, // * 2^784   <  10^255
        0xAA7EEBFB9DF9DE8E, // * 2^787   >  10^256
        0xD51EA6FA85785631, // * 2^790   <  10^257
        0x8533285C936B35DF, // * 2^794   >  10^258
        0xA67FF273B8460357, // * 2^797   >  10^259
        0xD01FEF10A657842C, // * 2^800   <  10^260
        0x8213F56A67F6B29C, // * 2^804   >  10^261
        0xA298F2C501F45F43, // * 2^807   >  10^262
        0xCB3F2F7642717713, // * 2^810   <  10^263
        0xFE0EFB53D30DD4D8, // * 2^813   >  10^264
        0x9EC95D1463E8A507, // * 2^817   >  10^265
        0xC67BB4597CE2CE49, // * 2^820   >  10^266
        0xF81AA16FDC1B81DB, // * 2^823   >  10^267
        0x9B10A4E5E9913129, // * 2^827   >  10^268
        0xC1D4CE1F63F57D73, // * 2^830   >  10^269
        0xF24A01A73CF2DCD0, // * 2^833   >  10^270
        0x976E41088617CA02, // * 2^837   >  10^271
        0xBD49D14AA79DBC82, // * 2^840   <  10^272
        0xEC9C459D51852BA3, // * 2^843   >  10^273
        0x93E1AB8252F33B46, // * 2^847   >  10^274
        0xB8DA1662E7B00A17, // * 2^850   <  10^275
        0xE7109BFBA19C0C9D, // * 2^853   <  10^276
        0x906A617D450187E2, // * 2^857   <  10^277
        0xB484F9DC9641E9DB, // * 2^860   >  10^278
        0xE1A63853BBD26451, // * 2^863   <  10^279
        0x8D07E33455637EB3, // * 2^867   >  10^280
        0xB049DC016ABC5E60, // * 2^870   >  10^281
        0xDC5C5301C56B75F7, // * 2^873   <  10^282
        0x89B9B3E11B6329BB, // * 2^877   >  10^283
        0xAC2820D9623BF429, // * 2^880   <  10^284
        0xD732290FBACAF134, // * 2^883   >  10^285
        0x867F59A9D4BED6C0, // * 2^887   <  10^286
        0xA81F301449EE8C70, // * 2^890   <  10^287
        0xD226FC195C6A2F8C, // * 2^893   <  10^288
        0x83585D8FD9C25DB8, // * 2^897   >  10^289
        0xA42E74F3D032F526, // * 2^900   >  10^290
        0xCD3A1230C43FB26F, // * 2^903   <  10^291
        0x80444B5E7AA7CF85, // * 2^907   <  10^292
        0xA0555E361951C367, // * 2^910   >  10^293
        0xC86AB5C39FA63441, // * 2^913   >  10^294
        0xFA856334878FC151, // * 2^916   >  10^295
        0x9C935E00D4B9D8D2, // * 2^920   <  10^296
        0xC3B8358109E84F07, // * 2^923   <  10^297
        0xF4A642E14C6262C9, // * 2^926   >  10^298
        0x98E7E9CCCFBD7DBE, // * 2^930   >  10^299
        0xBF21E44003ACDD2D, // * 2^933   >  10^300
        0xEEEA5D5004981478, // * 2^936   <  10^301
        0x95527A5202DF0CCB, // * 2^940   <  10^302
        0xBAA718E68396CFFE, // * 2^943   >  10^303
        0xE950DF20247C83FD, // * 2^946   <  10^304
        0x91D28B7416CDD27E, // * 2^950   <  10^305
        0xB6472E511C81471E, // * 2^953   >  10^306
        0xE3D8F9E563A198E5, // * 2^956   <  10^307
        0x8E679C2F5E44FF8F, // * 2^960   <  10^308
        0xB201833B35D63F73, // * 2^963   <  10^309
        0xDE81E40A034BCF50, // * 2^966   >  10^310
        0x8B112E86420F6192, // * 2^970   >  10^311
        0xADD57A27D29339F6, // * 2^973   <  10^312
        0xD94AD8B1C7380874, // * 2^976   <  10^313
        0x87CEC76F1C830549, // * 2^980   >  10^314
        0xA9C2794AE3A3C69B, // * 2^983   >  10^315
        0xD433179D9C8CB841, // * 2^986   <  10^316
        0x849FEEC281D7F329, // * 2^990   >  10^317
        0xA5C7EA73224DEFF3, // * 2^993   <  10^318
        0xCF39E50FEAE16BF0, // * 2^996   >  10^319
        0x81842F29F2CCE376, // * 2^1000  >  10^320
        0xA1E53AF46F801C53, // * 2^1003  <  10^321
        0xCA5E89B18B602368, // * 2^1006  <  10^322
        0xFCF62C1DEE382C42, // * 2^1009  <  10^323
        0x9E19DB92B4E31BA9, // * 2^1013  <  10^324
    };

    STRTOD_ASSERT(index >= 0);
    STRTOD_ASSERT(index < kCachedPowersSize);

    int const k = kCachedPowersMinDecExp + index * kCachedPowersDecExpStep;
    int const e = BinaryExponentFromDecimalExponent(k);

    return {kSignificands[index], e, k};
}

#endif // ^^^ !STRTOD_OPTIMIZE_SIZE

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

#if STRTOD_OPTIMIZE_SIZE
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
#endif

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

#if 1
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
