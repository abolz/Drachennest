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

#include "dtoa.h"

#if DTOA_UNNAMED_NAMESPACE
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

namespace impl {

inline constexpr int Min(int x, int y) { return y < x ? y : x; }
inline constexpr int Max(int x, int y) { return y < x ? x : y; }

inline bool IsDigit(char ch)
{
#if 0
    return static_cast<unsigned>(ch - '0') < 10;
#else
    return '0' <= ch && ch <= '9';
#endif
}

inline int DigitValue(char ch)
{
    DTOA_ASSERT(IsDigit(ch));
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
#define DTOA_CORRECT_DOUBLE_OPERATIONS 1
#elif defined(_M_IX86) || defined(__i386__) || defined(__i386)
#ifdef _WIN32
// Windows uses a 64bit wide floating point stack.
#define DTOA_CORRECT_DOUBLE_OPERATIONS 1
#endif
#endif

// 2^53 = 9007199254740992.
// Any integer with at most 15 decimal digits will hence fit into a double
// (which has a 53bit significand) without loss of precision.
constexpr int kMaxExactDoubleIntegerDecimalDigits = 15;

#if DTOA_CORRECT_DOUBLE_OPERATIONS

inline bool FastPath(double& result, uint64_t digits, int num_digits, int exponent)
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

    DTOA_ASSERT(num_digits <= kMaxExactDoubleIntegerDecimalDigits);

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

#else // ^^^ DTOA_CORRECT_DOUBLE_OPERATIONS

inline bool FastPath(double& /*result*/, uint64_t /*digits*/, int /*num_digits*/, int /*exponent*/)
{
    return false;
}

#endif // ^^^ !DTOA_CORRECT_DOUBLE_OPERATIONS

//--------------------------------------------------------------------------------------------------
// StrtodApprox
//--------------------------------------------------------------------------------------------------

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

// Normalize x
// and scale the error, so that the error is in ULP(x)
inline void Normalize(DiyFpWithError& num)
{
    int const s = CountLeadingZeros64(num.x.f);

    DTOA_ASSERT(((num.error << s) >> s) == num.error);

    num.x.f   <<= s;
    num.x.e    -= s;
    num.error <<= s;
}

// 2^64 = 18446744073709551616 > 10^19
// Any integer with at most 19 decimal digits will hence fit into an uint64_t.
constexpr int kMaxUint64DecimalDigits = 19;

template <typename Int>
inline Int ReadInt(char const* f, char const* l)
{
    DTOA_ASSERT(l - f <= std::numeric_limits<Int>::digits10);

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
        DTOA_ASSERT(IsDigit(*f));
        value = 10 * value + static_cast<unsigned char>(*f) - '0';
#else
        value = 10 * value + static_cast<uint32_t>(DigitValue(*f));
#endif
    }

    return value;
}

// Returns a cached power of ten x ~= 10^k such that
//  k <= e < k + kCachedPowersDecExpStep.
//
// PRE: e >= kCachedPowersMinDecExp
// PRE: e <  kCachedPowersMaxDecExp + kCachedPowersDecExpStep
inline CachedPower GetCachedPowerForDecimalExponent(int e)
{
    DTOA_ASSERT(e >= kCachedPowersMinDecExp);
    DTOA_ASSERT(e <  kCachedPowersMaxDecExp + kCachedPowersDecExpStep);

    int const index = static_cast<int>( static_cast<unsigned>(-kCachedPowersMinDecExp + e) / kCachedPowersDecExpStep );
    DTOA_ASSERT(index >= 0);
    DTOA_ASSERT(index < kCachedPowersSize);

    auto const cached = GetCachedPower(index);
    DTOA_ASSERT(e >= cached.k);
    DTOA_ASSERT(e <  cached.k + kCachedPowersDecExpStep);

    return cached;
}

// Returns 10^k as an exact DiyFp.
// PRE: 1 <= k < kCachedPowersDecExpStep
inline DiyFp GetAdjustmentPowerOfTen(int k)
{
    static_assert(kCachedPowersDecExpStep <= 8, "internal error");

    static constexpr uint64_t kSignificands[] = {
        0x8000000000000000, // e = -63, == 10^0 (unused)
        0xA000000000000000, // e = -60, == 10^1
        0xC800000000000000, // e = -57, == 10^2
        0xFA00000000000000, // e = -54, == 10^3
        0x9C40000000000000, // e = -50, == 10^4
        0xC350000000000000, // e = -47, == 10^5
        0xF424000000000000, // e = -44, == 10^6
        0x9896800000000000, // e = -40, == 10^7
    };

    DTOA_ASSERT(k > 0);
    DTOA_ASSERT(k < kCachedPowersDecExpStep);

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
inline int EffectiveSignificandSize(int order)
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
inline double LoadDouble(uint64_t f, int e)
{
    using Double = IEEE<double>;

    DTOA_ASSERT(f <= Double::HiddenBit + Double::SignificandMask);
    DTOA_ASSERT(e <= Double::MinExponent || (f & Double::HiddenBit) != 0);

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
inline bool StrtodApprox(double& result, char const* digits, int num_digits, int exponent)
{
    using Double = IEEE<double>;

    static_assert(DiyFp::SignificandSize == 64,
        "We use uint64's. This only works if the DiyFp uses uint64's too.");

    DTOA_ASSERT(num_digits > 0);
    DTOA_ASSERT(DigitValue(digits[0]) > 0);
//  DTOA_ASSERT(DigitValue(digits[num_digits - 1]) > 0);
    DTOA_ASSERT(num_digits + exponent <= kMaxDecimalPower);
    DTOA_ASSERT(num_digits + exponent >  kMinDecimalPower);

    // Compute an approximation 'input' for B = digits * 10^exponent using DiyFp's.
    // And keep track of the error.
    //
    //                       <-- error -->
    //                               B = digits * 10^exponent
    //  ---------(-----------|-------+---)------------------------------------
    //                       x
    //                       ~= (f * 2^e) * 10^exponent

    constexpr int kLogULP = DiyFpWithError::kDenominatorLog;
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
    DTOA_ASSERT(input.error <= 16 * (kULP / 2));

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

        DTOA_ASSERT(IsNormalized(input.x));
        DTOA_ASSERT(IsNormalized(adjustment_power));

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

            DTOA_ASSERT(input.error <= 17 * (kULP / 2));
        }

        // The result of the multiplication might not be normalized.
        // Normalize 'x' again and scale the error.
        Normalize(input);

        // Since both factors are normalized, input.f >= 2^(q-2), and the scaling
        // factor in the normalization step above is bounded by 2^1.
        DTOA_ASSERT(input.error <= 34 * (kULP / 2));
    }

    DTOA_ASSERT(IsNormalized(input.x));
    DTOA_ASSERT(IsNormalized(cached_power));

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

    DTOA_ASSERT(input.error <= 36 * (kULP / 2));

    // The result of the multiplication might not be normalized.
    // Normalize 'x' again and scale the error.
    Normalize(input);

    // Since both factors were normalized, the scaling factor in the
    // normalization step above is again bounded by 2^1.
    DTOA_ASSERT(input.error <= 72 * (kULP / 2));

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
    DTOA_ASSERT(prec >= 0);
    DTOA_ASSERT(prec <= 53);

    int excess_bits = DiyFp::SignificandSize - prec;
    if (excess_bits > DiyFp::SignificandSize - kLogULP - 1)
    {
        // In this case 'half' (see below) multiplied by kULP exceeds the range of an uint64_t.
        // This can only happen for very small subnormals (when excess_bits is large).

        int const s = excess_bits - (DiyFp::SignificandSize - kLogULP - 1);
        DTOA_ASSERT(s > 0);

#if 1
        uint64_t const discarded_bits = input.x.f & ((uint64_t{1} << s) - 1);

        // Move the discarded bits into the error: (f + err) * 2^e = (f - d + err + d) * 2^e
        DTOA_ASSERT(discarded_bits <= UINT32_MAX);
        input.error += static_cast<uint32_t>(discarded_bits);
        // Scale the error such that input.error is in ULP(input.x) again.
        input.error >>= s;
        // And add 1 so that input.error is still an upper bound.
        input.error += 1;

        // error = ((error + discarded_bits) div 2^s) + 1
        //       < ((error + 2^s           ) div 2^s) + 1
        //       = ((error div 2^s) + 1    )          + 1
        //      <= (error div 2) + 2
#else
        input.error = (input.error >> s) + 2;
#endif

        // x = f * 2^e ~= floor(f / 2^s) * 2^(e + s)
        input.x.f >>= s;
        input.x.e  += s;

        excess_bits = DiyFp::SignificandSize - kLogULP - 1;
    }

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

    DTOA_ASSERT(excess_bits >= 11);
    DTOA_ASSERT(excess_bits <  64);
    DTOA_ASSERT(excess_bits <= DiyFp::SignificandSize - kLogULP - 1);

    uint64_t const two_n = uint64_t{1} << excess_bits;

    uint64_t p2   = input.x.f & (two_n - 1);
    uint64_t half = two_n / 2;

    // error is scaled by kULP.
    // In order to compare p2 and half with error, these values need to be scaled, too.
    DTOA_ASSERT(p2   <= UINT64_MAX / kULP);
    DTOA_ASSERT(half <= UINT64_MAX / kULP);
    p2   *= kULP;
    half *= kULP;

    // Truncate the significand to p = q - n bits and move the discarded bits into the (binary) exponent.
    input.x.f >>= excess_bits;
    input.x.e  += excess_bits;

    DTOA_ASSERT(input.error > 0);
    DTOA_ASSERT(half >= input.error);

    // Note:
    // Since error is non-zero, we can safely use '<=' and '>=' in the comparisons below.

    bool success;
    if (p2 >= half + input.error)
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
            DTOA_ASSERT(input.x.f == (Double::HiddenBit << 1));

            input.x.f >>= 1;
            input.x.e  += 1;
        }
    }
    else if (p2 <= half - input.error)
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

inline bool ComputeGuess(double& result, char const* digits, int num_digits, int exponent)
{
    DTOA_ASSERT(num_digits > 0);
    DTOA_ASSERT(num_digits <= kMaxSignificantDigits);
    DTOA_ASSERT(DigitValue(digits[0]) > 0);
//  DTOA_ASSERT(DigitValue(digits[num_digits - 1]) > 0);

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

inline void AssignZero(DiyInt& x)
{
    x.size = 0;
    x.exponent = 0;
}

inline void AssignU32(DiyInt& x, uint32_t value)
{
    AssignZero(x);

    if (value == 0)
        return;

    x.bigits[0] = value;
    x.size = 1;
}

inline void AssignU64(DiyInt& x, uint64_t value)
{
    AssignZero(x);

    if (value == 0)
        return;

    x.bigits[0] = static_cast<uint32_t>(value);
    x.bigits[1] = static_cast<uint32_t>(value >> DiyInt::BigitSize);
    x.size = (x.bigits[1] == 0) ? 1 : 2;
}

// x := A * x + B
inline void MulAddU32(DiyInt& x, uint32_t A, uint32_t B = 0)
{
    DTOA_ASSERT(B == 0 || x.exponent == 0);

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
        DTOA_ASSERT(x.size < DiyInt::Capacity);
        x.bigits[x.size++] = carry;
    }
}

inline void AssignDecimalDigits(DiyInt& x, char const* digits, int num_digits)
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

inline void MulPow2(DiyInt& x, int exp) // aka left-shift
{
    DTOA_ASSERT(exp >= 0);

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
            DTOA_ASSERT(x.size < DiyInt::Capacity);
            x.bigits[x.size++] = carry;
        }
    }

    x.exponent += bigit_shift;
}

inline void MulPow5(DiyInt& x, int exp)
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

    DTOA_ASSERT(exp >= 0);
    if (exp == 0)
        return;

    while (exp > 0)
    {
        int const n = Min(exp, 13);
        MulAddU32(x, kPow5[n]);
        exp -= n;
    }
}

inline int Compare(DiyInt const& lhs, DiyInt const& rhs)
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
inline int CompareBufferWithDiyFp(char const* digits, int num_digits, int exponent, bool nonzero_tail, DiyFp v)
{
    DTOA_ASSERT(num_digits > 0);
    DTOA_ASSERT(num_digits + exponent <= kMaxDecimalPower);
    DTOA_ASSERT(num_digits + exponent >  kMinDecimalPower);
    DTOA_ASSERT(num_digits            <= kMaxSignificantDigits);

    DiyInt lhs;
    DiyInt rhs;

    AssignDecimalDigits(lhs, digits, num_digits);
    if (nonzero_tail)
    {
        MulAddU32(lhs, 10, 1);
        exponent--;
    }
    AssignU64(rhs, v.f);

    DTOA_ASSERT(lhs.size <= (2555 + 31) / 32); // bits <= log_2(10^769) = 2555
    DTOA_ASSERT(rhs.size <= (  64 + 31) / 32); // bits <= 64

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
        DTOA_ASSERT(lhs.size <= (1030 + 31) / 32);  // 1030 = log_2(10^(309 + 1)))
        DTOA_ASSERT(rhs.size <= (  64 + 31) / 32);
    }
    else if (rhs_exp5 > 0)
    {
        MulPow5(rhs, rhs_exp5);

        // kMinDecimalPower + 1 <= num_digits + exponent <= kMaxDecimalPower + 1
        // rhs_exp5 = -exponent <= -kMinDecimalPower - 1 + num_digits = 324 - 1 + num_digits <= 324 - 1 + 769
        // rhs_exp5 = -exponent >= -kMaxDecimalPower - 1 + num_digits
        DTOA_ASSERT(lhs.size <= (2555        + 31) / 32);
        DTOA_ASSERT(rhs.size <= (  64 + 2536 + 31) / 32); // 2536 = log_2(5^(324 - 1 + 769)) ---- XXX: 2504
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

    DTOA_ASSERT(lhs.size <= (2555        + 32 + 31) / 32);
    DTOA_ASSERT(rhs.size <= (  64 + 2536 + 32 + 31) / 32);

    return Compare(lhs, rhs);
}

//--------------------------------------------------------------------------------------------------
// DecimalToDouble
//--------------------------------------------------------------------------------------------------

// Returns whether the significand f of v = f * 2^e is even.
inline bool SignificandIsEven(double v)
{
    return (IEEE<double>(v).PhysicalSignificand() & 1) == 0;
}

// Returns the next larger double-precision value.
// If v is +Infinity returns v.
inline double NextFloat(double v)
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
    DTOA_ASSERT(num_digits >= 0);
    DTOA_ASSERT(exponent <= INT_MAX - num_digits);

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
        DTOA_ASSERT(DigitValue(digits[num_digits - 1]) > 0); // since trailing zeros have been trimmed above.

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

} // namespace impl

// Convert the decimal representation 'digits * 10^exponent' into an IEEE
// double-precision number.
//
// PRE: digits must contain only ASCII characters in the range '0'...'9'.
// PRE: num_digits >= 0
// PRE: num_digits + exponent must not overflow.
inline double DecimalToDouble(char const* digits, int num_digits, int exponent, bool nonzero_tail = false)
{
    return base_conv::impl::DecimalToDouble(digits, num_digits, exponent, nonzero_tail);
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
    using base_conv::impl::IsDigit;
    using base_conv::impl::DigitValue;

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
    value = base_conv::impl::DecimalToDouble(digits, num_digits, exponent, nonzero_tail);

L_done:
    result = is_neg ? -value : value;
    next = curr;

    return status;
}

inline double Strtod(char const* first, char const* last)
{
    double d = 0.0;
    base_conv::Strtod(d, first, last);
    return d;
}

} // namespace base_conv
#if DTOA_UNNAMED_NAMESPACE
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
