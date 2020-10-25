// Copyright 2020 Alexander Bolz
//
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#include "dragonbox.h"

#if 0

#include "jkj_dragonbox_to_chars.h"

char* dragonbox::Dtoa(char* buffer, double value)
{
    return jkj::dragonbox::to_chars_n(value, buffer);
}

#else

#include "jkj_dragonbox.h"

#include <cassert>
#include <cstdint>
#include <cstring>

#ifndef DRAGONBOX_ASSERT
#define DRAGONBOX_ASSERT(X) assert(X)
#endif

//==================================================================================================
//
//==================================================================================================

template <typename Dest, typename Source>
static inline Dest ReinterpretBits(Source source)
{
    static_assert(sizeof(Dest) == sizeof(Source), "size mismatch");

    Dest dest;
    std::memcpy(&dest, &source, sizeof(Source));
    return dest;
}

namespace {
struct Double
{
    static_assert(std::numeric_limits<double>::is_iec559
               && std::numeric_limits<double>::digits == 53
               && std::numeric_limits<double>::max_exponent == 1024,
        "IEEE-754 double-precision implementation required");

    using value_type = double;
    using bits_type = uint64_t;

//  static constexpr int32_t   MaxDigits10     = std::numeric_limits<value_type>::max_digits10;
    static constexpr int32_t   SignificandSize = std::numeric_limits<value_type>::digits; // = p   (includes the hidden bit)
    static constexpr int32_t   ExponentBias    = std::numeric_limits<value_type>::max_exponent - 1 + (SignificandSize - 1);
//  static constexpr int32_t   MaxExponent     = std::numeric_limits<value_type>::max_exponent - 1 - (SignificandSize - 1);
//  static constexpr int32_t   MinExponent     = std::numeric_limits<value_type>::min_exponent - 1 - (SignificandSize - 1);
    static constexpr bits_type MaxIeeeExponent = bits_type{2 * std::numeric_limits<value_type>::max_exponent - 1};
    static constexpr bits_type HiddenBit       = bits_type{1} << (SignificandSize - 1);   // = 2^(p-1)
    static constexpr bits_type SignificandMask = HiddenBit - 1;                           // = 2^(p-1) - 1
    static constexpr bits_type ExponentMask    = MaxIeeeExponent << (SignificandSize - 1);
    static constexpr bits_type SignMask        = ~(~bits_type{0} >> 1);

    bits_type bits;

    explicit Double(bits_type bits_) : bits(bits_) {}
    explicit Double(value_type value) : bits(ReinterpretBits<bits_type>(value)) {}

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
};
} // namespace

//==================================================================================================
// ToChars
//==================================================================================================

static inline void Utoa_2Digits(char* buf, uint32_t digits)
{
    static constexpr char Digits100[200] = {
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

    DRAGONBOX_ASSERT(digits <= 99);
    std::memcpy(buf, &Digits100[2 * digits], 2 * sizeof(char));
}

static inline int TrailingZeros_2Digits(uint32_t digits)
{
    static constexpr int8_t TrailingZeros100[100] = {
        2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };

    DRAGONBOX_ASSERT(digits <= 99);
    return TrailingZeros100[digits];
}

static inline int Utoa_8Digits_skip_trailing_zeros(char* buf, uint32_t digits)
{
    DRAGONBOX_ASSERT(digits >= 1);
    DRAGONBOX_ASSERT(digits <= 99999999);

    const uint32_t q = digits / 10000;
    const uint32_t r = digits % 10000;

    const uint32_t qH = q / 100;
    const uint32_t qL = q % 100;
    Utoa_2Digits(buf + 0, qH);
    Utoa_2Digits(buf + 2, qL);

    if (r == 0)
    {
        return TrailingZeros_2Digits(qL == 0 ? qH : qL) + (qL == 0 ? 6 : 4);
    }
    else
    {
        const uint32_t rH = r / 100;
        const uint32_t rL = r % 100;
        Utoa_2Digits(buf + 4, rH);
        Utoa_2Digits(buf + 6, rL);

        return TrailingZeros_2Digits(rL == 0 ? rH : rL) + (rL == 0 ? 2 : 0);
    }
}

static inline int PrintDecimalDigitsBackwards(char* buf, uint64_t output64)
{
    int tz = 0; // number of trailing zeros removed.
    int nd = 0; // number of decimal digits processed.

    // At most 17 digits remaining

    if (output64 >= 100000000)
    {
        const uint64_t q = output64 / 100000000;
        const uint32_t r = static_cast<uint32_t>(output64 % 100000000);
        output64 = q;
        buf -= 8;
        if (r != 0)
        {
            tz = Utoa_8Digits_skip_trailing_zeros(buf, r);
            DRAGONBOX_ASSERT(tz >= 0);
            DRAGONBOX_ASSERT(tz <= 7);
        }
        else
        {
            tz = 8;
        }
        nd = 8;
    }

    // At most 9 digits remaining.
    DRAGONBOX_ASSERT(output64 <= UINT32_MAX);
    uint32_t output = static_cast<uint32_t>(output64);

    if (output >= 10000)
    {
        const uint32_t q = output / 10000;
        const uint32_t r = output % 10000;
        output = q;
        buf -= 4;
        if (r != 0)
        {
            const uint32_t rH = r / 100;
            const uint32_t rL = r % 100;
            Utoa_2Digits(buf + 0, rH);
            Utoa_2Digits(buf + 2, rL);
            if (tz == nd)
            {
                tz += TrailingZeros_2Digits(rL == 0 ? rH : rL) + (rL == 0 ? 2 : 0);
            }
        }
        else
        {
            if (tz == nd)
                tz += 4;
            else
                std::memset(buf, '0', 4); // (actually not required...)
        }
        nd += 4;
    }

    // At most 5 digits remaining.

    if (output >= 100)
    {
        const uint32_t q = output / 100;
        const uint32_t r = output % 100;
        output = q;
        buf -= 2;
        Utoa_2Digits(buf, r);
        if (tz == nd)
        {
            tz += TrailingZeros_2Digits(r);
        }
        nd += 2;

        if (output >= 100)
        {
            const uint32_t q = output / 100;
            const uint32_t r = output % 100;
            output = q;
            buf -= 2;
            Utoa_2Digits(buf, r);
            if (tz == nd)
            {
                tz += TrailingZeros_2Digits(r);
            }
            nd += 2;
        }
    }

    // At most 2 digits remaining.

    DRAGONBOX_ASSERT(output >= 1);
    DRAGONBOX_ASSERT(output <= 99);

    if (output >= 10)
    {
        const uint32_t q = output;
        buf -= 2;
        Utoa_2Digits(buf, q);
        if (tz == nd)
        {
            tz += TrailingZeros_2Digits(q);
        }
//      nd += 2;
    }
    else
    {
        const uint32_t q = output;
        DRAGONBOX_ASSERT(q >= 1);
        DRAGONBOX_ASSERT(q <= 9);
        *--buf = static_cast<char>('0' + q);
    }

    return tz;
}

static inline int32_t DecimalLength(uint64_t v)
{
    DRAGONBOX_ASSERT(v >= 1);
    DRAGONBOX_ASSERT(v <= 99999999999999999ull);

    if (static_cast<uint32_t>(v >> 32) != 0)
    {
        if (v >= 10000000000000000ull) { return 17; }
        if (v >= 1000000000000000ull) { return 16; }
        if (v >= 100000000000000ull) { return 15; }
        if (v >= 10000000000000ull) { return 14; }
        if (v >= 1000000000000ull) { return 13; }
        if (v >= 100000000000ull) { return 12; }
        if (v >= 10000000000ull) { return 11; }
        return 10;
    }

    const uint32_t v32 = static_cast<uint32_t>(v);
    if (v32 >= 1000000000u) { return 10; }
    if (v32 >= 100000000u) { return 9; }
    if (v32 >= 10000000u) { return 8; }
    if (v32 >= 1000000u) { return 7; }
    if (v32 >= 100000u) { return 6; }
    if (v32 >= 10000u) { return 5; }
    if (v32 >= 1000u) { return 4; }
    if (v32 >= 100u) { return 3; }
    if (v32 >= 10u) { return 2; }
    return 1;
}

static inline char* FormatDigits(char* buffer, uint64_t digits, int32_t decimal_exponent, bool force_trailing_dot_zero = false)
{
    static constexpr int32_t MinFixedDecimalPoint = -6;
    static constexpr int32_t MaxFixedDecimalPoint =  17;
    static_assert(MinFixedDecimalPoint <= -1, "internal error");
    static_assert(MaxFixedDecimalPoint >= 17, "internal error");

    DRAGONBOX_ASSERT(digits >= 1);
    DRAGONBOX_ASSERT(digits <= 99999999999999999ull);
    DRAGONBOX_ASSERT(decimal_exponent >= -999);
    DRAGONBOX_ASSERT(decimal_exponent <=  999);

    int32_t num_digits = DecimalLength(digits);
    const int32_t decimal_point = num_digits + decimal_exponent;

    const bool use_fixed = MinFixedDecimalPoint <= decimal_point && decimal_point <= MaxFixedDecimalPoint;

    // Prepare the buffer.
    // Avoid calling memset/memcpy with variable arguments below...

    std::memset(buffer +  0, '0', 16);
    std::memset(buffer + 16, '0', 16);
    static_assert(MinFixedDecimalPoint >= -30, "internal error");
    static_assert(MaxFixedDecimalPoint <=  32, "internal error");

    int32_t decimal_digits_position;
    if (use_fixed)
    {
        if (decimal_point <= 0)
        {
            // 0.[000]digits
            decimal_digits_position = 2 - decimal_point;
        }
        else
        {
            // dig.its
            // digits[000]
            decimal_digits_position = 0;
        }
    }
    else
    {
        // dE+123 or d.igitsE+123
        decimal_digits_position = 1;
    }

    char* digits_end = buffer + decimal_digits_position + num_digits;

    const int tz = PrintDecimalDigitsBackwards(digits_end, digits);
    digits_end -= tz;
    num_digits -= tz;
//  decimal_exponent += tz; // => decimal_point unchanged.

    if (use_fixed)
    {
        if (decimal_point <= 0)
        {
            // 0.[000]digits
            buffer[1] = '.';
            buffer = digits_end;
        }
        else if (decimal_point < num_digits)
        {
            // dig.its
#if defined(_MSC_VER) && !defined(__clang__)
            // VC does not inline the memmove call below. (Even if compiled with /arch:AVX2.)
            // However, memcpy will be inlined.
            uint8_t tmp[16];
            char* const src = buffer + decimal_point;
            char* const dst = src + 1;
            std::memcpy(tmp, src, 16);
            std::memcpy(dst, tmp, 16);
#else
            std::memmove(buffer + decimal_point + 1, buffer + decimal_point, 16);
#endif
            buffer[decimal_point] = '.';
            buffer = digits_end + 1;
        }
        else
        {
            // digits[000]
            buffer += decimal_point;
            if (force_trailing_dot_zero)
            {
                std::memcpy(buffer, ".0", 2);
                buffer += 2;
            }
        }
    }
    else
    {
        // Copy the first digit one place to the left.
        buffer[0] = buffer[1];
        if (num_digits == 1)
        {
            // dE+123
            ++buffer;
        }
        else
        {
            // d.igitsE+123
            buffer[1] = '.';
            buffer = digits_end;
        }

        const int32_t scientific_exponent = decimal_point - 1;
//      SF_ASSERT(scientific_exponent != 0);

        std::memcpy(buffer, scientific_exponent < 0 ? "e-" : "e+", 2);
        buffer += 2;

        const uint32_t k = static_cast<uint32_t>(scientific_exponent < 0 ? -scientific_exponent : scientific_exponent);
        if (k < 10)
        {
            *buffer++ = static_cast<char>('0' + k);
        }
        else if (k < 100)
        {
            Utoa_2Digits(buffer, k);
            buffer += 2;
        }
        else
        {
            const uint32_t q = k / 100;
            const uint32_t r = k % 100;
            *buffer++ = static_cast<char>('0' + q);
            Utoa_2Digits(buffer, r);
            buffer += 2;
        }
    }

    return buffer;
}

static inline char* ToChars(char* buffer, double value, bool force_trailing_dot_zero = false)
{
    const Double v(value);

    if (!v.IsFinite())
    {
        if (v.IsNaN())
        {
            std::memcpy(buffer, "nan ", 4);
            return buffer + 3;
        }
        if (v.SignBit())
        {
            *buffer++ = '-';
        }
        std::memcpy(buffer, "inf ", 4);
        return buffer + 3;
    }

    if (v.SignBit())
    {
        value = -value;
        *buffer++ = '-';
    }

    if (v.IsZero())
    {
        std::memcpy(buffer, "0.0 ", 4);
        buffer += 1 + (force_trailing_dot_zero ? 2 : 0);
        return buffer;
    }

    const auto dec = jkj::dragonbox::to_decimal(
        value,
        jkj::dragonbox::policy::input_validation::do_nothing,
        jkj::dragonbox::policy::cache::normal,
        jkj::dragonbox::policy::sign::ignore,
        jkj::dragonbox::policy::rounding_mode::nearest_to_even,
        jkj::dragonbox::policy::trailing_zero::ignore);

    return FormatDigits(buffer, dec.significand, dec.exponent, force_trailing_dot_zero);
}

//==================================================================================================
//
//==================================================================================================

char* dragonbox::Dtoa(char* buffer, double value)
{
    return ToChars(buffer, value);
}

#endif
