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
#if _MSC_VER
#include <intrin.h>
#endif

#ifndef DTOA_ASSERT
#define DTOA_ASSERT(X) assert(X)
#endif

namespace dtoa {
namespace impl {

inline char* Utoa_2Digits(char* buf, uint32_t digits)
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

    DTOA_ASSERT(digits <= 99);
    std::memcpy(buf, &Digits100[2*digits], 2*sizeof(char));
    return buf + 2;
}

inline char* Utoa_4Digits(char* buf, uint32_t digits)
{
    DTOA_ASSERT(digits <= 9999);
    const uint32_t q = digits / 100;
    const uint32_t r = digits % 100;
    Utoa_2Digits(buf + 0, q);
    Utoa_2Digits(buf + 2, r);
    return buf + 4;
}

inline char* Utoa_8Digits(char* buf, uint32_t digits)
{
    DTOA_ASSERT(digits <= 99999999);
    const uint32_t q = digits / 10000;
    const uint32_t r = digits % 10000;
    Utoa_4Digits(buf + 0, q);
    Utoa_4Digits(buf + 4, r);
    return buf + 8;
}

#if 0 // Use intrinsics for DecimalLength

// Returns the number of leading 0-bits in x, starting at the most significant bit position.
// If x is 0, the result is undefined.
inline int CountLeadingZeros32(uint32_t x)
{
    DTOA_ASSERT(x != 0);

#if defined(__GNUC__)
    return __builtin_clz(x);
#elif defined(_MSC_VER) && (defined(_M_ARM) || defined(_M_ARM64))
    return static_cast<int>(_CountLeadingZeros(x));
#elif defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
    return static_cast<int>(__lzcnt(x));
#else
    int lz = 0;
    while ((x >> 31) == 0) {
        x <<= 1;
        ++lz;
    }
    return lz;
#endif
}

// Returns the number of leading 0-bits in x, starting at the most significant bit position.
// If x is 0, the result is undefined.
inline int CountLeadingZeros64(uint64_t x)
{
    DTOA_ASSERT(x != 0);

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

inline int DecimalLength(uint32_t v)
{
    DTOA_ASSERT(v >= 1);
    DTOA_ASSERT(v <= 999999999);

    static constexpr uint32_t Table[] = {
        0,
        9,
        99,
        999,
        9999,
        99999,
        999999,
        9999999,
        99999999,
        999999999,
    };

    constexpr int Bits = 32;
    const int y = ((19 * (Bits - 1) + (1 << 6)) - 19 * CountLeadingZeros32(v)) >> 6;
    return y + (Table[y] < v);
}

inline int DecimalLength(uint64_t v)
{
    DTOA_ASSERT(v >= 1);
    DTOA_ASSERT(v <= 99999999999999999ull); // 10^17 = 0x0163'4578'5D8A'0000

    static constexpr uint64_t Table[] = {
        0,
        9,
        99,
        999,
        9999,
        99999,
        999999,
        9999999,
        99999999,
        999999999,
        9999999999,
        99999999999,
        999999999999,
        9999999999999,
        99999999999999,
        999999999999999,
        9999999999999999,
        99999999999999999,
    };

    constexpr int Bits = 64;
    const int y = ((19 * (Bits - 1) + (1 << 6)) - 19 * CountLeadingZeros64(v)) >> 6;
    return y + (Table[y] < v);
}

inline void PrintDecimalDigits(char* buf, uint32_t digits, int num_digits)
{
    DTOA_ASSERT(digits >= 1);
    DTOA_ASSERT(digits <= 999999999);
    DTOA_ASSERT(num_digits >= 1);
    DTOA_ASSERT(num_digits <= 9);

    uint32_t n = digits;
    if (num_digits & 1)
    {
        const uint32_t q = n / 10;
        const uint32_t r = n % 10;
        n = q;
        num_digits--;
        buf[num_digits] = static_cast<char>('0' + r);
    }

    uint32_t q;
    uint32_t r;
    switch (num_digits)
    {
    case 8:
        q = n / 100;
        r = n % 100;
        n = q;
        num_digits -= 2;
        Utoa_2Digits(buf + num_digits, r);
        // fall through
    case 6:
        q = n / 100;
        r = n % 100;
        n = q;
        num_digits -= 2;
        Utoa_2Digits(buf + num_digits, r);
        // fall through
    case 4:
        q = n / 100;
        r = n % 100;
        n = q;
        num_digits -= 2;
        Utoa_2Digits(buf + num_digits, r);
        // fall through
    case 2:
        q = n / 100;
        r = n % 100;
//      n = q;
        num_digits -= 2;
        Utoa_2Digits(buf + num_digits, r);
    case 0:
        break;
    default:
#if __GNUC__
        __builtin_unreachable();
#elif _MSC_VER
        __assume(0);
#endif
    }
}

inline void PrintDecimalDigits(char* buf, uint64_t digits, int num_digits)
{
    DTOA_ASSERT(digits >= 1);
    DTOA_ASSERT(digits <= 99999999999999999ull);
    DTOA_ASSERT(num_digits >= 1);
    DTOA_ASSERT(num_digits <= 17);

    // We prefer 32-bit operations, even on 64-bit platforms.
    // We have at most 17 digits, and uint32_t can store 9 digits.
    // If output doesn't fit into uint32_t, we cut off 8 digits,
    // so the rest will fit into uint32_t.
    if (static_cast<uint32_t>(digits >> 32) != 0)
    {
        DTOA_ASSERT(num_digits >= 8);
        const uint64_t q = digits / 100000000;
        const uint64_t r = digits % 100000000;
        digits = q;
        num_digits -= 8;
        Utoa_8Digits(buf + num_digits, static_cast<uint32_t>(r));
    }

    DTOA_ASSERT(digits <= UINT32_MAX);
    DTOA_ASSERT(num_digits <= 10);

    uint32_t n = static_cast<uint32_t>(digits);
    if (num_digits & 1)
    {
        const uint32_t q = n / 10;
        const uint32_t r = n % 10;
        n = q;
        num_digits--;
        buf[num_digits] = static_cast<char>('0' + r);
    }

    uint32_t q;
    uint32_t r;
    switch (num_digits)
    {
    case 10:
        q = n / 100;
        r = n % 100;
        n = q;
        num_digits -= 2;
        Utoa_2Digits(buf + num_digits, r);
        // fall through
    case 8:
        q = n / 100;
        r = n % 100;
        n = q;
        num_digits -= 2;
        Utoa_2Digits(buf + num_digits, r);
        // fall through
    case 6:
        q = n / 100;
        r = n % 100;
        n = q;
        num_digits -= 2;
        Utoa_2Digits(buf + num_digits, r);
        // fall through
    case 4:
        q = n / 100;
        r = n % 100;
        n = q;
        num_digits -= 2;
        Utoa_2Digits(buf + num_digits, r);
        // fall through
    case 2:
        q = n / 100;
        r = n % 100;
        n = q;
        num_digits -= 2;
        Utoa_2Digits(buf + num_digits, r);
        // fall through
    case 0:
        break;
    default:
#if __GNUC__
        __builtin_unreachable();
#elif _MSC_VER
        __assume(0);
#endif
    }
}

#else

inline int DecimalLength(uint32_t v)
{
    DTOA_ASSERT(v >= 1);
    DTOA_ASSERT(v <= 999999999);

    if (v >= 100000000) { return 9; }
    if (v >= 10000000) { return 8; }
    if (v >= 1000000) { return 7; }
    if (v >= 100000) { return 6; }
    if (v >= 10000) { return 5; }
    if (v >= 1000) { return 4; }
    if (v >= 100) { return 3; }
    if (v >= 10) { return 2; }
    return 1;
}

inline int DecimalLength(uint64_t v)
{
    DTOA_ASSERT(v >= 1);
    DTOA_ASSERT(v <= 99999999999999999ull);

    if (v >= 10000000000000000ull) { return 17; }
    if (v >= 1000000000000000ull) { return 16; }
    if (v >= 100000000000000ull) { return 15; }
    if (v >= 10000000000000ull) { return 14; }
    if (v >= 1000000000000ull) { return 13; }
    if (v >= 100000000000ull) { return 12; }
    if (v >= 10000000000ull) { return 11; }
    if (v >= 1000000000ull) { return 10; }
    if (v >= 100000000ull) { return 9; }
    if (v >= 10000000ull) { return 8; }
    if (v >= 1000000ull) { return 7; }
    if (v >= 100000ull) { return 6; }
    if (v >= 10000ull) { return 5; }
    if (v >= 1000ull) { return 4; }
    if (v >= 100ull) { return 3; }
    if (v >= 10ull) { return 2; }
    return 1;
}

inline void PrintDecimalDigits(char* buf, uint32_t output, int output_length)
{
    while (output >= 10000)
    {
        DTOA_ASSERT(output_length > 4);
        const uint32_t q = output / 10000;
        const uint32_t r = output % 10000;
        output = q;
        output_length -= 4;
        Utoa_4Digits(buf + output_length, r);
    }

    if (output >= 100)
    {
        DTOA_ASSERT(output_length > 2);
        const uint32_t q = output / 100;
        const uint32_t r = output % 100;
        output = q;
        output_length -= 2;
        Utoa_2Digits(buf + output_length, r);
    }

    if (output >= 10)
    {
        DTOA_ASSERT(output_length == 2);
        Utoa_2Digits(buf, output);
    }
    else
    {
        DTOA_ASSERT(output_length == 1);
        buf[0] = static_cast<char>('0' + output);
    }
}

inline void PrintDecimalDigits(char* buf, uint64_t output, int output_length)
{
    // We prefer 32-bit operations, even on 64-bit platforms.
    // We have at most 17 digits, and uint32_t can store 9 digits.
    // If output doesn't fit into uint32_t, we cut off 8 digits,
    // so the rest will fit into uint32_t.
    if (static_cast<uint32_t>(output >> 32) != 0)
    {
        DTOA_ASSERT(output_length > 8);
        const uint64_t q = output / 100000000;
        const uint32_t r = static_cast<uint32_t>(output % 100000000);
        output = q;
        output_length -= 8;
        Utoa_8Digits(buf + output_length, r);
    }

    DTOA_ASSERT(output <= UINT32_MAX);
    uint32_t output2 = static_cast<uint32_t>(output);

    while (output2 >= 10000)
    {
        DTOA_ASSERT(output_length > 4);
        const uint32_t q = output2 / 10000;
        const uint32_t r = output2 % 10000;
        output2 = q;
        output_length -= 4;
        Utoa_4Digits(buf + output_length, r);
    }

    if (output2 >= 100)
    {
        DTOA_ASSERT(output_length > 2);
        const uint32_t q = output2 / 100;
        const uint32_t r = output2 % 100;
        output2 = q;
        output_length -= 2;
        Utoa_2Digits(buf + output_length, r);
    }

    if (output2 >= 10)
    {
        DTOA_ASSERT(output_length == 2);
        Utoa_2Digits(buf, output2);
    }
    else
    {
        DTOA_ASSERT(output_length == 1);
        buf[0] = static_cast<char>('0' + output2);
    }
}

#endif

} // namespace impl

// Print digits * 10^decimal_exponent in a form similar to printf("%g").
// PRE: sizeof(buffer) >= 32
template <typename UnsignedInt>
inline char* FormatDigits(char* buffer, UnsignedInt digits, int decimal_exponent, bool force_trailing_dot_zero = false)
{
    //
    // TODO:
    //
    // If the buffer were large enough (say > 64 bytes), the memset and memmove
    // call could use a compile-time constant, effectively removing them...
    //

    const int num_digits = dtoa::impl::DecimalLength(digits);
    const int decimal_point = num_digits + decimal_exponent;

    // NB:
    // These are the values used by JavaScript's ToString applied to Number
    // type. Printf uses the values -4 and max_digits10 resp. (sort of).
    constexpr int MinExp = -6; // -4;
    constexpr int MaxExp = 21; //  6;

    const bool use_fixed = MinExp < decimal_point && decimal_point <= MaxExp;

    char* first = buffer;
    if (use_fixed)
    {
        // Prepare the buffer.
        // Avoid calling memset with variable arguments below...
        // Need 21 '0's. Round up to a multiple of 8.
//      std::memset(buffer, '0', 24);

        if (decimal_point <= 0)
        {
            first += 2 + (-decimal_point);
        }
    }
    else
    {
        first += 1;
    }

    dtoa::impl::PrintDecimalDigits(first, digits, num_digits);

    if (use_fixed)
    {
        if (num_digits <= decimal_point)
        {
            // digits[000]
            // DTOA_ASSERT(buffer_capacity >= decimal_point + (force_trailing_dot_zero ? 2 : 0));

            //
            // TODO:
            //
            // decimal_point <= 21.
            // Pre-filling buffer with 21 '0's would therefore be sufficient.
            //

            std::memset(buffer + num_digits, '0', static_cast<unsigned>(decimal_point - num_digits));
            buffer += decimal_point;
            if (force_trailing_dot_zero)
            {
                *buffer++ = '.';
                *buffer++ = '0';
            }
        }
        else if (0 < decimal_point)
        {
            // dig.its
            // DTOA_ASSERT(buffer_capacity >= length + 1);

            //
            // TODO:
            //
            // 0 < decimal_point < num_digits
            //  ==> 1 <= num_digits - decimal_point <= 17 - 1 = 16
            // So we need to move at most 16 bytes one place to the right.
            // If we always copy 16 bytes, the buffer would need to be at least
            // 16 + 1 + 16 = 33 bytes large.
            //

            std::memmove(buffer + (decimal_point + 1), buffer + decimal_point, static_cast<unsigned>(num_digits - decimal_point));
            buffer[decimal_point] = '.';
            buffer += num_digits + 1;
        }
        else // decimal_point <= 0
        {
            // 0.[000]digits
            // DTOA_ASSERT(buffer_capacity >= 2 + (-decimal_point) + length);

            //
            // TODO:
            //
            // -5 <= decimal_point <= 0
            //  ==> 2 <= 2 + -decimal_point <= 7
            // Pre-filling buffer with 7 '0's would therefore be sufficient.
            //

            std::memset(buffer, '0', static_cast<unsigned>(2 + (-decimal_point)));
            buffer[1] = '.';
            buffer += (2 + (-decimal_point) + num_digits);
        }
    }
    else
    {
        // buffer = ?ddddd ==> d?dddd
        buffer[0] = buffer[1];

        if (num_digits == 1)
        {
            // dE+123
            // DTOA_ASSERT(buffer_capacity >= num_digits + 5);

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
            // DTOA_ASSERT(buffer_capacity >= num_digits + 1 + 5);

            buffer[1] = '.';
            buffer += 1 + num_digits;
        }

        int scientific_exponent = decimal_point - 1;
        *buffer++ = 'e';

        if (scientific_exponent < 0)
        {
            scientific_exponent = -scientific_exponent;
            *buffer++ = '-';
        }
        else
        {
             *buffer++ = '+';
        }

        const uint32_t k = static_cast<uint32_t>(scientific_exponent);
        if (k < 10)
        {
            *buffer++ = static_cast<char>('0' + k);
        }
        else if (k < 100)
        {
            buffer = dtoa::impl::Utoa_2Digits(buffer, k);
        }
        else
        {
            const uint32_t r = k % 10;
            const uint32_t q = k / 10;
            buffer = dtoa::impl::Utoa_2Digits(buffer, q);
            *buffer++ = static_cast<char>('0' + r);
        }
    }

    return buffer;
}

} // namespace dtoa
