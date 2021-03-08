// Copyright 2020 Ulf Adams
// Copyright 2020 Alexander Bolz
//
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#define RYU_STRTOD_FALLBACK() 1

namespace ryu {

// char* output_end = Ftoa(buffer, value);
//
// Converts the given single-precision number into decimal form and stores the result in the given
// buffer.
//
// The buffer must be large enough, i.e. >= FtoaMinBufferLength.
// The output format is similar to printf("%g").
// The output is _not_ null-terminted.
//
// The output is optimal, i.e. the output string
//  1. rounds back to the input number when read in (using round-to-nearest-even)
//  2. is as short as possible,
//  3. is as close to the input number as possible.
//
// Note:
// This function may temporarily write up to FtoaMinBufferLength characters into the buffer.

constexpr int FtoaMinBufferLength = 32;

char* Ftoa(char* buffer, float value);

// StrtofResult conversion_result = Strtof(first, last, value);
//
// Converts the given decimal floating-point number into a single-precision binary floating-point
// number.
// The function accepts the same inputs as std::strtof.
//
// If the input has more than 9 significant digits, the function may return
// StrtofStatus::input_too_long. In this case another algorithm must be used to convert the input
// string (e.g. std::strtof).
//
// Note:
// This function always succeeds to convert the output of Ftoa back into the correct binary
// floating-point number.

enum class StrtofStatus {
    invalid,
    integer,    // Add StrtofFormat ?
    fixed,      // Add StrtofFormat ?
    scientific, // Add StrtofFormat ?
    inf,
    nan,
#if !RYU_STRTOD_FALLBACK()
    input_too_long,
#endif
};

struct StrtofResult
{
    const char* next;
    StrtofStatus status;

    // Test for success.
    explicit operator bool() const noexcept
    {
#if !RYU_STRTOD_FALLBACK()
        return status != StrtofStatus::invalid && status != StrtofStatus::input_too_long;
#else
        return status != StrtofStatus::invalid;
#endif
    }
};

StrtofResult Strtof(const char* first, const char* last, float& value);

// Round10(x, n) returns: round(x * 10^-n) / 10^-n
//
// Use this function to round the given value to a specific number of decimal places.
// E.g.: Round10(1.005f, -2) == 1.01f
//       Round10(55.0f, 1) == 60.0f
float Round10(float value, int n);

} // namespace ryu
