// Copyright 2020 Ulf Adams
// Copyright 2020 Alexander Bolz
//
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#define RYU_STRTOD_FALLBACK() 1

namespace ryu {

// char* output_end = Dtoa(buffer, value);
//
// Converts the given double-precision number into decimal form and stores the result in the given
// buffer.
//
// The buffer must be large enough, i.e. >= DtoaMinBufferLength.
// The output format is similar to printf("%g").
// The output is _not_ null-terminted.
//
// The output is optimal, i.e. the output string
//  1. rounds back to the input number when read in (using round-to-nearest-even)
//  2. is as short as possible,
//  3. is as close to the input number as possible.
//
// Note:
// This function may temporarily write up to DtoaMinBufferLength characters into the buffer.

constexpr int DtoaMinBufferLength = 64;

char* Dtoa(char* buffer, double value);

// StrtodResult conversion_result = Strtod(first, last, value);
//
// Converts the given decimal floating-point number into a double-precision binary floating-point
// number.
// The function accepts the same inputs as std::strtod.
//
// If the input has more than 17 significant digits, the function may return
// StrtodStatus::input_too_long. In this case another algorithm must be used to convert the input
// string (e.g. std::strtod).
//
// Note:
// This function always succeeds to convert the output of Dtoa back into the correct binary
// floating-point number.

enum class StrtodStatus {
    invalid,
    integer,
    floating_point,
    inf,
    nan,
#if !RYU_STRTOD_FALLBACK()
    input_too_long,
#endif
};

struct StrtodResult
{
    const char* next;
    StrtodStatus status;

    // Test for success.
    explicit operator bool() const noexcept
    {
#if !RYU_STRTOD_FALLBACK()
        return status != StrtodStatus::invalid && status != StrtodStatus::input_too_long;
#else
        return status != StrtodStatus::invalid;
#endif
    }
};

StrtodResult Strtod(const char* next, const char* last, double& value);

} // namespace ryu
