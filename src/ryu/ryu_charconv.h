// Copyright 2019 Ulf Adams
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

#include <cstdint>

//--------------------------------------------------------------------------------------------------
// RyuDtoa and RyuFtoa compute a decimal representation of the floating-point number 'value'
// in a format similar to printf %g.
//
// The result is optimal, i.e.
//  1. rounds back to the input number when read in,
//  2. is as short as possible,
//  3. is as close to the input number as possible.
//
// Note: The result is not null-terminated.
// Note: NaN's are formatted as "NaN".
// Note: +/-Infinity is formatted as "Infinity" and "-Infinity", resp.
//
// PRE: The buffer must be large enough, i.e.,
//      >= RyuDtoaMinBufferLength or RyuFtoaMinBufferLength, resp.
//--------------------------------------------------------------------------------------------------

static constexpr int RyuDtoaMinBufferLength = 64;
static constexpr int RyuFtoaMinBufferLength = 64;

char* RyuDtoa(char* buffer, double value);
char* RyuFtoa(char* buffer, float value);

//--------------------------------------------------------------------------------------------------
// RyuToBinary64 and RyuToBinary32 compute the closest binary representation of the decimal
// floating-point number m10 * 10^e10.
//
// PRE: m10 != 0
// PRE: m10len = DigitLength(m10) <= 17 and 9, resp.
// PRE: m10len + e10 must not overflow
//--------------------------------------------------------------------------------------------------

double RyuToBinary64(uint64_t m10, int m10len, int e10);
float  RyuToBinary32(uint32_t m10, int m10len, int e10);

enum class StrtodStatus {
    invalid, // TODO: more detailed error code...
    zero,
    integer,
    decimal,
    nan,
    inf,
};

struct StrtodResult {
    const char* next;
    StrtodStatus status;
};

StrtodResult RyuStrtod(const char* next, const char* last, double& value);
