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

//--------------------------------------------------------------------------------------------------
// RyuFtoa and RyuDtoa compute a decimal representation of the floating-point number 'value'
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
//      >= RyuFtoaMinBufferLength or RyuDtoaMinBufferLength, resp.
//--------------------------------------------------------------------------------------------------

static constexpr int RyuFtoaMinBufferLength = 64;
static constexpr int RyuDtoaMinBufferLength = 64;

char* RyuFtoa(char* buffer, float value);
char* RyuDtoa(char* buffer, double value);
