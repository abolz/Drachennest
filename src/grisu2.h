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

namespace grisu2 {

// char* output_end = Dtoa(buffer, value);
//
// Converts the given double-precision number into decimal form and stores the result in the given
// buffer.
//
// The buffer must be large enough, i.e. >= DtoaMinBufferLength.
// The output format is similar to printf("%g").
// The output is _not_ null-terminted.
//
// Note:
// This function may temporarily write up to DtoaMinBufferLength characters into the buffer.

constexpr int DtoaMinBufferLength = 64;

char* Dtoa(char* buffer, double value);

} // namespace grisu2
