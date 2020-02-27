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

namespace charconv {

static constexpr int FtoaMinBufferLength = 32;

char* Ftoa(char* buffer, float value);

enum class StrtofStatus {
    invalid, // TODO: more detailed error code...
    zero,
    integer,
    decimal,
    nan,
    inf,
};

struct StrtofResult {
    const char* next;
    StrtofStatus status;
};

StrtofResult Strtof(const char* next, const char* last, float& value);

} // namespace charconv
