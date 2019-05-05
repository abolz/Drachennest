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
#include <limits>

namespace dtoa {

namespace impl {

template <typename Dest, typename Source>
inline Dest ReinterpretBits(Source source)
{
    static_assert(sizeof(Dest) == sizeof(Source), "size mismatch");

    Dest dest;
    std::memcpy(&dest, &source, sizeof(Source));
    return dest;
}

template <int Precision> struct BitsType;
template <> struct BitsType<24> { using type = uint32_t; };
template <> struct BitsType<53> { using type = uint64_t; };

} // namespace impl

template <typename Float>
struct IEEE
{
    // NB:
    // Works for double == long double.
    static_assert(std::numeric_limits<Float>::is_iec559 &&
                  ((std::numeric_limits<Float>::digits == 24 && std::numeric_limits<Float>::max_exponent == 128) ||
                   (std::numeric_limits<Float>::digits == 53 && std::numeric_limits<Float>::max_exponent == 1024)),
        "IEEE-754 single- or double-precision implementation required");

    using value_type = Float;
    using bits_type = typename dtoa::impl::BitsType<std::numeric_limits<Float>::digits>::type;

    static constexpr int       SignificandSize = std::numeric_limits<value_type>::digits; // = p   (includes the hidden bit)
    static constexpr int       ExponentBias    = std::numeric_limits<value_type>::max_exponent - 1 + (SignificandSize - 1);
    static constexpr int       MaxExponent     = std::numeric_limits<value_type>::max_exponent - 1 - (SignificandSize - 1);
    static constexpr int       MinExponent     = std::numeric_limits<value_type>::min_exponent - 1 - (SignificandSize - 1);
    static constexpr bits_type HiddenBit       = bits_type{1} << (SignificandSize - 1);   // = 2^(p-1)
    static constexpr bits_type SignificandMask = HiddenBit - 1;                           // = 2^(p-1) - 1
    static constexpr bits_type ExponentMask    = bits_type{2 * std::numeric_limits<value_type>::max_exponent - 1} << (SignificandSize - 1);
    static constexpr bits_type SignMask        = ~(~bits_type{0} >> 1);

    bits_type bits;

    explicit IEEE(bits_type bits_) : bits(bits_) {}
    explicit IEEE(value_type value) : bits(dtoa::impl::ReinterpretBits<bits_type>(value)) {}

    bits_type PhysicalSignificand() const {
        return bits & SignificandMask;
    }

    bits_type PhysicalExponent() const {
        return (bits & ExponentMask) >> (SignificandSize - 1);
    }

    // Returns the significand for a normalized double.
    bits_type NormalizedSignificand() const {
        return HiddenBit | PhysicalSignificand();
    }

    // Returns the exponent for a normalized double.
    int NormalizedExponent() const {
        return static_cast<int>(PhysicalExponent()) - ExponentBias;
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

    value_type Value() const {
        return dtoa::impl::ReinterpretBits<value_type>(bits);
    }

    value_type AbsValue() const {
        return dtoa::impl::ReinterpretBits<value_type>(bits & ~SignMask);
    }
};

} // namespace dtoa
