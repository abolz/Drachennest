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

#include "charconv_f32.h"

//#undef NDEBUG
#include <cassert>
#include <climits>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#if _MSC_VER
#include <intrin.h>
#endif

#ifndef RYU_ASSERT
#define RYU_ASSERT(X) assert(X)
#endif

#ifndef RYU_NEVER_INLINE
#if _MSC_VER
#define RYU_NEVER_INLINE __declspec(noinline) inline
#elif __GNUC__
#define RYU_NEVER_INLINE __attribute__((noinline)) inline
#else
#define RYU_NEVER_INLINE inline
#endif
#endif

#ifndef RYU_FORCE_INLINE
#if _MSC_VER
#define RYU_FORCE_INLINE __forceinline
#elif __GNUC__
#define RYU_FORCE_INLINE __attribute__((always_inline)) inline
#else
#define RYU_FORCE_INLINE inline
#endif
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
struct Single
{
    static_assert(std::numeric_limits<float>::is_iec559
               && std::numeric_limits<float>::digits == 24
               && std::numeric_limits<float>::max_exponent == 128,
        "IEEE-754 single-precision implementation required");

    using value_type = float;
    using bits_type = uint32_t;

//  static constexpr int       MaxDigits10     = std::numeric_limits<value_type>::max_digits10;
    static constexpr int       SignificandSize = std::numeric_limits<value_type>::digits; // = p   (includes the hidden bit)
    static constexpr int       ExponentBias    = std::numeric_limits<value_type>::max_exponent - 1 + (SignificandSize - 1);
//  static constexpr int       MaxExponent     = std::numeric_limits<value_type>::max_exponent - 1 - (SignificandSize - 1);
//  static constexpr int       MinExponent     = std::numeric_limits<value_type>::min_exponent - 1 - (SignificandSize - 1);
    static constexpr bits_type HiddenBit       = bits_type{1} << (SignificandSize - 1);   // = 2^(p-1)
    static constexpr bits_type SignificandMask = HiddenBit - 1;                           // = 2^(p-1) - 1
    static constexpr bits_type ExponentMask    = (bits_type{2 * std::numeric_limits<value_type>::max_exponent - 1}) << (SignificandSize - 1);
    static constexpr bits_type SignMask        = ~(~bits_type{0} >> 1);

    bits_type bits;

    explicit Single(bits_type bits_) : bits(bits_) {}
    explicit Single(value_type value) : bits(ReinterpretBits<bits_type>(value)) {}

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

    value_type Value() const {
        return ReinterpretBits<value_type>(bits);
    }

    value_type AbsValue() const {
        return ReinterpretBits<value_type>(bits & ~SignMask);
    }
};
} // namespace

//==================================================================================================
//
//==================================================================================================

static inline int Max(int x, int y)
{
    return y < x ? x : y;
}

// Returns floor(x / 2^n).
static inline int FloorDivPow2(int x, int n)
{
#if 1
    // Technically, right-shift of negative integers is implementation defined...
    // Should easily be optimized into SAR (or equivalent) instruction.
    return x < 0 ? ~(~x >> n) : (x >> n);
#else
    return x >> n;
#endif
}

static inline int FloorLog2Pow5(int e)
{
    RYU_ASSERT(e >= -1764);
    RYU_ASSERT(e <=  1763);
    return FloorDivPow2(e * 1217359, 19);
}

static inline int FloorLog10Pow2(int e)
{
    RYU_ASSERT(e >= -2620);
    RYU_ASSERT(e <=  2620);
    return FloorDivPow2(e * 315653, 20);
}

static inline int FloorLog10Pow5(int e)
{
    RYU_ASSERT(e >= -2620);
    RYU_ASSERT(e <=  2620);
    return FloorDivPow2(e * 732923, 20);
}

static inline uint32_t Lo32(uint64_t x)
{
    return static_cast<uint32_t>(x);
}

static inline uint32_t Hi32(uint64_t x)
{
    return static_cast<uint32_t>(x >> 32);
}

//==================================================================================================
// ToDecimal
//
// Single-precision implementation
//==================================================================================================
// Constant data: 808 (+ 104) bytes

static constexpr int BitsPerPow5_Single = 64;

static inline uint64_t ComputePow5_Single(int k)
{
    // Let e = FloorLog2Pow5(k) + 1 - 64
    // For k >= 0, stores 5^k in the form: ceil( 5^k / 2^e )
    // For k <= 0, stores 5^k in the form: ceil(2^-e / 5^-k)
    static constexpr int MinDecExp = -53;
    static constexpr int MaxDecExp =  47;
    static constexpr uint64_t Pow5[MaxDecExp - MinDecExp + 1] = {
        0xF53304714D9265E0, // e =  -187, k =  -53
        0x993FE2C6D07B7FAC, // e =  -184, k =  -52
        0xBF8FDB78849A5F97, // e =  -182, k =  -51
        0xEF73D256A5C0F77D, // e =  -180, k =  -50
        0x95A8637627989AAE, // e =  -177, k =  -49
        0xBB127C53B17EC15A, // e =  -175, k =  -48
        0xE9D71B689DDE71B0, // e =  -173, k =  -47
        0x9226712162AB070E, // e =  -170, k =  -46
        0xB6B00D69BB55C8D2, // e =  -168, k =  -45
        0xE45C10C42A2B3B06, // e =  -166, k =  -44
        0x8EB98A7A9A5B04E4, // e =  -163, k =  -43
        0xB267ED1940F1C61D, // e =  -161, k =  -42
        0xDF01E85F912E37A4, // e =  -159, k =  -41
        0x8B61313BBABCE2C7, // e =  -156, k =  -40
        0xAE397D8AA96C1B78, // e =  -154, k =  -39
        0xD9C7DCED53C72256, // e =  -152, k =  -38
        0x881CEA14545C7576, // e =  -149, k =  -37
        0xAA242499697392D3, // e =  -147, k =  -36
        0xD4AD2DBFC3D07788, // e =  -145, k =  -35
        0x84EC3C97DA624AB5, // e =  -142, k =  -34
        0xA6274BBDD0FADD62, // e =  -140, k =  -33
        0xCFB11EAD453994BB, // e =  -138, k =  -32
        0x81CEB32C4B43FCF5, // e =  -135, k =  -31
        0xA2425FF75E14FC32, // e =  -133, k =  -30
        0xCAD2F7F5359A3B3F, // e =  -131, k =  -29
        0xFD87B5F28300CA0E, // e =  -129, k =  -28
        0x9E74D1B791E07E49, // e =  -126, k =  -27
        0xC612062576589DDB, // e =  -124, k =  -26
        0xF79687AED3EEC552, // e =  -122, k =  -25
        0x9ABE14CD44753B53, // e =  -119, k =  -24
        0xC16D9A0095928A28, // e =  -117, k =  -23
        0xF1C90080BAF72CB2, // e =  -115, k =  -22
        0x971DA05074DA7BEF, // e =  -112, k =  -21
        0xBCE5086492111AEB, // e =  -110, k =  -20
        0xEC1E4A7DB69561A6, // e =  -108, k =  -19
        0x9392EE8E921D5D08, // e =  -105, k =  -18
        0xB877AA3236A4B44A, // e =  -103, k =  -17
        0xE69594BEC44DE15C, // e =  -101, k =  -16
        0x901D7CF73AB0ACDA, // e =   -98, k =  -15
        0xB424DC35095CD810, // e =   -96, k =  -14
        0xE12E13424BB40E14, // e =   -94, k =  -13
        0x8CBCCC096F5088CC, // e =   -91, k =  -12
        0xAFEBFF0BCB24AAFF, // e =   -89, k =  -11
        0xDBE6FECEBDEDD5BF, // e =   -87, k =  -10
        0x89705F4136B4A598, // e =   -84, k =   -9
        0xABCC77118461CEFD, // e =   -82, k =   -8
        0xD6BF94D5E57A42BD, // e =   -80, k =   -7
        0x8637BD05AF6C69B6, // e =   -77, k =   -6
        0xA7C5AC471B478424, // e =   -75, k =   -5
        0xD1B71758E219652C, // e =   -73, k =   -4
        0x83126E978D4FDF3C, // e =   -70, k =   -3
        0xA3D70A3D70A3D70B, // e =   -68, k =   -2
        0xCCCCCCCCCCCCCCCD, // e =   -66, k =   -1
        0x8000000000000000, // e =   -63, k =    0
        0xA000000000000000, // e =   -61, k =    1
        0xC800000000000000, // e =   -59, k =    2
        0xFA00000000000000, // e =   -57, k =    3
        0x9C40000000000000, // e =   -54, k =    4
        0xC350000000000000, // e =   -52, k =    5
        0xF424000000000000, // e =   -50, k =    6
        0x9896800000000000, // e =   -47, k =    7
        0xBEBC200000000000, // e =   -45, k =    8
        0xEE6B280000000000, // e =   -43, k =    9
        0x9502F90000000000, // e =   -40, k =   10
        0xBA43B74000000000, // e =   -38, k =   11
        0xE8D4A51000000000, // e =   -36, k =   12
        0x9184E72A00000000, // e =   -33, k =   13
        0xB5E620F480000000, // e =   -31, k =   14
        0xE35FA931A0000000, // e =   -29, k =   15
        0x8E1BC9BF04000000, // e =   -26, k =   16
        0xB1A2BC2EC5000000, // e =   -24, k =   17
        0xDE0B6B3A76400000, // e =   -22, k =   18
        0x8AC7230489E80000, // e =   -19, k =   19
        0xAD78EBC5AC620000, // e =   -17, k =   20
        0xD8D726B7177A8000, // e =   -15, k =   21
        0x878678326EAC9000, // e =   -12, k =   22
        0xA968163F0A57B400, // e =   -10, k =   23
        0xD3C21BCECCEDA100, // e =    -8, k =   24
        0x84595161401484A0, // e =    -5, k =   25
        0xA56FA5B99019A5C8, // e =    -3, k =   26
        0xCECB8F27F4200F3A, // e =    -1, k =   27
        0x813F3978F8940985, // e =     2, k =   28
        0xA18F07D736B90BE6, // e =     4, k =   29
        0xC9F2C9CD04674EDF, // e =     6, k =   30
        0xFC6F7C4045812297, // e =     8, k =   31
        0x9DC5ADA82B70B59E, // e =    11, k =   32
        0xC5371912364CE306, // e =    13, k =   33
        0xF684DF56C3E01BC7, // e =    15, k =   34
        0x9A130B963A6C115D, // e =    18, k =   35
        0xC097CE7BC90715B4, // e =    20, k =   36
        0xF0BDC21ABB48DB21, // e =    22, k =   37
        0x96769950B50D88F5, // e =    25, k =   38
        0xBC143FA4E250EB32, // e =    27, k =   39
        0xEB194F8E1AE525FE, // e =    29, k =   40
        0x92EFD1B8D0CF37BF, // e =    32, k =   41
        0xB7ABC627050305AE, // e =    34, k =   42
        0xE596B7B0C643C71A, // e =    36, k =   43
        0x8F7E32CE7BEA5C70, // e =    39, k =   44
        0xB35DBF821AE4F38C, // e =    41, k =   45
        0xE0352F62A19E306F, // e =    43, k =   46
        0x8C213D9DA502DE46, // e =    46, k =   47
    };

    RYU_ASSERT(k >= MinDecExp);
    RYU_ASSERT(k <= MaxDecExp);
    return Pow5[static_cast<unsigned>(k - MinDecExp)];
}

static inline uint64_t MulShift(uint32_t m, uint64_t mul, int j)
{
    RYU_ASSERT(j >= 32);
    RYU_ASSERT(j <= 95);

#if defined(__SIZEOF_INT128__)
    __extension__ using uint128_t = unsigned __int128;
    const uint64_t shifted_sum = static_cast<uint64_t>((uint128_t{mul} * m) >> j);
#elif defined(_MSC_VER) && defined(_M_X64)
    uint64_t hi;
    uint64_t lo = _umul128(m, mul, &hi);
    const uint64_t l = __shiftright128(lo, hi, static_cast<unsigned char>(j));
    const uint64_t h = __ull_rshift(hi, j); // Assume j >= 64: j - 64 == j % 64
    const uint64_t shifted_sum = (j & 64) ? h : l;
#else
    const uint64_t bits0 = uint64_t{m} * Lo32(mul);
    const uint64_t bits1 = uint64_t{m} * Hi32(mul);
    const uint64_t sum = bits1 + Hi32(bits0);

    const int shift = j - 32;
    const uint64_t shifted_sum = sum >> shift;
#endif

    return shifted_sum;
}

static inline void MulPow5DivPow2_Single(uint32_t u, uint32_t v, uint32_t w, int e5, int e2, uint64_t& a, uint64_t& b, uint64_t& c)
{
    // j >= 57 and m has at most 24 + 2 = 26 bits.
    // The product along with the subsequent shift therefore requires
    // 26 + 64 - 57 = 33 bits.

    const auto k = FloorLog2Pow5(e5) + 1 - BitsPerPow5_Single;
    const auto j = e2 - k;
    RYU_ASSERT(j >= BitsPerPow5_Single - 7); // 57
    RYU_ASSERT(j <= BitsPerPow5_Single - 1); // 63

    const auto pow5 = ComputePow5_Single(e5);

    a = MulShift(u, pow5, j);
    b = MulShift(v, pow5, j);
    c = MulShift(w, pow5, j);
}

// Returns whether value is divisible by 5^e5
static inline bool MultipleOfPow5(uint32_t value, int e5)
{
    RYU_ASSERT(e5 >= 0);
    RYU_ASSERT(e5 <= 12);

    struct MulCmp {
        uint32_t mul;
        uint32_t cmp;
    };

    static constexpr MulCmp Mod5[] = {
        {0x00000001u, 0xFFFFFFFFu}, // 5^0
        {0xCCCCCCCDu, 0x33333333u}, // 5^1
        {0xC28F5C29u, 0x0A3D70A3u}, // 5^2
        {0x26E978D5u, 0x020C49BAu}, // 5^3
        {0x3AFB7E91u, 0x0068DB8Bu}, // 5^4
        {0x0BCBE61Du, 0x0014F8B5u}, // 5^5
        {0x68C26139u, 0x000431BDu}, // 5^6
        {0xAE8D46A5u, 0x0000D6BFu}, // 5^7
        {0x22E90E21u, 0x00002AF3u}, // 5^8
        {0x3A2E9C6Du, 0x00000897u}, // 5^9
        {0x3ED61F49u, 0x000001B7u}, // 5^10
        {0x0C913975u, 0x00000057u}, // 5^11
        {0xCF503EB1u, 0x00000011u}, // 5^12
    };

    return value * Mod5[e5].mul <= Mod5[e5].cmp;
}

// Returns whether value is divisible by 2^e2
static inline bool MultipleOfPow2(uint32_t value, int e2)
{
    RYU_ASSERT(e2 >= 0);
    RYU_ASSERT(e2 <= 31);

    return (value & ((uint32_t{1} << e2) - 1)) == 0;
}

namespace {
struct ToDecimalResultSingle {
    uint32_t digits; // num_digits <= 9
    int exponent;
};
}

static inline ToDecimalResultSingle ToDecimal(float value)
{
    RYU_ASSERT(Single(value).IsFinite());
    RYU_ASSERT(value > 0);

    //
    // Step 1:
    // Decode the floating point number, and unify normalized and subnormal cases.
    //

    const Single ieee_value(value);

    // Decode bits into mantissa, and exponent.
    const uint32_t ieee_mantissa = ieee_value.PhysicalSignificand();
    const uint32_t ieee_exponent = ieee_value.PhysicalExponent();

    uint32_t m2;
    int e2;
    if (ieee_exponent == 0)
    {
        m2 = ieee_mantissa;
        e2 = 1 - Single::ExponentBias;
    }
    else
    {
        m2 = Single::HiddenBit | ieee_mantissa;
        e2 = static_cast<int>(ieee_exponent) - Single::ExponentBias;

        if /*unlikely*/ ((0 <= -e2 && -e2 < Single::SignificandSize) && MultipleOfPow2(m2, -e2))
        {
            // Since 2^23 <= m2 < 2^24 and 0 <= -e2 <= 23:
            //  1 <= value = m2 / 2^-e2 < 2^24.
            // Since m2 is divisible by 2^-e2, value is an integer.
            return {m2 >> -e2, 0};
        }
    }

    const bool is_even = (m2 % 2) == 0;
    const bool accept_lower = is_even;
    const bool accept_upper = is_even;

    //
    // Step 2:
    // Determine the interval of valid decimal representations.
    //

    const uint32_t lower_boundary_is_closer = (ieee_mantissa == 0 && ieee_exponent > 1);

    e2 -= 2;
    const uint32_t u = 4 * m2 - 2 + lower_boundary_is_closer;
    const uint32_t v = 4 * m2;
    const uint32_t w = 4 * m2 + 2;

    //
    // Step 3:
    // Convert to a decimal power base.
    //

    int e10;

    bool za = false; // a[0, ..., i-1] == 0
    bool zb = false; // b[0, ..., i-1] == 0
    bool zc = false; // c[0, ..., i-1] == 0

    if (e2 >= 0)
    {
        // We need
        //  (a,b,c) = (u,v,w) * 2^e2
        // and we need to remove at least q' = log_10(2^e2) digits from the
        // scaled values a,b,c, i.e. we want to compute
        //  (a,b,c) = (u,v,w) * 2^e2 / 10^(q')
        //          = (u,v,w) * 2^e2 / 10^(e10)
        //          = (u,v,w) * 5^(-e10) / 2^(e10 - e2)
        //
        // However, to correctly round the result we need to know the value of
        // the last removed digit. We therefore remove only q = q' - 1 digits in
        // the first step and make sure that we execute the loop below at least
        // once and determine the correct value of the last removed digit.

        const int q = FloorLog10Pow2(e2) - (e2 > 3); // == max(0, q' - 1)
        RYU_ASSERT(q >= 0);

        e10 = q;
        RYU_ASSERT(e10 >= 0);
        RYU_ASSERT(e10 - e2 <= 0);

        // Determine whether all the removed digits are 0.
        //
        // Z(x,e2,q) = (x * 2^e2) % 10^q == 0
        //           = p10(x * 2^e2) >= q
        //           = min(p2(x) + p2(e2), p5(x)) >= q
        //           = p2(x) + e2 >= q and p5(x) >= q
        //           = p5(x) >= q
        //           = x % 5^q == 0

        if (q <= 10) // 10 = floor(log_5(2^24))
        {
            za = MultipleOfPow5(u, q);
            zb = MultipleOfPow5(v, q);
            zc = MultipleOfPow5(w, q);
        }
    }
    else
    {
        // We need
        //  (a,b,c) = (u,v,w) * 2^e2 / 10^e2
        // and we need to remove at least q' = log_10(5^-e2) digits from the
        // scaled values a,b,c, i.e. we want to compute
        //  (a,b,c) = (u,v,w) * 2^e2 / 10^(e2 + q')
        //          = (u,v,w) * 2^e2 / 10^(e10),
        //          = (u,v,w) * 5^(-e10) / 2^(e10 - e2)

        const int q = FloorLog10Pow5(-e2) - (-e2 > 1); // == max(0, q' - 1)
        RYU_ASSERT(q >= 0);

        e10 = q + e2;
        RYU_ASSERT(e10 < 0);
        RYU_ASSERT(e10 - e2 >= 0);

        // Determine whether all the removed digits are 0.
        //
        // Z(x,e2,q) = (x * 5^-e2) % 10^q == 0
        //           = min(p2(x), p5(x) - e2) >= q
        //           = p2(x) >= q and p5(x) - e2 >= q
        //           = p2(x) >= q
        //           = x % 2^q == 0

        if (q <= Single::SignificandSize + 2)
        {
            za = MultipleOfPow2(u, q);
            zb = MultipleOfPow2(v, q);
            zc = MultipleOfPow2(w, q);
        }
    }

    uint64_t aq;
    uint64_t bq;
    uint64_t cq;
    MulPow5DivPow2_Single(u, v, w, -e10, e10 - e2, aq, bq, cq);

    //
    // Step 4:
    // Find the shortest decimal representation in the interval of valid representations.
    //

    cq -= !accept_upper && zc;

    // mask = 10^(number of digits removed),
    // i.e., (bq % mask) contains the actual digits removed from bq.
    // cq < 2^33 = 8'589'934'592,
    // and we will therefore remove at most 9 decimal digits, i.e. mask fits into an uint32_t.
    uint32_t mask = 1;

    // aq,bq,cq sometimes have 33 bits and we want to use 32-bit operations as much as
    // possible. In this case, we remove the first decimal digit and then use 32-bit
    // integers.

    uint32_t a = Lo32(aq);
    uint32_t b = Lo32(bq);
    uint32_t c = Lo32(cq);

    if (Hi32(cq) != 0)
    {
        RYU_ASSERT(aq / 10 < cq / 10);
        RYU_ASSERT(Hi32(aq / 2) == 0);
        RYU_ASSERT(Hi32(bq / 2) == 0);
        RYU_ASSERT(Hi32(cq / 2) == 0);

        mask = 10;
        a = Lo32(aq / 2) / 5; // = aq / 10
        b = Lo32(bq / 2) / 5; // = bq / 10
        c = Lo32(cq / 2) / 5; // = cq / 10
        ++e10;
    }

#if 0
    while (a / 100 < c / 100)
    {
        mask *= 100;
        a /= 100;
        b /= 100;
        c /= 100;
        e10 += 2;
    }
#else
    if (a / 100 < c / 100) // 2
    {
        mask *= 100;
        a /= 100;
        b /= 100;
        c /= 100;
        e10 += 2;
        if (a / 100 < c / 100) // 4
        {
            mask *= 100;
            a /= 100;
            b /= 100;
            c /= 100;
            e10 += 2;
            if (a / 100 < c / 100) // 6
            {
                mask *= 100;
                a /= 100;
                b /= 100;
                c /= 100;
                e10 += 2;
                if (a / 100 < c / 100) // 8
                {
                    mask *= 100;
                    a /= 100;
                    b /= 100;
                    c /= 100;
                    e10 += 2;
                }
            }
        }
    }
#endif

    if (a / 10 < c / 10)
    {
        mask *= 10;
        a /= 10;
        b /= 10;
//      c /= 10;
        ++e10;
    }

    if /*likely*/ (!za && !zb)
    {
        const uint32_t br = Lo32(bq) - b * mask; // Digits removed from bq
        const uint32_t half = mask / 2;

        b += (a == b || br >= half);
    }
    else
    {
        // za currently determines whether the first q removed digits were all
        // 0's. Still need to check whether the digits removed in the loop above
        // are all 0's.
        const bool can_use_lower = accept_lower && za && (Lo32(aq) - a * mask == 0);
        if (can_use_lower)
        {
            // If the loop is executed at least once, we have a == b == c when
            // the loop terminates.
            // We only remove 0's from a, so ar and za don't change.
            RYU_ASSERT(a != 0);
            for (;;)
            {
                const uint32_t q = a / 10;
                const uint32_t r = a - 10 * q; // = a % 10
                if (r != 0)
                    break;
                mask *= 10;
                a = q;
                b = q;
//              c = q;
                ++e10;
            }
        }

        const uint32_t br = Lo32(bq) - b * mask; // Digits removed from bq
        const uint32_t half = mask / 2;

        // A return value of b is valid if and only if a != b or za == true.
        // A return value of b + 1 is valid if and only if b + 1 <= c.
        const bool round_up = (a == b && !can_use_lower) // out of range
            || (br > half)
            || (br == half && (!zb || b % 2 != 0));

//      RYU_ASSERT(!round_up || b < c);
        b += round_up;
    }

    return {b, e10};
}

//==================================================================================================
// ToChars
//==================================================================================================

static inline char* Utoa_2Digits(char* buf, uint32_t digits)
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

    RYU_ASSERT(digits <= 99);
    std::memcpy(buf, &Digits100[2 * digits], 2 * sizeof(char));
    return buf + 2;
}

static inline char* Utoa_4Digits(char* buf, uint32_t digits)
{
    RYU_ASSERT(digits <= 9999);
    const uint32_t q = digits / 100;
    const uint32_t r = digits % 100;
    Utoa_2Digits(buf + 0, q);
    Utoa_2Digits(buf + 2, r);
    return buf + 4;
}

static inline int DecimalLength(uint32_t v)
{
    RYU_ASSERT(v >= 1);
    RYU_ASSERT(v <= 999999999);

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

static inline void PrintDecimalDigits(char* buf, uint32_t output, int output_length)
{
    while (output >= 10000)
    {
        RYU_ASSERT(output_length > 4);
        const uint32_t q = output / 10000;
        const uint32_t r = output % 10000;
        output = q;
        output_length -= 4;
        Utoa_4Digits(buf + output_length, r);
    }

    if (output >= 100)
    {
        RYU_ASSERT(output_length > 2);
        const uint32_t q = output / 100;
        const uint32_t r = output % 100;
        output = q;
        output_length -= 2;
        Utoa_2Digits(buf + output_length, r);
    }

    if (output >= 10)
    {
        RYU_ASSERT(output_length == 2);
        Utoa_2Digits(buf, output);
    }
    else
    {
        RYU_ASSERT(output_length == 1);
        buf[0] = static_cast<char>('0' + output);
    }
}

static inline char* FormatDigits(char* buffer, uint32_t digits, int decimal_exponent, bool force_trailing_dot_zero = false)
{
    RYU_ASSERT(digits >= 1);
    RYU_ASSERT(digits <= 999999999u);
    RYU_ASSERT(decimal_exponent >= -99);
    RYU_ASSERT(decimal_exponent <=  99);

    const int num_digits = DecimalLength(digits);
    const int decimal_point = num_digits + decimal_exponent;

    // In order to successfully parse all numbers output by Dtoa using the Strtod implementation
    // below, we have to make sure to never emit more than 9 (significant) digits.
    constexpr int MaxFixedDecimalPoint =  9;
    constexpr int MinFixedDecimalPoint = -4;

    const bool use_fixed = MinFixedDecimalPoint <= decimal_point && decimal_point <= MaxFixedDecimalPoint;

    // Prepare the buffer.
    // Avoid calling memset/memcpy with variable arguments below...

    int decimal_digits_position;
    if (use_fixed)
    {
        if (decimal_point <= 0)
        {
            // 0.[000]digits
            // -4 <= decimal_point <= 0
            //  ==> 2 <= 2 + -decimal_point <= 6
            // Pre-filling the buffer with 6 '0's is therefore sufficient.
            std::memset(buffer, '0', 8);
            decimal_digits_position = 2 + (-decimal_point);
        }
        else if (decimal_point < num_digits)
        {
            // dig.its
            // 0 < decimal_point <= Min(9 - 1, MaxExp)
            // We need to move at most 8 bytes to the right.
            decimal_digits_position = 0;
        }
        else
        {
            // digits[000]
            // 1 <= num_digits <= 9 <= decimal_point <= 13.
            // Pre-filling buffer with 13 '0's is therefore sufficient.
            static_assert(MaxFixedDecimalPoint <= 16, "invalid parameter");
            std::memset(buffer, '0', 16);
            decimal_digits_position = 0;
        }
    }
    else
    {
        // dE+123 or d.igitsE+123
        // We only need to copy the first digit one position to the left.
        decimal_digits_position = 1;
    }

    PrintDecimalDigits(buffer + decimal_digits_position, digits, num_digits);

    if (use_fixed)
    {
        if (decimal_point <= 0)
        {
            // 0.[000]digits
            buffer[1] = '.';
            buffer += 2 + (-decimal_point) + num_digits;
        }
        else if (decimal_point < num_digits)
        {
            // dig.its
            // We need to move at most 8 bytes one place to the right.
            std::memmove(buffer + (decimal_point + 1), buffer + decimal_point, 8);
            buffer[decimal_point] = '.';
            buffer += num_digits + 1;
        }
        else // 0 < num_digits <= decimal_point
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
            buffer += 1;
        }
        else
        {
            // d.igitsE+123
            buffer[1] = '.';
            buffer += 1 + num_digits;
        }

        const auto scientific_exponent = decimal_point - 1;
//      RYU_ASSERT(scientific_exponent != 0);

        std::memcpy(buffer, scientific_exponent < 0 ? "e-" : "e+", 2);
        buffer += 2;
        const uint32_t k = static_cast<uint32_t>(scientific_exponent < 0 ? -scientific_exponent : scientific_exponent);
        if (k < 10)
        {
            *buffer++ = static_cast<char>('0' + k);
        }
        else
        {
            buffer = Utoa_2Digits(buffer, k);
        }
    }

    return buffer;
}

static inline char* ToChars(char* buffer, float value, bool force_trailing_dot_zero = false)
{
    const Single v(value);

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
        value = v.AbsValue();
        *buffer++ = '-';
    }

    if (v.IsZero())
    {
        std::memcpy(buffer, "0.0 ", 4);
        buffer += 1 + (force_trailing_dot_zero ? 2 : 0);
        return buffer;
    }

    const auto dec = ToDecimal(value);
    return FormatDigits(buffer, dec.digits, dec.exponent, force_trailing_dot_zero);
}

//==================================================================================================
//
//==================================================================================================

char* charconv::Ftoa(char* buffer, float value)
{
    return ToChars(buffer, value);
}

//==================================================================================================
// ToBinary32
//==================================================================================================

// Maximum number of decimal digits in the significand the fast ToBinary method can handle.
// Inputs with more significant digits must be processed using another algorithm.
static constexpr int ToBinaryMaxDecimalDigits = 9;

// Any input <= 10^MinDecimalExponent is interpreted as 0.
// Any input >  10^MaxDecimalExponent is interpreted as +Infinity.
static constexpr int MinDecimalExponent = -46; // denorm_min / 2 =  7.00649232e-46 >=  1 * 10^-46
static constexpr int MaxDecimalExponent =  39; //            max = 3.402823466e+38 <= 10 * 10^+38

static inline int FloorLog2(uint32_t x)
{
    RYU_ASSERT(x != 0);

#if defined(__GNUC__) || defined(__clang__)
    return 31 - __builtin_clz(x);
#elif defined(_MSC_VER) && defined(_M_X64)
    unsigned long index;
    _BitScanReverse(&index, x);
    return static_cast<int>(index);
#else
    int l2 = 0;
    for (;;)
    {
        x >>= 1;
        if (x == 0)
            break;
        ++l2;
    }
    return l2;
#endif
}

static inline int FloorLog2Pow10(int e)
{
    RYU_ASSERT(e >= -1233);
    RYU_ASSERT(e <= 1233);
    return FloorDivPow2(e * 1741647, 19);
}

static inline int ExtractBit(uint32_t x, int n)
{
    RYU_ASSERT(n >= 0);
    RYU_ASSERT(n <= 31);
    return (x & (uint32_t{1} << n)) != 0;
}

static inline float ToBinary32(uint32_t m10, int m10_digits, int e10)
{
    static constexpr int MantissaBits = Single::SignificandSize - 1;
    static constexpr int ExponentBias = Single::ExponentBias - (Single::SignificandSize - 1);

    RYU_ASSERT(m10 > 0);
    RYU_ASSERT(m10_digits == DecimalLength(m10));
    RYU_ASSERT(m10_digits <= ToBinaryMaxDecimalDigits);
    RYU_ASSERT(e10 >  MinDecimalExponent - m10_digits);
    RYU_ASSERT(e10 <= MaxDecimalExponent - m10_digits);
    static_cast<void>(m10_digits);

    // Convert to binary float m2 * 2^e2, while retaining information about whether the conversion
    // was exact.

    const auto log2_m10 = FloorLog2(m10);
    RYU_ASSERT(log2_m10 >= 0);
    RYU_ASSERT(log2_m10 <= 29); // 29 = floor(log_2(10^9))

    // The length of m10 * 10^e10 in bits is: log2(m10 * 10^e10) = log2(m10) + log2(10^e10).
    // We want to compute the (MantissaBits + 1) top-most bits (+1 for the implicit leading
    // one in IEEE format). We therefore choose a binary output exponent of
    //   e2 = log2(m10 * 10^e10) - (MantissaBits + 1).
    //
    // We use floor(log2(5^e10)) so that we get at least this many bits; better to have an
    // additional bit than to not have enough bits.

    // We compute [m10 * 10^e10 / 2^e2] == [m10 * 5^e10 / 2^(e2 - e10)]
    //
    // Let b = floor(log_2(m10))
    // Let n = floor(log_2(5^e10))
    // Then
    //  j = ( e2 - e10 ) - ( n + 1 - BitsPerPow5 )
    //    = ( ( b + e10 + n - (MantissaBits + 1) ) - e10 ) - ( n + 1 - BitsPerPow5 )
    //    = b + BitsPerPow5 - MantissaBits - 2
    //    = b + 64 - 23 - 2
    //    = b + 39
    // Since 0 <= b <= 29, we have
    //    39 <= j <= 68
    // The product along with the subsequent shift therefore has (at most)
    //  b + 64 - (64 - 25 + b) = 25
    // bits.

    const auto log2_10_e10 = FloorLog2Pow10(e10);
    const auto e2 = log2_m10 + log2_10_e10 - (MantissaBits + 1);

    const auto pow5 = ComputePow5_Single(e10);
    const auto j = log2_m10 + (BitsPerPow5_Single - MantissaBits - 2);
    RYU_ASSERT(j >= 39);
    RYU_ASSERT(j <= 68);
    const auto product = MulShift(m10, pow5, j);
    RYU_ASSERT(product <= UINT32_MAX);
    const auto m2 = static_cast<uint32_t>(product);

    const auto log2_m2 = FloorLog2(m2);
    RYU_ASSERT(log2_m2 >= 24);
    RYU_ASSERT(log2_m2 <= 25);

    // We also compute if the result is exact, i.e., [m10 * 10^e10 / 2^e2] == m10 * 10^e10 / 2^e2.
    //  (See: Ryu Revisited, Section 4.3)

    bool is_exact = (e2 <= e10) || (e2 - e10 < 32 && MultipleOfPow2(m10, e2 - e10));
    if (e10 >= 0)
    {
        // 2^(e2 - e10) | m10 5^e10
        //  <==> p2(m10 5^e10)       >= e2 - e10
        //  <==> p2(m10) + e10 p2(5) >= e2 - e10
        //  <==> p2(m10)             >= e2 - e10

        // is_exact
        //  <==>   (e2 <= e10   OR   p2(m10) >= e2 - e10)

        // is_exact = (e2 <= e10) || (e2 - e10 < 32 && MultipleOfPow2(m10, e2 - e10));
    }
    else
    {
        // e2 <= e10:
        //
        // m10 10^e10 / 2^e2
        //  == m10 2^e10 5^e10 / 2^e2
        //  == m10 2^(e10 - e2) / 5^(-e10)
        //
        // 5^(-e10) | m10 2^(e10 - e2)
        //  <==> p5(m10 2^(e10 - e2))       >= -e10
        //  <==> p5(m10) + (e10 - e2) p5(2) >= -e10
        //  <==> p5(m10)                    >= -e10

        // e2 > e10:
        //
        // m10 10^e10 / 2^e2
        //  == m10 (2^e10 5^e10) / 2^e2
        //  == m10 / (5^(-e10) 2^(e2 - e10))
        //  == m10 / (10^(-e10) 2^e2)
        //
        // 5^(-e10) 2^(e2 - e10) | m10
        //  <==> 5^(-e10) | m10   AND   2^(e2 - e10) | m10
        //  <==> p5(m10) >= -e10   AND   p2(m10) >= e2 - e10

        // is_exact
        //  <==>   (e2 <= e10   OR   p2(m10) >= e2 - e10)   AND   p5(m10) >= -e10

        // e2 <= e10 ==> is_exact = true
        // In this case we need to check p5(m10) >= -e10.
        // Check that the test below works.
        RYU_ASSERT(e2 > e10 || is_exact);

        // 30 = ceil(log_2(10^9))
        // 12 = floor(log_5(2^30))
        is_exact = is_exact && (-e10 <= 12 && MultipleOfPow5(m10, -e10));
    }

    // Compute the final IEEE exponent.
    int ieee_e2 = Max(0, log2_m2 + e2 + ExponentBias);
    if (ieee_e2 >= 2 * std::numeric_limits<float>::max_exponent - 1)
    {
        // Overflow:
        // Final IEEE exponent is larger than the maximum representable.
        return std::numeric_limits<float>::infinity();
    }

    // We need to figure out how much we need to shift m2.
    // The tricky part is that we need to take the final IEEE exponent into account, so we need to
    // reverse the bias and also special-case the value 0.
    const int shift = (ieee_e2 == 0 ? 1 : ieee_e2) - e2 - (ExponentBias + MantissaBits);
    RYU_ASSERT(shift > 0);

    // We need to round up if the exact value is more than 0.5 above the value we computed. That's
    // equivalent to checking if the last removed bit was 1 and either the value was not just
    // trailing zeros or the result would otherwise be odd.
    const auto trailing_zeros
        = is_exact && MultipleOfPow2(m2, shift - 1);
    const auto last_removed_bit
        = ExtractBit(m2, shift - 1);
    const auto round_up
        = last_removed_bit != 0 && (!trailing_zeros || ExtractBit(m2, shift) != 0);

    uint32_t significand = (m2 >> shift) + round_up;
    RYU_ASSERT(significand <= 2 * Single::HiddenBit); // significand <= 2^(p+1) = 2^25

    significand &= Single::SignificandMask;

    // Rounding up may cause overflow...
    if (significand == 0 && round_up)
    {
        // Rounding up did overflow the p-bit significand.
        // Move a trailing zero of the significand into the exponent.
        // Due to how the IEEE represents +/-Infinity, we don't need to check for overflow here.
        ++ieee_e2;
    }

    RYU_ASSERT(ieee_e2 <= 2 * std::numeric_limits<float>::max_exponent - 1);
    const uint32_t ieee = static_cast<uint32_t>(ieee_e2) << MantissaBits | significand;
    return ReinterpretBits<float>(ieee);
}

//==================================================================================================
// Strtof
//==================================================================================================

#define RYU_ASSUME_NULL_TERMINATED_INPUT() 0

using charconv::StrtofStatus;
using charconv::StrtofResult;

static inline bool IsDigit(char ch)
{
    return static_cast<unsigned>(ch - '0') <= 9u;
}

static inline int DigitValue(char ch)
{
    RYU_ASSERT(IsDigit(ch));
    return ch - '0';
}

static inline bool IsLowerASCII(char ch)
{
    return 'a' <= ch && ch <= 'z';
}

static inline bool IsUpperASCII(char ch)
{
    return 'A' <= ch && ch <= 'Z';
}

static inline char ToLowerASCII(char ch)
{
//  return IsUpperASCII(ch) ? static_cast<char>(ch - 'A' + 'a') : ch;
    return static_cast<char>(static_cast<unsigned char>(ch) | 0x20);
}

static inline bool StartsWith(const char* next, const char* last, const char* lower_case_prefix)
{
    for ( ; next != last && *lower_case_prefix != '\0'; ++next, ++lower_case_prefix)
    {
        RYU_ASSERT(IsLowerASCII(*lower_case_prefix));
        if (ToLowerASCII(*next) != *lower_case_prefix)
            return false;
    }

    return *lower_case_prefix == '\0';
}

static inline StrtofResult ParseInfinity(const char* next, const char* last)
{
    RYU_ASSERT(*next == 'i' || *next == 'I');

    if (!StartsWith(next + 1, last, "nf"))
        return {next, StrtofStatus::invalid};

    next += 3;
    if (StartsWith(next, last, "inity"))
        next += 5;

    return {next, StrtofStatus::ok};
}

static inline bool IsNaNSequenceChar(char ch)
{
    return ch == '_' || IsDigit(ch) || IsUpperASCII(ch) || IsLowerASCII(ch);
}

// FIXME:
// Don't ignore the nan-sequence!!!
static inline StrtofResult ParseNaN(const char* next, const char* last)
{
    RYU_ASSERT(*next == 'n' || *next == 'N');

    if (!StartsWith(next + 1, last, "an"))
        return {next, StrtofStatus::invalid};

    next += 3;
    if (next != last && *next == '(')
    {
        for (const char* p = next + 1; p != last; ++p)
        {
            if (*p == ')')
                return {p + 1, StrtofStatus::ok};

            if (!IsNaNSequenceChar(*p))
                break; // invalid/incomplete nan-sequence
        }
    }

    return {next, StrtofStatus::ok};
}

static RYU_NEVER_INLINE StrtofResult ParseSpecial(bool is_negative, const char* next, const char* last, float& value)
{
    if (*next == 'i' || *next == 'I')
    {
        const auto res = ParseInfinity(next, last);
        if (res.status != StrtofStatus::invalid)
        {
            value = is_negative ? -std::numeric_limits<float>::infinity() : std::numeric_limits<float>::infinity();
        }
        return res;
    }

    if (*next == 'n' || *next == 'N')
    {
        const auto res = ParseNaN(next, last);
        if (res.status != StrtofStatus::invalid)
        {
            value = std::numeric_limits<float>::quiet_NaN();
        }
        return res;
    }

    return {next, StrtofStatus::invalid};
}

static RYU_NEVER_INLINE float ToBinarySlow(const char* next, const char* last)
{
    //
    // FIXME:
    // _strtof_l( ..., C_LOCALE )
    //

#if RYU_ASSUME_NULL_TERMINATED_INPUT()
    const char* const ptr = next;
#else
    // std::strtod expects null-terminated inputs. So we need to make a copy and null-terminate the input.
    // This function is actually almost never going to be called, so that should be ok.
    const std::string inp(next, last);

    const char* const ptr = inp.c_str();
#endif
    char* end;
    const auto flt = ::strtof(ptr, &end);

    // std::strtod should have consumed all of the input.
    RYU_ASSERT(ptr != end);
    RYU_ASSERT(last - next == end - ptr);

    return flt;
}

StrtofResult charconv::Strtof(const char* next, const char* last, float& value)
{
    if (next == last)
        return {next, StrtofStatus::invalid};

    // Decompose the input into the form significand * 10^exponent,
    // where significand has num_digits decimal digits.

    uint32_t significand = 0; // only valid iff num_digits <= 9
    int64_t  num_digits  = 0; // 64-bit to avoid overflow...
    int64_t  exponent    = 0; // 64-bit to avoid overflow...

// [-]

    const bool is_negative = (*next == '-');
    if (is_negative || *next == '+')
    {
        ++next;
        if (next == last)
            return {next, StrtofStatus::invalid};
    }

// int

    const char* const start = next;

    const bool has_leading_zero = (*next == '0');
    const bool has_leading_dot  = (*next == '.');

    if (has_leading_zero)
    {
        for (;;)
        {
            ++next;
            if (next == last || *next != '0')
                break;
        }
    }

    if (next != last && IsDigit(*next)) // non-0
    {
        const char* const p = next;

        significand = DigitValue(*next);
        ++next;
        while (next != last && IsDigit(*next))
        {
            significand = 10 * significand + DigitValue(*next);
            ++next;
        }

        num_digits = next - p;
    }
    else if (!has_leading_zero && !has_leading_dot)
    {
        return ParseSpecial(is_negative, next, last, value);
    }

// frac

    if (has_leading_dot || (next != last && *next == '.'))
    {
        ++next; // skip '.'
        if (next != last && IsDigit(*next))
        {
            // TODO:
            // Ignore trailing zeros...

            const char* const p = next;

            significand = 10 * significand + DigitValue(*next);
            ++next;
            while (next != last && IsDigit(*next))
            {
                significand = 10 * significand + DigitValue(*next);
                ++next;
            }

            const char* nz = p;
            if (num_digits == 0)
            {
                // Number is of the form "0.xxx...".
                // Move the leading zeros in the fractional part into the exponent.
                while (nz != next && *nz == '0')
                    ++nz;
            }

            num_digits += next - nz;
            exponent = -(next - p);
        }
        else
        {
            // No digits in the fractional part.
            // But at least one digit must appear in either the integral or the fractional part.
            if (has_leading_dot)
                return {next, StrtofStatus::invalid};
        }
    }

// exp

    // Exponents larger than this limit will be treated as +Infinity.
    // But we must still scan all the digits if this happens to be the case.
    static constexpr int MaxExp = 999999;
    static_assert(MaxExp >= 999, "invalid parameter");
    static_assert(MaxExp <= (INT_MAX - 9) / 10, "invalid parameter");

    int parsed_exponent = 0;
    if (next != last && (*next == 'e' || *next == 'E'))
    {
        // Possibly the start of an exponent...
        // We accept (and ignore!) invalid or incomplete exponents.
        // The 'next' pointer is updated if and only if a valid exponent has been found.
        const char* p = next;

        ++p; // skip 'e' or 'E'
        if (p != last)
        {
            const bool parsed_exponent_is_negative = (*p == '-');
            if (parsed_exponent_is_negative || *p == '+')
                ++p;

            if (p != last && IsDigit(*p))
            {
                next = p; // Found a valid exponent.

                parsed_exponent = DigitValue(*next);
                ++next;
                if (next != last && IsDigit(*next))
                {
                    parsed_exponent = 10 * parsed_exponent + DigitValue(*next);
                    ++next;
                }
                while (next != last && IsDigit(*next))
                {
                    if (parsed_exponent <= MaxExp)
                        parsed_exponent = 10 * parsed_exponent + DigitValue(*next);
                    ++next;
                }

                parsed_exponent = parsed_exponent_is_negative ? -parsed_exponent : parsed_exponent;

                // (Assume overflow does not happen here...)
                exponent += parsed_exponent;
            }
        }
    }

    RYU_ASSERT(num_digits >= 0);

    float flt;
    if (num_digits == 0)
    {
        flt = 0;
    }
    else if (parsed_exponent < -MaxExp || exponent + num_digits <= MinDecimalExponent)
    {
        // input = x * 10^-inf = 0
        // or
        // input < 10^MinDecimalExponent, which rounds to +-0.
        flt = 0;
    }
    else if (parsed_exponent > +MaxExp || exponent + num_digits > MaxDecimalExponent)
    {
        // input = x * 10^+inf = +inf
        // or
        // input >= 10^MaxDecimalExponent, which rounds to +-infinity.
        flt = std::numeric_limits<float>::infinity();
    }
    else if (num_digits <= ToBinaryMaxDecimalDigits)
    {
        RYU_ASSERT(exponent >= INT_MIN);
        RYU_ASSERT(exponent <= INT_MAX);
        flt = ToBinary32(significand, static_cast<int>(num_digits), static_cast<int>(exponent));
    }
    else
    {
        // We need to fall back to another algorithm if the input is too long.
        flt = ToBinarySlow(start, next);
    }

    value = is_negative ? -flt : flt;
    return {next, StrtofStatus::ok};
}
