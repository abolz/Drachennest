// Copyright 2019 Alexander Bolz
//
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#include "grisu2.h"

#define GRISU_SMALL_INT_OPTIMIZATION() 1
#define GRISU_ROUND() 0

#include <cassert>
#include <cstdint>
#include <cstring>
#include <limits>
#ifdef _MSC_VER
#include <intrin.h>
#endif

#ifndef GRISU_ASSERT
#define GRISU_ASSERT(X) assert(X)
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
    static constexpr bits_type HiddenBit       = bits_type{1} << (SignificandSize - 1);   // = 2^(p-1)
    static constexpr bits_type SignificandMask = HiddenBit - 1;                           // = 2^(p-1) - 1
    static constexpr bits_type ExponentMask    = (bits_type{2 * std::numeric_limits<value_type>::max_exponent - 1}) << (SignificandSize - 1);
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

// Returns: floor(x / 2^n)
static inline int32_t SAR(int32_t x, int32_t n)
{
    // Technically, right-shift of negative integers is implementation defined...
    // Should easily get optimized into SAR (or equivalent) instruction.
#if 0
    return x < 0 ? ~(~x >> n) : (x >> n);
#else
    return x >> n;
#endif
}

// Returns: floor(log_2(10^e))
static inline int32_t FloorLog2Pow10(int32_t e)
{
    GRISU_ASSERT(e >= -1233);
    GRISU_ASSERT(e <=  1232);
    return SAR(e * 1741647, 19);
}

// Returns: ceil(log_10(2^e))
static inline int32_t CeilLog10Pow2(int32_t e)
{
    GRISU_ASSERT(e >= -2620);
    GRISU_ASSERT(e <=  2620);
    return SAR(e * 315653 + ((1 << 20) - 1), 20);
}

//==================================================================================================
// Grisu2
//
// Implements the Grisu2 algorithm for (IEEE) binary to decimal floating-point conversion.
//
// References:
//
// [1]  Loitsch, "Printing Floating-Point Numbers Quickly and Accurately with Integers",
//      Proceedings of the ACM SIGPLAN 2010 Conference on Programming Language Design and Implementation, PLDI 2010
// [2]  Burger, Dybvig, "Printing Floating-Point Numbers Quickly and Accurately",
//      Proceedings of the ACM SIGPLAN 1996 Conference on Programming Language Design and Implementation, PLDI 1996
//==================================================================================================
// Constant data: 159 * 64 bits = 1272 bytes

//
// TODO:
// Clean up comments...
//

namespace {
struct DiyFp // f * 2^e
{
    static constexpr int32_t SignificandSize = 64; // = q

    uint64_t f = 0;
    int32_t e = 0;

    constexpr DiyFp() = default;
    constexpr DiyFp(uint64_t f_, int32_t e_) : f(f_), e(e_) {}
};
}

// Returns x * y.
// The result is rounded (ties up). (Only the upper q bits are returned.)
static inline uint64_t MultiplyHighRoundUp(uint64_t x, uint64_t y)
{
    // Computes:
    //  f = round((x.f * y.f) / 2^q)

#if defined(__SIZEOF_INT128__)
    __extension__ using uint128_t = unsigned __int128;

    const uint128_t p = uint128_t{x} * y;

    uint64_t h = static_cast<uint64_t>(p >> 64);
    uint64_t l = static_cast<uint64_t>(p);
    h += l >> 63; // round, ties up: [h, l] += 2^q / 2

    return h;
#elif defined(_MSC_VER) && defined(_M_X64)
    uint64_t h = 0;
    uint64_t l = _umul128(x, y, &h);
    h += l >> 63; // round, ties up: [h, l] += 2^q / 2

    return h;
#else
    const uint32_t x_lo = static_cast<uint32_t>(x);
    const uint32_t x_hi = static_cast<uint32_t>(x >> 32);
    const uint32_t y_lo = static_cast<uint32_t>(y);
    const uint32_t y_hi = static_cast<uint32_t>(y >> 32);

    const uint64_t b00 = uint64_t{x_lo} * y_lo;
    const uint64_t b01 = uint64_t{x_lo} * y_hi;
    const uint64_t b10 = uint64_t{x_hi} * y_lo;
    const uint64_t b11 = uint64_t{x_hi} * y_hi;

    const uint32_t b00_hi = static_cast<uint32_t>(b00 >> 32);

    const uint64_t mid1 = b10 + b00_hi;
    const uint32_t mid1_lo = static_cast<uint32_t>(mid1);
    const uint32_t mid1_hi = static_cast<uint32_t>(mid1 >> 32);

    const uint64_t mid2 = b01 + mid1_lo;
    const uint32_t mid2_lo = static_cast<uint32_t>(mid2);
    const uint32_t mid2_hi = static_cast<uint32_t>(mid2 >> 32);

    // NB: mid2_lo has the upper 32 bits of the low part of the product.
    const uint32_t r = mid2_lo >> 31;
    const uint64_t h = b11 + mid1_hi + mid2_hi + r;

    return h;
#endif
}

// Returns the number of leading 0-bits in x, starting at the most significant
// bit position.
// If x is 0, the result is undefined.
static inline int32_t CountLeadingZeros64(uint64_t x)
{
    GRISU_ASSERT(x != 0);

#if defined(__GNUC__) || defined(__clang__)
    return __builtin_clzll(x);
#elif defined(_MSC_VER) && (defined(_M_ARM) || defined(_M_ARM64))
    return static_cast<int32_t>(_CountLeadingZeros64(x));
#elif defined(_MSC_VER) && defined(_M_X64)
    return static_cast<int32_t>(__lzcnt64(x));
#elif defined(_MSC_VER) && defined(_M_IX86)
    int32_t lz = static_cast<int32_t>( __lzcnt(static_cast<uint32_t>(x >> 32)) );
    if (lz == 32) {
        lz += static_cast<int32_t>( __lzcnt(static_cast<uint32_t>(x)) );
    }
    return lz;
#else
    int32_t lz = 0;
    while ((x >> 63) == 0) {
        x <<= 1;
        ++lz;
    }
    return lz;
#endif
}

// Given normalized DiyFp w, Grisu needs to find a (normalized) cached
// power-of-ten c, such that the exponent of the product c * w = f * 2^e lies
// within a certain range [alpha, gamma] (Definition 3.2 from [1])
//
//      alpha <= e = e_c + e_w + q <= gamma
//
// or
//
//      f_c * f_w * 2^alpha <= f_c 2^(e_c) * f_w 2^(e_w) * 2^q
//                          <= f_c * f_w * 2^gamma
//
// Since c and w are normalized, i.e. 2^(q-1) <= f < 2^q, this implies
//
//      2^(q-1) * 2^(q-1) * 2^alpha <= c * w * 2^q < 2^q * 2^q * 2^gamma
//
// or
//
//      2^(q - 2 + alpha) <= c * w < 2^(q + gamma)
//
// The choice of (alpha,gamma) determines the size of the table and the form of
// the digit generation procedure.
//
// Now
//
//      alpha <= e_c + e + q <= gamma                                        (1)
//      ==> f_c * 2^alpha <= c * 2^e * 2^q
//
// and since the c's are normalized, 2^(q-1) <= f_c,
//
//      ==> 2^(q - 1 + alpha) <= c * 2^(e + q)
//      ==> 2^(alpha - e - 1) <= c
//
// If c were an exakt power of ten, i.e. c = 10^k, one may determine k as
//
//      k = ceil( log_10( 2^(alpha - e - 1) ) )
//        = ceil( (alpha - e - 1) * log_10(2) )
//
// From the paper:
// "In theory the result of the procedure could be wrong since c is rounded, and
//  the computation itself is approximated [...]. In practice, however, this
//  simple function is sufficient."
//
// For IEEE double precision floating-point numbers converted into normalized
// DiyFp's w = f * 2^e, with q = 64,
//
//      e >= -1022      (min IEEE exponent)
//           -52        (p - 1)
//           -52        (p - 1, possibly normalize denormal IEEE numbers)
//           -11        (normalize the DiyFp)
//         = -1137
//
// and
//
//      e <= +1023      (max IEEE exponent)
//           -52        (p - 1)
//           -11        (normalize the DiyFp)
//         = 960
//
// For IEEE single-precision the range is [-180, 96].
//
// One does not need to store a cached power for each k in this range. For each
// such k it suffices to find a cached power such that the exponent of the
// product lies in [alpha,gamma].
// This implies that the difference of the decimal exponents of adjacent table
// entries must be less than or equal to
//
//      floor( (gamma - alpha) * log_10(2) )
//
// (A smaller distance gamma-alpha would require a larger table.)

namespace {
struct CachedPower { // c = f * 2^e ~= 10^k
    uint64_t f;
    int32_t e; // binary exponent
    int32_t k; // decimal exponent
};
}

static constexpr int32_t kAlpha = -50;
static constexpr int32_t kGamma = -36;
// k_min = -304
// k_max =  327

static constexpr int32_t kCachedPowersSize       =  159;
static constexpr int32_t kCachedPowersMinDecExp  = -304;
static constexpr int32_t kCachedPowersMaxDecExp  =  328;
static constexpr int32_t kCachedPowersDecExpStep =    4;

// For a normalized DiyFp w = f * 2^e, this function returns a (normalized)
// cached power-of-ten c = f_c * 2^e_c, such that the exponent of the product
// w * c satisfies
//
//      kAlpha <= e_c + e + q <= kGamma.
//
static inline CachedPower GetCachedPowerForBinaryExponent(int32_t e)
{
    static constexpr uint64_t kSignificands[] = {
        0x8C71DCD9BA0B4926, // e = -1073, k = -304
        0xAB70FE17C79AC6CA, // e = -1060, k = -300
        0xD1476E2C07286FAA, // e = -1047, k = -296
        0xFF77B1FCBEBCDC4F, // e = -1034, k = -292
        0x9BECCE62836AC577, // e = -1020, k = -288
        0xBE5691EF416BD60C, // e = -1007, k = -284
        0xE858AD248F5C22CA, // e =  -994, k = -280
        0x8DD01FAD907FFC3C, // e =  -980, k = -276
        0xAD1C8EAB5EE43B67, // e =  -967, k = -272
        0xD3515C2831559A83, // e =  -954, k = -268
        0x80FA687F881C7F8E, // e =  -940, k = -264
        0x9D71AC8FADA6C9B5, // e =  -927, k = -260
        0xC0314325637A193A, // e =  -914, k = -256
        0xEA9C227723EE8BCB, // e =  -901, k = -252
        0x8F31CC0937AE58D3, // e =  -887, k = -248
        0xAECC49914078536D, // e =  -874, k = -244
        0xD5605FCDCF32E1D7, // e =  -861, k = -240
        0x823C12795DB6CE57, // e =  -847, k = -236
        0x9EFA548D26E5A6E2, // e =  -834, k = -232
        0xC21094364DFB5637, // e =  -821, k = -228
        0xECE53CEC4A314EBE, // e =  -808, k = -224
        0x9096EA6F3848984F, // e =  -794, k = -220
        0xB080392CC4349DED, // e =  -781, k = -216
        0xD77485CB25823AC7, // e =  -768, k = -212
        0x8380DEA93DA4BC60, // e =  -754, k = -208
        0xA086CFCD97BF97F4, // e =  -741, k = -204
        0xC3F490AA77BD60FD, // e =  -728, k = -200
        0xEF340A98172AACE5, // e =  -715, k = -196
        0x91FF83775423CC06, // e =  -701, k = -192
        0xB23867FB2A35B28E, // e =  -688, k = -188
        0xD98DDAEE19068C76, // e =  -675, k = -184
        0x84C8D4DFD2C63F3B, // e =  -661, k = -180
        0xA21727DB38CB0030, // e =  -648, k = -176
        0xC5DD44271AD3CDBA, // e =  -635, k = -172
        0xF18899B1BC3F8CA2, // e =  -622, k = -168
        0x936B9FCEBB25C996, // e =  -608, k = -164
        0xB3F4E093DB73A093, // e =  -595, k = -160
        0xDBAC6C247D62A584, // e =  -582, k = -156
        0x8613FD0145877586, // e =  -568, k = -152
        0xA3AB66580D5FDAF6, // e =  -555, k = -148
        0xC7CABA6E7C5382C9, // e =  -542, k = -144
        0xF3E2F893DEC3F126, // e =  -529, k = -140
        0x94DB483840B717F0, // e =  -515, k = -136
        0xB5B5ADA8AAFF80B8, // e =  -502, k = -132
        0xDDD0467C64BCE4A1, // e =  -489, k = -128
        0x87625F056C7C4A8B, // e =  -475, k = -124
        0xA54394FE1EEDB8FF, // e =  -462, k = -120
        0xC9BCFF6034C13053, // e =  -449, k = -116
        0xF64335BCF065D37D, // e =  -436, k = -112
        0x964E858C91BA2655, // e =  -422, k = -108
        0xB77ADA0617E3BBCB, // e =  -409, k = -104
        0xDFF9772470297EBD, // e =  -396, k = -100
        0x88B402F7FD75539B, // e =  -382, k =  -96
        0xA6DFBD9FB8E5B88F, // e =  -369, k =  -92
        0xCBB41EF979346BCA, // e =  -356, k =  -88
        0xF8A95FCF88747D94, // e =  -343, k =  -84
        0x97C560BA6B0919A6, // e =  -329, k =  -80
        0xB94470938FA89BCF, // e =  -316, k =  -76
        0xE2280B6C20DD5232, // e =  -303, k =  -72
        0x8A08F0F8BF0F156B, // e =  -289, k =  -68
        0xA87FEA27A539E9A5, // e =  -276, k =  -64
        0xCDB02555653131B6, // e =  -263, k =  -60
        0xFB158592BE068D2F, // e =  -250, k =  -56
        0x993FE2C6D07B7FAC, // e =  -236, k =  -52
        0xBB127C53B17EC159, // e =  -223, k =  -48
        0xE45C10C42A2B3B06, // e =  -210, k =  -44
        0x8B61313BBABCE2C6, // e =  -196, k =  -40
        0xAA242499697392D3, // e =  -183, k =  -36
        0xCFB11EAD453994BA, // e =  -170, k =  -32
        0xFD87B5F28300CA0E, // e =  -157, k =  -28
        0x9ABE14CD44753B53, // e =  -143, k =  -24
        0xBCE5086492111AEB, // e =  -130, k =  -20
        0xE69594BEC44DE15B, // e =  -117, k =  -16
        0x8CBCCC096F5088CC, // e =  -103, k =  -12
        0xABCC77118461CEFD, // e =   -90, k =   -8
        0xD1B71758E219652C, // e =   -77, k =   -4
        0x8000000000000000, // e =   -63, k =    0
        0x9C40000000000000, // e =   -50, k =    4
        0xBEBC200000000000, // e =   -37, k =    8
        0xE8D4A51000000000, // e =   -24, k =   12
        0x8E1BC9BF04000000, // e =   -10, k =   16
        0xAD78EBC5AC620000, // e =     3, k =   20
        0xD3C21BCECCEDA100, // e =    16, k =   24
        0x813F3978F8940984, // e =    30, k =   28
        0x9DC5ADA82B70B59E, // e =    43, k =   32
        0xC097CE7BC90715B3, // e =    56, k =   36
        0xEB194F8E1AE525FD, // e =    69, k =   40
        0x8F7E32CE7BEA5C70, // e =    83, k =   44
        0xAF298D050E4395D7, // e =    96, k =   48
        0xD5D238A4ABE98068, // e =   109, k =   52
        0x82818F1281ED44A0, // e =   123, k =   56
        0x9F4F2726179A2245, // e =   136, k =   60
        0xC2781F49FFCFA6D5, // e =   149, k =   64
        0xED63A231D4C4FB27, // e =   162, k =   68
        0x90E40FBEEA1D3A4B, // e =   176, k =   72
        0xB0DE65388CC8ADA8, // e =   189, k =   76
        0xD7E77A8F87DAF7FC, // e =   202, k =   80
        0x83C7088E1AAB65DB, // e =   216, k =   84
        0xA0DC75F1778E39D6, // e =   229, k =   88
        0xC45D1DF942711D9A, // e =   242, k =   92
        0xEFB3AB16C59B14A3, // e =   255, k =   96
        0x924D692CA61BE758, // e =   269, k =  100
        0xB2977EE300C50FE7, // e =   282, k =  104
        0xDA01EE641A708DEA, // e =   295, k =  108
        0x850FADC09923329E, // e =   309, k =  112
        0xA26DA3999AEF774A, // e =   322, k =  116
        0xC646D63501A1511E, // e =   335, k =  120
        0xF209787BB47D6B85, // e =   348, k =  124
        0x93BA47C980E98CE0, // e =   362, k =  128
        0xB454E4A179DD1877, // e =   375, k =  132
        0xDC21A1171D42645D, // e =   388, k =  136
        0x865B86925B9BC5C2, // e =   402, k =  140
        0xA402B9C5A8D3A6E7, // e =   415, k =  144
        0xC83553C5C8965D3D, // e =   428, k =  148
        0xF46518C2EF5B8CD1, // e =   441, k =  152
        0x952AB45CFA97A0B3, // e =   455, k =  156
        0xB616A12B7FE617AA, // e =   468, k =  160
        0xDE469FBD99A05FE3, // e =   481, k =  164
        0x87AA9AFF79042287, // e =   495, k =  168
        0xA59BC234DB398C25, // e =   508, k =  172
        0xCA28A291859BBF93, // e =   521, k =  176
        0xF6C69A72A3989F5C, // e =   534, k =  180
        0x969EB7C47859E744, // e =   548, k =  184
        0xB7DCBF5354E9BECE, // e =   561, k =  188
        0xE070F78D3927556B, // e =   574, k =  192
        0x88FCF317F22241E2, // e =   588, k =  196
        0xA738C6BEBB12D16D, // e =   601, k =  200
        0xCC20CE9BD35C78A5, // e =   614, k =  204
        0xF92E0C3537826146, // e =   627, k =  208
        0x98165AF37B2153DF, // e =   641, k =  212
        0xB9A74A0637CE2EE1, // e =   654, k =  216
        0xE2A0B5DC971F303A, // e =   667, k =  220
        0x8A5296FFE33CC930, // e =   681, k =  224
        0xA8D9D1535CE3B396, // e =   694, k =  228
        0xCE1DE40642E3F4B9, // e =   707, k =  232
        0xFB9B7CD9A4A7443C, // e =   720, k =  236
        0x9991A6F3D6BF1766, // e =   734, k =  240
        0xBB764C4CA7A44410, // e =   747, k =  244
        0xE4D5E82392A40515, // e =   760, k =  248
        0x8BAB8EEFB6409C1A, // e =   774, k =  252
        0xAA7EEBFB9DF9DE8E, // e =   787, k =  256
        0xD01FEF10A657842C, // e =   800, k =  260
        0xFE0EFB53D30DD4D8, // e =   813, k =  264
        0x9B10A4E5E9913129, // e =   827, k =  268
        0xBD49D14AA79DBC82, // e =   840, k =  272
        0xE7109BFBA19C0C9D, // e =   853, k =  276
        0x8D07E33455637EB3, // e =   867, k =  280
        0xAC2820D9623BF429, // e =   880, k =  284
        0xD226FC195C6A2F8C, // e =   893, k =  288
        0x80444B5E7AA7CF85, // e =   907, k =  292
        0x9C935E00D4B9D8D2, // e =   920, k =  296
        0xBF21E44003ACDD2D, // e =   933, k =  300
        0xE950DF20247C83FD, // e =   946, k =  304
        0x8E679C2F5E44FF8F, // e =   960, k =  308
        0xADD57A27D29339F6, // e =   973, k =  312
        0xD433179D9C8CB841, // e =   986, k =  316
        0x81842F29F2CCE376, // e =  1000, k =  320
        0x9E19DB92B4E31BA9, // e =  1013, k =  324
        0xC0FE908895CF3B44, // e =  1026, k =  328
    };

    GRISU_ASSERT(e >= -1137);
    GRISU_ASSERT(e <=   960);

    const int32_t k = CeilLog10Pow2(kAlpha - e - 1);
    GRISU_ASSERT(k >= kCachedPowersMinDecExp - (kCachedPowersDecExpStep - 1));
    GRISU_ASSERT(k <= kCachedPowersMaxDecExp);

    const unsigned index = static_cast<unsigned>(k - (kCachedPowersMinDecExp - (kCachedPowersDecExpStep - 1))) / kCachedPowersDecExpStep;
    GRISU_ASSERT(index < kCachedPowersSize);

    const int32_t k_cached = kCachedPowersMinDecExp + static_cast<int32_t>(index) * kCachedPowersDecExpStep;
    const int32_t e_cached = FloorLog2Pow10(k_cached) + 1 - 64;

    const CachedPower cached = {kSignificands[index], e_cached, k_cached};
    GRISU_ASSERT(kAlpha <= cached.e + e + 64);
    GRISU_ASSERT(kGamma >= cached.e + e + 64);

    return cached;
}

namespace {
struct FloatingDecimal64 {
    uint64_t digits;
    int32_t exponent;
};
}

static inline FloatingDecimal64 ToDecimal64(double value)
{
    static_assert(DiyFp::SignificandSize >= std::numeric_limits<double>::digits + 3,
        "Grisu2 requires q >= p + 3");
    static_assert(DiyFp::SignificandSize == 64,
        "This implementation requires q = 64");

    // Compute the boundaries m- and m+ of the floating-point value
    // v = f * 2^e.
    //
    // Determine v- and v+, the floating-point predecessor and successor if v,
    // respectively.
    //
    //      v- = v - 2^e        if f != 2^(p-1) or e == e_min                (A)
    //         = v - 2^(e-1)    if f == 2^(p-1) and e > e_min                (B)
    //
    //      v+ = v + 2^e
    //
    // Let m- = (v- + v) / 2 and m+ = (v + v+) / 2. All real numbers _strictly_
    // between m- and m+ round to v, regardless of how the input rounding
    // algorithm breaks ties.
    //
    //      ---+-------------+-------------+-------------+-------------+---  (A)
    //         v-            m-            v             m+            v+
    //
    //      -----------------+------+------+-------------+-------------+---  (B)
    //                       v-     m-     v             m+            v+

    GRISU_ASSERT(Double(value).IsFinite());
    GRISU_ASSERT(value > 0);

    const auto ieee_value = Double(value);
    const auto ieee_significand = ieee_value.PhysicalSignificand();
    const auto ieee_exponent    = ieee_value.PhysicalExponent();

    int32_t shared_exponent;

#if GRISU_ROUND()
    uint64_t m_minus;
    uint64_t v;
    uint64_t m_plus;
    if (ieee_exponent != 0) // normalized floating-point number
    {
        const bool lower_boundary_is_closer = (ieee_significand == 0 && ieee_exponent > 1);

        const auto f2 = ieee_significand | Double::HiddenBit;
        const auto e2 = static_cast<int32_t>(ieee_exponent) - Double::ExponentBias;

#if GRISU_SMALL_INT_OPTIMIZATION()
        if (0 <= -e2 && -e2 < 53)
        {
            const uint64_t d2 = f2 >> -e2;
            if (d2 << -e2 == f2)
                return {d2, 0};
        }
#endif

        const auto fm = 4 * f2 - 2 + (lower_boundary_is_closer ? 1 : 0);
        const auto fv = 4 * f2;
        const auto fp = 4 * f2 + 2;

        const auto shift = DiyFp::SignificandSize - Double::SignificandSize - 2;

        shared_exponent = e2 - 2 - shift;
        m_minus = uint64_t{fm} << shift;
        v       = uint64_t{fv} << shift;
        m_plus  = uint64_t{fp} << shift;
    }
    else
    {
        const auto f2 = ieee_significand;
        const auto e2 = 1 - Double::ExponentBias;

        const auto fm = 4 * f2 - 2;
        const auto fv = 4 * f2;
        const auto fp = 4 * f2 + 2;

        const int32_t shift = CountLeadingZeros64(fv);

        shared_exponent = e2 - 2 - shift;
        m_minus = uint64_t{fm} << shift;
        v       = uint64_t{fv} << shift;
        m_plus  = uint64_t{fp} << shift;
    }
#else // ^^^ GRISU_ROUND() ^^^
    uint64_t m_minus;
    uint64_t m_plus;
    if (ieee_exponent != 0) // normalized floating-point number
    {
        const bool lower_boundary_is_closer = (ieee_significand == 0 && ieee_exponent > 1);

        const auto f2 = ieee_significand | Double::HiddenBit;
        const auto e2 = static_cast<int32_t>(ieee_exponent) - Double::ExponentBias;

#if GRISU_SMALL_INT_OPTIMIZATION()
        if (0 <= -e2 && -e2 < 53)
        {
            const uint64_t d2 = f2 >> -e2;
            if (d2 << -e2 == f2)
                return {d2, 0};
        }
#endif

        const auto fm = 4 * f2 - 2 + (lower_boundary_is_closer ? 1 : 0);
        const auto fp = 4 * f2 + 2;

        const auto shift = DiyFp::SignificandSize - Double::SignificandSize - 2;

        shared_exponent = e2 - 2 - shift;
        m_minus = uint64_t{fm} << shift;
        m_plus  = uint64_t{fp} << shift;
    }
    else
    {
        const auto f2 = ieee_significand;
        const auto e2 = 1 - Double::ExponentBias;

        const auto fm = 4 * f2 - 2;
        const auto fp = 4 * f2 + 2;

        const auto shift = CountLeadingZeros64(fp);

        shared_exponent = e2 - 2 - shift;
        m_minus = uint64_t{fm} << shift;
        m_plus  = uint64_t{fp} << shift;
    }
#endif // ^^^ not GRISU_ROUND() ^^^

    //
    // Step 1:
    // Compute rounding interval
    //

    //  --------+-----------------------+-----------------------+--------    (A)
    //          m-                      v                       m+
    //
    //  --------------------+-----------+-----------------------+--------    (B)
    //                      m-          v                       m+
    //
    // First scale v (and m- and m+) such that the exponent is in the range
    // [alpha, gamma].

    const auto cached = GetCachedPowerForBinaryExponent(shared_exponent);

    const uint64_t w_minus = MultiplyHighRoundUp(m_minus, cached.f); // XXX: round down?
#if GRISU_ROUND()
    const uint64_t w       = MultiplyHighRoundUp(v,       cached.f); // XXX: compute from w_minus/w_plus?
#endif
    const uint64_t w_plus  = MultiplyHighRoundUp(m_plus,  cached.f);

    // The exponent of the products is = v.e + cached.e + q and is in the
    // range [alpha, gamma].
    const int32_t e = shared_exponent + cached.e + 64;
    GRISU_ASSERT(e >= kAlpha);
    GRISU_ASSERT(e <= kGamma);

    // Note:
    // The result of Multiply() is **NOT** neccessarily normalized.
    // But since m+ and c are normalized, w+ >= 2^(q - 2).
    GRISU_ASSERT(w_plus >= (uint64_t{1} << (64 - 2)));

    //  ----(---+---)---------------(---+---)---------------(---+---)----
    //          w-                      w                       w+
    //          = c*m-                  = c*v                   = c*m+
    //
    // Multiply rounds its result and c_minus_k is approximated too. w, w- and
    // w+ are now off by a small amount.
    // In fact:
    //
    //      w - v * 10^-k < 1 ulp
    //
    // To account for this inaccuracy, add resp. subtract 1 ulp.
    // Note: ulp(w-) = ulp(w) = ulp(w+).
    //
    //  ----(---+---[---------------(---+---)---------------]---+---)----
    //          w-  L                   w                   H   w+
    //
    // Now any number in [L, H] (bounds included) will round to w when input,
    // regardless of how the input rounding algorithm breaks ties.
    //
    // And DigitGen generates the shortest possible such number in [L, H].
    // Note that this does not mean that Grisu2 always generates the shortest
    // possible number in the interval (m-, m+).

    //GRISU_ASSERT(w_plus - w_minus >= (3 * Pow2(64 - Double::SignificandSize - 2)) / 2); // >= delta_m / 2
    // w_plus - w_minus >= 768          (double precision)
    // w_plus - w_minus >= 412316860416 (single precision)
    const uint64_t L = w_minus + 1;
    const uint64_t H = w_plus  - 1;
    //GRISU_ASSERT(H - L >= (3 * Pow2(64 - Double::SignificandSize - 2)) / 2 - 2);

    //
    // Step 2:
    // Generate digits
    //

    static_assert(kAlpha >= -50, "internal error");
    static_assert(kGamma <= -32, "internal error");

    // Generates the digits (and the exponent) of a decimal floating-point
    // number V = digits * 10^exponent in the range [L, H].
    // The DiyFp's w, L and H share the same exponent e, which satisfies
    // alpha <= e <= gamma.
    //
    //                                  <---- distance ----->
    //              <---------------------------- delta ---->
    //  ----(---+---[---------------(---+---)---------------]---+---)----
    //              L                   w                   H
    //
    // This routine generates the digits of H from left to right and stops as
    // soon as V is in [L, H].

#if GRISU_ROUND()
    GRISU_ASSERT(H >= w);
    uint64_t distance = H - w; // (significand of (H - w), implicit exponent is H.e)
#endif
    GRISU_ASSERT(H >= L);
    uint64_t delta    = H - L; // (significand of (H - L), implicit exponent is H.e)
    uint64_t rest;
    uint64_t ten_kappa;

    // Split H = f * 2^e into two parts p1 and p2 (note: e < 0):
    //
    //      H = f * 2^e
    //           = ((f div 2^-e) * 2^-e + (f mod 2^-e)) * 2^e
    //           = ((p1        ) * 2^-e + (p2        )) * 2^e
    //           = p1 + p2 * 2^e

    const DiyFp one(uint64_t{1} << -e, e); // one = 2^-e * 2^e

    uint32_t p1 = static_cast<uint32_t>(H >> -one.e); // p1 = f div 2^-e (Since -e >= 32, p1 fits into a 32-bit int32_t.)
    uint64_t p2 = H & (one.f - 1);                    // p2 = f mod 2^-e

    GRISU_ASSERT(p1 >= 4); // (2^(64-2) - 1) >> 60

    uint64_t digits = p1;
    int32_t exponent = 0;

    if (p2 > delta)
    {
        // We have
        //
        //      H = d[k-1]...d[1]d[0] + p2 * 2^e
        //        = digits            + p2 * 2^e
        //
        // Now generate the digits of the fractional part p2 * 2^e.
        // p2 actually represents the fraction
        //
        //      p2 * 2^e
        //          = p2 / 2^-e
        //          = d[-1] / 10^1 + d[-2] / 10^2 + ...
        //
        // Now generate the digits d[-m] of p1 from left to right (m = 1,2,...)
        //
        //      p2 * 2^e = d[-1]d[-2]...d[-m] * 10^-m
        //                      + 10^-m * (d[-m-1] / 10^1 + d[-m-2] / 10^2 + ...)
        //
        // using
        //
        //      10^m * p2 = ((10^m * p2) div 2^-e) * 2^-e + ((10^m * p2) mod 2^-e)
        //                = (                   d) * 2^-e + (                   r)
        //
        // or
        //      10^m * p2 * 2^e = d + r * 2^e
        //
        // i.e.
        //
        //      H = digits + p2 * 2^e
        //        = digits + 10^-m * (d + r * 2^e)
        //        = (digits * 10^m + d) * 10^-m + 10^-m * r * 2^e
        //
        // and stop as soon as 10^-m * r * 2^e <= delta * 2^e

        // unit = 1
        // m = 0
        for (;;)
        {
            GRISU_ASSERT(digits <= 99999999999999999ull);

            GRISU_ASSERT(p2 <= UINT64_MAX / 10000);
            const uint64_t s4 = 10000 * p2;
            const uint64_t r4 = s4 & (one.f - 1);
            const uint32_t q4 = static_cast<uint32_t>(s4 >> -one.e);
            GRISU_ASSERT(q4 <= 9999);

            if (r4 <= 10000 * delta)
            {
                const uint64_t s2 = 100 * p2;
                const uint64_t r2 = s2 & (one.f - 1);
                const uint32_t q2 = static_cast<uint32_t>(s2 >> -one.e);
                GRISU_ASSERT(q2 <= 99);

                if (r2 <= 100 * delta)
                {
                    const uint64_t s1 = 10 * p2;
                    const uint64_t r1 = s1 & (one.f - 1);
                    const uint32_t q1 = static_cast<uint32_t>(s1 >> -one.e);
                    GRISU_ASSERT(q1 <= 9);

                    if (r1 <= 10 * delta)
                    {
                        digits = 10 * digits + q1;
                        exponent -= 1; // m += 1
#if GRISU_ROUND()
                        p2 = r1;
                        delta *= 10;
                        distance *= 10;
#endif
                    }
                    else
                    {
                        digits = 100 * digits + q2;
                        exponent -= 2; // m += 2
#if GRISU_ROUND()
                        p2 = r2;
                        delta *= 100;
                        distance *= 100;
#endif
                    }
                }
                else
                {
                    const uint64_t s3 = 1000 * p2;
                    const uint64_t r3 = s3 & (one.f - 1);
                    const uint64_t q3 = static_cast<uint32_t>(s3 >> -one.e);
                    GRISU_ASSERT(q3 <= 999);

                    if (r3 <= 1000 * delta)
                    {
                        digits = 1000 * digits + q3;
                        exponent -= 3; // m += 3
#if GRISU_ROUND()
                        p2 = r3;
                        delta *= 1000;
                        distance *= 1000;
#endif
                    }
                    else
                    {
                        digits = 10000 * digits + q4;
                        exponent -= 4; // m += 4
#if GRISU_ROUND()
                        p2 = r4;
                        delta *= 10000;
                        distance *= 10000;
#endif
                    }
                }

                // V = digits * 10^-m, with L <= V <= H.
                // exponent = -m

#if GRISU_ROUND()
                rest = p2;
                ten_kappa = one.f; // one.f == 2^-e
#endif
                break;
            }

            digits = 10000 * digits + q4;
            exponent -= 4; // m += 4
            p2 = r4;
            delta *= 10000;
#if GRISU_ROUND()
            distance *= 10000;
#endif
        }
    }
    else // p2 <= delta
    {
        GRISU_ASSERT((uint64_t{p1} << -one.e) + p2 > delta); // Loop terminates.

        // In this case: p1 contains too many digits.
        //
        // Find the largest 0 <= n < k = length, such that
        //
        //      H = (p1 div 10^n) * 10^n + ((p1 mod 10^n) * 2^-e + p2) * 2^e
        //        = (p1 div 10^n) * 10^n + (                     rest) * 2^e
        //
        // and rest <= delta.
        //
        // Compute rest * 2^e = H mod 10^n = p1 + p2 * 2^e = (p1 * 2^-e + p2) * 2^e
        // and check if enough digits have been generated:
        //
        //      rest * 2^e <= delta * 2^e
        //

        rest = p2;

        // 10^n is now 1 ulp in the decimal representation V. The rounding
        // procedure works with DiyFp's with an implicit exponent of e.
        //
        //      10^n = (10^n * 2^-e) * 2^e = ten_kappa * 2^e
        //
        ten_kappa = one.f; // Start with 2^-e

        // n = 0
        for (;;)
        {
            GRISU_ASSERT(rest <= delta);

            // rn = d[n]...d[0] * 2^-e + p2
            const uint32_t q = p1 / 10;
            const uint32_t r = p1 % 10;
            const uint64_t r_next = ten_kappa * r + rest;

            if (r_next > delta)
            {
                digits = p1;
                break;
            }

            p1 = q;
            exponent += 1; // n += 1
            rest = r_next;
            ten_kappa *= 10;
        }
    }

    //
    // Step 3 (optional):
    // Round towards w.
    //

#if GRISU_ROUND()
    GRISU_ASSERT(digits >= 1);
    GRISU_ASSERT(distance <= delta);
    GRISU_ASSERT(rest <= delta);
    GRISU_ASSERT(ten_kappa > 0);

    // By generating the digits of H we got the largest (closest to H) value
    // that is still in the interval [L, H]. In the case where w < B <= H we
    // try to decrement this value.
    //
    //                                  <---- distance ----->
    //              <---------------------------- delta ---->
    //                                         <--- rest --->
    //                       <--- ten_kappa --->
    //  ----(---+---[---------------(---+---)--+------------]---+---)----
    //              L                   w      B            H
    //                                         = digits * 10^kappa
    //
    // ten_kappa represents a unit-in-the-last-place in the decimal
    // representation stored in 'digits'.
    //
    // There are three stopping conditions:
    // (The position of the numbers is measured relative to H.)
    //
    //  1)  B is already <= w
    //          rest >= distance
    //
    //  2)  Decrementing B would yield a number B' < L
    //          rest + ten_kappa > delta
    //
    //  3)  Decrementing B would yield a number B' < w and farther away from
    //      w than the current number B: w - B' > B - w
    //          rest + ten_kappa > distance &&
    //          rest + ten_kappa - distance >= distance - rest

    // The tests are written in this order to avoid overflow in unsigned
    // integer arithmetic.

    while (rest < distance
        && delta - rest >= ten_kappa
        && (rest + ten_kappa <= distance || rest + ten_kappa - distance < distance - rest))
    {
        GRISU_ASSERT(digits % 10 != 0);
        digits--;
        rest += ten_kappa;
    }
#endif

    //
    // Done.
    //

    // v = decimal_digits * 10^decimal_exponent

    return {digits, exponent - cached.k};
}

//==================================================================================================
// ToChars
//==================================================================================================

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

    GRISU_ASSERT(digits <= 99);
    std::memcpy(buf, &Digits100[2 * digits], 2);
    return buf + 2;
}

static inline char* Utoa_4Digits(char* buf, uint32_t digits)
{
    GRISU_ASSERT(digits <= 9999);
    const uint32_t q = digits / 100;
    const uint32_t r = digits % 100;
    Utoa_2Digits(buf + 0, q);
    Utoa_2Digits(buf + 2, r);
    return buf + 4;
}

static inline char* Utoa_8Digits(char* buf, uint32_t digits)
{
    GRISU_ASSERT(digits <= 99999999);
    const uint32_t q = digits / 10000;
    const uint32_t r = digits % 10000;
    Utoa_4Digits(buf + 0, q);
    Utoa_4Digits(buf + 4, r);
    return buf + 8;
}

static inline int32_t DecimalLength(uint64_t v)
{
    GRISU_ASSERT(v >= 1);
    GRISU_ASSERT(v <= 99999999999999999ull);

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

static inline void PrintDecimalDigits(char* buf, uint64_t output, int32_t output_length)
{
    // We prefer 32-bit operations, even on 64-bit platforms.
    // We have at most 17 digits, and uint32_t can store 9 digits.
    // If output doesn't fit into uint32_t, we cut off 8 digits,
    // so the rest will fit into uint32_t.
    if (static_cast<uint32_t>(output >> 32) != 0)
    {
        GRISU_ASSERT(output_length > 8);
        const uint64_t q = output / 100000000;
        const uint32_t r = static_cast<uint32_t>(output % 100000000);
        output = q;
        output_length -= 8;
        Utoa_8Digits(buf + output_length, r);
    }

    GRISU_ASSERT(output <= UINT32_MAX);
    uint32_t output2 = static_cast<uint32_t>(output);

    while (output2 >= 10000)
    {
        GRISU_ASSERT(output_length > 4);
        const uint32_t q = output2 / 10000;
        const uint32_t r = output2 % 10000;
        output2 = q;
        output_length -= 4;
        Utoa_4Digits(buf + output_length, r);
    }

    if (output2 >= 100)
    {
        GRISU_ASSERT(output_length > 2);
        const uint32_t q = output2 / 100;
        const uint32_t r = output2 % 100;
        output2 = q;
        output_length -= 2;
        Utoa_2Digits(buf + output_length, r);
    }

    if (output2 >= 10)
    {
        GRISU_ASSERT(output_length == 2);
        Utoa_2Digits(buf, output2);
    }
    else
    {
        GRISU_ASSERT(output_length == 1);
        buf[0] = static_cast<char>('0' + output2);
    }
}

static inline char* FormatDigits(char* buffer, uint64_t digits, int32_t decimal_exponent, bool force_trailing_dot_zero = false)
{
    GRISU_ASSERT(digits >= 1);
    GRISU_ASSERT(digits <= 99999999999999999ull);
    GRISU_ASSERT(decimal_exponent >= -999);
    GRISU_ASSERT(decimal_exponent <=  999);

    const int32_t num_digits = DecimalLength(digits);
    const int32_t decimal_point = num_digits + decimal_exponent;

    // In order to successfully parse all numbers output by Dtoa using the Strtod implementation
    // below, we have to make sure to never emit more than 17 (significant) digits.
    static constexpr int32_t MaxFixedDecimalPoint =  17;
    static constexpr int32_t MinFixedDecimalPoint = -6;
#if !GRISU_SMALL_INT_OPTIMIZATION()
    static_assert(MaxFixedDecimalPoint >= 17, "internal error");
#endif

    const bool use_fixed = MinFixedDecimalPoint <= decimal_point && decimal_point <= MaxFixedDecimalPoint;

    // Prepare the buffer.
    // Avoid calling memset/memcpy with variable arguments below...

    int32_t decimal_digits_position;
    if (use_fixed)
    {
        if (decimal_point <= 0)
        {
            // 0.[000]digits
            // -6 <= decimal_point <= 0
            //  ==> 2 <= 2 + -decimal_point <= 8
            // Pre-filling the buffer with 8 '0's is therefore sufficient.
            std::memset(buffer, '0', 8);
            decimal_digits_position = 2 + (-decimal_point);
        }
        else if (decimal_point < num_digits)
        {
            // dig.its
            // 0 < decimal_point <= Min(17 - 1, MaxExp)
            // We need to move at most 16 bytes to the right.
            decimal_digits_position = 0;
        }
        else
        {
            // digits[000]
            // 1 <= num_digits <= 17 decimal_point <= 21.
            // Pre-filling buffer with 21 '0's is therefore sufficient.
            static_assert(MaxFixedDecimalPoint <= 24, "invalid parameter");
            std::memset(buffer, '0', 24);
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
            // We need to move at most 16 bytes one place to the right.
            std::memmove(buffer + (decimal_point + 1), buffer + decimal_point, 16);
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
//      GRISU_ASSERT(scientific_exponent != 0);

        std::memcpy(buffer, scientific_exponent < 0 ? "e-" : "e+", 2);
        buffer += 2;
        const uint32_t k = static_cast<uint32_t>(scientific_exponent < 0 ? -scientific_exponent : scientific_exponent);
        if (k < 10)
        {
            *buffer++ = static_cast<char>('0' + k);
        }
        else if (k < 100)
        {
            buffer = Utoa_2Digits(buffer, k);
        }
        else
        {
            const uint32_t r = k % 10;
            const uint32_t q = k / 10;
            buffer = Utoa_2Digits(buffer, q);
            *buffer++ = static_cast<char>('0' + r);
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
        value = v.AbsValue();
        *buffer++ = '-';
    }

    if (v.IsZero())
    {
        std::memcpy(buffer, "0.0 ", 4);
        buffer += 1 + (force_trailing_dot_zero ? 2 : 0);
        return buffer;
    }

    const auto dec = ToDecimal64(value);
    return FormatDigits(buffer, dec.digits, dec.exponent, force_trailing_dot_zero);
}

//==================================================================================================
//
//==================================================================================================

char* grisu2::Dtoa(char* buffer, double value)
{
    return ToChars(buffer, value);
}
