#include "catch.hpp"

#include "double-conversion/double-conversion.h"

#include "grisu2.h"
#include "grisu2b.h"
#include "grisu3.h"
#include "ryu_32.h"
#include "ryu_64.h"
#include "schubfach_32.h"
#include "schubfach_64.h"

#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>

#include "scan_number.h"

#define TEST_OPTIMAL 0

constexpr int BufSize = 64;

//==================================================================================================
//
//==================================================================================================

static char* ftoa_double_conversion(char* buf, int buflen, float value)
{
    using namespace double_conversion;

    const auto& conv = DoubleToStringConverter::EcmaScriptConverter();
    StringBuilder builder(buf, buflen);
    conv.ToShortestSingle(value, &builder);
    return buf + builder.position();
}

static char* dtoa_double_conversion(char* buf, int buflen, double value)
{
    using namespace double_conversion;

    const auto& conv = DoubleToStringConverter::EcmaScriptConverter();
    StringBuilder builder(buf, buflen);
    conv.ToShortest(value, &builder);
    return buf + builder.position();
}

static float strtof_double_conversion(const char* buf, int len)
{
    double_conversion::StringToDoubleConverter s2d(0, 0.0, 1.0, "inf", "nan");
    int processed_characters_count = 0;
    return s2d.StringToFloat(buf, len, &processed_characters_count);
}

static double strtod_double_conversion(const char* buf, int len)
{
    double_conversion::StringToDoubleConverter s2d(0, 0.0, 1.0, "inf", "nan");
    int processed_characters_count = 0;
    return s2d.StringToDouble(buf, len, &processed_characters_count);
}

//==================================================================================================
//
//==================================================================================================

//struct D2SBase
//{
//    virtual ~D2SBase() {}
//
//    virtual bool Optimal() = 0;
//    virtual const char* Name() = 0;
//    virtual char* Ftoa(char* buf, int buflen, float f) = 0;
//    virtual char* Dtoa(char* buf, int buflen, double f) = 0;
//};

struct D2S_DoubleConversion
{
    bool Optimal() const { return true; }
    const char* Name() const { return "double-conversion"; }
    char* operator()(char* buf, int buflen, float f) { return ftoa_double_conversion(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) { return dtoa_double_conversion(buf, buflen, f); }
};

struct D2S_Grisu2
{
    bool Optimal() const { return false; }
    const char* Name() const { return "grisu2"; }
    char* operator()(char* buf, int buflen, double f) { return grisu2::Dtoa(buf, f); }
};

struct D2S_Grisu2b
{
    bool Optimal() const { return false; }
    const char* Name() const { return "grisu2b"; }
    char* operator()(char* buf, int buflen, double f) { return grisu2b::Dtoa(buf, f); }
};

struct D2S_Grisu3
{
    bool Optimal() const { return true; }
    const char* Name() const { return "grisu3"; }
    char* operator()(char* buf, int buflen, double f) { return grisu3::Dtoa(buf, f); }
};

struct D2S_Ryu
{
    bool Optimal() const { return true; }
    const char* Name() const { return "ryu"; }
    char* operator()(char* buf, int buflen, float f) { return ryu::Ftoa(buf, f); }
    char* operator()(char* buf, int buflen, double f) { return ryu::Dtoa(buf, f); }
};

struct D2S_Schubfach
{
    bool Optimal() const { return true; }
    const char* Name() const { return "schubfach"; }
    char* operator()(char* buf, int buflen, float f) { return schubfach::Ftoa(buf, f); }
    char* operator()(char* buf, int buflen, double f) { return schubfach::Dtoa(buf, f); }
};

//==================================================================================================
//
//==================================================================================================

template <typename Dest, typename Source>
inline Dest ReinterpretBits(Source source)
{
    static_assert(sizeof(Dest) == sizeof(Source), "size mismatch");

    Dest dest;
    std::memcpy(&dest, &source, sizeof(Source));
    return dest;
}

static float MakeSingle(uint32_t sign_bit, uint32_t biased_exponent, uint32_t significand)
{
    assert(sign_bit == 0 || sign_bit == 1);
    assert(biased_exponent <= 0xFF);
    assert(significand <= 0x007FFFFF);

    uint32_t bits = 0;

    bits |= sign_bit << 31;
    bits |= biased_exponent << 23;
    bits |= significand;

    return ReinterpretBits<float>(bits);
}

static float MakeSingle(uint64_t f, int e)
{
    constexpr uint64_t kHiddenBit = 0x00800000;
    constexpr uint64_t kSignificandMask = 0x007FFFFF;
    constexpr int kPhysicalSignificandSize = 23;  // Excludes the hidden bit.
    constexpr int kExponentBias = 0x7F + kPhysicalSignificandSize;
    constexpr int kDenormalExponent = -kExponentBias + 1;
    constexpr int kMaxExponent = 0xFF - kExponentBias;

    assert(f <= kHiddenBit + kSignificandMask);
    if (e >= kMaxExponent) {
        return std::numeric_limits<float>::infinity();
    }
    if (e < kDenormalExponent) {
        return 0.0;
    }
    while (e > kDenormalExponent && (f & kHiddenBit) == 0) {
        f <<= 1;
        e--;
    }

    uint64_t biased_exponent;
    if (e == kDenormalExponent && (f & kHiddenBit) == 0)
        biased_exponent = 0;
    else
        biased_exponent = static_cast<uint64_t>(e + kExponentBias);

    uint64_t bits = (f & kSignificandMask) | (biased_exponent << kPhysicalSignificandSize);

    return ReinterpretBits<float>(static_cast<uint32_t>(bits));
}

static double MakeDouble(uint64_t sign_bit, uint64_t biased_exponent, uint64_t significand)
{
    assert(sign_bit == 0 || sign_bit == 1);
    assert(biased_exponent <= 0x7FF);
    assert(significand <= 0x000FFFFFFFFFFFFF);

    uint64_t bits = 0;

    bits |= sign_bit << 63;
    bits |= biased_exponent << 52;
    bits |= significand;

    return ReinterpretBits<double>(bits);
}

static double MakeDouble(uint64_t f, int e)
{
    constexpr uint64_t kHiddenBit = 0x0010000000000000;
    constexpr uint64_t kSignificandMask = 0x000FFFFFFFFFFFFF;
    constexpr int kPhysicalSignificandSize = 52;  // Excludes the hidden bit.
    constexpr int kExponentBias = 0x3FF + kPhysicalSignificandSize;
    constexpr int kDenormalExponent = -kExponentBias + 1;
    constexpr int kMaxExponent = 0x7FF - kExponentBias;

    assert(f <= kHiddenBit + kSignificandMask);
    if (e >= kMaxExponent) {
        return std::numeric_limits<double>::infinity();
    }
    if (e < kDenormalExponent) {
        return 0.0;
    }
    while (e > kDenormalExponent && (f & kHiddenBit) == 0) {
        f <<= 1;
        e--;
    }

    uint64_t biased_exponent;
    if (e == kDenormalExponent && (f & kHiddenBit) == 0)
        biased_exponent = 0;
    else
        biased_exponent = static_cast<uint64_t>(e + kExponentBias);

    uint64_t bits = (f & kSignificandMask) | (biased_exponent << kPhysicalSignificandSize);

    return ReinterpretBits<double>(bits);
}

//==================================================================================================
//
//==================================================================================================

template <typename Converter>
static void CheckSingle(Converter d2s, float f0)
{
    CAPTURE(d2s.Name());
    CAPTURE(f0);

    // Dtoa
    char buf0[BufSize];
    char* end0 = d2s(buf0, BufSize, f0);
    *end0 = '\0';
    const int length0 = static_cast<int>(end0 - buf0);

    CAPTURE(buf0);

    // Strtod
    const float f1 = strtof_double_conversion(buf0, length0);

    const uint32_t bits0 = ReinterpretBits<uint32_t>(f0);
    const uint32_t bits1 = ReinterpretBits<uint32_t>(f1);
    CAPTURE(bits0);
    CAPTURE(bits1);
    CHECK(bits0 == bits1);

    char buf1[BufSize];
    char* end1 = ftoa_double_conversion(buf1, BufSize, f0);
    *end1 = '\0';

    CAPTURE(buf1);

    if (d2s.Optimal())
    {
        const auto num0 = ScanNumber(buf0, end0);
        const auto num1 = ScanNumber(buf1, end1);
        CHECK(num0.digits == num1.digits);
    }
    else
    {
#if TEST_OPTIMAL
        const auto num0 = ScanNumber(buf0, end0);
        const auto num1 = ScanNumber(buf1, end1);
        if (num0.digits.size() != num1.digits.size())
        {
            printf("%s: not short [0x%08X]\n  actual:   %s\n  expected: %s\n", Converter::Name(), bits0, num0.digits.c_str(), num1.digits.c_str());
        }
        else if (num0.digits != num1.digits)
        {
            printf("%s: not optimal [0x%08X]\n  actual:   %s\n  expected: %s\n", Converter::Name(), bits0, num0.digits.c_str(), num1.digits.c_str());
        }
#endif
    }
}

static void CheckSingle(float f)
{
    CheckSingle(D2S_DoubleConversion{}, f);
    CheckSingle(D2S_Ryu{}, f);
    CheckSingle(D2S_Schubfach{}, f);
}

inline void CheckSingleBits(uint32_t bits)
{
    CheckSingle(ReinterpretBits<float>(bits));
}

template <typename Converter>
static void CheckSingle(Converter d2s, float value, const std::string& expected)
{
    if (!d2s.Optimal())
    {
        CheckSingle(value);
        return;
    }

    char buf[BufSize];
    char* end = d2s(buf, BufSize, value);

    const auto num_actual = ScanNumber(buf, end);
    const auto num_expected = ScanNumber(expected);

    CAPTURE(d2s.Name());
    CAPTURE(value);
    CHECK(num_actual.digits == num_expected.digits);
    CHECK(num_actual.exponent == num_expected.exponent);
}

static void CheckSingle(float value, const std::string& expected)
{
    CheckSingle(D2S_DoubleConversion{}, value, expected);
    CheckSingle(D2S_Ryu{}, value, expected);
    CheckSingle(D2S_Schubfach{}, value, expected);
}

static void CheckSingleBits(uint32_t bits, const std::string& expected)
{
    CheckSingle(ReinterpretBits<float>(bits), expected);
}

template <typename Converter>
static void CheckDouble(Converter d2s, double f0)
{
    CAPTURE(d2s.Name());
    CAPTURE(f0);

    // Dtoa
    char buf0[BufSize];
    char* end0 = d2s(buf0, BufSize, f0);
    *end0 = '\0';
    const int length0 = static_cast<int>(end0 - buf0);

    CAPTURE(buf0);

    // Strtod
    const double f1 = strtod_double_conversion(buf0, length0);

    const uint64_t bits0 = ReinterpretBits<uint64_t>(f0);
    const uint64_t bits1 = ReinterpretBits<uint64_t>(f1);
    CAPTURE(bits0);
    CAPTURE(bits1);
    CHECK(bits0 == bits1);

    char buf1[BufSize];
    char* end1 = dtoa_double_conversion(buf1, BufSize, f0);
    *end1 = '\0';

    CAPTURE(buf1);

    if (d2s.Optimal())
    {
        const auto num0 = ScanNumber(buf0, end0);
        const auto num1 = ScanNumber(buf1, end1);
        CHECK(num0.digits == num1.digits);
    }
    else
    {
#if TEST_OPTIMAL
        const auto num0 = ScanNumber(buf0, end0);
        const auto num1 = ScanNumber(buf1, end1);
        if (num0.digits.size() != num1.digits.size())
        {
            printf("%s: not short [0x%016llX]\n  actual:   %s\n  expected: %s\n", Converter::Name(), bits0, num0.digits.c_str(), num1.digits.c_str());
        }
        else if (num0.digits != num1.digits)
        {
            printf("%s: not optimal [0x%016llX]\n  actual:   %s\n  expected: %s\n", Converter::Name(), bits0, num0.digits.c_str(), num1.digits.c_str());
        }
#endif
    }
}

static void CheckDouble(double f)
{
    CheckDouble(D2S_DoubleConversion{}, f);
    CheckDouble(D2S_Grisu2{}, f);
    CheckDouble(D2S_Grisu2b{}, f);
    CheckDouble(D2S_Grisu3{}, f);
    CheckDouble(D2S_Ryu{}, f);
    CheckDouble(D2S_Schubfach{}, f);
}

inline void CheckDoubleBits(uint64_t bits)
{
    CheckDouble(ReinterpretBits<double>(bits));
}

template <typename Converter>
static void CheckDouble(Converter d2s, double value, const std::string& expected)
{
    if (!d2s.Optimal())
    {
        CheckDouble(value);
        return;
    }

    char buf[BufSize];
    char* end = d2s(buf, BufSize, value);

    const auto num_actual = ScanNumber(buf, end);
    const auto num_expected = ScanNumber(expected);

    CAPTURE(d2s.Name());
    CAPTURE(value);
    CAPTURE(expected);
    CHECK(num_actual.digits == num_expected.digits);
    CHECK(num_actual.exponent == num_expected.exponent);
}

static void CheckDouble(double value, const std::string& expected)
{
    CheckDouble(D2S_DoubleConversion{}, value, expected);
    CheckDouble(D2S_Grisu2{}, value, expected);
    CheckDouble(D2S_Grisu2b{}, value, expected);
    CheckDouble(D2S_Grisu3{}, value, expected);
    CheckDouble(D2S_Ryu{}, value, expected);
    CheckDouble(D2S_Schubfach{}, value, expected);
}

static void CheckDoubleBits(uint64_t bits, const std::string& expected)
{
    CheckDouble(ReinterpretBits<double>(bits), expected);
}

//==================================================================================================
//
//==================================================================================================
TEST_CASE("000")
{
    //schubfach::DumpG();
    //std::exit(0);
}

//==================================================================================================
//
//==================================================================================================

TEST_CASE("Single")
{
    CheckSingle(MakeSingle(0,   0, 0x00000000), "0"            ); // +0
    CheckSingle(MakeSingle(0,   0, 0x00000001), "1e-45"        ); // min denormal
    CheckSingle(MakeSingle(0,   0, 0x007FFFFF), "1.1754942e-38"); // max denormal
    CheckSingle(MakeSingle(0,   1, 0x00000000), "1.1754944e-38"); // min normal
    CheckSingle(MakeSingle(0,   1, 0x00000001), "1.1754945e-38");
    CheckSingle(MakeSingle(0,   1, 0x007FFFFF), "2.3509886e-38");
    CheckSingle(MakeSingle(0,   2, 0x00000000), "2.3509887e-38");
    CheckSingle(MakeSingle(0,   2, 0x00000001), "2.350989e-38" );
    CheckSingle(MakeSingle(0,  24, 0x00000000), "9.8607613e-32"); // fail if no special case in normalized boundaries
    CheckSingle(MakeSingle(0,  30, 0x00000000), "6.3108872e-30"); // fail if no special case in normalized boundaries
    CheckSingle(MakeSingle(0,  31, 0x00000000), "1.2621775e-29"); // fail if no special case in normalized boundaries
    CheckSingle(MakeSingle(0,  57, 0x00000000), "8.4703295e-22"); // fail if no special case in normalized boundaries
    CheckSingle(MakeSingle(0, 254, 0x007FFFFE), "3.4028233e+38");
    CheckSingle(MakeSingle(0, 254, 0x007FFFFF), "3.4028235e+38"); // max normal
}

TEST_CASE("Single - Boundaries")
{
    for (uint32_t e = 2; e < 254; ++e)
    {
        CAPTURE(e);
        //CheckSingle(MakeSingle(0, e-1, 0x007FFFFE));
        CheckSingle(MakeSingle(0, e-1, 0x007FFFFF));
        CheckSingle(MakeSingle(0, e,   0x00000000));
        //CheckSingle(MakeSingle(0, e,   0x00000001));
        //CheckSingle(MakeSingle(0, e,   0x00000002));
        //CheckSingle(MakeSingle(0, e,   0x002AAAAA));
        //CheckSingle(MakeSingle(0, e,   0x00555555));
        //CheckSingle(MakeSingle(0, e,   0x00400000));
    }
}

TEST_CASE("Single - Paxson, Kahan")
{
    // V. Paxson and W. Kahan, "A Program for Testing IEEE Binary-Decimal Conversion", manuscript, May 1991,
    // ftp://ftp.ee.lbl.gov/testbase-report.ps.Z    (report)
    // ftp://ftp.ee.lbl.gov/testbase.tar.Z          (program)

    // Table 16: Stress Inputs for Converting 24-bit Binary to Decimal, < 1/2 ULP
    CheckSingle(MakeSingle(12676506, -102), "2.5e-24"        ); // digits  1, bits 32
    CheckSingle(MakeSingle(12676506, -103), "1.25e-24"       ); // digits  2, bits 29
    CheckSingle(MakeSingle(15445013,  +86), "1.195e+33"      ); // digits  3, bits 34
    CheckSingle(MakeSingle(13734123, -138), "3.9415e-35"     ); // digits  4, bits 32
    CheckSingle(MakeSingle(12428269, -130), "9.13085e-33"    ); // digits  5, bits 30
    CheckSingle(MakeSingle(15334037, -146), "1.719005e-37"   ); // digits  6, bits 31
    CheckSingle(MakeSingle(11518287,  -41), "0.0000052379105"); // digits  7, bits 30
    CheckSingle(MakeSingle(12584953, -145), "2.821644e-37"   ); // digits  8, bits 31
    CheckSingle(MakeSingle(15961084, -125), "3.7524328e-31"  ); // digits  9, bits 32
    CheckSingle(MakeSingle(14915817, -146), "1.6721209e-37"  ); // digits 10, bits 31
    CheckSingle(MakeSingle(10845484, -102), "2.1388946e-24"  ); // digits 11, bits 30
    CheckSingle(MakeSingle(16431059,  -61), "7.125836e-12"   ); // digits 12, bits 29

    // Table 17: Stress Inputs for Converting 24-bit Binary to Decimal, > 1/2 ULP
    CheckSingle(MakeSingle(16093626,  +69), "9.5e+27"              ); // digits  1, bits 30
    CheckSingle(MakeSingle( 9983778,  +25), "335000000000000"      ); // digits  2, bits 31
    CheckSingle(MakeSingle(12745034, +104), "2.585e+38"            ); // digits  3, bits 31
    CheckSingle(MakeSingle(12706553,  +72), "6.0005e+28"           ); // digits  4, bits 31
    CheckSingle(MakeSingle(11005028,  +45), "387205000000000000000"); // digits  5, bits 30
    CheckSingle(MakeSingle(15059547,  +71), "3.555835e+28"         ); // digits  6, bits 31
    CheckSingle(MakeSingle(16015691,  -99), "2.5268305e-23"        ); // digits  7, bits 29
    CheckSingle(MakeSingle( 8667859,  +56), "6.245851e+23"         ); // digits  8, bits 33
    CheckSingle(MakeSingle(14855922,  -82), "3.0721327e-18"        ); // digits  9, bits 35
    CheckSingle(MakeSingle(14855922,  -83), "1.5360663e-18"        ); // digits 10, bits 33
    CheckSingle(MakeSingle(10144164, -110), "7.81478e-27"          ); // digits 11, bits 32
    CheckSingle(MakeSingle(13248074,  +95), "5.2481028e+35"        ); // digits 12, bits 33
}

TEST_CASE("Single - Regression")
{
    CheckSingle(7.0385307e-26f);

    CheckSingleBits(0x4C000009, "33554468");
    CheckSingleBits(0x4C800009, "67108936");
    CheckSingleBits(0x4D00001D, "134218190");
    CheckSingleBits(0x4D80001D, "268436380");
    CheckSingleBits(0x4E00001D, "536872770");
    CheckSingleBits(0x4E80004F, "1073751900");
    CheckSingleBits(0x4F00004F, "2147503900");
    CheckSingleBits(0x4F80004F, "4295007700");
    CheckSingleBits(0x50000437, "8591039000");
    CheckSingleBits(0x50800437, "17182079000");
    CheckSingleBits(0x51000437, "34364158000");
    CheckSingleBits(0x51800437, "68728316000");
    CheckSingleBits(0x52000DFB, "137497590000");
    CheckSingleBits(0x52800DFB, "274995180000");
    CheckSingleBits(0x53000DFB, "549990370000");
    CheckSingleBits(0x53802665, "1100799900000");
    CheckSingleBits(0x54002665, "2201599900000");
    CheckSingleBits(0x54802665, "4403199700000");
    CheckSingleBits(0x55002665, "8806399000000");
    CheckSingleBits(0x55802665, "17612799000000");
    CheckSingleBits(0x56002665, "35225598000000");
    CheckSingleBits(0x56802665, "70451196000000");
    CheckSingleBits(0x57002665, "140902390000000");
    CheckSingleBits(0x57802665, "281804780000000");
    CheckSingleBits(0x58002665, "563609570000000");
    CheckSingleBits(0x58A3E9AB, "1441791900000000");
    CheckSingleBits(0x5923E9AB, "2883583900000000");
    CheckSingleBits(0x59A3E9AB, "5767167700000000");
    CheckSingleBits(0x5A5F8475, "15728639000000000");
    CheckSingleBits(0x5ADF8475, "31457279000000000");
    CheckSingleBits(0x5B5F8475, "62914558000000000");
    CheckSingleBits(0x5BDF8475, "125829116000000000");

    CheckSingleBits(0x4D00001E, "134218200");
    CheckSingleBits(0x4D80001E, "268436400");
    CheckSingleBits(0x4E800050, "1073752000");
    CheckSingleBits(0x4F800050, "4295008000");
    CheckSingleBits(0x50000438, "8591040000");
    CheckSingleBits(0x52000DFC, "137497600000");
    CheckSingleBits(0x52800DFC, "274995200000");
    CheckSingleBits(0x53802666, "1100800000000");
    CheckSingleBits(0x54802666, "4403200000000");
    CheckSingleBits(0x55002666, "8806400000000");
    CheckSingleBits(0x57002666, "140902400000000");
    CheckSingleBits(0x57802666, "281804800000000");
    CheckSingleBits(0x58A3E9AC, "1441792000000000");
    CheckSingleBits(0x59A3E9AC, "5767168000000000");
    CheckSingleBits(0x5A5F8476, "15728640000000000");
}

TEST_CASE("Single - Ryu")
{
    CheckSingleBits(0x3800000A, "0.000030517615");
    CheckSingleBits(0x3880001E, "0.000061035375");
    CheckSingleBits(0x390000FB, "0.00012207397");
    CheckSingleBits(0x39800091, "0.00024414485");
    CheckSingleBits(0x3A000024, "0.00048828335");
    CheckSingleBits(0x3A80020F, "0.0009766239");
    CheckSingleBits(0x3B000020, "0.0019531325");
    CheckSingleBits(0x3B800007, "0.0039062533");
    CheckSingleBits(0x3C800028, "0.015625075");
    CheckSingleBits(0x3D000014, "0.031250075");
    CheckSingleBits(0x3D80000A, "0.062500075");
    CheckSingleBits(0x3E000032, "0.12500075");
    CheckSingleBits(0x3E800019, "0.25000075");
    CheckSingleBits(0x3F000024, "0.50000215");
    CheckSingleBits(0x3F80008A, "1.0000165");
    CheckSingleBits(0x40000045, "2.0000165");
    CheckSingleBits(0x40800020, "4.0000153");
    CheckSingleBits(0x418000E8, "16.000443");
    CheckSingleBits(0x42000074, "32.000443");
    CheckSingleBits(0x4280004A, "64.000565");
    CheckSingleBits(0x4300003A, "128.00089");
    CheckSingleBits(0x43800091, "256.00443");
    CheckSingleBits(0x4400003F, "512.00385");
    CheckSingleBits(0x448000B3, "1024.0219");
    CheckSingleBits(0x45000111, "2048.0667");
    CheckSingleBits(0x45800015, "4096.0103");

    CheckSingleBits(0x39800000, "0.00024414062");
    CheckSingleBits(0x3B200000, "0.0024414062");
    CheckSingleBits(0x3B900000, "0.0043945312");
    CheckSingleBits(0x3C880000, "0.016601562");
    CheckSingleBits(0x3D040000, "0.032226562");
    CheckSingleBits(0x3E020000, "0.12695312");
    CheckSingleBits(0x3E810000, "0.25195312");
    CheckSingleBits(0x3F808000, "1.0039062");
    CheckSingleBits(0x40004000, "2.0039062");
    CheckSingleBits(0x40802000, "4.0039062");
    CheckSingleBits(0x41801000, "16.007812");
    CheckSingleBits(0x42000800, "32.007812");
    CheckSingleBits(0x43000400, "128.01562");
    CheckSingleBits(0x43800200, "256.01562");
    CheckSingleBits(0x44800100, "1024.0312");
    CheckSingleBits(0x45000080, "2048.0312");
    CheckSingleBits(0x45800040, "4096.0312");
    CheckSingleBits(0x46800020, "16384.062");
    CheckSingleBits(0x47000010, "32768.062");
    CheckSingleBits(0x48000008, "131072.12");
    CheckSingleBits(0x48800004, "262144.12");
    CheckSingleBits(0x49800002, "1048576.2");

    CheckSingleBits(0x4F80001E, "4294982700");
    CheckSingleBits(0x51000002, "34359747000");
    CheckSingleBits(0x51800142, "68722115000");
    CheckSingleBits(0x5300016F, "549779870000");
    CheckSingleBits(0x54801E57, "4402118700000");
    CheckSingleBits(0x56004279, "35255747000000");
    CheckSingleBits(0x5680DC48, "70841795000000");
    CheckSingleBits(0x580214C8, "572103070000000");
    CheckSingleBits(0x5984CC95, "4672454700000000");
    CheckSingleBits(0x5B00BEFC, "36238787000000000");
    CheckSingleBits(0x5B99C7AD, "86570435000000000");
}

TEST_CASE("Double")
{
    CheckDouble(MakeDouble(20, -1074), "1e-322");

    CheckDouble(MakeDouble(0,    0, 0x0000000000000000), "0"                      ); // +0
    CheckDouble(MakeDouble(0,    0, 0x0000000000000001), "5e-324"                 ); // min denormal
    CheckDouble(MakeDouble(0,    0, 0x000FFFFFFFFFFFFF), "2.225073858507201e-308" ); // max denormal
    CheckDouble(MakeDouble(0,    1, 0x0000000000000000), "2.2250738585072014e-308"); // min normal
    CheckDouble(MakeDouble(0,    1, 0x0000000000000001), "2.225073858507202e-308" );
    CheckDouble(MakeDouble(0,    1, 0x000FFFFFFFFFFFFF), "4.4501477170144023e-308");
    CheckDouble(MakeDouble(0,    2, 0x0000000000000000), "4.450147717014403e-308" );
    CheckDouble(MakeDouble(0,    2, 0x0000000000000001), "4.450147717014404e-308" );
    CheckDouble(MakeDouble(0,    4, 0x0000000000000000), "1.7800590868057611e-307"); // fail if no special case in normalized boundaries
    CheckDouble(MakeDouble(0,    5, 0x0000000000000000), "3.5601181736115222e-307"); // fail if no special case in normalized boundaries
    CheckDouble(MakeDouble(0,    6, 0x0000000000000000), "7.120236347223045e-307" ); // fail if no special case in normalized boundaries
    CheckDouble(MakeDouble(0,   10, 0x0000000000000000), "1.1392378155556871e-305"); // fail if no special case in normalized boundaries
    CheckDouble(MakeDouble(0, 2046, 0x000FFFFFFFFFFFFE), "1.7976931348623155e+308");
    CheckDouble(MakeDouble(0, 2046, 0x000FFFFFFFFFFFFF), "1.7976931348623157e+308"); // max normal
}

//inline void CheckPermutations(uint64_t e, uint64_t f)
//{
//    uint64_t p = next_permutation(f);
//
////  while (p != f && p <= 0x000FFFFFFFFFFFFF)
//    for (int i = 0; i < 10 && p <= 0x000FFFFFFFFFFFFF; ++i)
//    {
//        CAPTURE(p);
//        CheckDouble(MakeDouble(0, e, p));
//        p = next_permutation(p);
//    }
//}

TEST_CASE("Double - Boundaries")
{
    //CheckDouble(MakeDouble(0, 1080-1, 0x000FFFFFFFFFFFFF));

    for (uint64_t e = 2; e < 2046; ++e)
    {
        CAPTURE(e);
        //CheckDouble(MakeDouble(0, e-1, 0x000FFFFFFFFFFFFE));
        CheckDouble(MakeDouble(0, e-1, 0x000FFFFFFFFFFFFF));
        CheckDouble(MakeDouble(0, e,   0x0000000000000000));
        //CheckDouble(MakeDouble(0, e,   0x0000000000000001));
        //CheckDouble(MakeDouble(0, e,   0x0000000000000002));
        //CheckDouble(MakeDouble(0, e,   0x0005555555555555));
        //CheckDouble(MakeDouble(0, e,   0x000AAAAAAAAAAAAA));
        //CheckDouble(MakeDouble(0, e,   0x0008000000000000));
        //CheckDouble(MakeDouble(0, e,   0x000878678326EAC9));
    }
}

TEST_CASE("Double - Paxson, Kahan")
{
    // V. Paxson and W. Kahan, "A Program for Testing IEEE Binary-Decimal Conversion", manuscript, May 1991,
    // ftp://ftp.ee.lbl.gov/testbase-report.ps.Z    (report)
    // ftp://ftp.ee.lbl.gov/testbase.tar.Z          (program)

    // Table 3: Stress Inputs for Converting 53-bit Binary to Decimal, < 1/2 ULP
    CheckDouble( MakeDouble(8511030020275656,  -342), "9.5e-88"                 ); // digits  1, bits 63
    CheckDouble( MakeDouble(5201988407066741,  -824), "4.65e-233"               ); // digits  2, bits 63
    CheckDouble( MakeDouble(6406892948269899,  +237), "1.415e+87"               ); // digits  3, bits 62 (D3. [Calculate q'.] One correction step)
    CheckDouble( MakeDouble(8431154198732492,   +72), "3.9815e+37"              ); // digits  4, bits 61 (D3. [Calculate q'.] One correction step)
    CheckDouble( MakeDouble(6475049196144587,   +99), "4.10405e+45"             ); // digits  5, bits 64 (D3. [Calculate q'.] One correction step)
    CheckDouble( MakeDouble(8274307542972842,  +726), "2.920845e+234"           ); // digits  6, bits 64
    CheckDouble( MakeDouble(5381065484265332,  -456), "2.8919465e-122"          ); // digits  7, bits 64
    CheckDouble( MakeDouble(6761728585499734, -1057), "4.37877185e-303"         ); // digits  8, bits 64
    CheckDouble( MakeDouble(7976538478610756,  +376), "1.227701635e+129"        ); // digits  9, bits 67 (D6. [Add back.])
    CheckDouble( MakeDouble(5982403858958067,  +377), "1.8415524525e+129"       ); // digits 10, bits 63
    CheckDouble( MakeDouble(5536995190630837,   +93), "5.48357443505e+43"       ); // digits 11, bits 63
    CheckDouble( MakeDouble(7225450889282194,  +710), "3.891901811465e+229"     ); // digits 12, bits 66 (D6. [Add back.])
    CheckDouble( MakeDouble(7225450889282194,  +709), "1.9459509057325e+229"    ); // digits 13, bits 64
    CheckDouble( MakeDouble(8703372741147379,  +117), "1.44609583816055e+51"    ); // digits 14, bits 66
    CheckDouble( MakeDouble(8944262675275217, -1001), "4.173677474585315e-286"  ); // digits 15, bits 63
    CheckDouble( MakeDouble(7459803696087692,  -707), "1.1079507728788885e-197" ); // digits 16, bits 63
    CheckDouble( MakeDouble(6080469016670379,  -381), "1.234550136632744e-99"   ); // digits 17, bits 62
    CheckDouble( MakeDouble(8385515147034757,  +721), "9.25031711960365e+232"   ); // digits 18, bits 64
    CheckDouble( MakeDouble(7514216811389786,  -828), "4.19804715028489e-234"   ); // digits 19, bits 64
    CheckDouble( MakeDouble(8397297803260511,  -345), "1.1716315319786511e-88"  ); // digits 20, bits 64
    CheckDouble( MakeDouble(6733459239310543,  +202), "4.328100728446125e+76"   ); // digits 21, bits 63
    CheckDouble( MakeDouble(8091450587292794,  -473), "3.317710118160031e-127"  ); // digits 22, bits 63

    // Table 4: Stress Inputs for Converting 53-bit Binary to Decimal, > 1/2 ULP
    CheckDouble( MakeDouble(6567258882077402, +952), "2.5e+302"                ); // digits  1, bits 62
    CheckDouble( MakeDouble(6712731423444934, +535), "7.55e+176"               ); // digits  2, bits 65
    CheckDouble( MakeDouble(6712731423444934, +534), "3.775e+176"              ); // digits  3, bits 63
    CheckDouble( MakeDouble(5298405411573037, -957), "4.3495e-273"             ); // digits  4, bits 62
    CheckDouble( MakeDouble(5137311167659507, -144), "2.30365e-28"             ); // digits  5, bits 61
    CheckDouble( MakeDouble(6722280709661868, +363), "1.263005e+125"           ); // digits  6, bits 64
    CheckDouble( MakeDouble(5344436398034927, -169), "7.1422105e-36"           ); // digits  7, bits 61
    CheckDouble( MakeDouble(8369123604277281, -853), "1.39345735e-241"         ); // digits  8, bits 65
    CheckDouble( MakeDouble(8995822108487663, -780), "1.414634485e-219"        ); // digits  9, bits 63
    CheckDouble( MakeDouble(8942832835564782, -383), "4.5392779195e-100"       ); // digits 10, bits 66
    CheckDouble( MakeDouble(8942832835564782, -384), "2.26963895975e-100"      ); // digits 11, bits 64
    CheckDouble( MakeDouble(8942832835564782, -385), "1.134819479875e-100"     ); // digits 12, bits 61
    CheckDouble( MakeDouble(6965949469487146, -249), "7.7003665618895e-60"     ); // digits 13, bits 67
    CheckDouble( MakeDouble(6965949469487146, -250), "3.85018328094475e-60"    ); // digits 14, bits 65
    CheckDouble( MakeDouble(6965949469487146, -251), "1.925091640472375e-60"   ); // digits 15, bits 63
    CheckDouble( MakeDouble(7487252720986826, +548), "6.8985865317742005e+180" ); // digits 16, bits 63
    CheckDouble( MakeDouble(5592117679628511, +164), "1.3076622631878654e+65"  ); // digits 17, bits 65
    CheckDouble( MakeDouble(8887055249355788, +665), "1.3605202075612124e+216" ); // digits 18, bits 67
    CheckDouble( MakeDouble(6994187472632449, +690), "3.5928102174759597e+223" ); // digits 19, bits 64
    CheckDouble( MakeDouble(8797576579012143, +588), "8.912519771248455e+192"  ); // digits 20, bits 62
    CheckDouble( MakeDouble(7363326733505337, +272), "5.5876975736230114e+97"  ); // digits 21, bits 61
    CheckDouble( MakeDouble(8549497411294502, -448), "1.1762578307285404e-119" ); // digits 22, bits 66
}

TEST_CASE("Double - Regression")
{
    CheckDouble(1.5745340942675811e+257, "1.574534094267581e+257");
    CheckDouble(1.6521200219181297e-180, "1.6521200219181297e-180");
    CheckDouble(4.6663180925160944e-302, "4.6663180925160944e-302");

    CheckDouble(18776091678571.0 / 64.0);

    CheckDouble(2.0919495182368195e+19, "2.0919495182368195e+19");
    CheckDouble(2.6760179287532483e+19, "2.6760179287532483e+19");
    CheckDouble(3.2942957306323907e+19, "3.2942957306323907e+19");
    CheckDouble(3.9702293349085635e+19, "3.9702293349085635e+19");
    CheckDouble(4.0647939013152195e+19, "4.0647939013152195e+19");

    CheckDouble(1.8014398509481984E16, "1.8014398509481984E16");
    CheckDouble(1.8014398509481985E16, "1.8014398509481984E16");
}

// Some numbers to check different code paths in grisu2::Dtoa
TEST_CASE("Double - Grisu2 code paths")
{
    CheckDoubleBits(0x40C3880000000000, "10000"                 );
    CheckDoubleBits(0x41324F8000000000, "1200000"               );
    CheckDoubleBits(0x0000000000000001, "5e-324"                ); // DigitGen: exit integral loop
    CheckDoubleBits(0x000FFFFFFFFFFFFF, "2.225073858507201e-308"); // DigitGen: exit fractional loop
    CheckDoubleBits(0x2B70000000000000, "1.82877982605164e-99"  );
    CheckDoubleBits(0x3E13C42855500898, "1.1505466208671903e-9" );
    CheckDoubleBits(0x443E2A6B41CE4B23, "556458931337667200000" );
    CheckDoubleBits(0x404A8475527A8B30, "53.034830388866226"    );
    CheckDoubleBits(0x3F6141F8CE9A7906, "0.0021066531670178605" );
}

TEST_CASE("Double - Round to even")
{
    CheckDouble(1.00000000000000005, "1");
    CheckDouble(1.00000000000000015, "1.0000000000000002"); // 1.000000000000000222...
    CheckDouble(1.99999999999999985, "1.9999999999999998"); // 1.999999999999999777...
    CheckDouble(1.99999999999999995, "2");
    CheckDouble(1125899906842623.75, "1125899906842623.8");
    CheckDouble(1125899906842624.25, "1125899906842624.2");
    CheckDouble(562949953421312.25, "562949953421312.2");

    CheckDouble(2.20781707763671875, "22078170776367188e-16");
    CheckDouble(1.81835174560546875, "18183517456054688e-16");
    CheckDouble(3.94171905517578125, "39417190551757812e-16");
    CheckDouble(3.73860931396484375, "37386093139648438e-16");
    CheckDouble(3.96773529052734375, "39677352905273438e-16");
    CheckDouble(1.32802581787109375, "13280258178710938e-16");
    CheckDouble(3.92096710205078125, "39209671020507812e-16");
    CheckDouble(1.01523590087890625, "10152359008789062e-16");
    CheckDouble(1.33522796630859375, "13352279663085938e-16");
    CheckDouble(1.34452056884765625, "13445205688476562e-16");
    CheckDouble(2.87912750244140625, "28791275024414062e-16");
    CheckDouble(3.69583892822265625, "36958389282226562e-16");
    CheckDouble(1.84534454345703125, "18453445434570312e-16");
    CheckDouble(3.79395294189453125, "37939529418945312e-16");
    CheckDouble(3.21140289306640625, "32114028930664062e-16");
    CheckDouble(2.56597137451171875, "25659713745117188e-16");
    CheckDouble(0.96515655517578125, "9651565551757812e-16");
    CheckDouble(2.70000457763671875, "27000045776367188e-16");
    CheckDouble(0.76709747314453125, "7670974731445312e-16");
    CheckDouble(1.78044891357421875, "17804489135742188e-16");
    CheckDouble(2.62483978271484375, "26248397827148438e-16");
    CheckDouble(1.30529022216796875, "13052902221679688e-16");
    CheckDouble(3.83492279052734375, "38349227905273438e-16");
}

TEST_CASE("Double - Integers")
{
    CheckDouble(1.0, "1");
    CheckDouble(10.0, "10");
    CheckDouble(100.0, "100");
    CheckDouble(1000.0, "1000");
    CheckDouble(10000.0, "10000");
    CheckDouble(100000.0, "100000");
    CheckDouble(1000000.0, "1000000");
    CheckDouble(10000000.0, "10000000");
    CheckDouble(100000000.0, "100000000");
    CheckDouble(1000000000.0, "1000000000");
    CheckDouble(10000000000.0, "10000000000");
    CheckDouble(100000000000.0, "100000000000");
    CheckDouble(1000000000000.0, "1000000000000");
    CheckDouble(10000000000000.0, "10000000000000");
    CheckDouble(100000000000000.0, "100000000000000");
    CheckDouble(1000000000000000.0, "1000000000000000");
    CheckDouble(9007199254740000.0, "9007199254740000");
    CheckDouble(9007199254740992.0, "9007199254740992");
    CheckDouble(1e+22, "1e+22");
    CheckDouble(1e+23, "1e+23");
}

TEST_CASE("Double - Looks like pow5")
{
    // From
    // https://github.com/ulfjack/ryu/blob/master/ryu/tests/d2s_test.cc

    // These numbers have a mantissa that is a multiple of the largest power of 5 that fits,
    // and an exponent that causes the computation for q to result in 22, which is a corner
    // case for Ryu.
    CheckDoubleBits(0x4830F0CF064DD592, "5.764607523034235e+39");
    CheckDoubleBits(0x4840F0CF064DD592, "1.152921504606847e+40");
    CheckDoubleBits(0x4850F0CF064DD592, "2.305843009213694e+40");
}
