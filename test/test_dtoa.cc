#include "drachennest.h"
#include "double-conversion/double-conversion.h"

#include "ryu_32.h"
#include "ryu_64.h"

#include "catch.hpp"

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
    char* operator()(char* buf, int /*buflen*/, float f) { return drachennest::ftoa_grisu2(buf, f); }
    char* operator()(char* buf, int /*buflen*/, double f) { return drachennest::dtoa_grisu2(buf, f); }
};

struct D2S_Grisu3
{
    bool Optimal() const { return true; }
    const char* Name() const { return "grisu3-dragon4"; }
    char* operator()(char* buf, int /*buflen*/, float f) { return drachennest::ftoa_grisu3(buf, f); }
    char* operator()(char* buf, int /*buflen*/, double f) { return drachennest::dtoa_grisu3(buf, f); }
};

struct D2S_Ryu
{
    static_assert(BufSize >= ryu::DtoaMinBufferLength, "");

    bool Optimal() const { return true; }
    const char* Name() const { return "ryu"; }
    char* operator()(char* buf, int buflen, float f) { return ryu::Ftoa(buf, f); }
    char* operator()(char* buf, int buflen, double f) { return ryu::Dtoa(buf, f); }
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
    //CheckSingle(D2S_DoubleConversion{}, f);
    CheckSingle(D2S_Grisu2{}, f);
    CheckSingle(D2S_Grisu3{}, f);
    CheckSingle(D2S_Ryu{}, f);
    //CheckSingle(D2S_Swift{}, f);
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
    CheckSingle(D2S_Grisu2{}, value, expected);
    CheckSingle(D2S_Grisu3{}, value, expected);
    CheckSingle(D2S_Ryu{}, value, expected);
    //CheckSingle(D2S_Swift{}, value, expected);
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
    //CheckDouble(D2S_DoubleConversion{}, f);
    CheckDouble(D2S_Grisu2{}, f);
    CheckDouble(D2S_Grisu3{}, f);
    CheckDouble(D2S_Ryu{}, f);
    //CheckDouble(D2S_Swift{}, f);
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
    CHECK(num_actual.digits == num_expected.digits);
    CHECK(num_actual.exponent == num_expected.exponent);
}

static void CheckDouble(double value, const std::string& expected)
{
    CheckDouble(D2S_DoubleConversion{}, value, expected);
    CheckDouble(D2S_Grisu2{}, value, expected);
    CheckDouble(D2S_Grisu3{}, value, expected);
    CheckDouble(D2S_Ryu{}, value, expected);
    //CheckDouble(D2S_Swift{}, value, expected);
}

static void CheckDoubleBits(uint64_t bits, const std::string& expected)
{
    CheckDouble(ReinterpretBits<double>(bits), expected);
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
        CheckSingle(MakeSingle(0, e-1, 0x007FFFFE));
        CheckSingle(MakeSingle(0, e-1, 0x007FFFFF));
        CheckSingle(MakeSingle(0, e,   0x00000000));
        CheckSingle(MakeSingle(0, e,   0x00000001));
        CheckSingle(MakeSingle(0, e,   0x00000002));
        CheckSingle(MakeSingle(0, e,   0x002AAAAA));
        CheckSingle(MakeSingle(0, e,   0x00555555));
        CheckSingle(MakeSingle(0, e,   0x00400000));
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
    for (uint64_t e = 2; e < 2046; ++e)
    {
        CAPTURE(e);
        CheckDouble(MakeDouble(0, e-1, 0x000FFFFFFFFFFFFE));
        CheckDouble(MakeDouble(0, e-1, 0x000FFFFFFFFFFFFF));
        CheckDouble(MakeDouble(0, e,   0x0000000000000000));
        CheckDouble(MakeDouble(0, e,   0x0000000000000001));
        CheckDouble(MakeDouble(0, e,   0x0000000000000002));
        CheckDouble(MakeDouble(0, e,   0x0005555555555555));
        CheckDouble(MakeDouble(0, e,   0x000AAAAAAAAAAAAA));
        CheckDouble(MakeDouble(0, e,   0x0008000000000000));
        CheckDouble(MakeDouble(0, e,   0x000878678326EAC9));
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

TEST_CASE("Double - Ryu")
{
    CheckDoubleBits(0x3E5000000E8D4A52, "1.4901162001641227e-8");
    CheckDoubleBits(0x3E5000B1A2BC2ED2, "1.4903685548744407e-8");
    CheckDoubleBits(0x3E58000000000400, "2.2351741790774873e-8");
    CheckDoubleBits(0x3E6000000000001A, "2.9802322387695485e-8");
    CheckDoubleBits(0x3E60000000000100, "2.9802322387697007e-8");
    CheckDoubleBits(0x3E600000001DCD87, "2.9802322400620235e-8");
    CheckDoubleBits(0x3E60000002E90EF9, "2.9802322710812925e-8");
    CheckDoubleBits(0x3E6000071AFD49A6, "2.9802524336087215e-8");
    CheckDoubleBits(0x3E61000000000080, "3.1664967536927117e-8");
    CheckDoubleBits(0x3E61000000001000, "3.1664967536953375e-8");
    CheckDoubleBits(0x3E62000000000000, "3.3527612686157227e-8");
    CheckDoubleBits(0x3E90000000000296, "2.3841857910159755e-7");
    CheckDoubleBits(0x3E90000000950337, "2.3841857961855367e-7");
    CheckDoubleBits(0x3E9000071AFD498E, "2.3842019468869645e-7");
    CheckDoubleBits(0x3E92000000000100, "2.6822090148927137e-7");
    CheckDoubleBits(0x3E94000000000800, "2.9802322387706155e-7");
    CheckDoubleBits(0x3E94000000001000, "2.9802322387716997e-7");
    CheckDoubleBits(0x3E9878678326EAD2, "3.6463632393692487e-7");
    CheckDoubleBits(0x3EA4000000000000, "5.960464477539062e-7");
    CheckDoubleBits(0x3EAAAAAAAAAAAACA, "7.947285970052117e-7");
    CheckDoubleBits(0x3EB2000000000000, "1.0728836059570312e-6");
    CheckDoubleBits(0x3EB5555555555565, "1.2715657552083367e-6");
    CheckDoubleBits(0x3EC0000000000C3C, "1.9073486328138265e-6");
    CheckDoubleBits(0x3EC0000000000C4D, "1.9073486328138337e-6");
    CheckDoubleBits(0x3EC0000000003D12, "1.9073486328191213e-6");
    CheckDoubleBits(0x3EC00000001DCD78, "1.9073486336396887e-6");
    CheckDoubleBits(0x3EC000016BCC41EA, "1.9073512177519147e-6");
    CheckDoubleBits(0x3EC000071AFD49A3, "1.9073615575095805e-6");
    CheckDoubleBits(0x3EC5555555555563, "2.5431315104166725e-6");
    CheckDoubleBits(0x3ED000071AFD49A6, "3.8147231150191635e-6");
    CheckDoubleBits(0x3ED056BC75E2D632, "3.8954766223196325e-6");
    CheckDoubleBits(0x3ED5555555555566, "5.0862630208333475e-6");
    CheckDoubleBits(0x3F00000000013136, "3.0517578125529457e-5");
    CheckDoubleBits(0x3F0000000005F61D, "3.0517578127647385e-5");
    CheckDoubleBits(0x3F000000001DCD66, "3.0517578138234897e-5");
    CheckDoubleBits(0x3F0000000095031F, "3.0517578191174707e-5");
    CheckDoubleBits(0x3F0000000095033E, "3.0517578191174917e-5");
    CheckDoubleBits(0x3F00002386F27077, "3.0518612100766925e-5");
    CheckDoubleBits(0x3F0000B1A2BC2ED2, "3.0522748003828545e-5");
    CheckDoubleBits(0x3F0AAAAAAAAAAAB4, "5.0862630208333397e-5");
    CheckDoubleBits(0x3F1000000E8D4A58, "6.103515955872255e-5");
    CheckDoubleBits(0x3F101158E460914F, "6.129365019142307e-5");
    CheckDoubleBits(0x3F300000000000BD, "2.4414062500001025e-4");
    CheckDoubleBits(0x3F3AAAAAAAAAAACA, "4.0690104166666837e-4");
    CheckDoubleBits(0x3F40000000000006, "4.882812500000007e-4");
    CheckDoubleBits(0x3F60000000000031, "1.9531250000000213e-3");
    CheckDoubleBits(0x3F600000001DCD74, "1.9531250008470395e-3");
    CheckDoubleBits(0x3F60000002E90EFA, "1.9531250211758363e-3");
    CheckDoubleBits(0x3F6000071AFD4993, "1.9531382348898035e-3");
    CheckDoubleBits(0x3F6000071AFD49A5, "1.9531382348898113e-3");
    CheckDoubleBits(0x3F6000B1A2BC2ED8, "1.9534558722450295e-3");
    CheckDoubleBits(0x3F6003782DACE9E5, "1.9547793612251113e-3");
    CheckDoubleBits(0x3F61B1AE4D6E2EF9, "2.1599201531382587e-3");
    CheckDoubleBits(0x3F65555555555559, "2.6041666666666683e-3");
    CheckDoubleBits(0x3F7000000095033D, "3.9062500084703885e-3");
    CheckDoubleBits(0x3F701158E460916B, "3.9227936122511005e-3");
    CheckDoubleBits(0x3F7056BC75E2D64A, "3.9889680612553245e-3");
    CheckDoubleBits(0x3F8000000000028E, "7.812500000001135e-3");
    CheckDoubleBits(0x3F8000000005F5E8, "7.812500000677639e-3");
    CheckDoubleBits(0x3F9878678326EAD1, "2.3896806125530305e-2");
    CheckDoubleBits(0x3FA000000095033E, "3.1250000067763115e-2");
    CheckDoubleBits(0x3FB000000095034F, "6.250000013552647e-2");
    CheckDoubleBits(0x3FC003782DACE9E3, "1.2510587911840707e-1");
    CheckDoubleBits(0x3FD000016BCC41F1, "2.5000033881317935e-1");
    CheckDoubleBits(0x3FD056BC75E2D638, "2.5529395592033977e-1");
    CheckDoubleBits(0x3FE5555555555569, "6.666666666666689e-1");
    CheckDoubleBits(0x4000000000000018, "2.0000000000000107e+0");
    CheckDoubleBits(0x4000000000000021, "2.0000000000000147e+0");
    CheckDoubleBits(0x400000000000002A, "2.0000000000000187e+0");
    CheckDoubleBits(0x40000000000000DB, "2.0000000000000973e+0");
    CheckDoubleBits(0x4000000000000293, "2.0000000000002927e+0");
    CheckDoubleBits(0x400000000005F5E8, "2.0000000001734755e+0");
    CheckDoubleBits(0x400000000E8D4A83, "2.0000001084202395e+0");
    CheckDoubleBits(0x4000000048C273AC, "2.0000005421010965e+0");
    CheckDoubleBits(0x400000016BCC41EA, "2.0000027105054317e+0");
    CheckDoubleBits(0x4000002386F26FC8, "2.0000677626357835e+0");
    CheckDoubleBits(0x400003782DACE9F1, "2.0016940658945193e+0");
    CheckDoubleBits(0x40001158E460915E, "2.0084703294725577e+0");
    CheckDoubleBits(0x40001158E4609167, "2.0084703294725617e+0");
    CheckDoubleBits(0x400AAAAAAAAAAAB8, "3.3333333333333393e+0");
    CheckDoubleBits(0x4010000000003D0C, "4.0000000000138805e+0");
    CheckDoubleBits(0x401000071AFD49A5, "4.0000271050543335e+0");
    CheckDoubleBits(0x401056BC75E2D651, "4.0847032947254585e+0");
    CheckDoubleBits(0x402056BC75E2D638, "8.169406589450873e+0");
    CheckDoubleBits(0x402555555555555C, "1.0666666666666679e+1");
    CheckDoubleBits(0x405000B1A2BC2ED8, "6.401084202172513e+1");
    CheckDoubleBits(0x40501158E460915E, "6.427105054312185e+1");
    CheckDoubleBits(0x4060000000003D12, "1.2800000000044435e+2");
    CheckDoubleBits(0x4060000000013147, "1.2800000000222119e+2");
    CheckDoubleBits(0x406000000005F617, "1.2800000001110377e+2");
    CheckDoubleBits(0x406000016BCC41F1, "1.2800017347234783e+2");
    CheckDoubleBits(0x407000016BCC41EC, "2.5600034694469537e+2");
    CheckDoubleBits(0x4070002386F2705F, "2.5600867361738887e+2");
    CheckDoubleBits(0x407003782DACE9F1, "2.5621684043449847e+2");
    CheckDoubleBits(0x407056BC75E2D643, "2.6142101086242855e+2");
    CheckDoubleBits(0x4072000000008000, "2.8800000000186265e+2");
    CheckDoubleBits(0x4074000000008000, "3.2000000000186265e+2");
    CheckDoubleBits(0x4078000000008000, "3.8400000000186265e+2");
    CheckDoubleBits(0x407878678326EAD2, "3.9152527156068857e+2");
    CheckDoubleBits(0x4080000000000278, "5.120000000000719e+2");
    CheckDoubleBits(0x40900000001DCD8B, "1.0240000004440979e+3");
    CheckDoubleBits(0x4091B1AE4D6E2F37, "1.1324202172485655e+3");
    CheckDoubleBits(0x40A000000000000C, "2.0480000000000055e+3");
    CheckDoubleBits(0x40A0000000000C4E, "2.0480000000014325e+3");
    CheckDoubleBits(0x40A000000095032F, "2.0480000044409167e+3");
    CheckDoubleBits(0x40A000016BCC41F2, "2.0480027755575657e+3");
    CheckDoubleBits(0x40A000071AFD4995, "2.0480138777878115e+3");
    CheckDoubleBits(0x40A0002386F26FC8, "2.0480693889390423e+3");
    CheckDoubleBits(0x40A0002386F2707B, "2.0480693889391237e+3");
    CheckDoubleBits(0x40A056BC75E2D64C, "2.0913680868994325e+3");
    CheckDoubleBits(0x40A878678326EAD8, "3.1322021724855113e+3");
    CheckDoubleBits(0x40B0000000000006, "4.0960000000000055e+3");
    CheckDoubleBits(0x40B0000000000C53, "4.0960000000028695e+3");
    CheckDoubleBits(0x40B0000000013136, "4.0960000000710625e+3");
    CheckDoubleBits(0x40B000071AFD49A5, "4.0960277555756375e+3");
    CheckDoubleBits(0x40B056BC75E2D651, "4.1827361737988695e+3");
    CheckDoubleBits(0x40CAAAAAAAAAAAC9, "1.3653333333333389e+4");
    CheckDoubleBits(0x40E0000000000002, "3.2768000000000015e+4");
    CheckDoubleBits(0x40F0000002E90EF5, "6.553600071054309e+4");
    CheckDoubleBits(0x40F0000048C273A9, "6.553601776356869e+4");
    CheckDoubleBits(0x4110000000000008, "2.6214400000000047e+5");
    CheckDoubleBits(0x4110000000000C53, "2.6214400000018365e+5");
    CheckDoubleBits(0x411000000005F60F, "2.6214400002274005e+5");
    CheckDoubleBits(0x41100000001DCD6C, "2.6214400011368725e+5");
    CheckDoubleBits(0x4110000048C273A6, "2.6214407105427457e+5");
    CheckDoubleBits(0x411000016BCC41F1, "2.6214435527136835e+5");
    CheckDoubleBits(0x411000071AFD4995, "2.6214577635683987e+5");
    CheckDoubleBits(0x411056BC75E2D651, "2.6769511512312765e+5");
    CheckDoubleBits(0x4115555555555565, "3.4952533333333425e+5");
    CheckDoubleBits(0x4120000000000040, "5.242880000000075e+5");
    CheckDoubleBits(0x4121000000000040, "5.570560000000075e+5");
    CheckDoubleBits(0x4122000000000040, "5.898240000000075e+5");
    CheckDoubleBits(0x4130000000000020, "1.0485760000000075e+6");
    CheckDoubleBits(0x413003782DACE9E9, "1.0494641784197039e+6");
    CheckDoubleBits(0x4131000000000020, "1.1141120000000075e+6");
    CheckDoubleBits(0x4131B1AE4D6E2F3E, "1.1595983024625327e+6");
    CheckDoubleBits(0x414000000000000A, "2.0971520000000047e+6");
    CheckDoubleBits(0x4140000000000010, "2.0971520000000075e+6");
    CheckDoubleBits(0x4140000000000C36, "2.0971520000014557e+6");
    CheckDoubleBits(0x4140000000000C3C, "2.0971520000014585e+6");
    CheckDoubleBits(0x414000000001312E, "2.0971520000363803e+6");
    CheckDoubleBits(0x414000000095032F, "2.0971520045474987e+6");
    CheckDoubleBits(0x414000071AFD499A, "2.0971662108547213e+6");
    CheckDoubleBits(0x414056BC75E2D634, "2.1415609209850077e+6");
    CheckDoubleBits(0x414AAAAAAAAAAAC9, "3.4952533333333475e+6");
    CheckDoubleBits(0x4170000000000014, "1.6777216000000075e+7");
    CheckDoubleBits(0x4170000002E90EEB, "1.6777216181898993e+7");
    CheckDoubleBits(0x418000000000000A, "3.3554432000000075e+7");
    CheckDoubleBits(0x418000000005F627, "3.3554432002910905e+7");
    CheckDoubleBits(0x418000000E8D4A70, "3.3554433818989635e+7");
    CheckDoubleBits(0x4190000000000032, "6.710886400000075e+7");
    CheckDoubleBits(0x419056BC75E2D634, "6.852994947152025e+7");
    CheckDoubleBits(0x4198000000020000, "1.0066329600195312e+8");
    CheckDoubleBits(0x41A4000000010000, "1.6777216000195312e+8");
    CheckDoubleBits(0x41A8000000010000, "2.0132659200195312e+8");
    CheckDoubleBits(0x41B00000000000BD, "2.6843545600001127e+8");
    CheckDoubleBits(0x41B000000E8D4A58, "2.6843547055191565e+8");
    CheckDoubleBits(0x41B000016BCC41F2, "2.6843581979788125e+8");
    CheckDoubleBits(0x41B003782DACE9EA, "2.6866282967544425e+8");
    CheckDoubleBits(0x41B1B1AE4D6E2EF9, "2.9685716543040425e+8");
    CheckDoubleBits(0x41B555555555556A, "3.5791394133333457e+8");
    CheckDoubleBits(0x41C2000000008000, "6.039797760039062e+8");
    CheckDoubleBits(0x41C4000000008000, "6.710886400039062e+8");
    CheckDoubleBits(0x41C8000000008000, "8.053063680039062e+8");
    CheckDoubleBits(0x41D0000000000272, "1.0737418240001493e+9");
    CheckDoubleBits(0x41D0000000000C3A, "1.0737418240007463e+9");
    CheckDoubleBits(0x41D0000000003D22, "1.0737418240037313e+9");
    CheckDoubleBits(0x41D0000000013147, "1.0737418240186327e+9");
    CheckDoubleBits(0x41D003782DACE9F1, "1.0746513187017787e+9");
    CheckDoubleBits(0x41D1000000004000, "1.1408506880039062e+9");
    CheckDoubleBits(0x41D2000000004000, "1.2079595520039062e+9");
    CheckDoubleBits(0x41D4000000004000, "1.3421772800039062e+9");
    CheckDoubleBits(0x41D8000000004000, "1.6106127360039062e+9");
    CheckDoubleBits(0x41E0000000000020, "2.1474836480000153e+9");
    CheckDoubleBits(0x41E0000000000274, "2.1474836480002995e+9");
    CheckDoubleBits(0x41E0000000000C3C, "2.1474836480014935e+9");
    CheckDoubleBits(0x41E0000000002000, "2.1474836480039062e+9");
    CheckDoubleBits(0x41E0000000003D24, "2.1474836480074635e+9");
    CheckDoubleBits(0x41E000000E8D4A54, "2.1474837644153233e+9");
    CheckDoubleBits(0x41E000000E8D4A61, "2.1474837644153295e+9");
    CheckDoubleBits(0x41E000071AFD4996, "2.1474981999152327e+9");
    CheckDoubleBits(0x41E1000000000020, "2.2817013760000153e+9");
    CheckDoubleBits(0x41E1000000002000, "2.2817013760039062e+9");
    CheckDoubleBits(0x41E2000000002000, "2.4159191040039062e+9");
    CheckDoubleBits(0x41E4000000002000, "2.6843545600039062e+9");
    CheckDoubleBits(0x41E8000000002000, "3.2212254720039062e+9");
    CheckDoubleBits(0x4200000000001000, "8.589934592007812e+9");
    CheckDoubleBits(0x4200000000013139, "8.589934592149035e+9");
    CheckDoubleBits(0x420000000005F61D, "8.589934592745173e+9");
    CheckDoubleBits(0x4201000000001000, "9.126805504007812e+9");
    CheckDoubleBits(0x4202000000001000, "9.663676416007812e+9");
    CheckDoubleBits(0x4204000000001000, "1.0737418240007812e+10");
    CheckDoubleBits(0x4208000000001000, "1.2884901888007812e+10");
    CheckDoubleBits(0x4210000000000800, "1.7179869184007812e+10");
    CheckDoubleBits(0x4211000000000800, "1.8253611008007812e+10");
    CheckDoubleBits(0x4212000000000800, "1.9327352832007812e+10");
    CheckDoubleBits(0x4214000000000800, "2.1474836480007812e+10");
    CheckDoubleBits(0x4218000000000800, "2.5769803776007812e+10");
    CheckDoubleBits(0x42200000001DCD6C, "3.4359738382901215e+10");
    CheckDoubleBits(0x422AAAAAAAAAAAB4, "5.7266230613333405e+10");
    CheckDoubleBits(0x4230000000000400, "6.871947673601562e+10");
    CheckDoubleBits(0x423000071AFD4996, "6.871994239728745e+10");
    CheckDoubleBits(0x4231000000000400, "7.301444403201562e+10");
    CheckDoubleBits(0x4232000000000400, "7.730941132801562e+10");
    CheckDoubleBits(0x4234000000000400, "8.589934592001562e+10");
    CheckDoubleBits(0x4238000000000400, "1.0307921510401562e+11");
    CheckDoubleBits(0x4240000000000200, "1.3743895347201562e+11");
    CheckDoubleBits(0x4240000000000274, "1.3743895347201917e+11");
    CheckDoubleBits(0x4241000000000200, "1.4602888806401562e+11");
    CheckDoubleBits(0x4242000000000200, "1.5461882265601562e+11");
    CheckDoubleBits(0x4244000000000200, "1.7179869184001562e+11");
    CheckDoubleBits(0x4248000000000200, "2.0615843020801562e+11");
    CheckDoubleBits(0x4250000000000C53, "2.7487790694419257e+11");
    CheckDoubleBits(0x4250000000003D14, "2.7487790694495435e+11");
    CheckDoubleBits(0x425000000001313A, "2.7487790694876917e+11");
    CheckDoubleBits(0x425000000005F60F, "2.7487790696784467e+11");
    CheckDoubleBits(0x425000000E8D4A83, "2.7487792184516425e+11");
    CheckDoubleBits(0x4250000048C273A5, "2.7487798144980695e+11");
    CheckDoubleBits(0x4250002386F2703F, "2.7488722016975385e+11");
    CheckDoubleBits(0x4260000000000100, "5.497558138880312e+11");
    CheckDoubleBits(0x4261000000000100, "5.841155522560312e+11");
    CheckDoubleBits(0x4262000000000100, "6.184752906240312e+11");
    CheckDoubleBits(0x4264000000000100, "6.871947673600312e+11");
    CheckDoubleBits(0x4268000000000100, "8.246337208320312e+11");
    CheckDoubleBits(0x4270000000000080, "1.0995116277760312e+12");
    CheckDoubleBits(0x4270000000003D11, "1.0995116277798167e+12");
    CheckDoubleBits(0x4271000000000080, "1.1682311045120312e+12");
    CheckDoubleBits(0x4272000000000080, "1.2369505812480312e+12");
    CheckDoubleBits(0x4274000000000080, "1.3743895347200312e+12");
    CheckDoubleBits(0x4280000000000040, "2.1990232555520312e+12");
    CheckDoubleBits(0x4280000002E90EF3, "2.1990232793938687e+12");
    CheckDoubleBits(0x428000016BCC41EA, "2.1990262357842393e+12");
    CheckDoubleBits(0x4280002386F26FC8, "2.1990977613579727e+12");
    CheckDoubleBits(0x428003782DACE9EA, "2.2008859007012393e+12");
    CheckDoubleBits(0x4281000000000040, "2.3364622090240312e+12");
    CheckDoubleBits(0x4282000000000040, "2.4739011624960312e+12");
    CheckDoubleBits(0x4285555555555559, "2.9320310074026685e+12");
    CheckDoubleBits(0x42A0000000000020, "8.796093022208062e+12");
    CheckDoubleBits(0x42A1000000000020, "9.345848836096062e+12");
    CheckDoubleBits(0x42B0000000000010, "1.7592186044416062e+13");
    CheckDoubleBits(0x42B0002386F26FD0, "1.7592782090863812e+13");
    CheckDoubleBits(0x42D0000000000008, "7.036874417766412e+13");
    CheckDoubleBits(0x42D0000000003D28, "7.036874417790862e+13");
    CheckDoubleBits(0x42D000000005F5E8, "7.036874418376762e+13");
    CheckDoubleBits(0x42D000000E8D4A68, "7.036874799236162e+13");
    CheckDoubleBits(0x42D0002386F26FC8, "7.037112836345512e+13");
    CheckDoubleBits(0x42D878678326EAE8, "1.0762164716228362e+14");
    CheckDoubleBits(0x42E0000000000004, "1.4073748835532812e+14");
    CheckDoubleBits(0x42E0000000000014, "1.4073748835532862e+14");
    CheckDoubleBits(0x42E0000000000274, "1.4073748835534762e+14");
    CheckDoubleBits(0x42E0000000003D14, "1.4073748835581662e+14");
    CheckDoubleBits(0x42E0000000003D24, "1.4073748835581712e+14");
    CheckDoubleBits(0x42E000000005F5E4, "1.4073748836753512e+14");
    CheckDoubleBits(0x42E00000001DCD74, "1.4073748841636362e+14");
    CheckDoubleBits(0x42E000000E8D4A54, "1.4073749598472262e+14");
    CheckDoubleBits(0x42E000000E8D4A64, "1.4073749598472312e+14");
    CheckDoubleBits(0x42E000016BCC41F4, "1.4073767909019162e+14");
    CheckDoubleBits(0x42E0002386F26FC4, "1.4074225672691012e+14");
    CheckDoubleBits(0x42E000B1A2BC2ED4, "1.4076133021323862e+14");
    CheckDoubleBits(0x42E056BC75E2D634, "1.4371772059409762e+14");
    CheckDoubleBits(0x42E878678326EAD4, "2.1524329432456662e+14");
    CheckDoubleBits(0x42E878678326EAE4, "2.1524329432456712e+14");
    CheckDoubleBits(0x42EAAAAAAAAAAAB4, "2.3456248059221362e+14");
    CheckDoubleBits(0x4300000000000002, "5.629499534213122e+14");
    CheckDoubleBits(0x430000000000000A, "5.629499534213132e+14");
    CheckDoubleBits(0x4300000000000012, "5.629499534213142e+14");
    CheckDoubleBits(0x430000000000001A, "5.629499534213152e+14");
    CheckDoubleBits(0x4300000000000022, "5.629499534213162e+14");
    CheckDoubleBits(0x430000000000002A, "5.629499534213172e+14");
    CheckDoubleBits(0x4300000000000032, "5.629499534213182e+14");
    CheckDoubleBits(0x4300000000000272, "5.629499534213902e+14");
    CheckDoubleBits(0x4300000000000C3A, "5.629499534217032e+14");
    CheckDoubleBits(0x4300000000003D0A, "5.629499534232652e+14");
    CheckDoubleBits(0x4300000000003D12, "5.629499534232662e+14");
    CheckDoubleBits(0x4300000000003D22, "5.629499534232682e+14");
    CheckDoubleBits(0x430000000001313A, "5.629499534310792e+14");
    CheckDoubleBits(0x430000000005F5E2, "5.629499534701402e+14");
    CheckDoubleBits(0x43000000001DCD6A, "5.629499536654532e+14");
    CheckDoubleBits(0x43000000001DCD72, "5.629499536654542e+14");
    CheckDoubleBits(0x43000000009502FA, "5.629499546420152e+14");
    CheckDoubleBits(0x4300000002E90EFA, "5.629499595248312e+14");
    CheckDoubleBits(0x430000000E8D4A52, "5.629499839388902e+14");
    CheckDoubleBits(0x430000000E8D4A62, "5.629499839388922e+14");
    CheckDoubleBits(0x4300000048C2739A, "5.629501060092032e+14");
    CheckDoubleBits(0x4300000048C273AA, "5.629501060092052e+14");
    CheckDoubleBits(0x430000016BCC41EA, "5.629507163607652e+14");
    CheckDoubleBits(0x430000016BCC41F2, "5.629507163607662e+14");
    CheckDoubleBits(0x430000071AFD499A, "5.629537681185792e+14");
    CheckDoubleBits(0x4300002386F26FC2, "5.629690269076402e+14");
    CheckDoubleBits(0x430000B1A2BC2ECA, "5.630453208529532e+14");
    CheckDoubleBits(0x430000B1A2BC2ED2, "5.630453208529542e+14");
    CheckDoubleBits(0x430000B1A2BC2EE2, "5.630453208529562e+14");
    CheckDoubleBits(0x430003782DACE9DA, "5.634267905795152e+14");
    CheckDoubleBits(0x430003782DACE9EA, "5.634267905795172e+14");
    CheckDoubleBits(0x430003782DACE9F2, "5.634267905795182e+14");
    CheckDoubleBits(0x430056BC75E2D632, "5.748708823763902e+14");
    CheckDoubleBits(0x430056BC75E2D64A, "5.748708823763932e+14");
    CheckDoubleBits(0x4301B1AE4D6E2EFA, "6.225545981967032e+14");
    CheckDoubleBits(0x430555555555555A, "7.505999378950832e+14");
    CheckDoubleBits(0x430555555555556A, "7.505999378950852e+14");
    CheckDoubleBits(0x430878678326EACA, "8.609731772982652e+14");
    CheckDoubleBits(0x430878678326EAD2, "8.609731772982662e+14");
    CheckDoubleBits(0x430878678326EAE2, "8.609731772982682e+14");
    CheckDoubleBits(0x430AAAAAAAAAAAB2, "9.382499223688542e+14");
    CheckDoubleBits(0x430AAAAAAAAAAACA, "9.382499223688572e+14");
    CheckDoubleBits(0x4310000000000001, "1.1258999068426242e+15");
    CheckDoubleBits(0x4310000000000005, "1.1258999068426252e+15");
    CheckDoubleBits(0x4310000000000019, "1.1258999068426302e+15");
    CheckDoubleBits(0x431000000000007D, "1.1258999068426552e+15");
    CheckDoubleBits(0x4310000000000271, "1.1258999068427802e+15");
    CheckDoubleBits(0x4310000000000C35, "1.1258999068434052e+15");
    CheckDoubleBits(0x4310000000003D09, "1.1258999068465302e+15");
    CheckDoubleBits(0x431000000001312D, "1.1258999068621552e+15");
    CheckDoubleBits(0x431000000005F5E1, "1.1258999069402802e+15");
    CheckDoubleBits(0x43100000001DCD65, "1.1258999073309052e+15");
    CheckDoubleBits(0x43100000009502F9, "1.1258999092840302e+15");
    CheckDoubleBits(0x4310000002E90EDD, "1.1258999190496552e+15");
    CheckDoubleBits(0x431000000E8D4A51, "1.1258999678777802e+15");
    CheckDoubleBits(0x4310000048C27395, "1.1259002120184052e+15");
    CheckDoubleBits(0x431000016BCC41E9, "1.1259014327215302e+15");
    CheckDoubleBits(0x431000071AFD498D, "1.1259075362371552e+15");
    CheckDoubleBits(0x4310002386F26FC1, "1.1259380538152802e+15");
    CheckDoubleBits(0x431000B1A2BC2EC5, "1.1260906417059052e+15");
    CheckDoubleBits(0x431003782DACE9D9, "1.1268535811590302e+15");
    CheckDoubleBits(0x43101158E460913D, "1.1306682784246552e+15");
    CheckDoubleBits(0x431056BC75E2D631, "1.1497417647527802e+15");
    CheckDoubleBits(0x4311B1AE4D6E2EF5, "1.2451091963934052e+15");
    CheckDoubleBits(0x4315555555555555, "1.5011998757901652e+15");
    CheckDoubleBits(0x431878678326EAC9, "1.7219463545965302e+15");
    CheckDoubleBits(0x4350000000000001, "1.8014398509481988e+16");
    CheckDoubleBits(0x4350000000000002, "1.801439850948199e+16");
    CheckDoubleBits(0x4351000000000010, "1.914029841632467e+16");
    CheckDoubleBits(0x4352000000000000, "2.026619832316723e+16");
    CheckDoubleBits(0x4360000000000001, "3.6028797018963976e+16");
    CheckDoubleBits(0x4360000000000002, "3.602879701896398e+16");
    CheckDoubleBits(0x4361000000000010, "3.828059683264934e+16");
    CheckDoubleBits(0x4362000000000000, "4.053239664633446e+16");
    CheckDoubleBits(0x43C0000000000019, "2.305843009213707e+18");
    CheckDoubleBits(0x43F2250FA08398D4, "2.0919495182368195e+19");
    CheckDoubleBits(0x43F735F4147CF23B, "2.6760179287532483e+19");
    CheckDoubleBits(0x43FC92CE5217E4BB, "3.2942957306323907e+19");
    CheckDoubleBits(0x440137D70B9804D3, "3.9702293349085635e+19");
    CheckDoubleBits(0x4401A0D3E95D156E, "4.0647939013152195e+19");
    CheckDoubleBits(0x4403BC001F85257A, "4.5504374967041475e+19");
    CheckDoubleBits(0x4490FFFFFFFFFFFF, "2.007005755219599e+22");
    CheckDoubleBits(0x4491000000000000, "2.0070057552195992e+22");
    CheckDoubleBits(0x4532000000000020, "2.1760664753063463e+25");
    CheckDoubleBits(0x457100000000000F, "3.288278229351802e+26");
    CheckDoubleBits(0x4571000000000010, "3.2882782293518024e+26");
    CheckDoubleBits(0x45A4DF86245B6D9D, "3.2299620609902747e+27");
    CheckDoubleBits(0x46E00000001DCD65, "2.596148430393314e+33");
}

TEST_CASE("Double - Ryu 2/")
{
    // 0 <= e2 <= 79

    // vr_has_trailing_zeros = false
    // ...

    // vr_has_trailing_zeros = true
    CheckDoubleBits(0x43F168819406A592, "20070319503829443000");
    CheckDoubleBits(0x44045696499C0681, "46896767092804035000");
    CheckDoubleBits(0x44390CEE12EBB2D6, "462100304819178050000");
    CheckDoubleBits(0x444C7773810B6764, "1.0502323985253637e+21");
    CheckDoubleBits(0x4455DDB16F92D0FB, "1.6134251903100899e+21");
    CheckDoubleBits(0x446BB1A802862B8B, "4.0869085976499035e+21");
    CheckDoubleBits(0x447E8BEC73507936, "9.015758122664487e+21");
    CheckDoubleBits(0x448383091015911D, "1.1517767461729773e+22");
    CheckDoubleBits(0x449AC390930AE059, "3.1597265332826465e+22");
    CheckDoubleBits(0x44AAB1563B668DC7, "6.3026407964263685e+22");
    CheckDoubleBits(0x44B9646B8F16CD80, "1.1991158691766769e+23");
    CheckDoubleBits(0x44C9C0281DAABDDC, "2.4320765517632325e+23");
    CheckDoubleBits(0x44D6A5D863A869B9, "4.2780547154720475e+23");
    CheckDoubleBits(0x44E1D4E04D6588B7, "6.736568211094877e+23");
    CheckDoubleBits(0x44FB0EE97FF19527, "2.0444635981898903e+24");
    CheckDoubleBits(0x4505BAAE2A4148AF, "3.2836268944980665e+24");
    CheckDoubleBits(0x451637C2E527B9EC, "6.714923342200239e+24");
    CheckDoubleBits(0x45277ADE480A826F, "1.4192761465167933e+25");
    CheckDoubleBits(0x453683FEF09EA6DB, "2.7219700852251563e+25");
    CheckDoubleBits(0x454038A2F2CDF950, "3.9220543010657915e+25");
    CheckDoubleBits(0x455D5EF1BF4B9245, "1.4202884267123067e+26");
    CheckDoubleBits(0x456073F43B2E7238, "1.5912312423511635e+26");
    CheckDoubleBits(0x45742D8A284EA2CF, "3.9029714302626735e+26");
    CheckDoubleBits(0x4582EB75D9C02B17, "7.319230347578495e+26");
    CheckDoubleBits(0x459CB1EEDC2AB4FA, "2.2201720324133575e+27");
    CheckDoubleBits(0x45AF0C466F8319F0, "4.8044375184931557e+27");
    CheckDoubleBits(0x45BA8E45A3AB1D00, "8.218606584176863e+27");
    CheckDoubleBits(0x45C79618ADA27AE2, "1.4599221277246727e+28");
    CheckDoubleBits(0x45D60D0A43989895, "2.7297738889267863e+28");
    CheckDoubleBits(0x45E8957722446280, "6.0866662212506105e+28");
    CheckDoubleBits(0x45FD9878D292C966, "1.4655028124434765e+29");
    CheckDoubleBits(0x460242ADD71CE83A, "1.8084288698922457e+29");
    CheckDoubleBits(0x4610715FBE284A2F, "3.2568453807093795e+29");
    CheckDoubleBits(0x462D33409CDE2D9E, "1.1567392802273227e+30");
    CheckDoubleBits(0x463873CA143BA1F4, "1.9373109750347185e+30");
    CheckDoubleBits(0x4641ECEADE08FFEE, "2.8404023244669917e+30");
    CheckDoubleBits(0x465E4EDDF47682DB, "9.605012132974787e+30");
    CheckDoubleBits(0x4662ECA0A799CD3D, "1.1994716857424619e+31");
    CheckDoubleBits(0x467A66C833DEF707, "3.3467867623843855e+31");
    CheckDoubleBits(0x468FF7DCFAE25DF8, "8.1049055482103615e+31");
    CheckDoubleBits(0x469444F0D7B3B5CA, "1.0277756107351719e+32");
    CheckDoubleBits(0x46A52B8E7C6DF9CA, "2.1469075498117023e+32");
    CheckDoubleBits(0x46BEE24BC011AD52, "6.2640129641074525e+32");
    CheckDoubleBits(0x46C140C5EDF67C64, "6.998656437786843e+32");
    CheckDoubleBits(0x46D480385DA52DB2, "1.6632273649821755e+33");
    CheckDoubleBits(0x46E313C918A79E89, "3.0954668307992495e+33");
    CheckDoubleBits(0x46F40E79CC1BD0E7, "6.508721292586339e+33");
    CheckDoubleBits(0x47085A1F4A5CC99E, "1.5805377569533309e+34");
    CheckDoubleBits(0x4715DACC19B1A929, "2.8368992454971077e+34");
    CheckDoubleBits(0x473D5D528164527D, "1.5246940974482417e+35");
    CheckDoubleBits(0x4740565E4A304EA4, "1.6965701474007633e+35");
    CheckDoubleBits(0x475BB9F8F68F13D4, "5.7585594339080755e+35");
    CheckDoubleBits(0x476D5743CC091FF8, "1.2187723997271157e+36");
    CheckDoubleBits(0x477CA113E80F021C, "2.3784217142226719e+36");
    CheckDoubleBits(0x478D68A1451C83D0, "4.8863602118205645e+36");
    CheckDoubleBits(0x47912211984EFA34, "5.693442722405007e+36");
    CheckDoubleBits(0x47A7E91D94C307D2, "1.5891324523548635e+37");
    CheckDoubleBits(0x47B74B1BC55C64D8, "3.0962229429045025e+37");
    CheckDoubleBits(0x47D093B227873241, "8.813811588564917e+37");
    CheckDoubleBits(0x47E062F7FB5170E1, "1.7425218153072567e+38");
    CheckDoubleBits(0x47F39470B9388B1F, "4.1641725055298525e+38");
    CheckDoubleBits(0x480355FBEE7024F4, "8.224571817186897e+38");
    CheckDoubleBits(0x4817F24F8DB46DE2, "2.0371451699320473e+39");
    CheckDoubleBits(0x4821F3CBFE1F8A18, "3.0544309155624107e+39");
    CheckDoubleBits(0x4837A5228D49B193, "8.045996462237653e+39");
    CheckDoubleBits(0x484E197D4F6333D3, "2.0484704709600159e+40");

    // vm_has_trailing_zeros = false
    CheckDoubleBits(0x435F52E7FF1816E4, "35267522594888590");
    CheckDoubleBits(0x436F7126F7BADAFA, "70801091655555020");
    CheckDoubleBits(0x4373E5B803774CDA, "89609648838331800");
    CheckDoubleBits(0x438599E2E8490DDC, "194565579189894000");
    CheckDoubleBits(0x4398CCC70212714E, "446755490968130400");
    CheckDoubleBits(0x43AB01E6DF42F470, "973045180288088000");
    CheckDoubleBits(0x43B92E7D2B1C419C, "1814525323988016000");
    CheckDoubleBits(0x43CC7C614726EE02, "4105244976795616000");
    CheckDoubleBits(0x43DB1E5C4B5FFA5E, "7816403068511680000");
    CheckDoubleBits(0x43E2B685F89AB552, "10787299529904000000");
    CheckDoubleBits(0x43F700C46F19675A, "26520650309103360000");
    CheckDoubleBits(0x440BA80BBC21FD26, "63771383608430080000");
    CheckDoubleBits(0x441AD9CCC56B28D4, "123827370445491200000");
    CheckDoubleBits(0x4427493331CAB69E, "214774864015308800000");
    CheckDoubleBits(0x443D00B8160B4878, "535007393771008000000");
    CheckDoubleBits(0x444E5B9A920036C6, "1.120006141889536e+21");
    CheckDoubleBits(0x445EA75BB67E1506, "2.261847021164544e+21");
    CheckDoubleBits(0x446B274DDBA22FFC, "4.00715400976384e+21");
    CheckDoubleBits(0x44775DC8C47D6EEC, "6.8965276956672e+21");
    CheckDoubleBits(0x448411F29318A1A2, "1.18473004560384e+22");
    CheckDoubleBits(0x4499E3846338BC24, "3.05640281268224e+22");
    CheckDoubleBits(0x44A9AFD0AE88C40E, "6.065118969561088e+22");
    CheckDoubleBits(0x44BE159D94C5EB14, "1.420697310298112e+23");
    CheckDoubleBits(0x44CA240DE7FB27B6, "2.468932267737088e+23");
    CheckDoubleBits(0x44DD877EC519103C, "5.577922927525888e+23");

    // vm_has_trailing_zeros = true
    CheckDoubleBits(0x43C3D80484E57438, "2859795701043851300");
    CheckDoubleBits(0x43F366ED8E325072, "22369054869286232000");
    CheckDoubleBits(0x440CDA08F48F8906, "66527488782658224000");
    CheckDoubleBits(0x441337CA79692444, "88627074127587510000");
    CheckDoubleBits(0x44283D7AA131CEC2, "223575944093897920000");
    CheckDoubleBits(0x44361D56992C7A70, "407942415113018930000");
    CheckDoubleBits(0x444106D3248A2D6C, "628172852426524100000");
    CheckDoubleBits(0x445BEA8732CC24A2, "2.0598464878770331e+21");
    CheckDoubleBits(0x4469A20EE7F7DECA, "3.7827690222267023e+21");
    CheckDoubleBits(0x447B3950C0602122, "8.035073637880991e+21");
    CheckDoubleBits(0x4485B2C746B58628, "1.2808446993679071e+22");
    CheckDoubleBits(0x449399BB8F5B5FEA, "2.3140207534818772e+22");
    CheckDoubleBits(0x44AFB00F0E87BCAE, "7.4820536439868284e+22");
    CheckDoubleBits(0x44B6289A7E29DDF8, "1.0464106476744281e+23");
    CheckDoubleBits(0x44CA336511F89E34, "2.4745919075560152e+23");
    CheckDoubleBits(0x44DF91A2A564A448, "5.9631943497563573e+23");
    CheckDoubleBits(0x44E1FFC0DDF9E71A, "6.799843798910081e+23");
    CheckDoubleBits(0x44FD454B3D9EEFE8, "2.2116300001373651e+24");
    CheckDoubleBits(0x4501E87805F70CE6, "2.7061927495737213e+24");
    CheckDoubleBits(0x45160F5D2A386BC0, "6.667230529565941e+24");
    CheckDoubleBits(0x45285E7CED4F6044, "1.4730213308199071e+25");
    CheckDoubleBits(0x45352FBEC2980C20, "2.5612912339946862e+25");
    CheckDoubleBits(0x454C176680A7B44E, "6.7920856433575484e+25");
    CheckDoubleBits(0x4557CE93E436EE06, "1.1512331785031701e+26");
    CheckDoubleBits(0x4564E0CD7762C026, "2.0192093335712002e+26");
    CheckDoubleBits(0x4578423DEB78A6F6, "4.6923260924006663e+26");
    CheckDoubleBits(0x4583CCF88585ACA8, "7.660012079679891e+26");
    CheckDoubleBits(0x45971E66763947B6, "1.7887267156759741e+27");
    CheckDoubleBits(0x45AD89820F24363E, "4.5706511545297433e+27");
    CheckDoubleBits(0x45B85BDEA3487B2C, "7.538703862698511e+27");
    CheckDoubleBits(0x45C7D369C21BED0A, "1.4747476005981901e+28");
    CheckDoubleBits(0x45D3778906B859A8, "2.4098897789248592e+28");
    CheckDoubleBits(0x45EEFC76962F6000, "7.6718076886828454e+28");
    CheckDoubleBits(0x45FF6557D747653C, "1.5546482606914151e+29");
    CheckDoubleBits(0x460263F21D73633C, "1.8212983004432152e+29");
    CheckDoubleBits(0x46100CD68E68B526, "3.1790595074475954e+29");
    CheckDoubleBits(0x4623F247D5598B10, "7.901586506984111e+29");
    CheckDoubleBits(0x463A80ED4A68AB02, "2.0998331734331571e+30");
    CheckDoubleBits(0x46491143D9969C88, "3.9720946671706063e+30");
    CheckDoubleBits(0x465DCBD481075D64, "9.442796286009161e+30");
    CheckDoubleBits(0x4663DF78E443E416, "1.2595971152091691e+31");
    CheckDoubleBits(0x467E7F11AD345EB2, "3.8658733461594082e+31");
    CheckDoubleBits(0x46829D2AA78F9AAA, "4.7191924414955675e+31");
    CheckDoubleBits(0x469F175BC7A678B0, "1.5765133748732521e+32");
    CheckDoubleBits(0x46A209D39476B32C, "1.8293095357349822e+32");
    CheckDoubleBits(0x46B67364F5AD13FC, "4.5535549547462314e+32");
    CheckDoubleBits(0x46CECEB89A003CA4, "1.2497008420099811e+33");
    CheckDoubleBits(0x46DCC003491F6AAA, "2.3324811718398431e+33");
    CheckDoubleBits(0x46E0229E6D55F192, "2.6180907359538813e+33");
    CheckDoubleBits(0x46F6F2CB3236BEB0, "7.447185804329581e+33");
    CheckDoubleBits(0x470CF4066D7415F0, "1.8791716153180791e+34");
    CheckDoubleBits(0x47166E0777464EB6, "2.9115546863754542e+34");
    CheckDoubleBits(0x47245C43827310B2, "5.2858633756638775e+34");
    CheckDoubleBits(0x473C8F5F711978C6, "1.4829225829033581e+35");
    CheckDoubleBits(0x4742009A4B5E7B7A, "1.8694713583250392e+35");
    CheckDoubleBits(0x4756605F555A989C, "4.6474078120362184e+35");
    CheckDoubleBits(0x4767B7F45E02D368, "9.852309557634391e+35");
    CheckDoubleBits(0x4779F8C47B1C5B1A, "2.1576484138538471e+36");
    CheckDoubleBits(0x47863143A446C492, "3.6873512987625053e+36");
    CheckDoubleBits(0x47953D6898F1E838, "7.058159876999991e+36");
    CheckDoubleBits(0x47ADC1DEFB85B5E4, "1.9777123897020991e+37");
    CheckDoubleBits(0x47BFC1D964D20D30, "4.2212590433737652e+37");
    CheckDoubleBits(0x47C4AB74E630829C, "5.4949627351024705e+37");
    CheckDoubleBits(0x47D190FE216CD94C, "9.339888422491341e+37");
    CheckDoubleBits(0x47EB2D5B043697DE, "2.8899724222327332e+38");
    CheckDoubleBits(0x47F0BF89886FCA90, "3.5619465811651724e+38");
    CheckDoubleBits(0x48092D5C5CDA43F0, "1.0709192509274021e+39");
    CheckDoubleBits(0x4814C84C4E64D1F8, "1.7679722855381872e+39");
    CheckDoubleBits(0x482EBD9B105064DA, "5.2302501178662973e+39");
    CheckDoubleBits(0x483BAD7FC26B69EC, "9.418243715134821e+39");
    CheckDoubleBits(0x484BE29828E8C502, "1.8977638985547031e+40");

    // vp_has_trailing_zeros = false
    CheckDoubleBits(0x435A1EC64EE4092B, "29408745881609388");
    CheckDoubleBits(0x436C40520D1C8509, "63616163994937416");
    CheckDoubleBits(0x4378CC37BB98B27F, "111679026938718190");
    CheckDoubleBits(0x438CCDE1A8FD9757, "259445389901621980");
    CheckDoubleBits(0x439E67855805DE29, "547716589210274370");
    CheckDoubleBits(0x43A9547EBF5118F1, "912611554578167900");
    CheckDoubleBits(0x43BC4E6C20F9655B, "2039686570124335900");
    CheckDoubleBits(0x43CA04770AD01711, "3749508466435039700");
    CheckDoubleBits(0x43DE2768CF05BAC9, "8691282334372799000");
    CheckDoubleBits(0x43E414033E5F3423, "11574279573703039000");
    CheckDoubleBits(0x43F9D4E0B1BCBE51, "29781753599860478000");
    CheckDoubleBits(0x4406876037D21AE9, "51947903473753596000");
    CheckDoubleBits(0x44133335447E5AC3, "88544517044915190000");
    CheckDoubleBits(0x442B277CD7CBCD1D, "250453738162534380000");
    CheckDoubleBits(0x4434EAB8DA461DED, "385848389869465570000");
    CheckDoubleBits(0x4448B700EED09591, "911817320115199900000");
    CheckDoubleBits(0x44518ABD13F0649F, "1.2943672716963839e+21");
    CheckDoubleBits(0x446305342E46D9C1, "2.8069049036103677e+21");
    CheckDoubleBits(0x4472A09925D69F1F, "5.497819450367999e+21");
    CheckDoubleBits(0x4486B1F810368435, "1.3396876396380159e+22");
    CheckDoubleBits(0x449BEFBC48A89671, "3.2981558537584638e+22");
    CheckDoubleBits(0x44AFBC564A6934B5, "7.4933783376363516e+22");
    CheckDoubleBits(0x44B73C44E490F21B, "1.0972619800248319e+23");
    CheckDoubleBits(0x44CC7FDC90247FEB, "2.6916978252185598e+23");
    CheckDoubleBits(0x44DFABE0E2DBE903, "5.9825583584706557e+23");
    CheckDoubleBits(0x44E59DEEE49A600D, "8.166643921059839e+23");

    // vp_has_trailing_zeros = true
    CheckDoubleBits(0x43C1645898F54C71, "2506448020579083000");
    CheckDoubleBits(0x43F398143B6E45D5, "22590411658661550000");
    CheckDoubleBits(0x4409B87B98BB8E4B, "59307748562395820000");
    CheckDoubleBits(0x44141731223DB8A5, "92651509014852100000");
    CheckDoubleBits(0x4421237DE376A455, "158076049757590700000");
    CheckDoubleBits(0x443FA49B7E4089FB, "583710279144264400000");
    CheckDoubleBits(0x4440C6AEF1F64163, "618929102969543000000");
    CheckDoubleBits(0x4450206B13A503F1, "1.189935550442125e+21");
    CheckDoubleBits(0x446610F50CEE7245, "3.25640213371016e+21");
    CheckDoubleBits(0x4470D55CF301D051, "4.96835736954247e+21");
    CheckDoubleBits(0x44866FD3FFD07229, "1.324436592162047e+22");
    CheckDoubleBits(0x44929253D45265DD, "2.192546546750322e+22");
    CheckDoubleBits(0x44A516D3C5EF777D, "4.979539218804653e+22");
    CheckDoubleBits(0x44B41BCDB0F0B82F, "9.49602133584306e+22");
    CheckDoubleBits(0x44C6E4D15E09FBD3, "2.162260135574338e+23");
    CheckDoubleBits(0x44DAE7211B1DD401, "5.081804478754468e+23");
    CheckDoubleBits(0x44ECF8589B040883, "1.094459510015691e+24");
    CheckDoubleBits(0x44F15C6A2E123C4F, "1.311759707782919e+24");
    CheckDoubleBits(0x45068FACC585B885, "3.409356688942245e+24");
    CheckDoubleBits(0x4511B9AC08D6C3B5, "5.35713755241796e+24");
    CheckDoubleBits(0x4528CAEC62CCE837, "1.49862491256063e+25");
    CheckDoubleBits(0x453DAD8A90A4356B, "3.587837424355704e+25");
    CheckDoubleBits(0x4545A71B64D0AF0D, "5.235316548227021e+25");
    CheckDoubleBits(0x45570B723E0FFAE3, "1.11437389133365e+26");
    CheckDoubleBits(0x456C0AC6EADEE163, "2.712065279486174e+26");
    CheckDoubleBits(0x457ABBE3F517C4DB, "5.17109742623722e+26");
    CheckDoubleBits(0x4584636963A6A903, "7.88735192410781e+26");
    CheckDoubleBits(0x4591F3614ACF3573, "1.388868397667022e+27");
    CheckDoubleBits(0x45A482A29B7E82F5, "3.173814222354239e+27");
    CheckDoubleBits(0x45B4957C97414F8F, "6.37041850715785e+27");
    CheckDoubleBits(0x45CD8F354D120D63, "1.829638676829518e+28");
    CheckDoubleBits(0x45DC68A8B5334857, "3.516842104145073e+28");
    CheckDoubleBits(0x45E809E126AB6149, "5.951667051098373e+28");
    CheckDoubleBits(0x45F73BBC6CA38305, "1.150459465308417e+29");
    CheckDoubleBits(0x4606B829FDD7E3F7, "2.25001947727593e+29");
    CheckDoubleBits(0x46140F19A37C87E7, "3.973091301552558e+29");
    CheckDoubleBits(0x462986D86137ADAD, "1.01121832062317e+30");
    CheckDoubleBits(0x4632CC4E9BC203D9, "1.489336899019993e+30");
    CheckDoubleBits(0x46436D53934769F5, "3.078339980379608e+30");
    CheckDoubleBits(0x465D4BE8B5F76845, "9.28443767501044e+30");
    CheckDoubleBits(0x4665FBF73F773E79, "1.393416843720725e+31");
    CheckDoubleBits(0x467AA168F04C9703, "3.375817880029023e+31");
    CheckDoubleBits(0x468810A916BF42B5, "6.101222644426234e+31");
    CheckDoubleBits(0x469DF1D56BB6D04D, "1.51837479049326e+32");
    CheckDoubleBits(0x46A5883D77C91E1F, "2.18362327588074e+32");
    CheckDoubleBits(0x46BBE8F59AE7064D, "5.660820040948698e+32");
    CheckDoubleBits(0x46CB863F6F0A1A83, "1.116522529739208e+33");
    CheckDoubleBits(0x46D3B9E366807481, "1.600373478195183e+33");
    CheckDoubleBits(0x46E6DD82BDEEAF51, "3.710103182891225e+33");
    CheckDoubleBits(0x46F211639EFA7B3D, "5.86337732540143e+33");
    CheckDoubleBits(0x470ACDED77E3DE17, "1.739705330867014e+34");
    CheckDoubleBits(0x4715185C2A9FA8D9, "2.738307851051433e+34");
    CheckDoubleBits(0x472760F12DB17C1F, "6.06945235984045e+34");
    CheckDoubleBits(0x473B7EE8CCF8659B, "1.427660431594392e+35");
    CheckDoubleBits(0x474D25BE84F8AA5D, "3.026843051127245e+35");
    CheckDoubleBits(0x475915C20E3D0DCB, "5.209949069406996e+35");
    CheckDoubleBits(0x4768FB662F54B4C3, "1.037712842689022e+36");
    CheckDoubleBits(0x4777741CB1706A19, "1.94844576901796e+36");
    CheckDoubleBits(0x478709B24B242EB9, "3.827823849624933e+36");
    CheckDoubleBits(0x479F7C46DC06118B, "1.046283747013421e+37");
    CheckDoubleBits(0x47AD635F58FB4075, "1.953179157275341e+37");
    CheckDoubleBits(0x47B9B588415F8DBD, "3.417326921312679e+37");
    CheckDoubleBits(0x47C9860ECF7BCCD7, "6.78535361318899e+37");
    CheckDoubleBits(0x47DABCCFDB9CE547, "1.42161182232499e+38");
    CheckDoubleBits(0x47E4DC3B9C64001D, "2.218245942182767e+38");
    CheckDoubleBits(0x47F1A0DC52321213, "3.74913793088438e+38");
    CheckDoubleBits(0x48021DEF1926DD9F, "7.70608960692918e+38");
    CheckDoubleBits(0x4818FC043E91C9FD, "2.125441074821937e+39");
    CheckDoubleBits(0x4821EEFE5E55D20B, "3.051238628700367e+39");
    CheckDoubleBits(0x483F180E2C35EBCB, "1.058072843530204e+40");
    CheckDoubleBits(0x484F4426217A69DB, "2.127867772912241e+40");

    // -4 <= e2 <= -1

    // vr_has_trailing_zeros = false
    CheckDoubleBits(0x431E7F07D5C31F79, "2145980232288222.2");

    // vr_has_trailing_zeros = true
    // ...

    // vm_has_trailing_zeros = false
    // ...

    // vm_has_trailing_zeros = true
    // ...

    // vm_has_trailing_zeros = mm_shift == 0
    // ...

    // vp_has_trailing_zeros = false
    // ...

    // vp_has_trailing_zeros = true
    // ...

    // -81 <= e2 <= -5

    // vr_has_trailing_zeros = false
    CheckDoubleBits(0x417768282DB40000, "24543874.856445312");
    CheckDoubleBits(0x41AF8E05C4F90000, "264700642.48632812");
    CheckDoubleBits(0x41CDA769E2A48000, "995021765.2851562");
    CheckDoubleBits(0x41DE596C435E4000, "2036707597.4726562");
    CheckDoubleBits(0x41EF10D74E132000, "4169579120.5976562");
    CheckDoubleBits(0x420302F5A99E5000, "10206754099.789062");
    CheckDoubleBits(0x4218CCBD48432800, "26628542992.789062");
    CheckDoubleBits(0x4238E58B89AAA400, "106930342314.64062");
    CheckDoubleBits(0x4243FC12B4BC8A00, "171666925945.07812");
    CheckDoubleBits(0x426810EDE81F0900, "826905936120.2812");
    CheckDoubleBits(0x42791211C6F9DA80, "1722837397405.6562");
    CheckDoubleBits(0x42803ECF10B0AB40, "2232743499285.4062");
    CheckDoubleBits(0x42A2E1710E785920, "10379736857644.562");
    CheckDoubleBits(0x42B777FA696C0010, "25804069760000.062");
    CheckDoubleBits(0x42D32BB53EAB3788, "84313781218526.12");
    CheckDoubleBits(0x42E4AA685645A1F4, "181777019841807.62");
    CheckDoubleBits(0x43036CEBE9E47AD2, "683473131835226.2");

    // vr_has_trailing_zeros = true
    CheckDoubleBits(0x3E4B07B6B7526048, "1.2586885959114437e-8");
    CheckDoubleBits(0x3E5C66CEFA1CCE2A, "2.6451047247662053e-8");
    CheckDoubleBits(0x3E621CB17D84661C, "3.3736384077130207e-8");
    CheckDoubleBits(0x3E7A9B862D2E229E, "9.912072167055713e-8");
    CheckDoubleBits(0x3E86EC890D565978, "1.7079685813983345e-7");
    CheckDoubleBits(0x3E9C26F62C0AE592, "4.1950037758426137e-7");
    CheckDoubleBits(0x3EA79C6CAC2FBB0A, "7.036636237004639e-7");
    CheckDoubleBits(0x3EB29B6547EA8DA8, "0.0000011090644701860129");
    CheckDoubleBits(0x3ECA1D05F52FD80F, "0.0000031129565427568945");
    CheckDoubleBits(0x3ED5E52211FF8966, "0.0000052201869777688655");
    CheckDoubleBits(0x3EE1175553315F7A, "0.000008149693348766025");
    CheckDoubleBits(0x3EF8411727032279, "0.000023130664374789647");
    CheckDoubleBits(0x3F0A715F9B358211, "0.000050435762569692987");
    CheckDoubleBits(0x3F1658202E206B0E, "0.00008523651516194003");
    CheckDoubleBits(0x3F2B1AACFDA7A0B2, "0.00020678865151003073");
    CheckDoubleBits(0x3F3DF016742F3580, "0.00045681522550185955");
    CheckDoubleBits(0x3F4E37BCD8636423, "0.0009221717926052347");
    CheckDoubleBits(0x3F5B9E6A60FA8D71, "0.0016857184272474869");
    CheckDoubleBits(0x3F6A486540BB05E9, "0.0032083489985275005");
    CheckDoubleBits(0x3F741892E5CC29FD, "0.0049062479199712935");
    CheckDoubleBits(0x3F8E819F567EA188, "0.014895672633273419");
    CheckDoubleBits(0x3F947E070343F2D8, "0.020012006353669815");
    CheckDoubleBits(0x3FA49702EC3695C6, "0.040214625677701885");
    CheckDoubleBits(0x3FB784FE87CC2A05, "0.09187308135384605");
    CheckDoubleBits(0x3FCC2A49BBED7594, "0.22004052806998387");
    CheckDoubleBits(0x3FD1B268EDEF60F0, "0.27651427493903125");
    CheckDoubleBits(0x3FE4E14CEFE5384C, "0.6525025067765085");
    CheckDoubleBits(0x3FF86034DB896B0F, "1.5234879089027265");
    CheckDoubleBits(0x400B33B5ED7F05C8, "3.4002493433369843");
    CheckDoubleBits(0x4015C1449E1FA22A, "5.4387383181388405");
    CheckDoubleBits(0x4020AB2B83761A4E, "8.334316356818047");
    CheckDoubleBits(0x403B6098EB875178, "27.377333374535255");
    CheckDoubleBits(0x4048EB6C8F704EAF, "49.839250497663095");
    CheckDoubleBits(0x405CB8791E1B6F2A, "114.88239243201375");
    CheckDoubleBits(0x4063E6B552D0F5BD, "159.20963421642765");
    CheckDoubleBits(0x4079F46FDA644501, "415.27730788390915");
    CheckDoubleBits(0x4085044EB087B24A, "672.5384226418767");
    CheckDoubleBits(0x409C71DCD891F6B0, "1820.4656698996369");
    CheckDoubleBits(0x40ABF73C3C03556C, "3579.6176453630233");
    CheckDoubleBits(0x40B233F8EBBE9569, "4659.9723471750995");
    CheckDoubleBits(0x40C199A45A5E9773, "9011.284007858229");
    CheckDoubleBits(0x40D2939E380C3493, "19022.472170878737");
    CheckDoubleBits(0x40E413C257D86BA6, "41118.073223314525");
    CheckDoubleBits(0x40FB9075DA902050, "112903.36586010573");
    CheckDoubleBits(0x4101AD6F180F6DB0, "144813.88674817747");
    CheckDoubleBits(0x41195397DA76D485, "414949.96334392607");
    CheckDoubleBits(0x412C601BA0B45FD0, "929805.8138761465");
    CheckDoubleBits(0x413197890857EF95, "1152905.0325917949");
    CheckDoubleBits(0x414A125F914642A4, "3417279.1349566747");
    CheckDoubleBits(0x416180F3C22AB0BD, "9176990.067711229");
    CheckDoubleBits(0x41747775F8FCEF4D, "21460831.561751653");
    CheckDoubleBits(0x41876BFF62430B66, "49119212.282736585");
    CheckDoubleBits(0x419941301BF535D9, "105925638.98946323");
    CheckDoubleBits(0x41A6E3E10ECE0375, "192016519.40237013");
    CheckDoubleBits(0x41B4B38D8737E96A, "347311495.21840537");
    CheckDoubleBits(0x41C53A037E7C8DD4, "712247036.9730783");
    CheckDoubleBits(0x41D201EFE1996720, "1208467334.3969193");
    CheckDoubleBits(0x41EF69085D7F8240, "4215816939.9846497");
    CheckDoubleBits(0x42043AA48E049EAC, "10860401088.577477");
    CheckDoubleBits(0x4215CD8120E5A939, "23410526265.415257");
    CheckDoubleBits(0x422B367F17F97C84, "58439207932.743195");
    CheckDoubleBits(0x4239C99CA254CC74, "110756667988.79865");
    CheckDoubleBits(0x4243F02E1BECDCCB, "171267864537.72495");
    CheckDoubleBits(0x425DB2669AB36E5C, "510188481229.72437");
    CheckDoubleBits(0x4265A8D64E8A2DD5, "744215442513.4323");
    CheckDoubleBits(0x427A210CDBC84111, "1795578248324.0667");
    CheckDoubleBits(0x42849A25A2450F59, "2831536113825.9185");

    // vr_has_trailing_zeros = MultipleOfPow2(mr, q);
    // ...

    // vr_has_trailing_zeros = MultipleOfPow2(mr, q + 1);
    CheckDoubleBits(0x41E386606E782000, "2620588915.7539062");
    CheckDoubleBits(0x4288DDF4D9E56E40, "3417696844973.7812");
}

TEST_CASE("Double - Ryu regression - e2 >= 0 - mul_of_pow5")
{
    // e2 = 0
    // e2 = 1
    // e2 = 2
    // e2 = 3
    // e2 = 4
    // e2 = 5
    // e2 = 6
    // e2 = 7
    // e2 = 8
    // e2 = 9
    // e2 = 10
    CheckDouble(19443281670686147e3, "19443281670686147e3");
    CheckDouble(22321261069452227e3, "22321261069452227e3");
    CheckDouble(24575199045336003e3, "24575199045336003e3");
    CheckDouble(24735555225031107e3, "24735555225031107e3");
    CheckDouble(24776122185823683e3, "24776122185823683e3");
    CheckDouble(28047471983399363e3, "28047471983399363e3");
    CheckDouble(29922511680421315e3, "29922511680421315e3");
    CheckDouble(29964708819498947e3, "29964708819498947e3");
    CheckDouble(33292078465391043e3, "33292078465391043e3");
    CheckDouble(34130286528779203e3, "34130286528779203e3");
    // e2 = 11
    CheckDouble(37630861566675395e3, "37630861566675395e3");
    CheckDouble(42394783477136835e3, "42394783477136835e3");
    CheckDouble(43884354252443075e3, "43884354252443075e3");
    CheckDouble(45537172375138755e3, "45537172375138755e3");
    CheckDouble(50966770324466115e3, "50966770324466115e3");
    CheckDouble(52737176802894275e3, "52737176802894275e3");
    CheckDouble(63054434550130115e3, "63054434550130115e3");
    CheckDouble(71205164718639555e3, "71205164718639555e3");
    CheckDouble(72164368469336515e3, "72164368469336515e3");
    CheckDouble(73702488673445315e3, "73702488673445315e3");
    // e2 = 12
    // e2 = 13
    // e2 = 14
    CheckDouble(29545414111397315e4, "29545414111397315e4");
    CheckDouble(33643205435749827e4, "33643205435749827e4");
    CheckDouble(38051273361274307e4, "38051273361274307e4");
    CheckDouble(39162053149193667e4, "39162053149193667e4");
    CheckDouble(43233699074479555e4, "43233699074479555e4");
    CheckDouble(43296844175545795e4, "43296844175545795e4");
    CheckDouble(46920882870650307e4, "46920882870650307e4");
    CheckDouble(55003400810075587e4, "55003400810075587e4");
    CheckDouble(55635504934593987e4, "55635504934593987e4");
    CheckDouble(58873941419193795e4, "58873941419193795e4");
    // e2 = 15
    // e2 = 16
    // e2 = 17
    CheckDouble(26173040727029187e5, "26173040727029187e5");
    CheckDouble(30526409591469507e5, "30526409591469507e5");
    CheckDouble(33456744693085635e5, "33456744693085635e5");
    CheckDouble(33636787655128515e5, "33636787655128515e5");
    CheckDouble(34601474124051907e5, "34601474124051907e5");
    CheckDouble(36697409008301507e5, "36697409008301507e5");
    CheckDouble(37226686896027075e5, "37226686896027075e5");
    CheckDouble(43695137082717635e5, "43695137082717635e5");
    CheckDouble(46565388728874435e5, "46565388728874435e5");
    CheckDouble(46838971190457795e5, "46838971190457795e5");
    // e2 = 18
    // e2 = 19
    // e2 = 20
    CheckDouble(20139970741597635e6, "20139970741597635e6");
    CheckDouble(22239298935256515e6, "22239298935256515e6");
    CheckDouble(22514250475697603e6, "22514250475697603e6");
    CheckDouble(24989729569109443e6, "24989729569109443e6");
    CheckDouble(26298545665013187e6, "26298545665013187e6");
    CheckDouble(30255977718937027e6, "30255977718937027e6");
    CheckDouble(31851710521210307e6, "31851710521210307e6");
    CheckDouble(34251176752838083e6, "34251176752838083e6");
    CheckDouble(35210964908307907e6, "35210964908307907e6");
    CheckDouble(37740834282796483e6, "37740834282796483e6");
    // e2 = 21
    CheckDouble(46820507530753475e6, "46820507530753475e6");
    CheckDouble(46971337615013315e6, "46971337615013315e6");
    CheckDouble(59594741781689795e6, "59594741781689795e6");
    CheckDouble(61154817608971715e6, "61154817608971715e6");
    CheckDouble(61567164876125635e6, "61567164876125635e6");
    CheckDouble(64376781695284675e6, "64376781695284675e6");
    CheckDouble(66546669552334275e6, "66546669552334275e6");
    CheckDouble(68319547938633155e6, "68319547938633155e6");
    CheckDouble(69838500466587075e6, "69838500466587075e6");
    CheckDouble(73312490263737795e6, "73312490263737795e6");
    // e2 = 22
    // e2 = 23
    // e2 = 24
    CheckDouble(31028531437303235e7, "31028531437303235e7");
    CheckDouble(34474460663838147e7, "34474460663838147e7");
    CheckDouble(35505731076158915e7, "35505731076158915e7");
    CheckDouble(46257708875183555e7, "46257708875183555e7");
    CheckDouble(55134592048887235e7, "55134592048887235e7");
    CheckDouble(56941280658781635e7, "56941280658781635e7");
    CheckDouble(58910370331424195e7, "58910370331424195e7");
    CheckDouble(59550121485465027e7, "59550121485465027e7");
    CheckDouble(59693791124714947e7, "59693791124714947e7");
    CheckDouble(60359390875678147e7, "60359390875678147e7");
    // e2 = 25
    // e2 = 26
    // e2 = 27
    CheckDouble(24563026695353795e8, "24563026695353795e8");
    CheckDouble(24857789097113027e8, "24857789097113027e8");
    CheckDouble(27285978931066307e8, "27285978931066307e8");
    CheckDouble(30676781097350595e8, "30676781097350595e8");
    CheckDouble(32082188212696515e8, "32082188212696515e8");
    CheckDouble(33551879569470915e8, "33551879569470915e8");
    CheckDouble(37224523114280387e8, "37224523114280387e8");
    CheckDouble(43292160692909507e8, "43292160692909507e8");
    CheckDouble(47121813358572995e8, "47121813358572995e8");
    CheckDouble(47225634925639107e8, "47225634925639107e8");
    // e2 = 28
    // e2 = 29
    // e2 = 30
    CheckDouble(20818939200730563e9, "20818939200730563e9");
    CheckDouble(25406835547239875e9, "25406835547239875e9");
    CheckDouble(28501953657501123e9, "28501953657501123e9");
    CheckDouble(30078426226947523e9, "30078426226947523e9");
    CheckDouble(31948014837364163e9, "31948014837364163e9");
    CheckDouble(32101237082617283e9, "32101237082617283e9");
    CheckDouble(33335966905726403e9, "33335966905726403e9");
    CheckDouble(34644823116412355e9, "34644823116412355e9");
    CheckDouble(34778461317428675e9, "34778461317428675e9");
    CheckDouble(37720662183835075e9, "37720662183835075e9");
    // e2 = 31
    CheckDouble(42146555182511555e9, "42146555182511555e9");
    CheckDouble(42863542930044355e9, "42863542930044355e9");
    CheckDouble(43045162081383875e9, "43045162081383875e9");
    CheckDouble(44703324551312835e9, "44703324551312835e9");
    CheckDouble(46965518387180995e9, "46965518387180995e9");
    CheckDouble(49190942269830595e9, "49190942269830595e9");
    CheckDouble(51368577594881475e9, "51368577594881475e9");
    CheckDouble(64823491771823555e9, "64823491771823555e9");
    CheckDouble(68077116145071555e9, "68077116145071555e9");
    CheckDouble(73801260652361155e9, "73801260652361155e9");
    // e2 = 32
    // e2 = 33
    // e2 = 34
    CheckDouble(35748296465642947e10, "35748296465642947e10");
    CheckDouble(35806407708702147e10, "35806407708702147e10");
    CheckDouble(37539095293523395e10, "37539095293523395e10");
    CheckDouble(43303054013429187e10, "43303054013429187e10");
    CheckDouble(44819973178717635e10, "44819973178717635e10");
    CheckDouble(49727923240367555e10, "49727923240367555e10");
    CheckDouble(54133111476778435e10, "54133111476778435e10");
    CheckDouble(55314288216700355e10, "55314288216700355e10");
    CheckDouble(57946848021247427e10, "57946848021247427e10");
    CheckDouble(60187985779946947e10, "60187985779946947e10");
    // e2 = 35
    // e2 = 36
    // e2 = 37
    CheckDouble(28100591463495107e11, "28100591463495107e11");
    CheckDouble(31923125778707907e11, "31923125778707907e11");
    CheckDouble(35404649891165635e11, "35404649891165635e11");
    CheckDouble(36115239345386947e11, "36115239345386947e11");
    CheckDouble(41514742038721987e11, "41514742038721987e11");
    CheckDouble(42517794875045315e11, "42517794875045315e11");
    CheckDouble(43247212392150467e11, "43247212392150467e11");
    CheckDouble(45535696595252675e11, "45535696595252675e11");
    CheckDouble(45729086892996035e11, "45729086892996035e11");
    CheckDouble(46203326981273027e11, "46203326981273027e11");
    // e2 = 38
    // e2 = 39
    // e2 = 40
    CheckDouble(20310141679760835e12, "20310141679760835e12");
    CheckDouble(22056584903980483e12, "22056584903980483e12");
    CheckDouble(22317263649043907e12, "22317263649043907e12");
    CheckDouble(22779579297494467e12, "22779579297494467e12");
    CheckDouble(27016475414427075e12, "27016475414427075e12");
    CheckDouble(27777561872889283e12, "27777561872889283e12");
    CheckDouble(28490454503323075e12, "28490454503323075e12");
    CheckDouble(30290053316539843e12, "30290053316539843e12");
    CheckDouble(33417633469101507e12, "33417633469101507e12");
    CheckDouble(35735167757252035e12, "35735167757252035e12");
    // e2 = 41
    CheckDouble(41833522027623875e12, "41833522027623875e12");
    CheckDouble(42791138673358275e12, "42791138673358275e12");
    CheckDouble(50747909186581955e12, "50747909186581955e12");
    CheckDouble(51663703588271555e12, "51663703588271555e12");
    CheckDouble(51684716715767235e12, "51684716715767235e12");
    CheckDouble(52023396361893315e12, "52023396361893315e12");
    CheckDouble(55735869455791555e12, "55735869455791555e12");
    CheckDouble(59332045572732355e12, "59332045572732355e12");
    CheckDouble(77388947078575555e12, "77388947078575555e12");
    CheckDouble(78668600372164035e12, "78668600372164035e12");
    // e2 = 42
    // e2 = 43
    // e2 = 44
    CheckDouble(36844348430153155e13, "36844348430153155e13");
    CheckDouble(39829651348583875e13, "39829651348583875e13");
    CheckDouble(45639926056285635e13, "45639926056285635e13");
    CheckDouble(46157469615453635e13, "46157469615453635e13");
    CheckDouble(46773582674064835e13, "46773582674064835e13");
    CheckDouble(46943620429313475e13, "46943620429313475e13");
    CheckDouble(47949716518401475e13, "47949716518401475e13");
    CheckDouble(56695944970302915e13, "56695944970302915e13");
    CheckDouble(61355683838752195e13, "61355683838752195e13");
    CheckDouble(62226746156053955e13, "62226746156053955e13");
    // e2 = 45
    // e2 = 46
    // e2 = 47
    CheckDouble(30566832820254147e14, "30566832820254147e14");
    CheckDouble(31841957070828995e14, "31841957070828995e14");
    CheckDouble(31931429829539267e14, "31931429829539267e14");
    CheckDouble(33036610814145987e14, "33036610814145987e14");
    CheckDouble(35886922910463427e14, "35886922910463427e14");
    CheckDouble(37494993025824195e14, "37494993025824195e14");
    CheckDouble(39268814519072195e14, "39268814519072195e14");
    CheckDouble(41019408829183427e14, "41019408829183427e14");
    CheckDouble(49524474867414467e14, "49524474867414467e14");
    CheckDouble(50318734579529155e14, "50318734579529155e14");
    // e2 = 48
    // e2 = 49
    // e2 = 50
    CheckDouble(20508431729882563e15, "20508431729882563e15");
    CheckDouble(21888456261694915e15, "21888456261694915e15");
    CheckDouble(25341060211865027e15, "25341060211865027e15");
    CheckDouble(29876408237487555e15, "29876408237487555e15");
    CheckDouble(32400337179047363e15, "32400337179047363e15");
    CheckDouble(33102100475475395e15, "33102100475475395e15");
    CheckDouble(33508782338799043e15, "33508782338799043e15");
    CheckDouble(35060743001404867e15, "35060743001404867e15");
    CheckDouble(39157111009637827e15, "39157111009637827e15");
    CheckDouble(39529570573546947e15, "39529570573546947e15");
    // e2 = 51
    CheckDouble(54212036534007235e15, "54212036534007235e15");
    CheckDouble(55934146621011395e15, "55934146621011395e15");
    CheckDouble(56780770574398915e15, "56780770574398915e15");
    CheckDouble(62182121445848515e15, "62182121445848515e15");
    CheckDouble(68086498887005635e15, "68086498887005635e15");
    CheckDouble(69635435892635075e15, "69635435892635075e15");
    CheckDouble(72976576851539395e15, "72976576851539395e15");
    CheckDouble(73283065717781955e15, "73283065717781955e15");
    CheckDouble(73401263217767875e15, "73401263217767875e15");
    CheckDouble(78191010746267075e15, "78191010746267075e15");
    // e2 = 52
    // e2 = 53
    // e2 = 54
    CheckDouble(32881510955152835e16, "32881510955152835e16");
    CheckDouble(35971138629203395e16, "35971138629203395e16");
    CheckDouble(39401614907864515e16, "39401614907864515e16");
    CheckDouble(42139398861026755e16, "42139398861026755e16");
    CheckDouble(44085534442190275e16, "44085534442190275e16");
    CheckDouble(49126795255543235e16, "49126795255543235e16");
    CheckDouble(52210925371454915e16, "52210925371454915e16");
    CheckDouble(54245021882840515e16, "54245021882840515e16");
    CheckDouble(59088370603193795e16, "59088370603193795e16");
    CheckDouble(63810773044491715e16, "63810773044491715e16");
    // e2 = 55
    // e2 = 56
    // e2 = 57
    CheckDouble(26399889909413315e17, "26399889909413315e17");
    CheckDouble(27402644513945027e17, "27402644513945027e17");
    CheckDouble(27917215955744195e17, "27917215955744195e17");
    CheckDouble(33480744792290755e17, "33480744792290755e17");
    CheckDouble(33735831489934787e17, "33735831489934787e17");
    CheckDouble(33999714280601027e17, "33999714280601027e17");
    CheckDouble(36383455489619395e17, "36383455489619395e17");
    CheckDouble(39035477535815107e17, "39035477535815107e17");
    CheckDouble(40838676605367747e17, "40838676605367747e17");
    CheckDouble(42716642465609155e17, "42716642465609155e17");
    // e2 = 58
    // e2 = 59
    // e2 = 60
    CheckDouble(21254175491421635e18, "21254175491421635e18");
    CheckDouble(22098600421553603e18, "22098600421553603e18");
    CheckDouble(22626366002886083e18, "22626366002886083e18");
    CheckDouble(23752265909728707e18, "23752265909728707e18");
    CheckDouble(24156886188750275e18, "24156886188750275e18");
    CheckDouble(25001311118882243e18, "25001311118882243e18");
    CheckDouble(25036495490971075e18, "25036495490971075e18");
    CheckDouble(26567015676835267e18, "26567015676835267e18");
    CheckDouble(35503846187398595e18, "35503846187398595e18");
    CheckDouble(36418639861708227e18, "36418639861708227e18");
    // e2 = 61
    // e2 = 62
    // e2 = 63
    // e2 = 64
    CheckDouble(33568705722512835e19, "33568705722512835e19");
    CheckDouble(34272393164289475e19, "34272393164289475e19");
    CheckDouble(39198205256725955e19, "39198205256725955e19");
    CheckDouble(39901892698502595e19, "39901892698502595e19");
    CheckDouble(42012955023832515e19, "42012955023832515e19");
    CheckDouble(46938767116268995e19, "46938767116268995e19");
    CheckDouble(50457204325152195e19, "50457204325152195e19");
    CheckDouble(55383016417588675e19, "55383016417588675e19");
    CheckDouble(63827265718908355e19, "63827265718908355e19");
    // e2 = 65
    // e2 = 66
    // e2 = 67
    CheckDouble(28642893630076355e20, "28642893630076355e20");
    CheckDouble(29205843583497667e20, "29205843583497667e20");
    CheckDouble(32020593350604227e20, "32020593350604227e20");
    CheckDouble(34272393164289475e20, "34272393164289475e20");
    CheckDouble(37650092884817347e20, "37650092884817347e20");
    CheckDouble(43279592419030467e20, "43279592419030467e20");
    CheckDouble(48346141999822275e20, "48346141999822275e20");
    CheckDouble(51160891766928835e20, "51160891766928835e20");
    // e2 = 68
    // e2 = 69
    // e2 = 70
    CheckDouble(26391093816391107e21, "26391093816391107e21");
    CheckDouble(28642893630076355e21, "28642893630076355e21");
    CheckDouble(30894693443761603e21, "30894693443761603e21");
    CheckDouble(37650092884817347e21, "37650092884817347e21");
    // e2 = 71
    // e2 = 72
    // e2 = 73
    // e2 = 74
    CheckDouble(51160891766928835e22, "51160891766928835e22");
    // e2 = 75
    // e2 = 76
    // e2 = 77
    // e2 = 78
    // e2 = 79
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

#if 0
static inline float FloatFromBits(uint32_t bits)
{
    float f;
    std::memcpy(&f, &bits, sizeof(uint32_t));
    return f;
}

static inline uint32_t BitsFromFloat(float f)
{
    uint32_t u;
    std::memcpy(&u, &f, sizeof(uint32_t));
    return u;
}

TEST_CASE("Ftoa - all of them")
{
    constexpr int P = 24;
    constexpr uint32_t MaxF = (1u << (P - 1)) - 1;
    constexpr int ExpBias     = std::numeric_limits<float>::max_exponent - 1 + (P - 1);
    constexpr int MaxExponent = std::numeric_limits<float>::max_exponent - 1 - (P - 1);
    constexpr int MinExponent = std::numeric_limits<float>::min_exponent - 1 - (P - 1);
    constexpr int MinExp = 0; // 2 + ExpBias; // 0;
    constexpr int MaxExp = 255 - 1;
    //constexpr int MinExp = -40 + 2 + ExpBias;
    //constexpr int MaxExp =  -1 + 2 + ExpBias;

    for (int e = MinExp; e <= MaxExp; ++e)
    {
        printf("e = %3d ...\n", e);

        for (uint32_t f = 0; f <= MaxF; ++f)
        {
            const uint32_t bits = ((uint32_t)e << (P-1)) | f;
            if (bits == 0)
                continue;

            const auto value = FloatFromBits(bits);

            char buf[256];
            char* end = RyuFtoa(buf, value);
            end[0] = '\0';

            double_conversion::StringToDoubleConverter s2f(0, 0.0, 0.0, "inf", "nan");
            int unused;
            const float value_out = s2f.StringToFloat(buf, (int)(end - buf), &unused);

            const uint32_t bits_out = BitsFromFloat(value_out);
            if (bits != bits_out)
            {
                CHECK(false);
            }
        }
    }
}
#endif
