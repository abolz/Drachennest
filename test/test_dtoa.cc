#include "double-conversion.h"
#include "grisu2.h"
#include "grisu3.h"
#include "ryu.h"
// #include "milo.h"
// #include "swift.h"
// #include "floaxie.h"

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

struct D2S_DoubleConversion
{
    static const bool optimal = true;
    static const char* Name() { return "double-conversion"; }
    char* operator()(char* buf, int buflen, float f) { return double_conversion_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) { return double_conversion_Dtoa(buf, buflen, f); }
};

struct D2S_Grisu2
{
    static const bool optimal = false;
    static const char* Name() { return "grisu2"; }
    char* operator()(char* buf, int buflen, float f) { return grisu2_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) { return grisu2_Dtoa(buf, buflen, f); }
};

struct D2S_Grisu3
{
    static const bool optimal = true;
    static const char* Name() { return "grisu3-dragon4"; }
    char* operator()(char* buf, int buflen, float f) { return grisu3_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) { return grisu3_Dtoa(buf, buflen, f); }
};

struct D2S_Ryu
{
    static const bool optimal = true;
    static const char* Name() { return "ryu"; }
    char* operator()(char* buf, int buflen, float f) { return ryu_Ftoa(buf, buflen, f); }
    char* operator()(char* buf, int buflen, double f) { return ryu_Dtoa(buf, buflen, f); }
};

// struct D2S_Milo
// {
//     static const bool optimal = false;
//     static const char* Name() { return "milo"; }
//     char* operator()(char* buf, int buflen, float f) = delete;
//     char* operator()(char* buf, int buflen, double f) { return milo_Dtoa(buf, buflen, f); }
// };

// struct D2S_Swift
// {
//     static const bool optimal = true;
//     static const char* Name() { return "swift"; }
//     char* operator()(char* buf, int buflen, float f) { return swift_Ftoa(buf, buflen, f); }
//     char* operator()(char* buf, int buflen, double f) { return swift_Dtoa(buf, buflen, f); }
// };

// struct D2S_Floaxie
// {
//     static const bool optimal = false;
//     static const char* Name() { return "floaxie"; }
//     char* operator()(char* buf, int buflen, float f) { return floaxie_Ftoa(buf, buflen, f); }
//     char* operator()(char* buf, int buflen, double f) { return floaxie_Dtoa(buf, buflen, f); }
// };

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
    CAPTURE(Converter::Name());
    CAPTURE(f0);

    // Dtoa
    char buf0[BufSize];
    char* end0 = d2s(buf0, BufSize, f0);
    *end0 = '\0';
    const int length0 = static_cast<int>(end0 - buf0);

    CAPTURE(buf0);

    // Strtod
    const float f1 = double_conversion_Strtof(buf0, length0);

    const uint32_t bits0 = ReinterpretBits<uint32_t>(f0);
    const uint32_t bits1 = ReinterpretBits<uint32_t>(f1);
    CAPTURE(bits0);
    CAPTURE(bits1);
    CHECK(bits0 == bits1);

    char buf1[BufSize];
    char* end1 = double_conversion_Ftoa(buf1, BufSize, f0);
    *end1 = '\0';

    if (Converter::optimal)
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
    // CheckSingle(D2S_DoubleConversion{}, f);
    CheckSingle(D2S_Grisu2{}, f);
    CheckSingle(D2S_Grisu3{}, f);
    CheckSingle(D2S_Ryu{}, f);
    // CheckSingle(D2S_Swift{}, f);
    // CheckSingle(D2S_Floaxie{}, f);
}

template <typename Converter>
static void CheckDouble(Converter d2s, double f0)
{
    CAPTURE(Converter::Name());
    CAPTURE(f0);

    // Dtoa
    char buf0[BufSize];
    char* end0 = d2s(buf0, BufSize, f0);
    *end0 = '\0';
    const int length0 = static_cast<int>(end0 - buf0);

    CAPTURE(buf0);

    // Strtod
    const double f1 = double_conversion_Strtod(buf0, length0);

    const uint64_t bits0 = ReinterpretBits<uint64_t>(f0);
    const uint64_t bits1 = ReinterpretBits<uint64_t>(f1);
    CAPTURE(bits0);
    CAPTURE(bits1);
    CHECK(bits0 == bits1);

    char buf1[BufSize];
    char* end1 = double_conversion_Dtoa(buf1, BufSize, f0);
    *end1 = '\0';

    if (Converter::optimal)
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
    // CheckDouble(D2S_DoubleConversion{}, f);
    CheckDouble(D2S_Grisu2{}, f);
    CheckDouble(D2S_Grisu3{}, f);
    CheckDouble(D2S_Ryu{}, f);
    // CheckDouble(D2S_Milo{}, f);
    // CheckDouble(D2S_Swift{}, f);
    // CheckDouble(D2S_Floaxie{}, f);
}

template <typename Converter>
static void CheckDoubleString(Converter d2s, double value, const std::string& expected)
{
    // static_assert(Converter::optimal, "xxx");
    if (!Converter::optimal)
        return;

    char buf[32];
    char* end = d2s(buf, 32, value);

    const auto num_actual = ScanNumber(buf, end);
    const auto num_expected = ScanNumber(expected);

    CAPTURE(Converter::Name());
    CAPTURE(value);
    CHECK(num_actual.digits == num_expected.digits);
    CHECK(num_actual.exponent == num_expected.exponent);
}

static void CheckDoubleString(double value, const std::string& expected)
{
    CheckDoubleString(D2S_DoubleConversion{}, value, expected);
    CheckDoubleString(D2S_Grisu3{}, value, expected);
    CheckDoubleString(D2S_Ryu{}, value, expected);
    // CheckDoubleString(D2S_Swift{}, value, expected);
}

//==================================================================================================
//
//==================================================================================================

TEST_CASE("Single")
{
    CheckSingle(MakeSingle(0,   0, 0x00000000)); // +0
    CheckSingle(MakeSingle(0,   0, 0x00000001)); // min denormal
    CheckSingle(MakeSingle(0,   0, 0x007FFFFF)); // max denormal
    CheckSingle(MakeSingle(0,   1, 0x00000000)); // min normal
    CheckSingle(MakeSingle(0,   1, 0x00000001));
    CheckSingle(MakeSingle(0,   1, 0x007FFFFF));
    CheckSingle(MakeSingle(0,   2, 0x00000000));
    CheckSingle(MakeSingle(0,   2, 0x00000001));
    CheckSingle(MakeSingle(0,  24, 0x00000000)); // fail if no special case in normalized boundaries
    CheckSingle(MakeSingle(0,  30, 0x00000000)); // fail if no special case in normalized boundaries
    CheckSingle(MakeSingle(0,  31, 0x00000000)); // fail if no special case in normalized boundaries
    CheckSingle(MakeSingle(0,  57, 0x00000000)); // fail if no special case in normalized boundaries
    CheckSingle(MakeSingle(0, 254, 0x007FFFFE));
    CheckSingle(MakeSingle(0, 254, 0x007FFFFF)); // max normal
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
    }
}

TEST_CASE("Single - Paxson, Kahan")
{
    // V. Paxson and W. Kahan, "A Program for Testing IEEE Binary-Decimal Conversion", manuscript, May 1991,
    // ftp://ftp.ee.lbl.gov/testbase-report.ps.Z    (report)
    // ftp://ftp.ee.lbl.gov/testbase.tar.Z          (program)

    // Table 16: Stress Inputs for Converting 24-bit Binary to Decimal, < 1/2 ULP
    CheckSingle(MakeSingle(12676506, -102)); // digits  1, bits 32
    CheckSingle(MakeSingle(12676506, -103)); // digits  2, bits 29
    CheckSingle(MakeSingle(15445013,  +86)); // digits  3, bits 34
    CheckSingle(MakeSingle(13734123, -138)); // digits  4, bits 32
    CheckSingle(MakeSingle(12428269, -130)); // digits  5, bits 30
    CheckSingle(MakeSingle(15334037, -146)); // digits  6, bits 31
    CheckSingle(MakeSingle(11518287,  -41)); // digits  7, bits 30
    CheckSingle(MakeSingle(12584953, -145)); // digits  8, bits 31
    CheckSingle(MakeSingle(15961084, -125)); // digits  9, bits 32
    CheckSingle(MakeSingle(14915817, -146)); // digits 10, bits 31
    CheckSingle(MakeSingle(10845484, -102)); // digits 11, bits 30
    CheckSingle(MakeSingle(16431059,  -61)); // digits 12, bits 29

    // Table 17: Stress Inputs for Converting 24-bit Binary to Decimal, > 1/2 ULP
    CheckSingle(MakeSingle(16093626,  +69)); // digits  1, bits 30
    CheckSingle(MakeSingle( 9983778,  +25)); // digits  2, bits 31
    CheckSingle(MakeSingle(12745034, +104)); // digits  3, bits 31
    CheckSingle(MakeSingle(12706553,  +72)); // digits  4, bits 31
    CheckSingle(MakeSingle(11005028,  +45)); // digits  5, bits 30
    CheckSingle(MakeSingle(15059547,  +71)); // digits  6, bits 31
    CheckSingle(MakeSingle(16015691,  -99)); // digits  7, bits 29
    CheckSingle(MakeSingle( 8667859,  +56)); // digits  8, bits 33
    CheckSingle(MakeSingle(14855922,  -82)); // digits  9, bits 35
    CheckSingle(MakeSingle(14855922,  -83)); // digits 10, bits 33
    CheckSingle(MakeSingle(10144164, -110)); // digits 11, bits 32
    CheckSingle(MakeSingle(13248074,  +95)); // digits 12, bits 33
}

TEST_CASE("Single - Regression")
{
    CheckSingle(7.0385307e-26f);
}

TEST_CASE("Double")
{
    CheckDouble(MakeDouble(0,    0, 0x0000000000000000)); // +0
    CheckDouble(MakeDouble(0,    0, 0x0000000000000001)); // min denormal
    CheckDouble(MakeDouble(0,    0, 0x000FFFFFFFFFFFFF)); // max denormal
    CheckDouble(MakeDouble(0,    1, 0x0000000000000000)); // min normal
    CheckDouble(MakeDouble(0,    1, 0x0000000000000001));
    CheckDouble(MakeDouble(0,    1, 0x000FFFFFFFFFFFFF));
    CheckDouble(MakeDouble(0,    2, 0x0000000000000000));
    CheckDouble(MakeDouble(0,    2, 0x0000000000000001));
    CheckDouble(MakeDouble(0,    4, 0x0000000000000000)); // fail if no special case in normalized boundaries
    CheckDouble(MakeDouble(0,    5, 0x0000000000000000)); // fail if no special case in normalized boundaries
    CheckDouble(MakeDouble(0,    6, 0x0000000000000000)); // fail if no special case in normalized boundaries
    CheckDouble(MakeDouble(0,   10, 0x0000000000000000)); // fail if no special case in normalized boundaries
    CheckDouble(MakeDouble(0, 2046, 0x000FFFFFFFFFFFFE));
    CheckDouble(MakeDouble(0, 2046, 0x000FFFFFFFFFFFFF)); // max normal
}

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
    }
}

TEST_CASE("Double - Paxson, Kahan")
{
    // V. Paxson and W. Kahan, "A Program for Testing IEEE Binary-Decimal Conversion", manuscript, May 1991,
    // ftp://ftp.ee.lbl.gov/testbase-report.ps.Z    (report)
    // ftp://ftp.ee.lbl.gov/testbase.tar.Z          (program)

    // Table 3: Stress Inputs for Converting 53-bit Binary to Decimal, < 1/2 ULP
    CheckDouble(MakeDouble(8511030020275656,  -342)); // digits  1, bits 63
    CheckDouble(MakeDouble(5201988407066741,  -824)); // digits  2, bits 63
    CheckDouble(MakeDouble(6406892948269899,  +237)); // digits  3, bits 62 (D3. [Calculate q'.] One correction step)
    CheckDouble(MakeDouble(8431154198732492,   +72)); // digits  4, bits 61 (D3. [Calculate q'.] One correction step)
    CheckDouble(MakeDouble(6475049196144587,   +99)); // digits  5, bits 64 (D3. [Calculate q'.] One correction step)
    CheckDouble(MakeDouble(8274307542972842,  +726)); // digits  6, bits 64
    CheckDouble(MakeDouble(5381065484265332,  -456)); // digits  7, bits 64
    CheckDouble(MakeDouble(6761728585499734, -1057)); // digits  8, bits 64
    CheckDouble(MakeDouble(7976538478610756,  +376)); // digits  9, bits 67 (D6. [Add back.])
    CheckDouble(MakeDouble(5982403858958067,  +377)); // digits 10, bits 63
    CheckDouble(MakeDouble(5536995190630837,   +93)); // digits 11, bits 63
    CheckDouble(MakeDouble(7225450889282194,  +710)); // digits 12, bits 66 (D6. [Add back.])
    CheckDouble(MakeDouble(7225450889282194,  +709)); // digits 13, bits 64
    CheckDouble(MakeDouble(8703372741147379,  +117)); // digits 14, bits 66
    CheckDouble(MakeDouble(8944262675275217, -1001)); // digits 15, bits 63
    CheckDouble(MakeDouble(7459803696087692,  -707)); // digits 16, bits 63
    CheckDouble(MakeDouble(6080469016670379,  -381)); // digits 17, bits 62
    CheckDouble(MakeDouble(8385515147034757,  +721)); // digits 18, bits 64
    CheckDouble(MakeDouble(7514216811389786,  -828)); // digits 19, bits 64
    CheckDouble(MakeDouble(8397297803260511,  -345)); // digits 20, bits 64
    CheckDouble(MakeDouble(6733459239310543,  +202)); // digits 21, bits 63
    CheckDouble(MakeDouble(8091450587292794,  -473)); // digits 22, bits 63

    // Table 4: Stress Inputs for Converting 53-bit Binary to Decimal, > 1/2 ULP
    CheckDouble(MakeDouble(6567258882077402, +952)); // digits  1, bits 62
    CheckDouble(MakeDouble(6712731423444934, +535)); // digits  2, bits 65
    CheckDouble(MakeDouble(6712731423444934, +534)); // digits  3, bits 63
    CheckDouble(MakeDouble(5298405411573037, -957)); // digits  4, bits 62
    CheckDouble(MakeDouble(5137311167659507, -144)); // digits  5, bits 61
    CheckDouble(MakeDouble(6722280709661868, +363)); // digits  6, bits 64
    CheckDouble(MakeDouble(5344436398034927, -169)); // digits  7, bits 61
    CheckDouble(MakeDouble(8369123604277281, -853)); // digits  8, bits 65
    CheckDouble(MakeDouble(8995822108487663, -780)); // digits  9, bits 63
    CheckDouble(MakeDouble(8942832835564782, -383)); // digits 10, bits 66
    CheckDouble(MakeDouble(8942832835564782, -384)); // digits 11, bits 64
    CheckDouble(MakeDouble(8942832835564782, -385)); // digits 12, bits 61
    CheckDouble(MakeDouble(6965949469487146, -249)); // digits 13, bits 67
    CheckDouble(MakeDouble(6965949469487146, -250)); // digits 14, bits 65
    CheckDouble(MakeDouble(6965949469487146, -251)); // digits 15, bits 63
    CheckDouble(MakeDouble(7487252720986826, +548)); // digits 16, bits 63
    CheckDouble(MakeDouble(5592117679628511, +164)); // digits 17, bits 65
    CheckDouble(MakeDouble(8887055249355788, +665)); // digits 18, bits 67
    CheckDouble(MakeDouble(6994187472632449, +690)); // digits 19, bits 64
    CheckDouble(MakeDouble(8797576579012143, +588)); // digits 20, bits 62
    CheckDouble(MakeDouble(7363326733505337, +272)); // digits 21, bits 61
    CheckDouble(MakeDouble(8549497411294502, -448)); // digits 22, bits 66
}

TEST_CASE("Double - Regression")
{
    CheckDouble(1.5745340942675811e+257);
    CheckDouble(1.6521200219181297e-180);
    CheckDouble(4.6663180925160944e-302);
}

// Some numbers to check different code paths in grisu2::Dtoa
TEST_CASE("Double - Grisu2 code paths")
{
    CheckDouble(1e+4);
    CheckDouble(1.2e+6);
    CheckDouble(4.9406564584124654e-324);    // DigitGen: exit integral loop
    CheckDouble(2.2250738585072009e-308);    // DigitGen: exit fractional loop
    CheckDouble(1.82877982605164e-99);
    CheckDouble(1.1505466208671903e-09);
    CheckDouble(5.5645893133766722e+20);
    CheckDouble(53.034830388866226);
    CheckDouble(0.0021066531670178605);
}

TEST_CASE("Double - Round to even")
{
    CheckDoubleString(1.00000000000000005, "1");
    CheckDoubleString(1.00000000000000015, "1.0000000000000002"); // 1.000000000000000222...
    CheckDoubleString(1.99999999999999985, "1.9999999999999998"); // 1.999999999999999777...
    CheckDoubleString(1.99999999999999995, "2");
    CheckDoubleString(1125899906842623.75, "1125899906842623.8");
    CheckDoubleString(1125899906842624.25, "1125899906842624.2");
    CheckDoubleString(562949953421312.25, "562949953421312.2");
}

TEST_CASE("Double - Integers")
{
    CheckDoubleString(1.0, "1");
    CheckDoubleString(10.0, "10");
    CheckDoubleString(100.0, "100");
    CheckDoubleString(1000.0, "1000");
    CheckDoubleString(10000.0, "10000");
    CheckDoubleString(100000.0, "100000");
    CheckDoubleString(1000000.0, "1000000");
    CheckDoubleString(10000000.0, "10000000");
    CheckDoubleString(100000000.0, "100000000");
    CheckDoubleString(1000000000.0, "1000000000");
    CheckDoubleString(10000000000.0, "10000000000");
    CheckDoubleString(100000000000.0, "100000000000");
    CheckDoubleString(1000000000000.0, "1000000000000");
    CheckDoubleString(10000000000000.0, "10000000000000");
    CheckDoubleString(100000000000000.0, "100000000000000");
    CheckDoubleString(1000000000000000.0, "1000000000000000");
    CheckDoubleString(9007199254740000.0, "9007199254740000");
    CheckDoubleString(9007199254740992.0, "9007199254740992");
}
