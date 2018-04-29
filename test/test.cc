#include "../src/dtoa.h"
#if 0
#include "../src/strtod.h"
#endif

#include <double-conversion/double-conversion.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>

#define TEST_ALL_SINGLE         0
#define TEST_P1_DIGITS          0
#define TEST_RANDOM_DOUBLES     1
#define TEST_DTOA               0

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------

// Use the double-conversion library.
// MinGW's strtod and printf(%a) are almost completely useless...

static float StringToSingle(char const* str, char const* end)
{
    double_conversion::StringToDoubleConverter conv(0, 0.0, 0.0, "inf", "nan");

    int processed_characters_count = 0;
    return conv.StringToFloat(str, static_cast<int>(end - str), &processed_characters_count);
}

static double StringToDouble(char const* str, char const* end)
{
    double_conversion::StringToDoubleConverter conv(0, 0.0, 0.0, "inf", "nan");

    int processed_characters_count = 0;
    return conv.StringToDouble(str, static_cast<int>(end - str), &processed_characters_count);
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------

template <typename Target, typename Source>
static Target ReinterpretBits(Source source)
{
    static_assert(sizeof(Target) == sizeof(Source), "ouch");

    Target target;
    std::memcpy(&target, &source, sizeof(Source));
    return target;
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

// ldexp -- convert f * 2^e to IEEE single precision
#if 1
static float MakeSingle(uint64_t f, int e)
{
    constexpr uint64_t kHiddenBit = 0x00800000;
    constexpr uint64_t kSignificandMask = 0x007FFFFF;
    constexpr int kPhysicalSignificandSize = 23;  // Excludes the hidden bit.
    // constexpr int kSignificandSize = 24;
    constexpr int kExponentBias = 0x7F + kPhysicalSignificandSize;
    constexpr int kDenormalExponent = -kExponentBias + 1;
    constexpr int kMaxExponent = 0xFF - kExponentBias;

    assert(f <= kHiddenBit + kSignificandMask);
    while (f > kHiddenBit + kSignificandMask) {
        f >>= 1;
        e++;
    }

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
#endif

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

// ldexp -- convert f * 2^e to IEEE double precision
#if 1
static double MakeDouble(uint64_t f, int e)
{
    constexpr uint64_t kHiddenBit = 0x0010000000000000;
    constexpr uint64_t kSignificandMask = 0x000FFFFFFFFFFFFF;
    constexpr int kPhysicalSignificandSize = 52;  // Excludes the hidden bit.
    // constexpr int kSignificandSize = 53;
    constexpr int kExponentBias = 0x3FF + kPhysicalSignificandSize;
    constexpr int kDenormalExponent = -kExponentBias + 1;
    constexpr int kMaxExponent = 0x7FF - kExponentBias;

    assert(f <= kHiddenBit + kSignificandMask);
    while (f > kHiddenBit + kSignificandMask) {
        f >>= 1;
        e++;
    }

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
#endif

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------

static bool CheckFloat(float d0)
{
    auto const b0 = ReinterpretBits<uint32_t>(d0);

#if TEST_DTOA
    char str[32];
    auto const end = base_conv::Dtoa(str, str + 32, d0);
    *end = '\0';
    assert(end - str <= 26);
#else
    char str[1024*4];
    char const* end = str + std::snprintf(str, 1024*2, "%.1500g", d0);
#endif

    // printf("check single: %08x = '%s'\n", ReinterpretBits<uint32_t>(d0), str);

    bool result = true;
    {
        auto const d1 = StringToSingle(str, end);
        auto const b1 = ReinterpretBits<uint32_t>(d1);
        if (b0 != b1)
        {
            printf("FAIL: single: StringToSingle expected[%08x] != actual[%08x] -- [%s] [%.17g] [%.17g]\n", b0, b1, str, d0, d1);
            result = false;
        }
    }
    {
        auto const x1 = StringToDouble(str, end);
        auto const d1 = static_cast<float>(x1);
        auto const b1 = ReinterpretBits<uint32_t>(d1);
        if (b0 != b1)
        {
            printf("FAIL: single: (float)StringToDouble expected[%08x] != actual[%08x] -- [%s] [%.17g] [%.17g]\n", b0, b1, str, d0, d1);
            result = false;
        }
    }
#if 0
    {
        float d1;
        auto const ok = base_conv::Strtod(d1, str, end);
        assert(ok);
        auto const b1 = ReinterpretBits<uint32_t>(d1);
        if (b0 != b1)
        {
            printf("FAIL: single: base_conv::Strtof expected[%08x] != actual[%08x] -- [%s] [%.17g] [%.17g]\n", b0, b1, str, d0, d1);
            result = false;
        }
    }
    {
        double d1_;
        auto const ok = base_conv::Strtod(d1_, str, end);
        assert(ok);
        auto const d1 = static_cast<float>(d1_);
        auto const b1 = ReinterpretBits<uint32_t>(d1);
        if (b0 != b1)
        {
            printf("FAIL: single: (float)base_conv::Strtod expected[%08x] != actual[%08x] -- [%s] [%.17g] [%.17g]\n", b0, b1, str, d0, d1);
            result = false;
        }
    }
#endif

    return result;
}

static bool CheckFloat(double d0)
{
    auto const b0 = ReinterpretBits<uint64_t>(d0);

#if TEST_DTOA
    char str[32];
    auto const end = base_conv::Dtoa(str, str + 32, d0);
    *end = '\0';
    assert(end - str <= 26);
#else
    char str[1024*4];
    char const* end = str + std::snprintf(str, 1024*2, "%.1500g", d0);
#endif

    bool result = true;
    {
        auto const d1 = StringToDouble(str, end);
        auto const b1 = ReinterpretBits<uint64_t>(d1);
        if (b0 != b1)
        {
            printf("FAIL: double: StringToDouble expected=[%016llx] != actual[%016llx] -- [%s] expected=[%.17g] actual=[%.17g]\n", b0, b1, str, d0, d1);
            result = false;
        }
    }
#if 0
    {
        double d1;
        auto const ok = base_conv::Strtod(d1, str, end);
        assert(ok);
        auto const b1 = ReinterpretBits<uint64_t>(d1);
        if (b0 != b1)
        {
            printf("FAIL: double: base_conv::Strtod expected[%016llx] != actual[%016llx] -- [%s] [%.17g] [%.17g]\n", b0, b1, str, d0, d1);
            result = false;
        }
    }
#endif
#if 0
    {
        auto const d1 = std::strtod(str, nullptr);
        auto const b1 = ReinterpretBits<uint64_t>(d1);
        if (b0 != b1)
        {
            printf("FAIL: double: std::strtod expected=[%016llx] != actual[%016llx] -- [%s] expected=[%.17g] actual=[%.17g]\n", b0, b1, str, d0, d1);
            result = false;
        }
    }
#endif

    return result;
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------

static void VerifySingle()
{
    printf("Check single precision...\n");

    CheckFloat(MakeSingle(0,   0, 0x00000000)); // +0
    CheckFloat(MakeSingle(0,   0, 0x00000001)); // min denormal
    CheckFloat(MakeSingle(0,   0, 0x007FFFFF)); // max denormal
    CheckFloat(MakeSingle(0,   1, 0x00000000)); // min normal
    CheckFloat(MakeSingle(0,   1, 0x00000001));
    CheckFloat(MakeSingle(0,   1, 0x007FFFFF));
    CheckFloat(MakeSingle(0,   2, 0x00000000));
    CheckFloat(MakeSingle(0,   2, 0x00000001));
    CheckFloat(MakeSingle(0,  24, 0x00000000)); // fail if no special case in normalized boundaries
    CheckFloat(MakeSingle(0,  30, 0x00000000)); // fail if no special case in normalized boundaries
    CheckFloat(MakeSingle(0,  31, 0x00000000)); // fail if no special case in normalized boundaries
    CheckFloat(MakeSingle(0,  57, 0x00000000)); // fail if no special case in normalized boundaries
    CheckFloat(MakeSingle(0, 254, 0x007FFFFE));
    CheckFloat(MakeSingle(0, 254, 0x007FFFFF)); // max normal

    for (int e = 2; e < 254; ++e)
    {
        CheckFloat(MakeSingle(0, e-1, 0x007FFFFF));
        CheckFloat(MakeSingle(0, e,   0x00000000));
        CheckFloat(MakeSingle(0, e,   0x00000001));
    }

    // V. Paxson and W. Kahan, "A Program for Testing IEEE Binary-Decimal Conversion", manuscript, May 1991,
    // ftp://ftp.ee.lbl.gov/testbase-report.ps.Z    (report)
    // ftp://ftp.ee.lbl.gov/testbase.tar.Z          (program)

    // Table 16: Stress Inputs for Converting 24-bit Binary to Decimal, < 1/2 ULP
    CheckFloat(MakeSingle(12676506, -102)); // digits  1, bits 32
    CheckFloat(MakeSingle(12676506, -103)); // digits  2, bits 29
    CheckFloat(MakeSingle(15445013,  +86)); // digits  3, bits 34
    CheckFloat(MakeSingle(13734123, -138)); // digits  4, bits 32
    CheckFloat(MakeSingle(12428269, -130)); // digits  5, bits 30
    CheckFloat(MakeSingle(15334037, -146)); // digits  6, bits 31
    CheckFloat(MakeSingle(11518287,  -41)); // digits  7, bits 30
    CheckFloat(MakeSingle(12584953, -145)); // digits  8, bits 31
    CheckFloat(MakeSingle(15961084, -125)); // digits  9, bits 32
    CheckFloat(MakeSingle(14915817, -146)); // digits 10, bits 31
    CheckFloat(MakeSingle(10845484, -102)); // digits 11, bits 30
    CheckFloat(MakeSingle(16431059,  -61)); // digits 12, bits 29

    // Table 17: Stress Inputs for Converting 24-bit Binary to Decimal, > 1/2 ULP
    CheckFloat(MakeSingle(16093626,  +69)); // digits  1, bits 30
    CheckFloat(MakeSingle( 9983778,  +25)); // digits  2, bits 31
    CheckFloat(MakeSingle(12745034, +104)); // digits  3, bits 31
    CheckFloat(MakeSingle(12706553,  +72)); // digits  4, bits 31
    CheckFloat(MakeSingle(11005028,  +45)); // digits  5, bits 30
    CheckFloat(MakeSingle(15059547,  +71)); // digits  6, bits 31
    CheckFloat(MakeSingle(16015691,  -99)); // digits  7, bits 29
    CheckFloat(MakeSingle( 8667859,  +56)); // digits  8, bits 33
    CheckFloat(MakeSingle(14855922,  -82)); // digits  9, bits 35
    CheckFloat(MakeSingle(14855922,  -83)); // digits 10, bits 33
    CheckFloat(MakeSingle(10144164, -110)); // digits 11, bits 32
    CheckFloat(MakeSingle(13248074,  +95)); // digits 12, bits 33
}

static void VerifyDouble()
{
    printf("Check double precision...\n");

    CheckFloat(MakeDouble(0,    0, 0x0000000000000000)); // +0
    CheckFloat(MakeDouble(0,    0, 0x0000000000000001)); // min denormal
    CheckFloat(MakeDouble(0,    0, 0x000FFFFFFFFFFFFF)); // max denormal
    CheckFloat(MakeDouble(0,    1, 0x0000000000000000)); // min normal
    CheckFloat(MakeDouble(0,    1, 0x0000000000000001));
    CheckFloat(MakeDouble(0,    1, 0x000FFFFFFFFFFFFF));
    CheckFloat(MakeDouble(0,    2, 0x0000000000000000));
    CheckFloat(MakeDouble(0,    2, 0x0000000000000001));
    CheckFloat(MakeDouble(0,    4, 0x0000000000000000)); // fail if no special case in normalized boundaries
    CheckFloat(MakeDouble(0,    5, 0x0000000000000000)); // fail if no special case in normalized boundaries
    CheckFloat(MakeDouble(0,    6, 0x0000000000000000)); // fail if no special case in normalized boundaries
    CheckFloat(MakeDouble(0,   10, 0x0000000000000000)); // fail if no special case in normalized boundaries
    CheckFloat(MakeDouble(0, 2046, 0x000FFFFFFFFFFFFE));
    CheckFloat(MakeDouble(0, 2046, 0x000FFFFFFFFFFFFF)); // max normal

    for (int e = 2; e < 2046; ++e)
    {
        CheckFloat(MakeDouble(0, e-1, 0x000FFFFFFFFFFFFF));
        CheckFloat(MakeDouble(0, e,   0x0000000000000000));
        CheckFloat(MakeDouble(0, e,   0x0000000000000001));
    }

    // Some numbers to check different code paths in fast_dtoa
    CheckFloat(-1.0);
    CheckFloat(1e+4);
    CheckFloat(1.2e+6);
    CheckFloat(4.9406564584124654e-324);    // DigitGen: exit integral loop
    CheckFloat(2.2250738585072009e-308);    // DigitGen: exit fractional loop
    CheckFloat(1.82877982605164e-99);
    CheckFloat(1.1505466208671903e-09);
    CheckFloat(5.5645893133766722e+20);
    CheckFloat(53.034830388866226);
    CheckFloat(0.0021066531670178605);

    // V. Paxson and W. Kahan, "A Program for Testing IEEE Binary-Decimal Conversion", manuscript, May 1991,
    // ftp://ftp.ee.lbl.gov/testbase-report.ps.Z    (report)
    // ftp://ftp.ee.lbl.gov/testbase.tar.Z          (program)

    // Table 3: Stress Inputs for Converting 53-bit Binary to Decimal, < 1/2 ULP
    CheckFloat(MakeDouble(8511030020275656,  -342)); // digits  1, bits 63
    CheckFloat(MakeDouble(5201988407066741,  -824)); // digits  2, bits 63
    CheckFloat(MakeDouble(6406892948269899,  +237)); // digits  3, bits 62
    CheckFloat(MakeDouble(8431154198732492,   +72)); // digits  4, bits 61
    CheckFloat(MakeDouble(6475049196144587,   +99)); // digits  5, bits 64
    CheckFloat(MakeDouble(8274307542972842,  +726)); // digits  6, bits 64
    CheckFloat(MakeDouble(5381065484265332,  -456)); // digits  7, bits 64
    CheckFloat(MakeDouble(6761728585499734, -1057)); // digits  8, bits 64
    CheckFloat(MakeDouble(7976538478610756,  +376)); // digits  9, bits 67
    CheckFloat(MakeDouble(5982403858958067,  +377)); // digits 10, bits 63
    CheckFloat(MakeDouble(5536995190630837,   +93)); // digits 11, bits 63
    CheckFloat(MakeDouble(7225450889282194,  +710)); // digits 12, bits 66
    CheckFloat(MakeDouble(7225450889282194,  +709)); // digits 13, bits 64
    CheckFloat(MakeDouble(8703372741147379,  +117)); // digits 14, bits 66
    CheckFloat(MakeDouble(8944262675275217, -1001)); // digits 15, bits 63
    CheckFloat(MakeDouble(7459803696087692,  -707)); // digits 16, bits 63
    CheckFloat(MakeDouble(6080469016670379,  -381)); // digits 17, bits 62
    CheckFloat(MakeDouble(8385515147034757,  +721)); // digits 18, bits 64
    CheckFloat(MakeDouble(7514216811389786,  -828)); // digits 19, bits 64
    CheckFloat(MakeDouble(8397297803260511,  -345)); // digits 20, bits 64
    CheckFloat(MakeDouble(6733459239310543,  +202)); // digits 21, bits 63
    CheckFloat(MakeDouble(8091450587292794,  -473)); // digits 22, bits 63

    // Table 4: Stress Inputs for Converting 53-bit Binary to Decimal, > 1/2 ULP
    CheckFloat(MakeDouble(6567258882077402, +952)); // digits  1, bits 62
    CheckFloat(MakeDouble(6712731423444934, +535)); // digits  2, bits 65
    CheckFloat(MakeDouble(6712731423444934, +534)); // digits  3, bits 63
    CheckFloat(MakeDouble(5298405411573037, -957)); // digits  4, bits 62
    CheckFloat(MakeDouble(5137311167659507, -144)); // digits  5, bits 61
    CheckFloat(MakeDouble(6722280709661868, +363)); // digits  6, bits 64
    CheckFloat(MakeDouble(5344436398034927, -169)); // digits  7, bits 61
    CheckFloat(MakeDouble(8369123604277281, -853)); // digits  8, bits 65
    CheckFloat(MakeDouble(8995822108487663, -780)); // digits  9, bits 63
    CheckFloat(MakeDouble(8942832835564782, -383)); // digits 10, bits 66
    CheckFloat(MakeDouble(8942832835564782, -384)); // digits 11, bits 64
    CheckFloat(MakeDouble(8942832835564782, -385)); // digits 12, bits 61
    CheckFloat(MakeDouble(6965949469487146, -249)); // digits 13, bits 67
    CheckFloat(MakeDouble(6965949469487146, -250)); // digits 14, bits 65
    CheckFloat(MakeDouble(6965949469487146, -251)); // digits 15, bits 63
    CheckFloat(MakeDouble(7487252720986826, +548)); // digits 16, bits 63
    CheckFloat(MakeDouble(5592117679628511, +164)); // digits 17, bits 65
    CheckFloat(MakeDouble(8887055249355788, +665)); // digits 18, bits 67
    CheckFloat(MakeDouble(6994187472632449, +690)); // digits 19, bits 64
    CheckFloat(MakeDouble(8797576579012143, +588)); // digits 20, bits 62
    CheckFloat(MakeDouble(7363326733505337, +272)); // digits 21, bits 61
    CheckFloat(MakeDouble(8549497411294502, -448)); // digits 22, bits 66

#if 0
    // Table 20: Stress Inputs for Converting 56-bit Binary to Decimal, < 1/2 ULP
    CheckFloat(MakeDouble(50883641005312716, -172)); // digits  1, bits 65
    CheckFloat(MakeDouble(38162730753984537, -170)); // digits  2, bits 64
    CheckFloat(MakeDouble(50832789069151999, -101)); // digits  3, bits 64
    CheckFloat(MakeDouble(51822367833714164, -109)); // digits  4, bits 62
    CheckFloat(MakeDouble(66840152193508133, -172)); // digits  5, bits 64
    CheckFloat(MakeDouble(55111239245584393, -138)); // digits  6, bits 64
    CheckFloat(MakeDouble(71704866733321482, -112)); // digits  7, bits 62
    CheckFloat(MakeDouble(67160949328233173, -142)); // digits  8, bits 61
    CheckFloat(MakeDouble(53237141308040189, -152)); // digits  9, bits 63
    CheckFloat(MakeDouble(62785329394975786, -112)); // digits 10, bits 62
    CheckFloat(MakeDouble(48367680154689523,  -77)); // digits 11, bits 61
    CheckFloat(MakeDouble(42552223180606797, -102)); // digits 12, bits 62
    CheckFloat(MakeDouble(63626356173011241, -112)); // digits 13, bits 62
    CheckFloat(MakeDouble(43566388595783643,  -99)); // digits 14, bits 64
    CheckFloat(MakeDouble(54512669636675272, -159)); // digits 15, bits 61
    CheckFloat(MakeDouble(52306490527514614, -167)); // digits 16, bits 67
    CheckFloat(MakeDouble(52306490527514614, -168)); // digits 17, bits 65
    CheckFloat(MakeDouble(41024721590449423,  -89)); // digits 18, bits 62
    CheckFloat(MakeDouble(37664020415894738, -132)); // digits 19, bits 60
    CheckFloat(MakeDouble(37549883692866294,  -93)); // digits 20, bits 62
    CheckFloat(MakeDouble(69124110374399839, -104)); // digits 21, bits 65
    CheckFloat(MakeDouble(69124110374399839, -105)); // digits 22, bits 62

    // Table 21: Stress Inputs for Converting 56-bit Binary to Decimal, > 1/2 ULP
    CheckFloat(MakeDouble(49517601571415211,  -94)); // digits  1, bits 63
    CheckFloat(MakeDouble(49517601571415211,  -95)); // digits  2, bits 60
    CheckFloat(MakeDouble(54390733528642804, -133)); // digits  3, bits 63
    CheckFloat(MakeDouble(71805402319113924, -157)); // digits  4, bits 62
    CheckFloat(MakeDouble(40435277969631694, -179)); // digits  5, bits 61
    CheckFloat(MakeDouble(57241991568619049, -165)); // digits  6, bits 61
    CheckFloat(MakeDouble(65224162876242886,  +58)); // digits  7, bits 65
    CheckFloat(MakeDouble(70173376848895368, -138)); // digits  8, bits 61
    CheckFloat(MakeDouble(37072848117383207,  -99)); // digits  9, bits 61
    CheckFloat(MakeDouble(56845051585389697, -176)); // digits 10, bits 64
    CheckFloat(MakeDouble(54791673366936431, -145)); // digits 11, bits 64
    CheckFloat(MakeDouble(66800318669106231, -169)); // digits 12, bits 64
    CheckFloat(MakeDouble(66800318669106231, -170)); // digits 13, bits 61
    CheckFloat(MakeDouble(66574323440112438, -119)); // digits 14, bits 65
    CheckFloat(MakeDouble(65645179969330963, -173)); // digits 15, bits 62
    CheckFloat(MakeDouble(61847254334681076, -109)); // digits 16, bits 63
    CheckFloat(MakeDouble(39990712921393606, -145)); // digits 17, bits 62
    CheckFloat(MakeDouble(59292318184400283, -149)); // digits 18, bits 62
    CheckFloat(MakeDouble(69116558615326153, -143)); // digits 19, bits 65
    CheckFloat(MakeDouble(69116558615326153, -144)); // digits 20, bits 62
    CheckFloat(MakeDouble(39462549494468513, -152)); // digits 21, bits 63
    CheckFloat(MakeDouble(39462549494468513, -153)); // digits 22, bits 61
#endif
}

#if 0
static void TestStrtod()
{
    printf("TestStrtod...\n");

    auto check_double = [](std::string const& number, double expected)
    {
        double d;
        auto const ok = base_conv::Strtod(d, number.data(), number.data() + number.size());
        if (d != expected)
        {
            printf("FAIL: Strtod: \"%s\" --- actual: %.17g --- expected: %.17g\n", number.c_str(), d, expected);
        }
    };

    check_double("1e-2147483649", 1e-2147483649);
    check_double("1e-2147483648", 1e-2147483648);
    check_double("1e-2147483647", 1e-2147483647);
    check_double("1e+2147483647", std::numeric_limits<double>::infinity());
    check_double("1e+2147483648", std::numeric_limits<double>::infinity());
    check_double("1.7976931348623159e+308", std::numeric_limits<double>::infinity());
    check_double("1e-100000", 1e-100000);
    check_double("1e-1000", 1e-1000);
    check_double("1e-325", 1e-325);
    check_double("4.9406564584124653e-324", 4.9406564584124653e-324);
    check_double("4.94065645841246539999999999999999999999999999999999999999999999999999999999e-324", 4.94065645841246539999999999999999999999999999999999999999999999999999999999e-324);
    check_double("4.9406564584124654e-324", 4.9406564584124654e-324); // min denormal
    check_double("4.94065645841246540000000000000000000000000000000000000000000000000000000001e-324", 4.94065645841246540000000000000000000000000000000000000000000000000000000001e-324);
    check_double("4.9406564584124655e-324", 4.9406564584124655e-324);
    check_double("1e-324", 1e-324);
    check_double("2e-324", 2e-324);
    check_double("2.4703282292062327e-324", 0.0);
    check_double("2.4703282292062328e-324", 2.4703282292062328e-324);
    check_double("2.48e-324", 2.48e-324);
    check_double("2.5e-324", 2.5e-324);
    check_double("2.500000000000000000000000000000000000000000000000000000000000000000000000001e-324", 2.500000000000000000000000000000000000000000000000000000000000000000000000001e-324);
    check_double("3e-324", 3e-324);
    check_double("4e-324", 4e-324);
    check_double("5e-324", 5e-324); // min denormal
    check_double("2.225073858507201e-308", 2.225073858507201e-308); // max denormal
    check_double("2.2250738585072014e-308", 2.2250738585072014e-308); // min normal
    check_double("1.7976931348623157e+308", 1.7976931348623157e+308); // max normal
    check_double("1.7976931348623156999999999999999999999999999999999999999999999999999e+308", 1.7976931348623156999999999999999999999999999999999999999999999999999e+308);
    check_double("1.7976931348623157000000000000000000000000000000000000000000000000001e+308", 1.7976931348623157000000000000000000000000000000000000000000000000001e+308);
    check_double("1.7976931348623158e+308", 1.7976931348623158e+308);
    check_double(
        "1797693134862315708145274237317043567980705675258449965989174768031572607800285387605"
        "8955863276687817154045895351438246423432132688946418276846754670353751698604991057655"
        "1282076245490090389328944075868508455133942304583236903222948165808559332123348274797"
        "826204144723168738177180919299881250404026184124858368",
        std::numeric_limits<double>::max());
    check_double(
        "0.0000000000000000000000000000000000000000000000000000000000"
        "000000000000000000000000000000000000000000000000000000000000"
        "000000000000000000000000000000000000000000000000000000000000"
        "000000000000000000000000000000000000000000000000000000000000"
        "000000000000000000000000000000000000000000000000000000000000"
        "000000000000000000000000049406564584124654417656879286822137"
        "236505980261432476442558568250067550727020875186529983636163"
        "599237979656469544571773092665671035593979639877479601078187"
        "812630071319031140452784581716784898210368871863605699873072"
        "305000638740915356498438731247339727316961514003171538539807"
        "412623856559117102665855668676818703956031062493194527159149"
        "245532930545654440112748012970999954193198940908041656332452"
        "475714786901472678015935523861155013480352649347201937902681"
        "071074917033322268447533357208324319360923828934583680601060"
        "115061698097530783422773183292479049825247307763759272478746"
        "560847782037344696995336470179726777175851256605511991315048"
        "911014510378627381672509558373897335989936648099411642057026"
        "37090279242767544565229087538682506419718265533447265625",
        std::numeric_limits<double>::denorm_min());
    check_double(
        "243546080556034731077856379609316893158278902575447060151047"
        "212703405344938119816206067372775299130836050315842578309818"
        "316450894337978612745889730079163798234256495613858256849283"
        "467066859489192118352020514036083287319232435355752493038825"
        "828481044358810649108367633313557305310641892225870327827273"
        "41408256.000000",
        2.4354608055603473e+307);
    check_double("2.2250738585072011e-308", 2.2250738585072011e-308);
    check_double(
        "2.4703282292062327208828439643411068618252990130716238221279"
        "284125033775363510437593264991818081799618989828234772285886"
        "546332835517796989819938739800539093906315035659515570226392"
        "290858392449105184435931802849936536152500319370457678249219"
        "365623669863658480757001585769269903706311928279558551332927"
        "834338409351978015531246597263579574622766465272827220056374"
        "006485499977096599470454020828166226237857393450736339007967"
        "761930577506740176324673600968951340535537458516661134223766"
        "678604162159680461914467291840300530057530849048765391711386"
        "591646239524912623653881879636239373280423891018672348497668"
        "235089863388587925628302755995657524455507255189313690836254"
        "779186948667994968324049705821028513185451396213837722826145"
        "437693412532098591327667236328125001e-324",
        5e-324);
    check_double(
        "0.0000000000000000000000000000000000000000000000000000000000000000000000000000"
        "0000000000000000000000000000000000000000000000000000000000000000000000000000000"
        "00000000000000000000000000000000000000000000000000000000000000000000000000000"
        "000000000000000000000000000000000000000000000000000000000000000000000000000000"
        "00000000000002470328229206232720882843964341106861825299013071623822127928412503"
        "3775363510437593264991818081799618989828234772285886546332835517796989819938739"
        "80053909390631503565951557022639229085839244910518443593180284993653615250"
        "0319370457678249219365623669863658480757001585769269903706311928279558551332"
        "9278343384093519780155312465972635795746227664652728272200563740064854999770"
        "9659947045402082816622623785739345073633900796776193057750674017632467360096"
        "8951340535537458516661134223766678604162159680461914467291840300530057530849"
        "0487653917113865916462395249126236538818796362393732804238910186723484976682"
        "3508986338858792562830275599565752445550725518931369083625477918694866799496"
        "8324049705821028513185451396213837722826145437693412532098591327667236328125", 0.0);
    check_double(
        "0.0000000000000000000000000000000000000000000000000000000000000000000000000000"
        "0000000000000000000000000000000000000000000000000000000000000000000000000000000"
        "00000000000000000000000000000000000000000000000000000000000000000000000000000"
        "000000000000000000000000000000000000000000000000000000000000000000000000000000"
        "00000000000002470328229206232720882843964341106861825299013071623822127928412503"
        "3775363510437593264991818081799618989828234772285886546332835517796989819938739"
        "80053909390631503565951557022639229085839244910518443593180284993653615250"
        "0319370457678249219365623669863658480757001585769269903706311928279558551332"
        "9278343384093519780155312465972635795746227664652728272200563740064854999770"
        "9659947045402082816622623785739345073633900796776193057750674017632467360096"
        "8951340535537458516661134223766678604162159680461914467291840300530057530849"
        "0487653917113865916462395249126236538818796362393732804238910186723484976682"
        "3508986338858792562830275599565752445550725518931369083625477918694866799496"
        "8324049705821028513185451396213837722826145437693412532098591327667236328126", 5e-324);
    check_double("0.500000000000000166533453693773481063544750213623046875", 0.500000000000000166533453693773481063544750213623046875);
    check_double("3.518437208883201171875e13", 3.518437208883201171875e13);
    check_double("62.5364939768271845828", 62.5364939768271845828);
    check_double("8.10109172351e-10", 8.10109172351e-10);
    check_double("1.50000000000000011102230246251565404236316680908203125", 1.50000000000000011102230246251565404236316680908203125);
    check_double("9007199254740991.4999999999999999999999999999999995", 9007199254740991.4999999999999999999999999999999995);
    check_double("1.2345678901234567e22", 1.2345678901234567e22);
    check_double("2.2250738585072011e-308", 2.2250738585072011e-308);
    check_double(
        "6.6312368714697582767853966302759672433990999473553031442499717"
        "587362866301392654396180682007880487441059604205526018528897150"
        "063763256665955396033303618005191075917832333584923372080578494"
        "993608994251286407188566165030934449228547591599881603044399098"
        "682919739314266256986631577498362522745234853124423586512070512"
        "924530832781161439325697279187097860044978723221938561502254152"
        "119972830784963194121246401117772161481107528151017752957198119"
        "743384519360959074196224175384736794951486324803914359317679811"
        "223967034438033355297560033532098300718322306892013830155987921"
        "841729099279241763393155074022348361207309147831684007154624400"
        "538175927027662135590421159867638194826541287705957668068727833"
        "49146967171293949598850675682115696218943412532098591327667236328125E-316",
        6.631236846766476e-316);
    check_double(
        "3.2378839133029012895883524125015321748630376694231080599012970"
        "495523019706706765657868357425877995578606157765598382834355143"
        "910841531692526891905643964595773946180389283653051434639551003"
        "566966656292020173313440317300443693602052583458034314716600326"
        "995807313009548483639755486900107515300188817581841745696521731"
        "104736960227499346384253806233697747365600089974040609674980283"
        "891918789639685754392222064169814626901133425240027243859416510"
        "512935526014211553334302252372915238433223313261384314778235911"
        "424088000307751706259156707286570031519536642607698224949379518"
        "458015308952384398197084033899378732414634842056080000272705311"
        "068273879077914449185347715987501628125488627684932015189916680"
        "28251730299953143924168545708663913273994694463908672332763671875E-319",
        3.2379086165851934e-319);
    check_double(
        "6.953355807847677105972805215521891690222119817145950754416205607980030"
        "13154963668880611572639944188006538639986402869127553953941465283158479"
        "56685600829998895513577849614468960421131982842131079351102171626549398"
        "02416034676213829409720583759540476786936413816541621287843248433202369"
        "20991661224967600557302270324479971462211654218883777037602237117207955"
        "91258533828013962195524188394697705149041926576270603193728475623010741"
        "40442660237844114174497210955449896389180395827191602886654488182452409"
        "58398138944278337700150546201574501784875457466834216175949666176602002"
        "87528887833870748507731929971029979366198762266880963149896457660004790"
        "09083731736585750335262099860150896718774401964796827166283225641992040"
        "747894382698751809812609536720628966577351093292236328125E-310",
        6.9533558078476524e-310);
    check_double(
        "3.339068557571188581835713701280943911923401916998521771655656997328440"
        "31455961531816884914907466260909999811300946556642680817037843406572299"
        "16596426194677060348844249897410807907667784563321682004646515939958173"
        "71782125010668346652995912233993254584461125868481633343674905074271064"
        "40976309070801785658401977687881242531200881232626036303547481153223685"
        "33599053346255754042160606228586332807443018924703005556787346899784768"
        "70369853549413277156622170245846166991655321535529623870646888786637528"
        "99559280043617790174628627227337447170145299143304725786386460142425202"
        "47915673681950560773208853293843223323915646452641434007986196650406080"
        "77549162173963649264049738362290606875883456826586710961041737908872035"
        "803481241600376705491726170293986797332763671875E-319",
        3.3390932608534806e-319);
    check_double("2.2250738585072012e-308", 2.2250738585072012e-308);
    check_double("2.2250738585072011e-308", 2.2250738585072011e-308);

    check_double("6114917000000003e-14", 6114917000000003e-14);
}
#endif

//------------------------------------------------------------------------------
// "7.038531e-26"
//
// is the only single-precision float, which does not round-trip with
// (float)strtod but with strtof
//------------------------------------------------------------------------------
// exp = 43
// FAIL: single strtod [15ae43fd] != [15ae43fe] -- [7.038531e-26] [7.0385306918512091e-26] [7.0385313081487913e-26]
//------------------------------------------------------------------------------
// strtof("7.038531e-26")
//  f   = 15AE'43FD                         (IEEE bits)
//      = 1010'1110'0100'0011'1111'1101     (IEEE bits)
//      = 11420669 * 2^-107
//      = 7.038530691851209120859188017140306974105991300039164570989669300615787506103515625 * 10^-26
//
//  f-  = 15AE43FC                          (IEEE bits)
//      = 1010'1110'0100'0011'1111'1100     (IEEE bits)
//      = 11420668 * 2^-107
//      = 7.0385300755536269169437150392273653469292493678466371420654468238353729248046875 * 10^-26
//
//  f+  = 15AE43FE                          (IEEE bits)
//      = 1010'1110'0100'0011'1111'1110     (IEEE bits)
//      = 11420670 * 2^-107
//      = 7.03853130814879132477466099505324860128273323223169199991389177739620208740234375 * 10^-26
//
// strtod("7.038531e-26")
//  d   = 3AB5C87FB0000000
//      = 6131425250115584 * 2^-136
//      = 7.0385310000000002228169245060967777876943622661354282854517805390059947967529296875 * 10^-26
//
//  d - f- =  3 / 324518553658426726783156020576256
//         =  9.244463733058732094668694124407651128982887911433863337151706218719482421875 * 10^-33
//  d - f  =  1 / 324518553658426726783156020576256
//         =  3.081487911019577364889564708135883709660962637144621112383902072906494140625 * 10^-33
//  d - f+ = -1 / 324518553658426726783156020576256
//         = -3.081487911019577364889564708135883709660962637144621112383902072906494140625 * 10^-33
//
// Cast d to single precision: (round to nearest, ties to even)
//  ==> f+
//------------------------------------------------------------------------------
// From:
// http://www.exploringbinary.com/floating-point-converter/
//
// strtof("7.0385307e-26") = 15AE43FD
//                         = 11420669 * 2^-107
// strtod("7.0385307e-26") = 3AB5C87FA06C50E6
//                         = 3065712494389363 * 2^-135
//                         = 6131424988778726 * 2^-136
//------------------------------------------------------------------------------
//   0 <= exp <= 114 ==> all optimal
// 149 <= exp <= 151 ==> all optimal
// 184 <= exp <= 255 ==> all optimal
//
//      XXX:      115 <= exp <= 183
//
//------------------------------------------------------------------------------
#if TEST_ALL_SINGLE
static void TestAllSingle()
{
    printf("Testing all finite single precision values...\n");

    using Clock = std::chrono::steady_clock;

    int const min_exp = 0;
    int const max_exp = (1 << 8) - 1; // exclusive!

    uint64_t num_tested = 0;
    uint64_t num_shortest = 0;
    uint64_t num_optimal = 0;

    uint32_t bits = min_exp << 23;

    auto const t_beg = Clock::now();

    auto t_lap = t_beg;
    auto curr_exp = min_exp;

    printf("exp = %d\n", curr_exp);
    for (;;)
    {
        float const f = ReinterpretBits<float>(bits);

#if 0
        if (f != 0.0)
        {
            ++num_tested;

            char buf1[32];
            char buf2[32];
            int len1 = 0;
            int len2 = 0;
            {
                auto const boundaries = base_conv::ComputeBoundaries(f);
                int k1 = 0;
                base_conv::Grisu2(buf1, len1, k1, boundaries.m_minus, boundaries.v, boundaries.m_plus);
            }
            {
                using double_conversion::DoubleToStringConverter;

                bool sign = false;
                int point = 0;
                DoubleToStringConverter::DoubleToAscii(f, DoubleToStringConverter::SHORTEST_SINGLE, -1 /*unused*/, buf2, 32, &sign, &len2, &point);
            }
            assert(len1 >= len2);
            if (len1 == len2)
            {
                ++num_shortest;
                if (memcmp(buf1, buf2, len1) == 0)
                {
                    ++num_optimal;
                }
            }
        }
#else
        CheckFloat(f);
#endif

        ++bits;

        int const next_exp = bits >> 23;
        if (next_exp == max_exp || curr_exp != next_exp)
        {
            auto const t_now = Clock::now();
            printf("   time: %f sec\n", std::chrono::duration<double>(t_now - t_lap).count());
            // printf("   num_tested   %u\n", num_tested);
            // printf("   num_shortest %u %.17g%%\n", num_shortest, 100.0 * (static_cast<double>(num_shortest) / static_cast<double>(num_tested)));
            // printf("   num_optimal  %u %.17g%%\n", num_optimal,  100.0 * (static_cast<double>(num_optimal)  / static_cast<double>(num_tested)));
            printf("exp = %d\n", next_exp);
            t_lap = t_now;
            if (next_exp == max_exp)
                break;
            // num_tested = 0;
            // num_shortest = 0;
            // num_optimal = 0;
        }

        curr_exp = next_exp;
    }

    auto const t_end = Clock::now();
    printf("all-floats time: %f sec\n", std::chrono::duration<double>(t_end - t_beg).count());
    printf("   num_tested   %llu\n", num_tested);
    printf("   num_shortest %llu %.17g%%\n", num_shortest, 100.0 * (static_cast<double>(num_shortest) / static_cast<double>(num_tested)));
    printf("   num_optimal  %llu %.17g%%\n", num_optimal,  100.0 * (static_cast<double>(num_optimal)  / static_cast<double>(num_tested)));
}
#endif

struct RandomDoubles
{
    // Test uniformly distributed bit patterns instead of uniformly distributed
    // floating-points...

    std::random_device rd_;
    std::mt19937_64 random_;
    std::uniform_int_distribution<uint64_t> gen_;

    RandomDoubles()
        : rd_()
        , random_(rd_())
        , gen_(0, (uint64_t{0x7FF} << 52) - 1)
    {
    }

    double operator()()
    {
        auto bits = gen_(random_);
        return ReinterpretBits<double>(bits);
    }
};

struct RandomUniformDoubles
{
    // Test uniformly distributed bit patterns instead of uniformly distributed
    // floating-points...

    std::random_device rd_;
    std::mt19937_64 random_;
    std::uniform_real_distribution<double> gen_;

    RandomUniformDoubles()
        : rd_()
        , random_(rd_())
        , gen_(0.0, std::numeric_limits<double>::max())
        //, gen_(0.0, 1.0)
        //, gen_(1.0e-10, 1.0e-1)
    {
    }

    double operator()()
    {
        return gen_(random_);
    }
};

#if TEST_RANDOM_DOUBLES
static void TestDoubles()
{
    printf("Testing random double precision values...\n");

    using Clock = std::chrono::steady_clock;

    //RandomDoubles rng;
    RandomUniformDoubles rng;

    auto t_start = Clock::now();

    uint64_t num_checked = 0;
    uint64_t num_shortest = 0;
    uint64_t num_optimal = 0;

    uint64_t const kNumDoubles = uint64_t{1} << 30;
    for (uint64_t i = 0; i < kNumDoubles; ++i)
    {
        double const value = rng();

        CheckFloat(value);
        num_checked++;

#if 1
        char buf1[32];
        char buf2[32];
        int len1 = 0;
        int len2 = 0;
        {
            auto const boundaries = base_conv::ComputeBoundaries(value);
            int k = 0;
            base_conv::Grisu2(buf1, len1, k, boundaries.m_minus, boundaries.v, boundaries.m_plus);
        }
        {
            using double_conversion::DoubleToStringConverter;

            bool sign = false;
            int point = 0;
            DoubleToStringConverter::DoubleToAscii(value, DoubleToStringConverter::SHORTEST, -1 /*unused*/, buf2, 32, &sign, &len2, &point);
        }
        assert(len1 >= len2);
        if (len1 == len2)
        {
            ++num_shortest;
            if (memcmp(buf1, buf2, len1) == 0)
                ++num_optimal;
        }
#endif

        auto const t_now = Clock::now();
        auto const t_sec = std::chrono::duration<double>(t_now - t_start).count();
        if (t_sec > 5.0)
        {
            fprintf(stderr, "%.2f%% [fp/sec %.3f] [shortest: %.17g%%] [optimal: %.17g%%]\n",
                100.0 * (static_cast<double>(i) / static_cast<double>(kNumDoubles)),
                num_checked / 1000.0 / t_sec,
                100.0 * (static_cast<double>(num_shortest) / static_cast<double>(num_checked)),
                100.0 * (static_cast<double>(num_optimal ) / static_cast<double>(num_checked))
                );
            t_start = t_now;
            num_checked = 0;
            num_shortest = 0;
            num_optimal = 0;
        }
    }
}
#endif

#if 0
static void TestDoubleRange()
{
    printf("Testing some finite double precision values...\n");

    using Clock = std::chrono::steady_clock;

    int const min_exp = 0;
    int const max_exp = (1 << 11) - 1; // exclusive!
    int curr_exp = min_exp;

    uint64_t num_tested = 0;
    uint64_t num_shortest = 0;
    uint64_t num_optimal = 0;

    auto const t_beg = Clock::now();
    auto t_lap = t_beg;

    std::mt19937 rng;

    printf("exp = %d\n", curr_exp);
    for (;;)
    {
        std::uniform_int_distribution<uint64_t> gen(0, (uint64_t{1} << 52) - 1);
        uint64_t significand = gen(rng);

        double const f = ReinterpretBits<double>((uint64_t(curr_exp) << 52) | significand);

        if (curr_exp != 0 || significand != 0)
        {
            ++num_tested;

            char buf1[32];
            char buf2[32];
            int len1 = 0;
            int len2 = 0;
            {
                auto const boundaries = base_conv::ComputeBoundaries(f);
                int k1 = 0;
                base_conv::Grisu2(buf1, len1, k1, boundaries.m_minus, boundaries.v, boundaries.m_plus);
            }
            {
                using double_conversion::DoubleToStringConverter;

                bool sign = false;
                int point = 0;
                DoubleToStringConverter::DoubleToAscii(f, DoubleToStringConverter::SHORTEST, -1 /*unused*/, buf2, 32, &sign, &len2, &point);
            }
            assert(len1 >= len2);
            if (len1 == len2)
            {
                ++num_shortest;
                if (memcmp(buf1, buf2, len1) == 0)
                {
                    ++num_optimal;
                }
            }
        }

        if (num_tested > (1 << 10))
        {        
            auto const t_now = Clock::now();

            printf("   time: %f sec\n", std::chrono::duration<double>(t_now - t_lap).count());
            printf("   num_shortest %.17g%%\n", 100.0 * (static_cast<double>(num_shortest) / static_cast<double>(num_tested)));
            printf("   num_optimal  %.17g%%\n", 100.0 * (static_cast<double>(num_optimal)  / static_cast<double>(num_tested)));
            if (++curr_exp == max_exp)
                break;
            printf("exp = %d\n", curr_exp);
            num_tested = 0;
            num_shortest = 0;
            num_optimal = 0;

            t_lap = t_now;
        }
    }

    auto const t_end = Clock::now();
    printf("some-doubles time: %f sec\n", std::chrono::duration<double>(t_end - t_beg).count());
}
#endif

static void FindMaxP1()
{
    constexpr int kExpMin = -1137;
    constexpr int kExpMax =   960;
    constexpr uint64_t kMaxF = UINT64_MAX; // ((uint64_t{1} << 53) - 1) << 11;

    uint64_t max_p1 = 0;
    uint64_t min_p1 = UINT64_MAX;
    for (int e = kExpMin; e <= kExpMax; ++e)
    {
        auto const v = base_conv::DiyFp(kMaxF, e);
        auto const cached = base_conv::GetCachedPowerForBinaryExponent(e);
        auto const c_minus_k = base_conv::DiyFp(cached.f, cached.e);
        auto const w = base_conv::Multiply(v, c_minus_k);
        auto const p1 = w.f >> -w.e;
        if (max_p1 < p1)
            max_p1 = p1;
        if (min_p1 > p1)
            min_p1 = p1;
    }

    printf("max_p1 = %llu [%llX]\n", max_p1, max_p1);
    printf("min_p1 = %llu [%llX]\n", min_p1, min_p1);
}

#if TEST_P1_DIGITS
static int CountDecimalDigits(uint32_t n)
{
    if (n >= 1000000000) { return 10; }
    if (n >=  100000000) { return  9; }
    if (n >=   10000000) { return  8; }
    if (n >=    1000000) { return  7; }
    if (n >=     100000) { return  6; }
    if (n >=      10000) { return  5; }
    if (n >=       1000) { return  4; }
    if (n >=        100) { return  3; }
    if (n >=         10) { return  2; }
    return 1;
}

static void TestP1Digits()
{
    printf("Testing P1 integral distribution...\n");

    using Clock = std::chrono::steady_clock;

    //RandomDoubles rng;
    RandomUniformDoubles rng;

    auto t_start = Clock::now();

    uint64_t hist[11] = {0,0,0,0,0, 0,0,0,0,0, 0};
    uint64_t num_checked = 0;

    uint64_t const kNumDoubles = uint64_t{1} << 30;
    for (uint64_t i = 0; i < kNumDoubles; ++i)
    {
        double const value = rng();
        ++num_checked;

        auto const boundaries = grisu::ComputeBoundaries(value);
        auto const cached = grisu::GetCachedPowerForBinaryExponent(boundaries.v.e);
        auto const c_minus_k = grisu::DiyFp(cached.f, cached.e);
        auto const w_plus = grisu::Multiply(boundaries.m_plus, c_minus_k);
        auto const p1 = static_cast<uint32_t>(w_plus.f >> -w_plus.e);

        hist[CountDecimalDigits(p1)]++;

        auto const t_now = Clock::now();
        auto const t_sec = std::chrono::duration<double>(t_now - t_start).count();
        if (t_sec > 5.0)
        {
            for (int k = 1; k <= 10; ++k) {
                fprintf(stderr, "hist[%2d] = %.3f%%\n", k, hist[k] / static_cast<double>(num_checked));
                //hist[k] = 0;
            }
            //num_checked = 0;

            t_start = t_now;
        }
    }
}
#endif

int main()
{
#if 0
    const float fff = ReinterpretBits<float>(0x15AE43FDu);
    printf("%.9g\n", fff);
    printf("%08x\n", ReinterpretBits<uint32_t>(fff));
    char buf[32];
    char* bufend = base_conv::Dtoa(buf, buf + 32, fff);
    *bufend = '\0';
    printf("|%s|\n", buf);
    //char const* inp = "7.0385307e-26";
    //char const* inp = "7.038531e-26";
    //char const* inp = "7.03853069e-26";
    char const* inp = buf;
    char const* end = inp + std::strlen(inp);
    uint32_t b0;
    uint32_t b1;
    {
        double d = StringToDouble(inp, end);
        auto const v = base_conv::DiyFpFromFloat(d);
        printf("%lld * 2^%d\n", v.f, v.e);
        b0 = ReinterpretBits<uint32_t>(static_cast<float>(d));
    }
    {
        double d;
        base_conv::Strtod(d, inp, end);
        auto const v = base_conv::DiyFpFromFloat(d);
        printf("%lld * 2^%d\n", v.f, v.e);
        b1 = ReinterpretBits<uint32_t>(static_cast<float>(d));
    }
    printf("b0 = %08x\n", b0);
    printf("b1 = %08x\n", b1);
#else
    FindMaxP1();

    VerifySingle();
    VerifyDouble();
#if 0
    TestStrtod();
#endif
#if 0
   TestDoubleRange();
#endif
#if TEST_ALL_SINGLE
    TestAllSingle();
#endif
#if TEST_P1_DIGITS
    TestP1Digits();
#endif
#if TEST_RANDOM_DOUBLES
    TestDoubles();
#endif
#endif
}
