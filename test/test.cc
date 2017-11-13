#include "../src/fast_dtoa.h"

#include <double-conversion/double-conversion.h>

#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>

#define TEST_ALL_SINGLE         0
#define TEST_RANDOM_DOUBLES     0
#define TEST_DOUBLE_CONVERSION  0

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

#if TEST_DOUBLE_CONVERSION
static char* SingleToString(char* next, char* last, float value)
{
    auto const& conv = double_conversion::DoubleToStringConverter::EcmaScriptConverter();

    double_conversion::StringBuilder builder(next, static_cast<int>(last - next));
    conv.ToShortestSingle(value, &builder);

    char* end = next + builder.position();
    *end = '\0';
    return end;
}

static char* DoubleToString(char* next, char* last, double value)
{
    auto const& conv = double_conversion::DoubleToStringConverter::EcmaScriptConverter();

    double_conversion::StringBuilder builder(next, static_cast<int>(last - next));
    conv.ToShortest(value, &builder);

    char* end = next + builder.position();
    *end = '\0';
    return end;
}
#endif

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
    char str[32];
#if TEST_DOUBLE_CONVERSION
    auto const end = SingleToString(str, str + 32, d0);
#else
    auto const end = fast_dtoa::ToString(str, str + 32, d0);
#endif
    assert(end - str <= 26);
    *end = '\0';

    // printf("check single: %08x = '%s'\n", ReinterpretBits<uint32_t>(d0), str);

    {
        auto const d1 = StringToSingle(str, end);
        auto const b0 = ReinterpretBits<uint32_t>(d0);
        auto const b1 = ReinterpretBits<uint32_t>(d1);
        if (b0 != b1)
        {
            printf("FAIL: single strtof [%08x] != [%08x] -- [%s] [%.17g] [%.17g]\n", b0, b1, str, d0, d1);
            return false;
        }
    }

#if 1
    {
        auto const d1 = static_cast<float>(StringToDouble(str, end));
        auto const b0 = ReinterpretBits<uint32_t>(d0);
        auto const b1 = ReinterpretBits<uint32_t>(d1);
        if (b0 != b1)
        {
            printf("FAIL: single strtod [%08x] != [%08x] -- [%s] [%.17g] [%.17g]\n", b0, b1, str, d0, d1);
            return false;
        }
    }
#endif

    return true;
}

static bool CheckFloat(double d0)
{
    char str[32];
#if TEST_DOUBLE_CONVERSION
    auto const end = DoubleToString(str, str + 32, d0);
#else
    auto const end = fast_dtoa::ToString(str, str + 32, d0);
#endif
    assert(end - str <= 26);
    *end = '\0';

    // printf("check double: %016llx = '%s'\n", ReinterpretBits<uint64_t>(d0), str);

    {
        auto const d1 = StringToDouble(str, end);
        auto const b0 = ReinterpretBits<uint64_t>(d0);
        auto const b1 = ReinterpretBits<uint64_t>(d1);
        if (b0 != b1)
        {
            printf("FAIL: double [%016llx] != [%016llx] -- [%s] [%.17g] [%.17g]\n", b0, b1, str, d0, d1);
            return false;
        }
    }

    return true;
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
}

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
//
// strtof("7.0385307e-26") = 15AE43FD
// strtod("7.0385307e-26") = 3AB5C87FA06C50E6
//                         = 6131424988778726 * 2^-136
//------------------------------------------------------------------------------
#if TEST_ALL_SINGLE
static void TestAllSingle()
{
    printf("Testing all finite single precision values...\n");

    using Clock = std::chrono::steady_clock;

    int const min_exp = 0;
    int const max_exp = (1 << 8) - 1; // exclusive!

    uint32_t bits = min_exp << 23;

    int curr_exp = min_exp;
    printf("exp = %d\n", curr_exp);

    auto const t_beg = Clock::now();

    auto t_lap = t_beg;
    for (;;)
    {
        float const f = ReinterpretBits<float>(bits);
        CheckFloat(f);

        ++bits;

        int next_exp = bits >> 23;
        if (next_exp == max_exp)
        {
            auto const t_now = Clock::now();
            printf("   time: %f sec\n", std::chrono::duration<double>(t_now - t_lap).count());
            break;
        }

        if (curr_exp != next_exp)
        {
            auto const t_now = Clock::now();
            printf("   time: %f sec\n", std::chrono::duration<double>(t_now - t_lap).count());
            printf("exp = %d\n", next_exp);
            t_lap = t_now;
        }

        curr_exp = next_exp;
    }

    auto const t_end = Clock::now();
    printf("all-floats time: %f sec\n", std::chrono::duration<double>(t_end - t_beg).count());
}
#endif

#if TEST_RANDOM_DOUBLES
struct RandomDoubles
{
    // Test uniformly distributed bit patterns instead of uniformly distributed
    // floating-points...

    std::random_device rd_;
    std::mt19937_64 random_;
    std::uniform_int_distribution<uint64_t> gen_;

#if 0
    RandomDoubles()
        : random_()
        , gen_(0, (uint64_t{0x7FF} << 52) - 1)
    {
    }
#else
    RandomDoubles()
        : rd_()
        , random_(rd_())
        , gen_(0, (uint64_t{0x7FF} << 52) - 1)
    {
    }
#endif

    double operator()()
    {
        auto const bits = gen_(random_);
        return ReinterpretBits<double>(bits);
    }
};

static void TestDoubles()
{
    printf("Testing random double precision values...\n");

    using Clock = std::chrono::steady_clock;

    RandomDoubles rng;

    uint64_t const kNumDoubles = uint64_t{1} << 30;

    auto t_start = Clock::now();

    uint64_t num_processed = 0;
    for (uint64_t i = 0; i < kNumDoubles; ++i)
    {
        CheckFloat(rng());
        num_processed++;

        auto const t_now = Clock::now();
        auto const t_sec = std::chrono::duration<double>(t_now - t_start).count();
        if (t_sec > 5.0)
        {
            fprintf(stderr, "%.2f%% [fp/sec %.2f]\n", 100.0 * (double)i / (double)kNumDoubles, num_processed / 1000.0 / t_sec);
            t_start = t_now;
            num_processed = 0;
        }
    }
}
#endif

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------

int main()
{
    VerifySingle();
    VerifyDouble();
#if TEST_ALL_SINGLE
    TestAllSingle();
#endif
#if TEST_RANDOM_DOUBLES
    TestDoubles();
#endif
}
