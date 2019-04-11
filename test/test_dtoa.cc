//#include "floaxie.h"
//#include "fmt.h"
#include "grisu2.h"
//#include "milo.h"

#include "catch.hpp"

#include <string>
#include <limits>
#include <iostream>

#include <double-conversion/double-conversion.h>

#define USE_DOUBLE_CONVERSION_DTOA 0
#define TEST_HEX 0
#define TEST_STRTOD_FOR_SINGLE 0

static inline char* Float32ToChars(char* buf, int buflen, float f)
{
    auto const& f2s = double_conversion::DoubleToStringConverter::EcmaScriptConverter();
    double_conversion::StringBuilder builder(buf, buflen);
    f2s.ToShortestSingle(f, &builder);
    int const length = builder.position();
    builder.Finalize();
    return buf + length;
}

static inline char* Float64ToChars(char* buf, int buflen, double f)
{
    auto const& f2s = double_conversion::DoubleToStringConverter::EcmaScriptConverter();
    double_conversion::StringBuilder builder(buf, buflen);
    f2s.ToShortest(f, &builder);
    int const length = builder.position();
    builder.Finalize();
    return buf + length;
}

static inline std::string Float32ToString(float f)
{
    char buf[32];
    Float32ToChars(buf, 32, f);
    return buf;
}

static inline std::string Float64ToString(double f)
{
    char buf[32];
    Float64ToChars(buf, 64, f);
    return buf;
}

static inline char* Dtoa(char* next, char* last, double value)
{
    //char* end = Float64ToChars(next, (int)(last - next), value);
    char* end = grisu2_Dtoa(next, last, value);
    //char* end = milo_Dtoa(next, last, value);
    //char* end = fmt_Dtoa(next, last, value);
    //char* end = floaxie_Dtoa(next, last, value);         // FAIL: incorrect (lower) boundary ?!!! TODO: https://github.com/aclex/floaxie/issues
    end[0] = '\0';
    return end;
}
static inline std::string Dtoa(double value)
{
    char buf[32];
    Dtoa(buf, buf + 32, value);
    return buf;
}

static inline char* Ftoa(char* next, char* last, float value)
{
    //char* end = Float32ToChars(next, (int)(last - next), value);
    char* end = grisu2_Ftoa(next, last, value);
    //char* end = fmt_Ftoa(next, last, value);
    //char* end = floaxie_Ftoa(next, last, value);         // FAIL: incorrect (lower) boundary ?!!! TODO: https://github.com/aclex/floaxie/issues
    end[0] = '\0';
    return end;
}
static inline std::string Ftoa(float value)
{
    char buf[32];
    Ftoa(buf, buf + 32, value);
    return buf;
}

template <typename Dest, typename Source>
inline Dest ReinterpretBits(Source source)
{
    static_assert(sizeof(Dest) == sizeof(Source), "size mismatch");

    Dest dest;
    std::memcpy(&dest, &source, sizeof(Source));
    return dest;
}

static bool CheckSingle(float f)
{
    // decimal
    {
#if USE_DOUBLE_CONVERSION_DTOA
        auto const& f2s = double_conversion::DoubleToStringConverter::EcmaScriptConverter();

        char buf[32];
        double_conversion::StringBuilder builder(buf, 32);
        f2s.ToShortestSingle(f, &builder);
        int const length = builder.position();
        builder.Finalize();
#else
        char buf[32];
        char* first = buf;
        char* last  = Ftoa(first, first + 32, f);
        *last = '\0';
        int const length = static_cast<int>(last - first);
#endif

        double_conversion::StringToDoubleConverter s2f(0, 0.0, 0.0, "inf", "nan");

        int processed_characters_count = 0;
        {
            auto const f_out = s2f.StringToFloat(buf, length, &processed_characters_count);

            const uint32_t f_bits = ReinterpretBits<uint32_t>(f);
            const uint32_t f_out_bits = ReinterpretBits<uint32_t>(f_out);
            CAPTURE(f_bits);
            CAPTURE(f_out_bits);
            CAPTURE(&buf[0]);
            CHECK(f == f_out);
        }
#if TEST_STRTOD_FOR_SINGLE
        {
            auto const d_out = s2f.StringToDouble(buf, length, &processed_characters_count);
            auto const f_out = static_cast<float>(d_out);

            const uint32_t f_bits = ReinterpretBits<uint32_t>(f);
            const uint32_t f_out_bits = ReinterpretBits<uint32_t>(f_out);
            CAPTURE(f_bits);
            CAPTURE(f_out_bits);
            CAPTURE(&buf[0]);
            CHECK(f == f_out);
        }
#endif
    }

#if TEST_HEX
    // hex
    {
        uint32_t f_bits = ReinterpretBits<uint32_t>(f);
        CAPTURE(f_bits);

        char buf[32] = {0};
        char* first = buf;
        char* last  = base_conv_HexDtoa(first, first + 32, f);

        std::string str = {first, last};
        if (str[0] == '-')
            str.insert(1, "0x");
        else
            str.insert(0, "0x");

        CAPTURE(str);

        auto const f_out = std::strtof(str.c_str(), nullptr);
        uint32_t f_out_bits = ReinterpretBits<uint32_t>(f_out);
        CAPTURE(f_out_bits);

        CHECK(f == f_out);
    }
#endif

    return true;
}

static bool CheckDouble(double f)
{
    // decimal
    {
#if USE_DOUBLE_CONVERSION_DTOA
        auto const& f2s = double_conversion::DoubleToStringConverter::EcmaScriptConverter();

        char buf[32];
        double_conversion::StringBuilder builder(buf, 32);
        f2s.ToShortest(f, &builder);
        int const length = builder.position();
        builder.Finalize();
#else
        char buf[32];
        char* first = buf;
        char* last  = Dtoa(first, first + 32, f);
        *last = '\0';
        int const length = static_cast<int>(last - first);
#if 0
        {
            auto const& f2s = double_conversion::DoubleToStringConverter::EcmaScriptConverter();

            char bufdc[32];
            double_conversion::StringBuilder builder(bufdc, 32);
            f2s.ToShortest(f, &builder);
            builder.Finalize();

            CHECK(std::string(buf) == std::string(bufdc));
        }
#endif
#endif

        double_conversion::StringToDoubleConverter s2f(0, 0.0, 0.0, "inf", "nan");

        int processed_characters_count = 0;
        auto const f_out = s2f.StringToDouble(buf, length, &processed_characters_count);

        const uint64_t f_bits = ReinterpretBits<uint64_t>(f);
        const uint64_t f_out_bits = ReinterpretBits<uint64_t>(f_out);
        CAPTURE(f_bits);
        CAPTURE(f_out_bits);
        CAPTURE(&buf[0]);
        CHECK(f == f_out);
    }

#if TEST_HEX
    // hex
    {
        char buf[32] = {0};
        char* first = buf;
        char* last  = base_conv_HexDtoa(first, first + 32, f);

        std::string str = {first, last};
        if (str[0] == '-')
            str.insert(1, "0x");
        else
            str.insert(0, "0x");

        auto const f_out = std::strtod(str.c_str(), nullptr);

        CHECK(f == f_out);
    }
#endif

    return true;
}

static std::string Dtostr(float value, bool force_trailing_dot_zero = false)
{
    char buf[32];
    char* first = buf;
    char* last  = grisu2_Dtoa(first, first + 32, value, force_trailing_dot_zero);
    return {first, last};
}

static std::string Dtostr(double value, bool force_trailing_dot_zero = false)
{
    char buf[32];
    char* first = buf;
    char* last  = grisu2_Dtoa(first, first + 32, value, force_trailing_dot_zero);
    return {first, last};
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

#define CHECK_SINGLE(VALUE) CHECK(CheckSingle(VALUE))
#define CHECK_DOUBLE(VALUE) CHECK(CheckDouble(VALUE))

TEST_CASE("Single - xxx")
{
    CHECK_SINGLE(7.0385307e-26f);
}

TEST_CASE("Dtoa - single 1")
{
    CHECK_SINGLE(MakeSingle(0,   0, 0x00000000)); // +0
    CHECK_SINGLE(MakeSingle(0,   0, 0x00000001)); // min denormal
    CHECK_SINGLE(MakeSingle(0,   0, 0x007FFFFF)); // max denormal
    CHECK_SINGLE(MakeSingle(0,   1, 0x00000000)); // min normal
    CHECK_SINGLE(MakeSingle(0,   1, 0x00000001));
    CHECK_SINGLE(MakeSingle(0,   1, 0x007FFFFF));
    CHECK_SINGLE(MakeSingle(0,   2, 0x00000000));
    CHECK_SINGLE(MakeSingle(0,   2, 0x00000001));
    CHECK_SINGLE(MakeSingle(0,  24, 0x00000000)); // fail if no special case in normalized boundaries
    CHECK_SINGLE(MakeSingle(0,  30, 0x00000000)); // fail if no special case in normalized boundaries
    CHECK_SINGLE(MakeSingle(0,  31, 0x00000000)); // fail if no special case in normalized boundaries
    CHECK_SINGLE(MakeSingle(0,  57, 0x00000000)); // fail if no special case in normalized boundaries
    CHECK_SINGLE(MakeSingle(0, 254, 0x007FFFFE));
    CHECK_SINGLE(MakeSingle(0, 254, 0x007FFFFF)); // max normal

    for (uint32_t e = 2; e < 254; ++e)
    {
        CAPTURE(e);
        CHECK_SINGLE(MakeSingle(0, e-1, 0x007FFFFF));
        CHECK_SINGLE(MakeSingle(0, e,   0x00000000));
        CHECK_SINGLE(MakeSingle(0, e,   0x00000001));
    }

    // V. Paxson and W. Kahan, "A Program for Testing IEEE Binary-Decimal Conversion", manuscript, May 1991,
    // ftp://ftp.ee.lbl.gov/testbase-report.ps.Z    (report)
    // ftp://ftp.ee.lbl.gov/testbase.tar.Z          (program)

    // Table 16: Stress Inputs for Converting 24-bit Binary to Decimal, < 1/2 ULP
    CHECK_SINGLE(MakeSingle(12676506, -102)); // digits  1, bits 32
    CHECK_SINGLE(MakeSingle(12676506, -103)); // digits  2, bits 29
    CHECK_SINGLE(MakeSingle(15445013,  +86)); // digits  3, bits 34
    CHECK_SINGLE(MakeSingle(13734123, -138)); // digits  4, bits 32
    CHECK_SINGLE(MakeSingle(12428269, -130)); // digits  5, bits 30
    CHECK_SINGLE(MakeSingle(15334037, -146)); // digits  6, bits 31
    CHECK_SINGLE(MakeSingle(11518287,  -41)); // digits  7, bits 30
    CHECK_SINGLE(MakeSingle(12584953, -145)); // digits  8, bits 31
    CHECK_SINGLE(MakeSingle(15961084, -125)); // digits  9, bits 32
    CHECK_SINGLE(MakeSingle(14915817, -146)); // digits 10, bits 31
    CHECK_SINGLE(MakeSingle(10845484, -102)); // digits 11, bits 30
    CHECK_SINGLE(MakeSingle(16431059,  -61)); // digits 12, bits 29

    // Table 17: Stress Inputs for Converting 24-bit Binary to Decimal, > 1/2 ULP
    CHECK_SINGLE(MakeSingle(16093626,  +69)); // digits  1, bits 30
    CHECK_SINGLE(MakeSingle( 9983778,  +25)); // digits  2, bits 31
    CHECK_SINGLE(MakeSingle(12745034, +104)); // digits  3, bits 31
    CHECK_SINGLE(MakeSingle(12706553,  +72)); // digits  4, bits 31
    CHECK_SINGLE(MakeSingle(11005028,  +45)); // digits  5, bits 30
    CHECK_SINGLE(MakeSingle(15059547,  +71)); // digits  6, bits 31
    CHECK_SINGLE(MakeSingle(16015691,  -99)); // digits  7, bits 29
    CHECK_SINGLE(MakeSingle( 8667859,  +56)); // digits  8, bits 33
    CHECK_SINGLE(MakeSingle(14855922,  -82)); // digits  9, bits 35
    CHECK_SINGLE(MakeSingle(14855922,  -83)); // digits 10, bits 33
    CHECK_SINGLE(MakeSingle(10144164, -110)); // digits 11, bits 32
    CHECK_SINGLE(MakeSingle(13248074,  +95)); // digits 12, bits 33
}

TEST_CASE("Dtoa - double 1")
{
    CHECK_DOUBLE(MakeDouble(0,    0, 0x0000000000000000)); // +0
    CHECK_DOUBLE(MakeDouble(0,    0, 0x0000000000000001)); // min denormal
    CHECK_DOUBLE(MakeDouble(0,    0, 0x000FFFFFFFFFFFFF)); // max denormal
    CHECK_DOUBLE(MakeDouble(0,    1, 0x0000000000000000)); // min normal
    CHECK_DOUBLE(MakeDouble(0,    1, 0x0000000000000001));
    CHECK_DOUBLE(MakeDouble(0,    1, 0x000FFFFFFFFFFFFF));
    CHECK_DOUBLE(MakeDouble(0,    2, 0x0000000000000000));
    CHECK_DOUBLE(MakeDouble(0,    2, 0x0000000000000001));
    CHECK_DOUBLE(MakeDouble(0,    4, 0x0000000000000000)); // fail if no special case in normalized boundaries
    CHECK_DOUBLE(MakeDouble(0,    5, 0x0000000000000000)); // fail if no special case in normalized boundaries
    CHECK_DOUBLE(MakeDouble(0,    6, 0x0000000000000000)); // fail if no special case in normalized boundaries
    CHECK_DOUBLE(MakeDouble(0,   10, 0x0000000000000000)); // fail if no special case in normalized boundaries
    CHECK_DOUBLE(MakeDouble(0, 2046, 0x000FFFFFFFFFFFFE));
    CHECK_DOUBLE(MakeDouble(0, 2046, 0x000FFFFFFFFFFFFF)); // max normal

    for (uint64_t e = 2; e < 2046; ++e)
    {
        CAPTURE(e);
        CHECK_DOUBLE(MakeDouble(0, e-1, 0x000FFFFFFFFFFFFF));
        CHECK_DOUBLE(MakeDouble(0, e,   0x0000000000000000));
        CHECK_DOUBLE(MakeDouble(0, e,   0x0000000000000001));
    }

    // Some numbers to check different code paths in fast_dtoa
    CHECK_DOUBLE(-1.0);
    CHECK_DOUBLE(1e+4);
    CHECK_DOUBLE(1.2e+6);
    CHECK_DOUBLE(4.9406564584124654e-324);    // DigitGen: exit integral loop
    CHECK_DOUBLE(2.2250738585072009e-308);    // DigitGen: exit fractional loop
    CHECK_DOUBLE(1.82877982605164e-99);
    CHECK_DOUBLE(1.1505466208671903e-09);
    CHECK_DOUBLE(5.5645893133766722e+20);
    CHECK_DOUBLE(53.034830388866226);
    CHECK_DOUBLE(0.0021066531670178605);

    // V. Paxson and W. Kahan, "A Program for Testing IEEE Binary-Decimal Conversion", manuscript, May 1991,
    // ftp://ftp.ee.lbl.gov/testbase-report.ps.Z    (report)
    // ftp://ftp.ee.lbl.gov/testbase.tar.Z          (program)

    // Table 3: Stress Inputs for Converting 53-bit Binary to Decimal, < 1/2 ULP
    CHECK_DOUBLE(MakeDouble(8511030020275656,  -342)); // digits  1, bits 63
    CHECK_DOUBLE(MakeDouble(5201988407066741,  -824)); // digits  2, bits 63
    CHECK_DOUBLE(MakeDouble(6406892948269899,  +237)); // digits  3, bits 62
    CHECK_DOUBLE(MakeDouble(8431154198732492,   +72)); // digits  4, bits 61
    CHECK_DOUBLE(MakeDouble(6475049196144587,   +99)); // digits  5, bits 64
    CHECK_DOUBLE(MakeDouble(8274307542972842,  +726)); // digits  6, bits 64
    CHECK_DOUBLE(MakeDouble(5381065484265332,  -456)); // digits  7, bits 64
    CHECK_DOUBLE(MakeDouble(6761728585499734, -1057)); // digits  8, bits 64
    CHECK_DOUBLE(MakeDouble(7976538478610756,  +376)); // digits  9, bits 67
    CHECK_DOUBLE(MakeDouble(5982403858958067,  +377)); // digits 10, bits 63
    CHECK_DOUBLE(MakeDouble(5536995190630837,   +93)); // digits 11, bits 63
    CHECK_DOUBLE(MakeDouble(7225450889282194,  +710)); // digits 12, bits 66
    CHECK_DOUBLE(MakeDouble(7225450889282194,  +709)); // digits 13, bits 64
    CHECK_DOUBLE(MakeDouble(8703372741147379,  +117)); // digits 14, bits 66
    CHECK_DOUBLE(MakeDouble(8944262675275217, -1001)); // digits 15, bits 63
    CHECK_DOUBLE(MakeDouble(7459803696087692,  -707)); // digits 16, bits 63
    CHECK_DOUBLE(MakeDouble(6080469016670379,  -381)); // digits 17, bits 62
    CHECK_DOUBLE(MakeDouble(8385515147034757,  +721)); // digits 18, bits 64
    CHECK_DOUBLE(MakeDouble(7514216811389786,  -828)); // digits 19, bits 64
    CHECK_DOUBLE(MakeDouble(8397297803260511,  -345)); // digits 20, bits 64
    CHECK_DOUBLE(MakeDouble(6733459239310543,  +202)); // digits 21, bits 63
    CHECK_DOUBLE(MakeDouble(8091450587292794,  -473)); // digits 22, bits 63

    // Table 4: Stress Inputs for Converting 53-bit Binary to Decimal, > 1/2 ULP
    CHECK_DOUBLE(MakeDouble(6567258882077402, +952)); // digits  1, bits 62
    CHECK_DOUBLE(MakeDouble(6712731423444934, +535)); // digits  2, bits 65
    CHECK_DOUBLE(MakeDouble(6712731423444934, +534)); // digits  3, bits 63
    CHECK_DOUBLE(MakeDouble(5298405411573037, -957)); // digits  4, bits 62
    CHECK_DOUBLE(MakeDouble(5137311167659507, -144)); // digits  5, bits 61
    CHECK_DOUBLE(MakeDouble(6722280709661868, +363)); // digits  6, bits 64
    CHECK_DOUBLE(MakeDouble(5344436398034927, -169)); // digits  7, bits 61
    CHECK_DOUBLE(MakeDouble(8369123604277281, -853)); // digits  8, bits 65
    CHECK_DOUBLE(MakeDouble(8995822108487663, -780)); // digits  9, bits 63
    CHECK_DOUBLE(MakeDouble(8942832835564782, -383)); // digits 10, bits 66
    CHECK_DOUBLE(MakeDouble(8942832835564782, -384)); // digits 11, bits 64
    CHECK_DOUBLE(MakeDouble(8942832835564782, -385)); // digits 12, bits 61
    CHECK_DOUBLE(MakeDouble(6965949469487146, -249)); // digits 13, bits 67
    CHECK_DOUBLE(MakeDouble(6965949469487146, -250)); // digits 14, bits 65
    CHECK_DOUBLE(MakeDouble(6965949469487146, -251)); // digits 15, bits 63
    CHECK_DOUBLE(MakeDouble(7487252720986826, +548)); // digits 16, bits 63
    CHECK_DOUBLE(MakeDouble(5592117679628511, +164)); // digits 17, bits 65
    CHECK_DOUBLE(MakeDouble(8887055249355788, +665)); // digits 18, bits 67
    CHECK_DOUBLE(MakeDouble(6994187472632449, +690)); // digits 19, bits 64
    CHECK_DOUBLE(MakeDouble(8797576579012143, +588)); // digits 20, bits 62
    CHECK_DOUBLE(MakeDouble(7363326733505337, +272)); // digits 21, bits 61
    CHECK_DOUBLE(MakeDouble(8549497411294502, -448)); // digits 22, bits 66
}

TEST_CASE("Dtoa format special")
{
    CHECK("NaN" == Dtostr(std::numeric_limits<float>::quiet_NaN()));
    CHECK("NaN" == Dtostr(std::numeric_limits<double>::quiet_NaN()));
    CHECK("Infinity" == Dtostr(std::numeric_limits<float>::infinity()));
    CHECK("Infinity" == Dtostr(std::numeric_limits<double>::infinity()));
    CHECK("-Infinity" == Dtostr(-std::numeric_limits<float>::infinity()));
    CHECK("-Infinity" == Dtostr(-std::numeric_limits<double>::infinity()));
    CHECK("-0" == Dtostr(-0.0f));
    CHECK("-0" == Dtostr(-0.0));
}

TEST_CASE("Dtoa format trailing dot-zero")
{
    CHECK("0" == Dtostr(0.0f, false));
    CHECK("0.0" == Dtostr(0.0f, true));
    CHECK("10" == Dtostr(10.0f, false));
    CHECK("10.0" == Dtostr(10.0f, true));

    CHECK("0" == Dtostr(0.0, false));
    CHECK("0.0" == Dtostr(0.0, true));
    CHECK("10" == Dtostr(10.0, false));
    CHECK("10.0" == Dtostr(10.0, true));
}
