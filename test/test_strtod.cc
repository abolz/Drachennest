#include "catch.hpp"

#include "ryu_32.h"
#include "ryu_64.h"

#include "double-conversion/double-conversion.h"

#include <cmath>
#include <cstring>
#include <limits>
#include <cfloat>

static float strtof_double_conversion(const char* buf, int len)
{
    double_conversion::StringToDoubleConverter s2d(0, 0.0, 1.0, "inf", "nan");
    int processed_characters_count = 0;
    return s2d.StringToFloat(buf, len, &processed_characters_count);
}

//==================================================================================================
// Strtof
//==================================================================================================
#if 1

static float FloatFromBits(uint32_t bits)
{
    float f;
    std::memcpy(&f, &bits, sizeof(uint32_t));
    return f;
}

static uint32_t BitsFromFloat(float f)
{
    uint32_t u;
    std::memcpy(&u, &f, sizeof(uint32_t));
    return u;
}

static float Strtof(const std::string& str)
{
    float flt;
    const auto res = ryu::Strtof(str.data(), str.data() + str.size(), flt);
    CHECK(res.status != ryu::StrtofStatus::invalid);
    return flt;
}

struct Ftoa1 {
    char* operator()(char* buf, int /*buflen*/, float value) const {
        return ryu::Ftoa(buf, value);
    }
};

struct Ftoa2 {
    char* operator()(char* buf, int buflen, float value) const {
        return buf + std::snprintf(buf, buflen, "%.9g", value);
    }
};

struct Ftoa3 {
    char* operator()(char* buf, int buflen, float value) const {
        return buf + std::snprintf(buf, buflen, "%.8e", value);
    }
};

template <typename FtoaFn>
static bool CheckStrtofImpl(float value)
{
    static constexpr int BUFLEN = 128;

    const auto bits = BitsFromFloat(value);
    char buf[BUFLEN];

    char* const end = FtoaFn{}(buf, BUFLEN, value);
    *end = '\0';

    float value2;
    const auto res = ryu::Strtof(buf, end, value2);
    CHECK(res.status != ryu::StrtofStatus::invalid);
    CHECK(res.next == end);

    if (std::isnan(value))
    {
        //CHECK(std::isnan(value2));
        if (!std::isnan(value2))
        {
            fprintf(stderr, "FAIL: nan 0x%08X\n", bits);
            return false;
        }
    }
    else
    {
        const auto bits2 = BitsFromFloat(value2);
        //CHECK(bits == bits2);
        if (bits != bits2)
        {
            fprintf(stderr, "FAIL: bits = 0x%08X != 0x%08X = bits2\n", bits, bits2);
            return false;
        }
    }

    return true;
}

static bool CheckStrtof(float value)
{
    if (!CheckStrtofImpl<Ftoa1>(value))
        return false;
    if (!CheckStrtofImpl<Ftoa2>(value))
        return false;
    if (!CheckStrtofImpl<Ftoa3>(value))
        return false;

    return true;
}

TEST_CASE("Strtof - Regression")
{
    CheckStrtof(FloatFromBits(0x00400001u));
    CheckStrtof(FloatFromBits(0x00800000u));
    CheckStrtof(FloatFromBits(0x00800001u));
    CheckStrtof(FloatFromBits(0x01000000u));

    CheckStrtof(16777215.0f);
    CheckStrtof(16777216.0f);
    CheckStrtof(16777217.0f); // == 16777216.0f
    CheckStrtof(16777218.0f);

    CheckStrtof(100000000.0f);
    CheckStrtof(10000000.0f);
    CheckStrtof(1000000.0f);

    CHECK(0.0f == Strtof("0.00000001e-45"));
    CHECK(0.0f == Strtof("0.00000001e-46"));
    CHECK(0.0f == Strtof("1.00000000e-46"));
    CHECK(0.0f == Strtof("1.00000000e-47"));

    CHECK(70064924e-53f == Strtof("70064924e-53"));

    CHECK(2.68435495e+07f == Strtof("2.68435495e+07"));
    CHECK(5.00000025e+07f == Strtof("5.00000025e+07"));
    CHECK(9.99999895e+07f == Strtof("9.99999895e+07"));

    CHECK(1.17549429e-38f == Strtof("1.17549429e-38"));
    CHECK(1.17549430e-38f == Strtof("1.17549430e-38"));
    CHECK(1.17549431e-38f == Strtof("1.17549431e-38"));
    CHECK(1.17549432e-38f == Strtof("1.17549432e-38"));
    CHECK(1.17549433e-38f == Strtof("1.17549433e-38"));
    CHECK(1.17549434e-38f == Strtof("1.17549434e-38"));
    CHECK(1.17549435e-38f == Strtof("1.17549435e-38"));
}

TEST_CASE("Strtof - 1")
{
    CheckStrtof(std::numeric_limits<float>::min());
    CheckStrtof(std::numeric_limits<float>::max());
    CheckStrtof(std::numeric_limits<float>::denorm_min());
    CheckStrtof(std::numeric_limits<float>::epsilon());

    CHECK(999999999.0f == Strtof("999999999"));
    CHECK(9999.00009f == Strtof("9999.00009"));
    CHECK(999999999.0f == Strtof("999999999e+00"));
    CHECK(999999999.0f == Strtof("99999999900000000e-8"));
    CHECK(0.00000000999999999f == Strtof("0.00000000999999999"));
    CHECK(9999.0009f == Strtof("9999.000900000000000000000000000"));
    CHECK(9999.0009f == Strtof("9999.000900000000000000000000000e+0"));
    CHECK(999999999.0f == Strtof("999999999.0"));
    CHECK(999999999.0f == Strtof("999999999.0000000000000000000000000000000000000000000000000000000000000000000000e+00"));
    CHECK(0.000999999999f == Strtof("0.000999999999"));
}

TEST_CASE("Strtof - Special")
{
    CHECK( 0.0f == Strtof("0"));
    CHECK( 0.0f == Strtof("0.0000000000000000000000000000000"));
    CHECK(-0.0f == Strtof("-0"));
    CHECK( 0.0f == Strtof("+0"));

    CheckStrtof(0.0f);
    CheckStrtof(-0.0f);
    CheckStrtof(std::numeric_limits<float>::infinity());
    CheckStrtof(-std::numeric_limits<float>::infinity());
    CheckStrtof(std::numeric_limits<float>::quiet_NaN());

    CHECK(std::isnan(Strtof("nan")));
    CHECK(std::isnan(Strtof("NaN")));
    CHECK(std::isnan(Strtof("nAn(_nananana123)")));
    CHECK(std::isnan(Strtof("nan(")));
    CHECK(std::isnan(Strtof("nan(123")));

    CHECK(std::isinf(Strtof("Inf")));
    CHECK(std::isinf(Strtof("Infinity")));
    CHECK(std::isinf(Strtof("-INF")));
}

TEST_CASE("Strtof - Long input")
{
    CHECK(1280.0f == Strtof("128.000000000000000000000000000000000000000000000000000000000000"
                            "0000000000000000000000000000000000000000000000000000000000000e+1"));
    CHECK(1280.0f == Strtof("128.000000000000000000000000000000000000000000000000000000000000"
                            "00000000000000000000000000000000000000000000000000000000000000000e+1"));
    //CHECK(1280.0f == Strtof("128.0000000000000000000000000000000000000000000000000000000000000000"
    //                            "0000000000000000000000000000000000000000000000000000000000000000e+1"));
}

static bool incbuf(char* buf, int n)
{
    for (int i = n - 1; i >= 0; --i) {
        if (buf[i] == '.')
            continue;
        if (buf[i] != '9') {
            ++buf[i];
            return true;
        }
        buf[i] = '0';
    }

    return false;
}

#if 0
TEST_CASE("Strtof - All")
{
    printf("Strtof - All...\n");

    char buf[1 + 1 + 8 + 1 + 1 + 2 + 1] = {0};
    buf[ 0] = '0';
    buf[ 1] = '.';
    buf[ 2] = '0';
    buf[ 3] = '0';
    buf[ 4] = '0';
    buf[ 5] = '0';
    buf[ 6] = '0';
    buf[ 7] = '0';
    buf[ 8] = '0';
    buf[ 9] = '0';
    buf[10] = 'e';
    buf[11] = '+';
    buf[12] = '0';
    buf[13] = '0';
    buf[14] = '\0';

    const int e_min = -31; // -38; // -45; // 8; // 7; // 0; // -35; // -44; // -45;
    const int e_max =  38;

    constexpr uint32_t Decimals   = 999999999;
    constexpr uint32_t OnePercent =  10000000;

    //uint32_t d = 268435495;
    //uint32_t d = 268435535;
    //uint32_t d = 268435575;
    //uint32_t d = 272598455;
    //uint32_t d = 277331735;
    //uint32_t d   = 500000000;
    // 2.68435495e+07
    // 5.00000025e+07
    // 9.99999895e+07
    //uint32_t d = 0;
    uint32_t d = 33 * OnePercent;
    for (int e = e_min; e <= e_max; ++e)
    {
        printf("Exp = %d\n", e);

        uint32_t percent_done = (uint32_t)-1;
        for ( ; d <= 999999999; ++d)
        {
            uint32_t percent_done_new = (uint32_t)(int32_t)(d / 10000000.0);
            if (percent_done != percent_done_new)
            {
                percent_done = percent_done_new;
                printf("%3u%%\n", percent_done);
            }

            sprintf(buf + 1, "%.09ue%+.02d", d, e);

            float flt1;
#if USE_RYU_UPSTREAM()
            const auto res1 = s2f_n(buf, 14, &flt1);
            assert(res1 == Status::SUCCESS);
#else
            const auto res1 = charconv::Strtof(buf, buf + 14, flt1);
            assert(res1.status == charconv::StrtofStatus::ok);

            //const auto res1 = charconv::Strtof(d, e, flt1);
            //if (res1 != charconv::StrtofStatus::ok)
            //    continue;
#endif

            // ref:
            //float flt2 = std::strtof(buf, nullptr);
            // ref:
            //float flt2;
            //std::from_chars(buf, buf + 14, flt2);
            // ref:
            float flt2 = strtof_double_conversion(buf, 14);

            if (flt1 != flt2)
            {
                printf("FAIL: %.9g != %.9g\n", flt1, flt2);
                printf("  buf = '%s'\n", buf);
                //return;
            }

            //incbuf(buf, 9);
        }
        d = 0;
    }

    printf("OK.\n");
}
#endif

#endif // 0

//==================================================================================================
// Strtod
//==================================================================================================
#if 1

static double FloatFromBits(uint64_t bits)
{
    double f;
    std::memcpy(&f, &bits, sizeof(uint64_t));
    return f;
}

static uint64_t BitsFromFloat(double f)
{
    uint64_t u;
    std::memcpy(&u, &f, sizeof(uint64_t));
    return u;
}

//static double MakeDouble(uint32_t e, uint64_t f)
//{
//    return FloatFromBits(uint64_t{e} << 52 | f);
//}

static double Strtod(const std::string& str)
{
    double flt;
    const auto res = ryu::Strtod(str.data(), str.data() + str.size(), flt);
    CHECK(res.status != ryu::StrtodStatus::invalid);
    //const auto res = absl::from_chars(str.data(), str.data() + str.size(), flt);
    //const auto res = std::from_chars(str.data(), str.data() + str.size(), flt);
    return flt;
}

struct Dtoa1 {
    char* operator()(char* buf, int /*buflen*/, double value) const {
        return ryu::Dtoa(buf, value);
    }
};

struct Dtoa2 {
    char* operator()(char* buf, int buflen, double value) const {
        return buf + std::snprintf(buf, buflen, "%.17g", value);
    }
};

struct Dtoa3 {
    char* operator()(char* buf, int buflen, double value) const {
        return buf + std::snprintf(buf, buflen, "%.16e", value);
    }
};

template <typename DtoaFn>
static bool CheckStrtodImpl(double value)
{
    static constexpr int BUFLEN = 128;

    const auto bits = BitsFromFloat(value);
    char buf[BUFLEN];

    char* const end = DtoaFn{}(buf, BUFLEN, value);
    *end = '\0';

    double value2;
    const auto res = ryu::Strtod(buf, end, value2);
    CHECK(res.status != ryu::StrtodStatus::invalid);
    CHECK(res.next == end);

    if (std::isnan(value))
    {
        CHECK(std::isnan(value2));
        //if (!std::isnan(value2))
        //{
        //    fprintf(stderr, "FAIL: nan 0x%016llX\n", bits);
        //    return false;
        //}
    }
    else
    {
        const auto bits2 = BitsFromFloat(value2);
        CHECK(bits == bits2);
        //if (bits != bits2)
        //{
        //    fprintf(stderr, "FAIL: bits = 0x%016llX != 0x%016llX = bits2\n", bits, bits2);
        //    return false;
        //}
    }

    return true;
}

static void CheckStrtod(double value)
{
    CheckStrtodImpl<Dtoa1>(value);
    CheckStrtodImpl<Dtoa2>(value);
    CheckStrtodImpl<Dtoa3>(value);
}

TEST_CASE("Strtod - 1")
{
    CheckStrtod(std::numeric_limits<double>::min());
    CheckStrtod(std::numeric_limits<double>::max());
    CheckStrtod(std::numeric_limits<double>::denorm_min());
    CheckStrtod(std::numeric_limits<double>::epsilon());

    CheckStrtod(9007199254740991.0);
    CheckStrtod(9007199254740992.0);
    CheckStrtod(9007199254740993.0); // == 9007199254740992.0
    CheckStrtod(9007199254740994.0);

    CheckStrtod(10000000000000000.0);
    CheckStrtod(1000000000000000.0);
    CheckStrtod(100000000000000.0);

    CheckStrtod(1e23);
    CHECK(7.2057594037927933e+16 == Strtod("7.2057594037927933e+16"));
}

TEST_CASE("Strtod - regression")
{
    CHECK(1.2999999999999999E+154 == Strtod("1.2999999999999999E+154"));
    CHECK(7.3177701707893310e+15 == Strtod("7.3177701707893310e+15"));
    CHECK(7.2057594037927933e+16 == Strtod("7.2057594037927933e+16"));

    for (int i = 0; i < 53; ++i)
    {
        CheckStrtod(FloatFromBits(uint64_t{1} << i));
    }

    CheckStrtod(FloatFromBits(uint64_t{0x1} << 51));
    CheckStrtod(FloatFromBits(uint64_t{0x1} << 52));
    CheckStrtod(FloatFromBits(uint64_t{0x1} << 53));
    CheckStrtod(FloatFromBits(uint64_t{0x3} << 51));
    CheckStrtod(FloatFromBits(uint64_t{0x3} << 52));
    CheckStrtod(FloatFromBits(uint64_t{0x3} << 53));

    double d = std::numeric_limits<double>::denorm_min();
    for (int i = 0; i < 100; ++i)
    {
        CheckStrtod(d);
        d *= 2;
    }

    d = std::numeric_limits<double>::denorm_min();
    for (int i = 0; i < 100; ++i)
    {
        CheckStrtod(d);
        d /= 2;
    }

    CHECK(0.0 == Strtod("0.0000000000000001e-325"));
    CHECK(0.0 == Strtod("1.0000000000000000e-325"));
    CHECK(0.0 == Strtod("0.0000000000000001e-324"));
    CHECK(0.0 == Strtod("0.0000000000000010e-324"));
    CHECK(0.0 == Strtod("0.0000000000000100e-324"));
    CHECK(0.0 == Strtod("0.0000000000001000e-324"));
    CHECK(0.0 == Strtod("0.0000000000010000e-324"));
    CHECK(0.0 == Strtod("0.0000000000100000e-324"));
    CHECK(0.0 == Strtod("0.0000000001000000e-324"));
    CHECK(0.0 == Strtod("0.0000000010000000e-324"));
    CHECK(0.0 == Strtod("0.0000000100000000e-324"));
    CHECK(0.0 == Strtod("0.0000001000000000e-324"));
    CHECK(0.0 == Strtod("0.0000010000000000e-324"));
    CHECK(0.0 == Strtod("0.0000100000000000e-324"));
    CHECK(0.0 == Strtod("0.0001000000000000e-324"));
    CHECK(0.0 == Strtod("0.0010000000000000e-324"));
    CHECK(0.0 == Strtod("0.0100000000000000e-324"));
    CHECK(0.0 == Strtod("0.1000000000000000e-324"));
    CHECK(0.0 == Strtod("1.0000000000000000e-324"));
    CHECK(0.0 == Strtod("1e-324"));
}

TEST_CASE("Strtod - Special")
{
    CHECK( 0.0 == Strtod("0"));
    CHECK( 0.0 == Strtod("0.0000000000000000000000000000000"));
    CHECK(-0.0 == Strtod("-0"));
    //CHECK( 0.0 == Strtod("+0"));

    CheckStrtod(0.0f);
    CheckStrtod(-0.0f);
    CheckStrtod(std::numeric_limits<double>::infinity());
    CheckStrtod(-std::numeric_limits<double>::infinity());
    CheckStrtod(std::numeric_limits<double>::quiet_NaN());

    CHECK(std::isnan(Strtod("nan")));
    CHECK(std::isnan(Strtod("NaN")));
    CHECK(std::isnan(Strtod("nAn(_nananana123)")));
    CHECK(std::isnan(Strtod("nan(")));
    CHECK(std::isnan(Strtod("nan(xxx")));
    CHECK(std::isnan(Strtod("nan(xxx)")));

    CHECK(std::isinf(Strtod("Inf")));
    CHECK(std::isinf(Strtod("Infinity")));
    CHECK(std::isinf(Strtod("-INF")));
}

TEST_CASE("Strtod - Long input")
{
    CHECK(1280.0 == Strtod("128.000000000000000000000000000000000000000000000000000000000000"
                           "0000000000000000000000000000000000000000000000000000000000000e+1"));
    CHECK(1280.0 == Strtod("128.000000000000000000000000000000000000000000000000000000000000"
                           "00000000000000000000000000000000000000000000000000000000000000000e+1"));
    //CHECK(1280.0 == Strtod("128.0000000000000000000000000000000000000000000000000000000000000000"
    //                           "0000000000000000000000000000000000000000000000000000000000000000e+1"));
}

TEST_CASE("Strtod - Paxson, Kahan")
{
    //
    // V. Paxson and W. Kahan, "A Program for Testing IEEE Binary-Decimal Conversion", manuscript, May 1991,
    // ftp://ftp.ee.lbl.gov/testbase-report.ps.Z    (report)
    // ftp://ftp.ee.lbl.gov/testbase.tar.Z          (program)
    //

    //
    // Table 1:
    // Stress Inputs for Conversion to 53-bit Binary, < 1/2 ULP
    //

    CheckStrtod(5e+125); //, Strtod("5", 125));
    CheckStrtod(69e+267); //, Strtod("69", 267));
    CheckStrtod(999e-26); //, Strtod("999", -26));
    CheckStrtod(7861e-34); //, Strtod("7861", -34));
    CheckStrtod(75569e-254); //, Strtod("75569", -254));
    CheckStrtod(928609e-261); //, Strtod("928609", -261));
    CheckStrtod(9210917e+80); //, Strtod("9210917", 80));
    CheckStrtod(84863171e+114); //, Strtod("84863171", 114));
    CheckStrtod(653777767e+273); //, Strtod("653777767", 273));
    CheckStrtod(5232604057e-298); //, Strtod("5232604057", -298));
    CheckStrtod(27235667517e-109); //, Strtod("27235667517", -109));
    CheckStrtod(653532977297e-123); //, Strtod("653532977297", -123));
    CheckStrtod(3142213164987e-294); //, Strtod("3142213164987", -294));
    CheckStrtod(46202199371337e-72); //, Strtod("46202199371337", -72));
    CheckStrtod(231010996856685e-73); //, Strtod("231010996856685", -73));
    CheckStrtod(9324754620109615e+212); //, Strtod("9324754620109615", 212));
    CheckStrtod(78459735791271921e+49); //, Strtod("78459735791271921", 49));
    CheckStrtod(272104041512242479e+200);
    CheckStrtod(6802601037806061975e+198);
    CheckStrtod(20505426358836677347e-221);
    CheckStrtod(836168422905420598437e-234);
    CheckStrtod(4891559871276714924261e+222);

    //
    // Table 2:
    // Stress Inputs for Conversion to 53-bit Binary, > 1/2 ULP
    //

    CheckStrtod(9e-265); //, Strtod("9", -265));
    CheckStrtod(85e-37); //, Strtod("85", -37));
    CheckStrtod(623e+100); //, Strtod("623", 100));
    CheckStrtod(3571e+263); //, Strtod("3571", 263));
    CheckStrtod(81661e+153); //, Strtod("81661", 153));
    CheckStrtod(920657e-23); //, Strtod("920657", -23));
    CheckStrtod(4603285e-24); //, Strtod("4603285", -24));
    CheckStrtod(87575437e-309); //, Strtod("87575437", -309));
    CheckStrtod(245540327e+122); //, Strtod("245540327", 122));
    CheckStrtod(6138508175e+120); //, Strtod("6138508175", 120));
    CheckStrtod(83356057653e+193); //, Strtod("83356057653", 193));
    CheckStrtod(619534293513e+124); //, Strtod("619534293513", 124));
    CheckStrtod(2335141086879e+218); //, Strtod("2335141086879", 218));
    CheckStrtod(36167929443327e-159); //, Strtod("36167929443327", -159));
    CheckStrtod(609610927149051e-255); //, Strtod("609610927149051", -255));
    CheckStrtod(3743626360493413e-165); //, Strtod("3743626360493413", -165));
    CheckStrtod(94080055902682397e-242); //, Strtod("94080055902682397", -242));
    CheckStrtod(899810892172646163e+283);
    CheckStrtod(7120190517612959703e+120);
    CheckStrtod(25188282901709339043e-252);
    CheckStrtod(308984926168550152811e-52);
    CheckStrtod(6372891218502368041059e+064);

    //
    // Table 18:
    // Stress Inputs for Conversion to 56-bit Binary, < 1/2 ULP
    //

    CheckStrtod(7e-27); //, Strtod("7", -27));
    CheckStrtod(37e-29); //, Strtod("37", -29));
    CheckStrtod(743e-18); //, Strtod("743", -18));
    CheckStrtod(7861e-33); //, Strtod("7861", -33));
    CheckStrtod(46073e-30); //, Strtod("46073", -30));
    CheckStrtod(774497e-34); //, Strtod("774497", -34));
    CheckStrtod(8184513e-33); //, Strtod("8184513", -33));
    CheckStrtod(89842219e-28); //, Strtod("89842219", -28));
    CheckStrtod(449211095e-29); //, Strtod("449211095", -29));
    CheckStrtod(8128913627e-40); //, Strtod("8128913627", -40));
    CheckStrtod(87365670181e-18); //, Strtod("87365670181", -18));
    CheckStrtod(436828350905e-19); //, Strtod("436828350905", -19));
    CheckStrtod(5569902441849e-49); //, Strtod("5569902441849", -49));
    CheckStrtod(60101945175297e-32); //, Strtod("60101945175297", -32));
    CheckStrtod(754205928904091e-51); //, Strtod("754205928904091", -51));
    CheckStrtod(5930988018823113e-37); //, Strtod("5930988018823113", -37));
    CheckStrtod(51417459976130695e-27); //, Strtod("51417459976130695", -27));
    CheckStrtod(826224659167966417e-41);
    CheckStrtod(9612793100620708287e-57);
    CheckStrtod(93219542812847969081e-39);
    CheckStrtod(544579064588249633923e-48);
    CheckStrtod(4985301935905831716201e-48);

    //
    // Table 19:
    // Stress Inputs for Conversion to 56-bit Binary, > 1/2 ULP
    //

    CheckStrtod(9e+26); //, Strtod("9", 26));
    CheckStrtod(79e-8); //, Strtod("79", -8));
    CheckStrtod(393e+26); //, Strtod("393", 26));
    CheckStrtod(9171e-40); //, Strtod("9171", -40));
    CheckStrtod(56257e-16); //, Strtod("56257", -16));
    CheckStrtod(281285e-17); //, Strtod("281285", -17));
    CheckStrtod(4691113e-43); //, Strtod("4691113", -43));
    CheckStrtod(29994057e-15); //, Strtod("29994057", -15));
    CheckStrtod(834548641e-46); //, Strtod("834548641", -46));
    CheckStrtod(1058695771e-47); //, Strtod("1058695771", -47));
    CheckStrtod(87365670181e-18); //, Strtod("87365670181", -18));
    CheckStrtod(872580695561e-36); //, Strtod("872580695561", -36));
    CheckStrtod(6638060417081e-51); //, Strtod("6638060417081", -51));
    CheckStrtod(88473759402752e-52); //, Strtod("88473759402752", -52));
    CheckStrtod(412413848938563e-27); //, Strtod("412413848938563", -27));
    CheckStrtod(5592117679628511e-48); //, Strtod("5592117679628511", -48));
    CheckStrtod(83881765194427665e-50); //, Strtod("83881765194427665", -50));
    CheckStrtod(638632866154697279e-35);
    CheckStrtod(3624461315401357483e-53);
    CheckStrtod(75831386216699428651e-30);
    CheckStrtod(356645068918103229683e-42);
    CheckStrtod(7022835002724438581513e-33);
}

TEST_CASE("Strtod - Boundaries")
{
    //for (uint32_t e = 0; e < 2047; ++e)
    //{
    //    CheckStrtod(MakeDouble(e, 0));
    //    CheckStrtod(MakeDouble(e, 1));
    //    CheckStrtod(MakeDouble(e, 1ull << 51));
    //    CheckStrtod(MakeDouble(e, (1ull << 52) - 1));
    //    CheckStrtod(MakeDouble(e, ((1ull << 51) - 1) << 1));
    //}

    // Boundary cases. Boundaries themselves should round to even.
    //
    // 0x1FFFFFFFFFFFF * 2^3 = 72057594037927928
    //                   next: 72057594037927936
    //               boundary: 72057594037927932  should round up.
    CheckStrtod(72057594037927928e0);
    CheckStrtod(72057594037927936e0);
    CheckStrtod(72057594037927932e0);
    CheckStrtod(7205759403792793199999e-5);
    CheckStrtod(7205759403792793200001e-5);

    // 0x1FFFFFFFFFFFF * 2^10 = 9223372036854774784
    //                    next: 9223372036854775808
    //                boundary: 9223372036854775296 should round up.
    CheckStrtod(9223372036854774784e0);
    CheckStrtod(9223372036854775808e0);
    CheckStrtod(9223372036854775296e0);
    CheckStrtod(922337203685477529599999e-5);
    CheckStrtod(922337203685477529600001e-5);

    // 0x1FFFFFFFFFFFF * 2^50 = 10141204801825834086073718800384
    //                    next: 10141204801825835211973625643008
    //                boundary: 10141204801825834649023672221696 should round up.
    CheckStrtod(10141204801825834086073718800384e0);
    CheckStrtod(10141204801825835211973625643008e0);
    CheckStrtod(10141204801825834649023672221696e0);
    CheckStrtod(1014120480182583464902367222169599999e-5);
    CheckStrtod(1014120480182583464902367222169600001e-5);

    // 0x1FFFFFFFFFFFF * 2^99 = 5708990770823838890407843763683279797179383808
    //                    next: 5708990770823839524233143877797980545530986496
    //                boundary: 5708990770823839207320493820740630171355185152
    // The boundary should round up.
    CheckStrtod(5708990770823838890407843763683279797179383808e0);
    CheckStrtod(5708990770823839524233143877797980545530986496e0);
    CheckStrtod(5708990770823839207320493820740630171355185152e0);
    CheckStrtod(5708990770823839207320493820740630171355185151999e-3);
    CheckStrtod(5708990770823839207320493820740630171355185152001e-3);
}

#endif // 0

TEST_CASE("Strtod - syntax")
{
    using ryu::StrtodStatus;

    auto check = [](const std::string& input, StrtodStatus ec, int consumed = -1)
    {
        if (consumed < 0)
            consumed = static_cast<int>(input.size());

        const char* next = input.data();
        const char* last = input.data() + input.size();

        double value;
        const auto res = ryu::Strtod(next, last, value);

        return (res.status == ec) && (res.next == next + consumed);
    };

    CHECK(check("0", StrtodStatus::ok));
    CHECK(check("-0", StrtodStatus::ok));
    CHECK(check("123e65", StrtodStatus::ok));
    CHECK(check("0e+1", StrtodStatus::ok));
    CHECK(check("0e1", StrtodStatus::ok));
    CHECK(check("4", StrtodStatus::ok));
    CHECK(check("-0.0000000000000000000000000000001", StrtodStatus::ok));
    CHECK(check("20e1", StrtodStatus::ok));
    CHECK(check("-123", StrtodStatus::ok));
    CHECK(check("-1", StrtodStatus::ok));
    CHECK(check("1E22", StrtodStatus::ok));
    CHECK(check("1E-2", StrtodStatus::ok));
    CHECK(check("1E+2", StrtodStatus::ok));
    CHECK(check("123e45", StrtodStatus::ok));
    CHECK(check("123.456e78", StrtodStatus::ok));
    CHECK(check("1e-2", StrtodStatus::ok));
    CHECK(check("1e+2", StrtodStatus::ok));
    CHECK(check("123", StrtodStatus::ok));
    CHECK(check("123.456789", StrtodStatus::ok));
    CHECK(check("123.456e-789", StrtodStatus::ok));
    CHECK(check("-1e+9999", StrtodStatus::ok));
    CHECK(check("1.5e+9999", StrtodStatus::ok));
    CHECK(check("-123123e999990", StrtodStatus::ok));
    CHECK(check("123123e999999", StrtodStatus::ok));
    CHECK(check("123123e-1000000", StrtodStatus::ok)); // 0
    CHECK(check("123123e+1000000", StrtodStatus::ok)); // +inf
    CHECK(check("-123123123123123123123123123123", StrtodStatus::ok));
    CHECK(check("100000000000000000000", StrtodStatus::ok));
    CHECK(check("-237462374673276894279832749832423479823246327846", StrtodStatus::ok));

    CHECK(check("Infinity", StrtodStatus::ok, 8));
    CHECK(check("-Infinity", StrtodStatus::ok, 9));
    CHECK(check("NaN", StrtodStatus::ok, 3));
    CHECK(check("-NaN", StrtodStatus::ok, 4));

    CHECK(check("-1.0.", StrtodStatus::ok, 4));
    CHECK(check("0.1.2", StrtodStatus::ok, 3));
    CHECK(check("1 000.0", StrtodStatus::ok, 1));
    CHECK(check("1+2", StrtodStatus::ok, 1));
    CHECK(check("0x1", StrtodStatus::ok, 1));
    CHECK(check("0x42", StrtodStatus::ok, 1));
    CHECK(check("-123.123foo", StrtodStatus::ok, 8));
    CHECK(check("123\345", StrtodStatus::ok, 3));
    CHECK(check("1e1\345", StrtodStatus::ok, 3));
    CHECK(check("1.1e1\345", StrtodStatus::ok, 5));
    CHECK(check("0\345", StrtodStatus::ok, 1));
    CHECK(check("-1x", StrtodStatus::ok, 2));
    CHECK(check("1.2a-3", StrtodStatus::ok, 3));
    CHECK(check("1.8011670033376514H-308", StrtodStatus::ok, 18));

    CHECK(check("Infinity1234", StrtodStatus::ok, 8));
    CHECK(check("-Infinity1234", StrtodStatus::ok, 9));
    CHECK(check("NaN1234", StrtodStatus::ok, 3));
    CHECK(check("-NaN1234", StrtodStatus::ok, 4));

    CHECK(check("", StrtodStatus::invalid, 0));
    CHECK(check("-", StrtodStatus::invalid, 1));
    CHECK(check("++1234", StrtodStatus::invalid, 1));
    CHECK(check("+1", StrtodStatus::ok));
    CHECK(check("+Inf", StrtodStatus::ok));
    CHECK(check("+Infinity", StrtodStatus::ok));
    CHECK(check("+NaN", StrtodStatus::ok));
    CHECK(check("-01", StrtodStatus::ok));
    CHECK(check("-2.", StrtodStatus::ok));
    CHECK(check(".-1", StrtodStatus::invalid, 1));
    CHECK(check(".2e-3", StrtodStatus::ok));
    CHECK(check("0.e1", StrtodStatus::ok));
    CHECK(check("2.e+3", StrtodStatus::ok));
    CHECK(check("2.e-3", StrtodStatus::ok));
    CHECK(check("2.e3", StrtodStatus::ok));

    CHECK(check("Inf", StrtodStatus::ok));

    CHECK(check("-foo", StrtodStatus::invalid, 1));
    CHECK(check("- 1", StrtodStatus::invalid, 1));
    CHECK(check("-012", StrtodStatus::ok));
    CHECK(check("-.123", StrtodStatus::ok));
    CHECK(check("1.", StrtodStatus::ok));

    CHECK(check(".123", StrtodStatus::ok));
    CHECK(check("\357\274\221", StrtodStatus::invalid, 0));
    CHECK(check("012", StrtodStatus::ok));

    CHECK(check("+Infinity", StrtodStatus::ok));
    CHECK(check("+NaN", StrtodStatus::ok));
    CHECK(check("+Infinity1234", StrtodStatus::ok, 9));
    CHECK(check("+NaN1234", StrtodStatus::ok, 4));

    CHECK(check("123.000000456", StrtodStatus::ok));
    CHECK(check("0123.000000456", StrtodStatus::ok));
    CHECK(check("00000123.000000456", StrtodStatus::ok));
    CHECK(check("123.000000456", StrtodStatus::ok));
    CHECK(check("0123.000000456", StrtodStatus::ok));
    CHECK(check("00000123.000000456", StrtodStatus::ok));

#if 1
    CHECK(Strtod("123123e-1000000") == +0.0);
    CHECK(Strtod("123123e+1000000") == +std::numeric_limits<double>::infinity());
    CHECK(Strtod("-123123e-00000000000000000000000000000999999") == -0.0);
    CHECK(Strtod("-123123e+00000000000000000000000000000999999") == -std::numeric_limits<double>::infinity());
    CHECK(Strtod(".000000456") == .000000456);
    CHECK(Strtod("0.000000456") == 0.000000456);
    CHECK(Strtod("00000.000000456") == 00000.000000456);
    CHECK(Strtod(".000000456") == .000000456);
    CHECK(Strtod("0.000000456") == 0.000000456);
    CHECK(Strtod("00000.000000456") == 00000.000000456);
#endif

#if 1
    CHECK(check("0.3e+", StrtodStatus::ok, 3));
    CHECK(check("0.3e", StrtodStatus::ok, 3));
    CHECK(check("0e+", StrtodStatus::ok, 1));
    CHECK(check("0e", StrtodStatus::ok, 1));
    CHECK(check("0E+", StrtodStatus::ok, 1));
    CHECK(check("0E", StrtodStatus::ok, 1));
    CHECK(check("1.0e+", StrtodStatus::ok, 3));
    CHECK(check("1.0e-", StrtodStatus::ok, 3));
    CHECK(check("1.0e", StrtodStatus::ok, 3));
    CHECK(check("1eE2", StrtodStatus::ok, 1));
    CHECK(check("9.e+", StrtodStatus::ok, 2));
    CHECK(check("0e+-1", StrtodStatus::ok, 1));
    CHECK(check("1ea", StrtodStatus::ok, 1));
    CHECK(check("1e\345", StrtodStatus::ok, 1));
#else
    CHECK(check("0.3e+", StrtodStatus::invalid, 5));
    CHECK(check("0.3e", StrtodStatus::invalid, 4));
    CHECK(check("0e+", StrtodStatus::invalid, 3));
    CHECK(check("0e", StrtodStatus::invalid, 2));
    CHECK(check("0E+", StrtodStatus::invalid, 3));
    CHECK(check("0E", StrtodStatus::invalid, 2));
    CHECK(check("1.0e+", StrtodStatus::invalid, 5));
    CHECK(check("1.0e-", StrtodStatus::invalid, 5));
    CHECK(check("1.0e", StrtodStatus::invalid, 4));
    CHECK(check("1eE2", StrtodStatus::invalid, 2));
    CHECK(check("9.e+", StrtodStatus::invalid, 4));
    CHECK(check("0e+-1", StrtodStatus::invalid, 3));
    CHECK(check("1ea", StrtodStatus::invalid, 2));
    CHECK(check("1e\345", StrtodStatus::invalid, 2));
#endif
}

static double Strtod(const std::string& digits, int exponent)
{
    std::string input = digits;
    input += 'e';
    input += std::to_string(exponent);

    return Strtod(input);
}

#define CHECK_EQ(EXPECTED, ...) CHECK((EXPECTED) == (__VA_ARGS__))

TEST_CASE("Strtod - double_conversion - part 2")
{
    constexpr double Inf = std::numeric_limits<double>::infinity();

#if 1
    CHECK_EQ(0.0,        Strtod("0", 12345));
    //CHECK_EQ(0.0,        Strtod("", 1324));
    CHECK_EQ(0.0,        Strtod("000000000", 123));
    CHECK_EQ(0.0,        Strtod("2", -324));
    CHECK_EQ(4e-324,     Strtod("3", -324));
    CHECK_EQ(0.0,        Strtod("1", -325));
    CHECK_EQ(0.0,        Strtod("1", -325));
    CHECK_EQ(0.0,        Strtod("20000", -328));
    CHECK_EQ(40000e-328, Strtod("30000", -328));
    CHECK_EQ(0.0,        Strtod("10000", -329));
    CHECK_EQ(0.0,        Strtod("90000", -329));
    CHECK_EQ(0.0,        Strtod("000000001", -325));
    CHECK_EQ(0.0,        Strtod("000000001", -325));
    CHECK_EQ(0.0,        Strtod("0000000020000", -328));
    CHECK_EQ(40000e-328, Strtod("00000030000", -328));
    CHECK_EQ(0.0,        Strtod("0000000010000", -329));
    CHECK_EQ(0.0,        Strtod("0000000090000", -329));

    CHECK_EQ(Inf,                     Strtod("1", 309));
    CHECK_EQ(1e308,                   Strtod("1", 308));
    CHECK_EQ(1234e305,                Strtod("1234", 305));
    CHECK_EQ(1234e304,                Strtod("1234", 304));
    CHECK_EQ(Inf,                     Strtod("18", 307));
    CHECK_EQ(17e307,                  Strtod("17", 307));
    CHECK_EQ(Inf,                     Strtod("0000001", 309));
    CHECK_EQ(1e308,                   Strtod("00000001", 308));
    CHECK_EQ(1234e305,                Strtod("00000001234", 305));
    CHECK_EQ(1234e304,                Strtod("000000001234", 304));
    CHECK_EQ(Inf,                     Strtod("0000000018", 307));
    CHECK_EQ(17e307,                  Strtod("0000000017", 307));
    CHECK_EQ(Inf,                     Strtod("1000000", 303));
    CHECK_EQ(1e308,                   Strtod("100000", 303));
    CHECK_EQ(1234e305,                Strtod("123400000", 300));
    CHECK_EQ(1234e304,                Strtod("123400000", 299));
    CHECK_EQ(Inf,                     Strtod("180000000", 300));
    CHECK_EQ(17e307,                  Strtod("170000000", 300));
    CHECK_EQ(Inf,                     Strtod("00000001000000", 303));
    CHECK_EQ(1e308,                   Strtod("000000000000100000", 303));
    CHECK_EQ(1234e305,                Strtod("00000000123400000", 300));
    CHECK_EQ(1234e304,                Strtod("0000000123400000", 299));
    CHECK_EQ(Inf,                     Strtod("00000000180000000", 300));
    CHECK_EQ(17e307,                  Strtod("00000000170000000", 300));
    CHECK_EQ(1.7976931348623157E+308, Strtod("17976931348623157", 292));
    CHECK_EQ(1.7976931348623158E+308, Strtod("17976931348623158", 292));
    CHECK_EQ(Inf,                     Strtod("17976931348623159", 292));
#endif

    // The following number is the result of 89255.0/1e-22. Both floating-point
    // numbers can be accurately represented with doubles. However on Linux,x86
    // the floating-point stack is set to 80bits and the double-rounding
    // introduces an error.
    CHECK_EQ(89255e-22, Strtod("89255", -22));

    // Some random values.
    CHECK_EQ(358416272e-33, Strtod("358416272", -33));
    CHECK_EQ(104110013277974872254e-225, Strtod("104110013277974872254", -225));

    CHECK_EQ(123456789e108, Strtod("123456789", 108));
    CHECK_EQ(123456789e109, Strtod("123456789", 109));
    CHECK_EQ(123456789e110, Strtod("123456789", 110));
    CHECK_EQ(123456789e111, Strtod("123456789", 111));
    CHECK_EQ(123456789e112, Strtod("123456789", 112));
    CHECK_EQ(123456789e113, Strtod("123456789", 113));
    CHECK_EQ(123456789e114, Strtod("123456789", 114));
    CHECK_EQ(123456789e115, Strtod("123456789", 115));

#if 1
    CHECK_EQ(1234567890123456789012345e108, Strtod("1234567890123456789012345", 108));
    CHECK_EQ(1234567890123456789012345e109, Strtod("1234567890123456789012345", 109));
    CHECK_EQ(1234567890123456789012345e110, Strtod("1234567890123456789012345", 110));
    CHECK_EQ(1234567890123456789012345e111, Strtod("1234567890123456789012345", 111));
    CHECK_EQ(1234567890123456789012345e112, Strtod("1234567890123456789012345", 112));
    CHECK_EQ(1234567890123456789012345e113, Strtod("1234567890123456789012345", 113));
    CHECK_EQ(1234567890123456789012345e114, Strtod("1234567890123456789012345", 114));
    CHECK_EQ(1234567890123456789012345e115, Strtod("1234567890123456789012345", 115));

    CHECK_EQ(1234567890123456789052345e108, Strtod("1234567890123456789052345", 108));
    CHECK_EQ(1234567890123456789052345e109, Strtod("1234567890123456789052345", 109));
    CHECK_EQ(1234567890123456789052345e110, Strtod("1234567890123456789052345", 110));
    CHECK_EQ(1234567890123456789052345e111, Strtod("1234567890123456789052345", 111));
    CHECK_EQ(1234567890123456789052345e112, Strtod("1234567890123456789052345", 112));
    CHECK_EQ(1234567890123456789052345e113, Strtod("1234567890123456789052345", 113));
    CHECK_EQ(1234567890123456789052345e114, Strtod("1234567890123456789052345", 114));
    CHECK_EQ(1234567890123456789052345e115, Strtod("1234567890123456789052345", 115));

    CHECK_EQ(5.445618932859895e-255,
        Strtod("5445618932859895362967233318697132813618813095743952975"
               "4392982234069699615600475529427176366709107287468930197"
               "8628345413991790019316974825934906752493984055268219809"
               "5012176093045431437495773903922425632551857520884625114"
               "6241265881735209066709685420744388526014389929047617597"
               "0302268848374508109029268898695825171158085457567481507"
               "4162979705098246243690189880319928315307816832576838178"
               "2563074014542859888710209237525873301724479666744537857"
               "9026553346649664045621387124193095870305991178772256504"
               "4368663670643970181259143319016472430928902201239474588"
               "1392338901353291306607057623202353588698746085415097902"
               "6640064319118728664842287477491068264828851624402189317"
               "2769161449825765517353755844373640588822904791244190695"
               "2998382932630754670573838138825217065450843010498555058"
               "88186560731", -1035));

    // Boundary cases. Boundaries themselves should round to even.
    //
    // 0x1FFFFFFFFFFFF * 2^3 = 72057594037927928
    //                   next: 72057594037927936
    //               boundary: 72057594037927932  should round up.
    CHECK_EQ(72057594037927928.0, Strtod("72057594037927928", 0));
    CHECK_EQ(72057594037927936.0, Strtod("72057594037927936", 0));
    CHECK_EQ(72057594037927936.0, Strtod("72057594037927932", 0));
    CHECK_EQ(72057594037927928.0, Strtod("7205759403792793199999", -5));
    CHECK_EQ(72057594037927936.0, Strtod("7205759403792793200001", -5));

    // 0x1FFFFFFFFFFFF * 2^10 = 9223372036854774784
    //                    next: 9223372036854775808
    //                boundary: 9223372036854775296 should round up.
    CHECK_EQ(9223372036854774784.0, Strtod("9223372036854774784", 0));
    CHECK_EQ(9223372036854775808.0, Strtod("9223372036854775808", 0));
    CHECK_EQ(9223372036854775808.0, Strtod("9223372036854775296", 0));
    CHECK_EQ(9223372036854774784.0, Strtod("922337203685477529599999", -5));
    CHECK_EQ(9223372036854775808.0, Strtod("922337203685477529600001", -5));

    // 0x1FFFFFFFFFFFF * 2^50 = 10141204801825834086073718800384
    //                    next: 10141204801825835211973625643008
    //                boundary: 10141204801825834649023672221696 should round up.
    CHECK_EQ(10141204801825834086073718800384.0, Strtod("10141204801825834086073718800384", 0));
    CHECK_EQ(10141204801825835211973625643008.0, Strtod("10141204801825835211973625643008", 0));
    CHECK_EQ(10141204801825835211973625643008.0, Strtod("10141204801825834649023672221696", 0));
    CHECK_EQ(10141204801825834086073718800384.0, Strtod("1014120480182583464902367222169599999", -5));
    CHECK_EQ(10141204801825835211973625643008.0, Strtod("1014120480182583464902367222169600001", -5));

    // 0x1FFFFFFFFFFFF * 2^99 = 5708990770823838890407843763683279797179383808
    //                    next: 5708990770823839524233143877797980545530986496
    //                boundary: 5708990770823839207320493820740630171355185152
    // The boundary should round up.
    CHECK_EQ(5708990770823838890407843763683279797179383808.0, Strtod("5708990770823838890407843763683279797179383808", 0));
    CHECK_EQ(5708990770823839524233143877797980545530986496.0, Strtod("5708990770823839524233143877797980545530986496", 0));
    CHECK_EQ(5708990770823839524233143877797980545530986496.0, Strtod("5708990770823839207320493820740630171355185152", 0));
    CHECK_EQ(5708990770823838890407843763683279797179383808.0, Strtod("5708990770823839207320493820740630171355185151999", -3));
    CHECK_EQ(5708990770823839524233143877797980545530986496.0, Strtod("5708990770823839207320493820740630171355185152001", -3));

    // The following test-cases got some public attention in early 2011 when they
    // sent Java and PHP into an infinite loop.
    CHECK_EQ(2.225073858507201e-308, Strtod("22250738585072011", -324));
    CHECK_EQ(2.22507385850720138309e-308,
        Strtod("22250738585072011360574097967091319759348195463516456480"
               "23426109724822222021076945516529523908135087914149158913"
               "03962110687008643869459464552765720740782062174337998814"
               "10632673292535522868813721490129811224514518898490572223"
               "07285255133155755015914397476397983411801999323962548289"
               "01710708185069063066665599493827577257201576306269066333"
               "26475653000092458883164330377797918696120494973903778297"
               "04905051080609940730262937128958950003583799967207254304"
               "36028407889577179615094551674824347103070260914462157228"
               "98802581825451803257070188608721131280795122334262883686"
               "22321503775666622503982534335974568884423900265498198385"
               "48794829220689472168983109969836584681402285424333066033"
               "98508864458040010349339704275671864433837704860378616227"
               "71738545623065874679014086723327636718751", -1076));
#endif
}

#if 1
TEST_CASE("Strtod - Exponents")
{
    constexpr double Inf = std::numeric_limits<double>::infinity();

    CHECK_EQ(0.0, Strtod("0e+0"));
    CHECK_EQ(0.0, Strtod("0e-0"));
    CHECK_EQ(0.0, Strtod("0e+100"));
    CHECK_EQ(0.0, Strtod("0e-100"));
    CHECK_EQ(0.0, Strtod("0e+2147483647"));
    CHECK_EQ(0.0, Strtod("0e-2147483647"));
    CHECK_EQ(0.0, Strtod("0.0e+2147483647"));
    CHECK_EQ(0.0, Strtod("0.0e-2147483647"));
    CHECK_EQ(0.0, Strtod("0.00000000000000000000000000000000000000000000000000000000000000000000e+2147483647"));
    CHECK_EQ(0.0, Strtod("0.00000000000000000000000000000000000000000000000000000000000000000000e-2147483647"));
    CHECK_EQ(0.0, Strtod("0.00000000000000000000000000000000000000000000000000000000000000000001e-2147483647"));
    CHECK_EQ(0.0, Strtod("1.00000000000000000000000000000000000000000000000000000000000000000000e-2147483647"));
    CHECK_EQ(0.0, Strtod("0e-2147483648"));
    CHECK_EQ(0.0, Strtod("1e-2147483649"));
    CHECK_EQ(0.0, Strtod("1e-2147483648"));
    CHECK_EQ(0.0, Strtod("1e-2147483647"));
    CHECK_EQ(0.0, Strtod("1e-1000"));
    CHECK_EQ(0.0, Strtod("1e-100000"));
    CHECK_EQ(0.0, Strtod("1e-99999999")); // 1e-99999999
    CHECK_EQ(0.0, Strtod("1e-100000000")); // 1e-Inf

    CHECK_EQ(Inf, Strtod("1e+2147483647"));
    CHECK_EQ(Inf, Strtod("1e+2147483648"));
    CHECK_EQ(Inf, Strtod("0.00000000000000000000000000000000000000000000000000000000000000000001e+2147483647"));
    CHECK_EQ(Inf, Strtod("1e+99999999")); // 1e+99999999
    CHECK_EQ(Inf, Strtod("1e+100000000")); // 1e+Inf

    CHECK_EQ(1.0, Strtod( "0.1e+0000000000000000000000000000000000000000000000000000000000000000000000000001"));
    CHECK_EQ(1.0, Strtod( "1.0e+0000000000000000000000000000000000000000000000000000000000000000000000000000"));
    CHECK_EQ(1.0, Strtod("10.0e-0000000000000000000000000000000000000000000000000000000000000000000000000001"));
}
#endif

#if 1
TEST_CASE("Strtod - Boundaries - part 2")
{
    constexpr double Inf = std::numeric_limits<double>::infinity();

    // 9007199254740991 * 2^-1074 = (2^53 - 1) * 2^-1074
    CHECK_EQ(4.450147717014402272e-308,
        Strtod("4.450147717014402272114819593418263951869639092703291296046852219449644444042153"
               "89103305904781627017582829831782607924221374017287738918929105531441481564124348"
               "67599762821265346585071045737627442980259622449029037796981144446145705102663115"
               "10031828794952795966823603998647925096578034214163701381261333311989876551545144"
               "03152612538132666529513060001849177663286607555958373922409899478075565940981010"
               "21612198814605258742579179000071675999344145086087205681577915435923018910334964"
               "86942061405218289243144579760516365090360651414037721744226256159024466852576737"
               "24464300755133324500796506867194913776884780053099639677097589658441378944337966"
               "21993967316936280457084866613206797017728916080020698679408551343728867675409720"
               "757232455434770912461317493580281734466552734375e-308"));
    // 9007199254740990 * 2^-1074
    CHECK_EQ(4.450147717014401778e-308,
        Strtod("4.450147717014401778049173752171719775300846224481918930987049605124880018456471"
               "39035755177760751831052846195619008686241717547743167145836439860405887584484471"
               "19639655002484083577939142623582164522087943959208000909794783876158397872163051"
               "22622675229968408654350206725478309956546318828765627255022767720818849892988457"
               "26333908582101604036318532842699932130356061901518261174396928478121372742040102"
               "17446565569357687263889031732270082446958029584739170416643195242132750803227473"
               "16608838720742955671061336566907126801014814608027120593609275183716632624844904"
               "31985250929886016737037234388448352929102742708402644340627409931664203093081360"
               "70794835812045179006047003875039546061891526346421705014598610179523165038319441"
               "51446491086954182492263498716056346893310546875e-308"));
    // half way between the two numbers above.
    // round to nearest even.
    CHECK_EQ(4.450147717014401778e-308,
        Strtod("4.450147717014402025081996672794991863585242658592605113516950912287262231249312"
               "64069530541271189424317838013700808305231545782515453032382772695923684574304409"
               "93619708911874715081505094180604803751173783204118519353387964161152051487413083"
               "16327252012460602310586905362063117526562176521464664318142050516404363222266800"
               "64743260560117135282915796422274554896821334728738317548403413978098469341510556"
               "19529382191981473003234105366170879223151087335413188049110555339027884856781219"
               "01775450062980622457102958163711745945687733011032421168917765671370549738710820"
               "78224775842509670618916870627821633352993761380751142008862499795052791018709663"
               "46394401564490729731565935244123171539810221213221201847003580761626016356864581"
               "1358486831521563686919762403704226016998291015625e-308"));
    CHECK_EQ(4.450147717014401778e-308,
        Strtod("4.450147717014402025081996672794991863585242658592605113516950912287262231249312"
               "64069530541271189424317838013700808305231545782515453032382772695923684574304409"
               "93619708911874715081505094180604803751173783204118519353387964161152051487413083"
               "16327252012460602310586905362063117526562176521464664318142050516404363222266800"
               "64743260560117135282915796422274554896821334728738317548403413978098469341510556"
               "19529382191981473003234105366170879223151087335413188049110555339027884856781219"
               "01775450062980622457102958163711745945687733011032421168917765671370549738710820"
               "78224775842509670618916870627821633352993761380751142008862499795052791018709663"
               "46394401564490729731565935244123171539810221213221201847003580761626016356864581"
               "13584868315215636869197624037042260169982910156250000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000e-308"));
    // ... round up
#if 0
    // BUG in MSVC's strtod...
    CHECK_EQ(4.450147717014402272e-308,
        Strtod("4.450147717014402025081996672794991863585242658592605113516950912287262231249312"
               "64069530541271189424317838013700808305231545782515453032382772695923684574304409"
               "93619708911874715081505094180604803751173783204118519353387964161152051487413083"
               "16327252012460602310586905362063117526562176521464664318142050516404363222266800"
               "64743260560117135282915796422274554896821334728738317548403413978098469341510556"
               "19529382191981473003234105366170879223151087335413188049110555339027884856781219"
               "01775450062980622457102958163711745945687733011032421168917765671370549738710820"
               "78224775842509670618916870627821633352993761380751142008862499795052791018709663"
               "46394401564490729731565935244123171539810221213221201847003580761626016356864581"
               "13584868315215636869197624037042260169982910156250000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000001e-308"));
#endif
    // ... round down
    CHECK_EQ(4.450147717014401778e-308,
        Strtod("4.450147717014402025081996672794991863585242658592605113516950912287262231249312"
               "64069530541271189424317838013700808305231545782515453032382772695923684574304409"
               "93619708911874715081505094180604803751173783204118519353387964161152051487413083"
               "16327252012460602310586905362063117526562176521464664318142050516404363222266800"
               "64743260560117135282915796422274554896821334728738317548403413978098469341510556"
               "19529382191981473003234105366170879223151087335413188049110555339027884856781219"
               "01775450062980622457102958163711745945687733011032421168917765671370549738710820"
               "78224775842509670618916870627821633352993761380751142008862499795052791018709663"
               "46394401564490729731565935244123171539810221213221201847003580761626016356864581"
               "13584868315215636869197624037042260169982910156249999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999e-308"));

    // 9007199254740989 * 2^-1074
    CHECK_EQ(4.450147717014401284e-308,
        Strtod("4.450147717014401283983527910925175598732053356260546565927246990800115592870788"
               "88968204450739876644522862559455409448262061078198595372743774189370293604844593"
               "71679547183702820570807239509536886063916265469386964022608423306171090641662987"
               "35213521664984021341876809452308694816514603443367553128784202129647823234431770"
               "49515204626070541543124005683550686597425516247078148426383957478167179543099194"
               "13280932324110115785198884464468488894571914083391135151708475048342482696119981"
               "46275616036267622098978093373297888511668977802016519442992294208408798397113071"
               "39506201104638708973277961909701792081320705363705649004157230204887027241824755"
               "19595704307154077555009141136872295106054136612822711349788669015317462401229162"
               "271697366304312737383952480740845203399658203125e-308"));

    // min denormal = 2^-1074
    CHECK_EQ(4.940656458412465442e-324,
        Strtod("4.940656458412465441765687928682213723650598026143247644255856825006755072702087"
               "51865299836361635992379796564695445717730926656710355939796398774796010781878126"
               "30071319031140452784581716784898210368871863605699873072305000638740915356498438"
               "73124733972731696151400317153853980741262385655911710266585566867681870395603106"
               "24931945271591492455329305456544401127480129709999541931989409080416563324524757"
               "14786901472678015935523861155013480352649347201937902681071074917033322268447533"
               "35720832431936092382893458368060106011506169809753078342277318329247904982524730"
               "77637592724787465608477820373446969953364701797267771758512566055119913150489110"
               "14510378627381672509558373897335989936648099411642057026370902792427675445652290"
               "87538682506419718265533447265625e-324"));
    // 2 * 2^-1074
    CHECK_EQ(9.881312916824930884e-324,
        Strtod("9.881312916824930883531375857364427447301196052286495288511713650013510145404175"
               "03730599672723271984759593129390891435461853313420711879592797549592021563756252"
               "60142638062280905569163433569796420737743727211399746144610001277481830712996877"
               "46249467945463392302800634307707961482524771311823420533171133735363740791206212"
               "49863890543182984910658610913088802254960259419999083863978818160833126649049514"
               "29573802945356031871047722310026960705298694403875805362142149834066644536895066"
               "71441664863872184765786916736120212023012339619506156684554636658495809965049461"
               "55275185449574931216955640746893939906729403594535543517025132110239826300978220"
               "29020757254763345019116747794671979873296198823284114052741805584855350891304581"
               "7507736501283943653106689453125e-324"));
    // half-way between the two smallest (subnormal) numbers: (1 * 2^-1074 + 2 * 2^-1074) / 2
    // round to nearest even
    CHECK_EQ(9.881312916824930884e-324,
        Strtod("7.410984687618698162648531893023320585475897039214871466383785237510132609053131"
               "27797949754542453988569694847043168576596389985065533909694598162194016172817189"
               "45106978546710679176872575177347315553307795408549809608457500958111373034747658"
               "09687100959097544227100475730780971111893578483867565399878350301522805593404659"
               "37397917907387238682993958184816601691220194564999312897984113620624844986787135"
               "72180352209017023903285791732520220528974020802906854021606612375549983402671300"
               "03581248647904138574340187552090159017259254714629617513415977493871857473787096"
               "16456389087181198412716730560170454930047052695901657637768849082679869725733665"
               "21765567941072508764337560846003984904972149117463085539556354188641513168478436"
               "313080237596295773983001708984375e-324"));
    // round up
    CHECK_EQ(9.881312916824930884e-324,
        Strtod("7.410984687618698162648531893023320585475897039214871466383785237510132609053131"
               "27797949754542453988569694847043168576596389985065533909694598162194016172817189"
               "45106978546710679176872575177347315553307795408549809608457500958111373034747658"
               "09687100959097544227100475730780971111893578483867565399878350301522805593404659"
               "37397917907387238682993958184816601691220194564999312897984113620624844986787135"
               "72180352209017023903285791732520220528974020802906854021606612375549983402671300"
               "03581248647904138574340187552090159017259254714629617513415977493871857473787096"
               "16456389087181198412716730560170454930047052695901657637768849082679869725733665"
               "21765567941072508764337560846003984904972149117463085539556354188641513168478436"
               "31308023759629577398300170898437500000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000001e-324"));
    // round down
#if 0
    // BUG in MSVC's strtod...
    CHECK_EQ(4.940656458412465442e-324,
        Strtod("7.410984687618698162648531893023320585475897039214871466383785237510132609053131"
               "27797949754542453988569694847043168576596389985065533909694598162194016172817189"
               "45106978546710679176872575177347315553307795408549809608457500958111373034747658"
               "09687100959097544227100475730780971111893578483867565399878350301522805593404659"
               "37397917907387238682993958184816601691220194564999312897984113620624844986787135"
               "72180352209017023903285791732520220528974020802906854021606612375549983402671300"
               "03581248647904138574340187552090159017259254714629617513415977493871857473787096"
               "16456389087181198412716730560170454930047052695901657637768849082679869725733665"
               "21765567941072508764337560846003984904972149117463085539556354188641513168478436"
               "31308023759629577398300170898437499999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999e-324"));
#endif

    // 9007199254740991 * 2^971 (max normal)
    CHECK_EQ(1.797693134862315708e+308,
        Strtod("1.797693134862315708145274237317043567980705675258449965989174768031572607800285"
               "38760589558632766878171540458953514382464234321326889464182768467546703537516986"
               "04991057655128207624549009038932894407586850845513394230458323690322294816580855"
               "9332123348274797826204144723168738177180919299881250404026184124858368e+308"));
    // 9007199254740992 * 2^971 ("infinity")
    CHECK_EQ(Inf,
        Strtod("1.797693134862315907729305190789024733617976978942306572734300811577326758055009"
               "63132708477322407536021120113879871393357658789768814416622492847430639474124377"
               "76789342486548527630221960124609411945308295208500576883815068234246288147391311"
               "0540827237163350510684586298239947245938479716304835356329624224137216e+308"));
    // half way between max-normal and infinity
    // should round to infinity in nearest-even mode.
    CHECK_EQ(Inf,
        Strtod("1.797693134862315807937289714053034150799341327100378269361737789804449682927647"
               "50946649017977587207096330286416692887910946555547851940402630657488671505820681"
               "90890200070838367627385484581771153176447573027006985557136695962284291481986083"
               "49364752927190741684443655107043427115596995080930428801779041744977920000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000e+308"));
    // ...round down
    CHECK_EQ(1.797693134862315708e+308,
        Strtod("1.797693134862315807937289714053034150799341327100378269361737789804449682927647"
               "50946649017977587207096330286416692887910946555547851940402630657488671505820681"
               "90890200070838367627385484581771153176447573027006985557136695962284291481986083"
               "49364752927190741684443655107043427115596995080930428801779041744977919999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999"
               "99999999999999999999999999999999999999999999999999999999999999999999999999999999e+308"));
    // ...round up
    CHECK_EQ(Inf,
        Strtod("1.797693134862315807937289714053034150799341327100378269361737789804449682927647"
               "50946649017977587207096330286416692887910946555547851940402630657488671505820681"
               "90890200070838367627385484581771153176447573027006985557136695962284291481986083"
               "49364752927190741684443655107043427115596995080930428801779041744977920000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000000000001e+308"));

    CHECK_EQ(2.225073858507202371e-308,
        Strtod("2.22507385850720212418870147920222032907240528279439037814303133837435107319244"
               "1946867544064325638818513821882185024380699999477330130056498841077919287413419"
               "2929720097048195199306799329096904278406473168204156592672863293363047467012331"
               "6852983422152744517260835859654566319282835244787787799894310779783833699159288"
               "5945552137141811284582511455843192230798975043950868594124572308917389461693683"
               "7232119137365897797772328669884035639025104444303545739673370658398105542045669"
               "3824658413747607155981176573877626747665912387199931904006317334709003012790188"
               "1752034471902500280612777779167983910905785840064647159438105114891542827750411"
               "7468219413395246668250343130618158782937900420539237507208336669324158000275839"
               "1118854188641513168478436313080237596295773983001708984375e-308"));
}
#endif

TEST_CASE("Strtod - Integers")
{
    constexpr double Inf = std::numeric_limits<double>::infinity();

    CHECK_EQ(  0.0, Strtod("0"));
    CHECK_EQ( -0.0, Strtod("-0"));
    CHECK_EQ(  1.0, Strtod("1"));
    CHECK_EQ( 12.0, Strtod("12"));
    CHECK_EQ( -1.0, Strtod("-1"));
    CHECK_EQ(-12.0, Strtod("-12"));

    CHECK_EQ(9.0, Strtod("9"));
    CHECK_EQ(99.0, Strtod("99"));
    CHECK_EQ(999.0, Strtod("999"));
    CHECK_EQ(9999.0, Strtod("9999"));
    CHECK_EQ(99999.0, Strtod("99999"));
    CHECK_EQ(999999.0, Strtod("999999"));
    CHECK_EQ(9999999.0, Strtod("9999999"));
    CHECK_EQ(99999999.0, Strtod("99999999"));
    CHECK_EQ(999999999.0, Strtod("999999999"));
    CHECK_EQ(9999999999.0, Strtod("9999999999"));
    CHECK_EQ(99999999999.0, Strtod("99999999999"));
    CHECK_EQ(999999999999.0, Strtod("999999999999"));
    CHECK_EQ(9999999999999.0, Strtod("9999999999999"));
    CHECK_EQ(99999999999999.0, Strtod("99999999999999"));
    CHECK_EQ(999999999999999.0, Strtod("999999999999999"));
    CHECK_EQ(9999999999999999.0, Strtod("9999999999999999"));
    CHECK_EQ(99999999999999999.0, Strtod("99999999999999999"));
    CHECK_EQ(999999999999999999.0, Strtod("999999999999999999"));
    CHECK_EQ(9999999999999999999.0, Strtod("9999999999999999999"));
    CHECK_EQ(99999999999999999999.0, Strtod("99999999999999999999"));

    CHECK_EQ(-9.0, Strtod("-9"));
    CHECK_EQ(-99.0, Strtod("-99"));
    CHECK_EQ(-999.0, Strtod("-999"));
    CHECK_EQ(-9999.0, Strtod("-9999"));
    CHECK_EQ(-99999.0, Strtod("-99999"));
    CHECK_EQ(-999999.0, Strtod("-999999"));
    CHECK_EQ(-9999999.0, Strtod("-9999999"));
    CHECK_EQ(-99999999.0, Strtod("-99999999"));
    CHECK_EQ(-999999999.0, Strtod("-999999999"));
    CHECK_EQ(-9999999999.0, Strtod("-9999999999"));
    CHECK_EQ(-99999999999.0, Strtod("-99999999999"));
    CHECK_EQ(-999999999999.0, Strtod("-999999999999"));
    CHECK_EQ(-9999999999999.0, Strtod("-9999999999999"));
    CHECK_EQ(-99999999999999.0, Strtod("-99999999999999"));
    CHECK_EQ(-999999999999999.0, Strtod("-999999999999999"));
    CHECK_EQ(-9999999999999999.0, Strtod("-9999999999999999"));
    CHECK_EQ(-99999999999999999.0, Strtod("-99999999999999999"));
    CHECK_EQ(-999999999999999999.0, Strtod("-999999999999999999"));
    CHECK_EQ(-9999999999999999999.0, Strtod("-9999999999999999999"));
    CHECK_EQ(-99999999999999999999.0, Strtod("-99999999999999999999"));

    CHECK_EQ( 2147483647.0, Strtod("2147483647")); // 2^31 - 1
    CHECK_EQ( 2147483648.0, Strtod("2147483648"));
    CHECK_EQ(-2147483647.0, Strtod("-2147483647"));
    CHECK_EQ(-2147483648.0, Strtod("-2147483648"));
    CHECK_EQ(-2147483649.0, Strtod("-2147483649"));
    CHECK_EQ( 4294967295.0, Strtod("4294967295")); // 2^32 - 1
    CHECK_EQ( 4294967296.0, Strtod("4294967296"));
    CHECK_EQ(-4294967295.0, Strtod("-4294967295"));
    CHECK_EQ(-4294967296.0, Strtod("-4294967296"));
    CHECK_EQ(-4294967297.0, Strtod("-4294967297"));
    CHECK_EQ( 9223372036854775807.0, Strtod("9223372036854775807")); // 2^63 - 1
    CHECK_EQ( 9223372036854775808.0, Strtod("9223372036854775808"));
    CHECK_EQ(-9223372036854775807.0, Strtod("-9223372036854775807"));
    CHECK_EQ(-9223372036854775808.0, Strtod("-9223372036854775808"));
    CHECK_EQ(-9223372036854775809.0, Strtod("-9223372036854775809"));
    CHECK_EQ( 18446744073709551615.0, Strtod("18446744073709551615")); // 2^64 - 1
    CHECK_EQ( 18446744073709551616.0, Strtod("18446744073709551616"));
    CHECK_EQ(-18446744073709551615.0, Strtod("-18446744073709551615"));
    CHECK_EQ(-18446744073709551616.0, Strtod("-18446744073709551616"));
    CHECK_EQ(-18446744073709551617.0, Strtod("-18446744073709551617"));

    // 10^799
    CHECK_EQ(Inf,
        Strtod("10000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000.0"));
}

TEST_CASE("Strtod - Regression")
{
    constexpr double Max = std::numeric_limits<double>::max();
    constexpr double Min = std::numeric_limits<double>::denorm_min();
    constexpr double MinNormal = std::numeric_limits<double>::min();

    CHECK_EQ(10000000000000000001e+19, Strtod("10000000000000000001", 19));

    CHECK_EQ(0.0, Strtod("0.0000"));
    CHECK_EQ(-0.0, Strtod("-0.0000"));

    CHECK_EQ(10000000000000000009e+0, Strtod("10000000000000000009e+0"));
    CHECK_EQ(10000000000000000009e+1, Strtod("10000000000000000009e+1"));
    CHECK_EQ(10000000000000000009e+2, Strtod("10000000000000000009e+2"));
    CHECK_EQ(10000000000000000009e+3, Strtod("10000000000000000009e+3"));
    CHECK_EQ(10000000000000000009e+4, Strtod("10000000000000000009e+4"));
    CHECK_EQ(10000000000000000009e+5, Strtod("10000000000000000009e+5"));
    CHECK_EQ(10000000000000000009e+6, Strtod("10000000000000000009e+6"));
    CHECK_EQ(10000000000000000009e+7, Strtod("10000000000000000009e+7"));
    CHECK_EQ(10000000000000000009e+8, Strtod("10000000000000000009e+8"));
    CHECK_EQ(10000000000000000009e+9, Strtod("10000000000000000009e+9"));
    CHECK_EQ(10000000000000000009e+10, Strtod("10000000000000000009e+10"));
    CHECK_EQ(10000000000000000009e+11, Strtod("10000000000000000009e+11"));
    CHECK_EQ(10000000000000000009e+12, Strtod("10000000000000000009e+12"));
    CHECK_EQ(10000000000000000009e+13, Strtod("10000000000000000009e+13"));
    CHECK_EQ(10000000000000000009e+14, Strtod("10000000000000000009e+14"));
    CHECK_EQ(10000000000000000009e+15, Strtod("10000000000000000009e+15"));
    CHECK_EQ(10000000000000000009e+16, Strtod("10000000000000000009e+16"));

    CHECK_EQ(1000000000000000000.0000000000000000001, Strtod("1000000000000000000.0000000000000000001"));

    CHECK_EQ(59.79470570797252226166574973080902316556696507444245101698,
        Strtod("59.79470570797252226166574973080902316556696507444245101698"));

    CHECK_EQ(0.0, Strtod("1e-324"));
    CHECK_EQ(0.0, Strtod("2e-324"));
    CHECK_EQ(3e-324, Strtod("3e-324"));
    CHECK_EQ(4e-324, Strtod("4e-324"));
    CHECK_EQ(5e-324, Strtod("5e-324")); // min denormal

    CHECK_EQ(4.9406564584124653e-324, Strtod("4.9406564584124653e-324"));
    CHECK_EQ(4.9406564584124654e-324, Strtod("4.9406564584124654e-324"));
    CHECK_EQ(4.9406564584124655e-324, Strtod("4.9406564584124655e-324"));
    CHECK_EQ(4.94065645841246539999999999999999999999999999999999999999999999999999999999e-324,
        Strtod("4.94065645841246539999999999999999999999999999999999999999999999999999999999e-324"));
    CHECK_EQ(4.94065645841246540000000000000000000000000000000000000000000000000000000001e-324,
        Strtod("4.94065645841246540000000000000000000000000000000000000000000000000000000001e-324"));

    CHECK_EQ(0.0, Strtod("2.4703282292062327e-324"));
    CHECK_EQ(2.4703282292062328e-324, Strtod("2.4703282292062328e-324"));
    CHECK_EQ(2.48e-324, Strtod("2.48e-324"));
    CHECK_EQ(2.5e-324, Strtod("2.5e-324"));
    CHECK_EQ(2.500000000000000000000000000000000000000000000000000000000000000000000000001e-324,
        Strtod("2.500000000000000000000000000000000000000000000000000000000000000000000000001e-324"));
    CHECK_EQ(2.225073858507201e-308, Strtod("2.225073858507201e-308")); // max denormal
    CHECK_EQ(2.2250738585072014e-308, Strtod("2.2250738585072014e-308")); // min normal
    CHECK_EQ(1.7976931348623157e+308, Strtod("1.7976931348623157e+308")); // max normal
    CHECK_EQ(1.7976931348623156999999999999999999999999999999999999999999999999999e+308,
        Strtod("1.7976931348623156999999999999999999999999999999999999999999999999999e+308"));
    CHECK_EQ(1.7976931348623157000000000000000000000000000000000000000000000000001e+308,
        Strtod("1.7976931348623157000000000000000000000000000000000000000000000000001e+308"));
    CHECK_EQ(1e-323, Strtod("1e-323"));
    CHECK_EQ(2e-323, Strtod("2e-323"));
    CHECK_EQ(3e-323, Strtod("3e-323"));
    CHECK_EQ(4e-323, Strtod("4e-323"));
    CHECK_EQ(1.7976931348623158e+308, Strtod("1.7976931348623158e+308"));
    CHECK_EQ(Max,
        Strtod("17976931348623157081452742373170435679807056752584499659891747680315726"
               "07800285387605895586327668781715404589535143824642343213268894641827684"
               "67546703537516986049910576551282076245490090389328944075868508455133942"
               "30458323690322294816580855933212334827479782620414472316873817718091929"
               "9881250404026184124858368"));
    CHECK_EQ(Min,
        Strtod("0.0000000000000000000000000000000000000000000000000000000000"
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
               "37090279242767544565229087538682506419718265533447265625"));
    CHECK_EQ(2.4354608055603473e+307,
        Strtod("243546080556034731077856379609316893158278902575447060151047"
               "212703405344938119816206067372775299130836050315842578309818"
               "316450894337978612745889730079163798234256495613858256849283"
               "467066859489192118352020514036083287319232435355752493038825"
               "828481044358810649108367633313557305310641892225870327827273"
               "41408256.000000"));
    CHECK_EQ(2.2250738585072011e-308, Strtod("2.2250738585072011e-308"));
    // 2^-1075
    CHECK_EQ(0.0,
        Strtod("2.4703282292062327208828439643411068618252990130716238221279"
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
               "437693412532098591327667236328125e-324"));
    CHECK_EQ(0.0,
        Strtod("0.00000000000000000000000000000000000000000000000000000000000"
               "0000000000000000000000000000000000000000000000000000000000000"
               "0000000000000000000000000000000000000000000000000000000000000"
               "0000000000000000000000000000000000000000000000000000000000000"
               "0000000000000000000000000000000000000000000000000000000000000"
               "0000000000000000000024703282292062327208828439643411068618252"
               "9901307162382212792841250337753635104375932649918180817996189"
               "8982823477228588654633283551779698981993873980053909390631503"
               "5659515570226392290858392449105184435931802849936536152500319"
               "3704576782492193656236698636584807570015857692699037063119282"
               "7955855133292783433840935197801553124659726357957462276646527"
               "2827220056374006485499977096599470454020828166226237857393450"
               "7363390079677619305775067401763246736009689513405355374585166"
               "6113422376667860416215968046191446729184030053005753084904876"
               "5391711386591646239524912623653881879636239373280423891018672"
               "3484976682350898633885879256283027559956575244555072551893136"
               "9083625477918694866799496832404970582102851318545139621383772"
               "2826145437693412532098591327667236328125"));
    CHECK_EQ(5e-324,
        Strtod("2.4703282292062327208828439643411068618252990130716238221279"
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
               "437693412532098591327667236328125001e-324"));
    CHECK_EQ(5e-324,
        Strtod("0.000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000247032822920623272088284396434"
               "11068618252990130716238221279284125033775363510437593264991818081799618"
               "98982823477228588654633283551779698981993873980053909390631503565951557"
               "02263922908583924491051844359318028499365361525003193704576782492193656"
               "23669863658480757001585769269903706311928279558551332927834338409351978"
               "01553124659726357957462276646527282722005637400648549997709659947045402"
               "08281662262378573934507363390079677619305775067401763246736009689513405"
               "35537458516661134223766678604162159680461914467291840300530057530849048"
               "76539171138659164623952491262365388187963623937328042389101867234849766"
               "82350898633885879256283027559956575244555072551893136908362547791869486"
               "67994968324049705821028513185451396213837722826145437693412532098591327"
               "6672363281255"));
    CHECK_EQ(5e-324,
        Strtod("0.000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000000000000000000000000000000000"
               "00000000000000000000000000000000000000000247032822920623272088284396434"
               "11068618252990130716238221279284125033775363510437593264991818081799618"
               "98982823477228588654633283551779698981993873980053909390631503565951557"
               "02263922908583924491051844359318028499365361525003193704576782492193656"
               "23669863658480757001585769269903706311928279558551332927834338409351978"
               "01553124659726357957462276646527282722005637400648549997709659947045402"
               "08281662262378573934507363390079677619305775067401763246736009689513405"
               "35537458516661134223766678604162159680461914467291840300530057530849048"
               "76539171138659164623952491262365388187963623937328042389101867234849766"
               "82350898633885879256283027559956575244555072551893136908362547791869486"
               "67994968324049705821028513185451396213837722826145437693412532098591327"
               "667236328126"));
    CHECK_EQ(0.500000000000000166533453693773481063544750213623046875,
        Strtod("0.500000000000000166533453693773481063544750213623046875"));
    CHECK_EQ(3.518437208883201171875e13, Strtod("3.518437208883201171875e13"));
    CHECK_EQ(62.5364939768271845828, Strtod("62.5364939768271845828"));
    CHECK_EQ(8.10109172351e-10, Strtod("8.10109172351e-10"));
    CHECK_EQ(1.50000000000000011102230246251565404236316680908203125,
        Strtod("1.50000000000000011102230246251565404236316680908203125"));
    CHECK_EQ(9007199254740991.4999999999999999999999999999999995,
        Strtod("9007199254740991.4999999999999999999999999999999995"));
    CHECK_EQ(1.2345678901234567e22, Strtod("1.2345678901234567e22"));
    CHECK_EQ(2.2250738585072011e-308, Strtod("2.2250738585072011e-308"));
    CHECK_EQ(6.631236846766476e-316,
        Strtod("6.6312368714697582767853966302759672433990999473553031442499717"
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
               "49146967171293949598850675682115696218943412532098591327667236328125E-316"));
    CHECK_EQ(3.2379086165851934e-319,
        Strtod("3.2378839133029012895883524125015321748630376694231080599012970"
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
               "28251730299953143924168545708663913273994694463908672332763671875E-319"));
    CHECK_EQ(6.9533558078476524e-310,
        Strtod("6.953355807847677105972805215521891690222119817145950754416205607980030"
               "13154963668880611572639944188006538639986402869127553953941465283158479"
               "56685600829998895513577849614468960421131982842131079351102171626549398"
               "02416034676213829409720583759540476786936413816541621287843248433202369"
               "20991661224967600557302270324479971462211654218883777037602237117207955"
               "91258533828013962195524188394697705149041926576270603193728475623010741"
               "40442660237844114174497210955449896389180395827191602886654488182452409"
               "58398138944278337700150546201574501784875457466834216175949666176602002"
               "87528887833870748507731929971029979366198762266880963149896457660004790"
               "09083731736585750335262099860150896718774401964796827166283225641992040"
               "747894382698751809812609536720628966577351093292236328125E-310"));
    CHECK_EQ(3.3390932608534806e-319,
        Strtod("3.339068557571188581835713701280943911923401916998521771655656997328440"
               "31455961531816884914907466260909999811300946556642680817037843406572299"
               "16596426194677060348844249897410807907667784563321682004646515939958173"
               "71782125010668346652995912233993254584461125868481633343674905074271064"
               "40976309070801785658401977687881242531200881232626036303547481153223685"
               "33599053346255754042160606228586332807443018924703005556787346899784768"
               "70369853549413277156622170245846166991655321535529623870646888786637528"
               "99559280043617790174628627227337447170145299143304725786386460142425202"
               "47915673681950560773208853293843223323915646452641434007986196650406080"
               "77549162173963649264049738362290606875883456826586710961041737908872035"
               "803481241600376705491726170293986797332763671875E-319"));
    CHECK_EQ(2.2250738585072012e-308, Strtod("2.2250738585072012e-308"));
    CHECK_EQ(2.2250738585072011e-308, Strtod("2.2250738585072011e-308"));

    CHECK_EQ(6114917000000003e-14, Strtod("6114917000000003e-14"));
}
