#include "charconv_f32.h"
#include "charconv_f64.h"

#include "catch.hpp"

#include <cmath>
#include <cstring>
#include <limits>

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
    const auto res = charconv::Strtof(str.data(), str.data() + str.size(), flt);
    CHECK(res.status != charconv::StrtofStatus::invalid);
    return flt;
}

struct Ftoa1 {
    char* operator()(char* buf, int /*buflen*/, float value) const {
        return charconv::Ftoa(buf, value);
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
    const auto bits = BitsFromFloat(value);
    char buf[charconv::FtoaMinBufferLength];

    char* const end = FtoaFn{}(buf, charconv::FtoaMinBufferLength, value);

    float value2;
    const auto res = charconv::Strtof(buf, end, value2);
    CHECK(res.status != charconv::StrtofStatus::invalid);
    CHECK(res.next == end);

    if (std::isnan(value))
    {
        CHECK(std::isnan(value2));
        //if (!std::isnan(value2))
        //{
        //    fprintf(stderr, "FAIL: nan 0x%08X\n", bits);
        //    return false;
        //}
    }
    else
    {
        const auto bits2 = BitsFromFloat(value2);
        CHECK(bits == bits2);
        //if (bits != bits2)
        //{
        //    fprintf(stderr, "FAIL: bits = 0x%08X != 0x%08X = bits2\n", bits, bits2);
        //    return false;
        //}
    }

    return true;
}

static void CheckStrtof(float value)
{
    CheckStrtofImpl<Ftoa1>(value);
    CheckStrtofImpl<Ftoa2>(value);
    CheckStrtofImpl<Ftoa3>(value);
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

#if 0
TEST_CASE("Strtof - All")
{
    printf("Strtof - All...\n");

    uint32_t curr_e = UINT32_MAX;
    uint32_t f = 0;
    for (;;)
    {
        uint32_t e = f >> 23;
        if (curr_e != e)
        {
            curr_e = e;
            printf("%u ...\n", e);
        }
        const auto value = FloatFromBits(f);
        if (!CheckStrtof(value))
            return;
        if (f == UINT32_MAX)
            break;
        ++f;
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
    const auto res = charconv::Strtod(str.data(), str.data() + str.size(), flt);
    CHECK(res.status != charconv::StrtodStatus::invalid);
    return flt;
}

struct Dtoa1 {
    char* operator()(char* buf, int /*buflen*/, double value) const {
        return charconv::Dtoa(buf, value);
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
    const auto bits = BitsFromFloat(value);
    char buf[charconv::DtoaMinBufferLength];

    char* const end = DtoaFn{}(buf, charconv::DtoaMinBufferLength, value);

    double value2;
    const auto res = charconv::Strtod(buf, end, value2);
    CHECK(res.status != charconv::StrtodStatus::invalid);
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
}

TEST_CASE("Strtod - Special")
{
    CHECK( 0.0 == Strtod("0"));
    CHECK( 0.0 == Strtod("0.0000000000000000000000000000000"));
    CHECK(-0.0 == Strtod("-0"));
    CHECK( 0.0 == Strtod("+0"));

    CheckStrtod(0.0f);
    CheckStrtod(-0.0f);
    CheckStrtod(std::numeric_limits<double>::infinity());
    CheckStrtod(-std::numeric_limits<double>::infinity());
    CheckStrtod(std::numeric_limits<double>::quiet_NaN());

    CHECK(std::isnan(Strtod("nan")));
    CHECK(std::isnan(Strtod("NaN")));
    CHECK(std::isnan(Strtod("nAn(_nananana123)")));

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
