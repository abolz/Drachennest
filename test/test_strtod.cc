#include "catch.hpp"
#include "../src/ryu.h"
#include "../src/ryu_charconv.h"

#define CHECK_EQ(X, Y) CHECK(X == Y)

static double Strtod(const std::string& digits, int exponent)
{
    uint64_t m10 = 0;
    for (char ch : digits)
        m10 = 10 * m10 + static_cast<uint32_t>(ch - '0');

    return RyuToBinary64(m10, static_cast<int>(digits.size()), exponent);
}

static void CheckStrtod(double value)
{
    const auto dec = ryu::ToDecimal(value);
    const auto m10 = dec.digits;
    const auto m10len = dtoa::impl::DecimalLength(m10);
    const auto e10 = dec.exponent;
    const auto value2 = RyuToBinary64(m10, m10len, e10);
    CHECK_EQ(value, value2);
}

static void CheckStrtof(float value)
{
    const auto dec = ryu::ToDecimal(value);
    const auto m10 = dec.digits;
    const auto m10len = dtoa::impl::DecimalLength(m10);
    const auto e10 = dec.exponent;
    const auto value2 = RyuToBinary32(m10, m10len, e10);
    CHECK_EQ(value, value2);
}

TEST_CASE("Strtof - 0")
{
    const auto value = RyuToBinary32(999999999, 9, 0);
    CHECK_EQ(999999999.0f, value);
}

TEST_CASE("Strtod - 0")
{
    constexpr auto Inf = std::numeric_limits<double>::infinity();

    CHECK_EQ(5e-324, Strtod("5", -324));
    CHECK_EQ(Inf, Strtod("1", 309));
    CheckStrtod(DBL_MIN);
    CheckStrtod(DBL_MAX);
    CheckStrtod(DBL_TRUE_MIN);
    CheckStrtod(DBL_EPSILON);
    CHECK_EQ(1e-324, Strtod("1", -324));
    CHECK_EQ(2e-324, Strtod("2", -324));
    CHECK_EQ(3e-324, Strtod("3", -324));
    CHECK_EQ(4e-324, Strtod("4", -324));
    CHECK_EQ(5e-324, Strtod("5", -324));

    const auto value = RyuToBinary64(99999999999999999ull, 17, 0);
    CHECK_EQ(99999999999999999.0, value);
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

    CHECK_EQ(5e+125, Strtod("5", 125));
    CHECK_EQ(69e+267, Strtod("69", 267));
    CHECK_EQ(999e-26, Strtod("999", -26));
    CHECK_EQ(7861e-34, Strtod("7861", -34));
    CHECK_EQ(75569e-254, Strtod("75569", -254));
    CHECK_EQ(928609e-261, Strtod("928609", -261));
    CHECK_EQ(9210917e+80, Strtod("9210917", 80));
    CHECK_EQ(84863171e+114, Strtod("84863171", 114));
    CHECK_EQ(653777767e+273, Strtod("653777767", 273));
    CHECK_EQ(5232604057e-298, Strtod("5232604057", -298));
    CHECK_EQ(27235667517e-109, Strtod("27235667517", -109));
    CHECK_EQ(653532977297e-123, Strtod("653532977297", -123));
    CHECK_EQ(3142213164987e-294, Strtod("3142213164987", -294));
    CHECK_EQ(46202199371337e-72, Strtod("46202199371337", -72));
    CHECK_EQ(231010996856685e-73, Strtod("231010996856685", -73));
    CHECK_EQ(9324754620109615e+212, Strtod("9324754620109615", 212));
    CHECK_EQ(78459735791271921e+49, Strtod("78459735791271921", 49));
    CheckStrtod(272104041512242479e+200);
    CheckStrtod(6802601037806061975e+198);
    CheckStrtod(20505426358836677347e-221);
    CheckStrtod(836168422905420598437e-234);
    CheckStrtod(4891559871276714924261e+222);

    //
    // Table 2:
    // Stress Inputs for Conversion to 53-bit Binary, > 1/2 ULP
    //

    CHECK_EQ(9e-265, Strtod("9", -265));
    CHECK_EQ(85e-37, Strtod("85", -37));
    CHECK_EQ(623e+100, Strtod("623", 100));
    CHECK_EQ(3571e+263, Strtod("3571", 263));
    CHECK_EQ(81661e+153, Strtod("81661", 153));
    CHECK_EQ(920657e-23, Strtod("920657", -23));
    CHECK_EQ(4603285e-24, Strtod("4603285", -24));
    CHECK_EQ(87575437e-309, Strtod("87575437", -309));
    CHECK_EQ(245540327e+122, Strtod("245540327", 122));
    CHECK_EQ(6138508175e+120, Strtod("6138508175", 120));
    CHECK_EQ(83356057653e+193, Strtod("83356057653", 193));
    CHECK_EQ(619534293513e+124, Strtod("619534293513", 124));
    CHECK_EQ(2335141086879e+218, Strtod("2335141086879", 218));
    CHECK_EQ(36167929443327e-159, Strtod("36167929443327", -159));
    CHECK_EQ(609610927149051e-255, Strtod("609610927149051", -255));
    CHECK_EQ(3743626360493413e-165, Strtod("3743626360493413", -165));
    CHECK_EQ(94080055902682397e-242, Strtod("94080055902682397", -242));
    CheckStrtod(899810892172646163e+283);
    CheckStrtod(7120190517612959703e+120);
    CheckStrtod(25188282901709339043e-252);
    CheckStrtod(308984926168550152811e-52);
    CheckStrtod(6372891218502368041059e+064);

    //
    // Table 18:
    // Stress Inputs for Conversion to 56-bit Binary, < 1/2 ULP
    //

    CHECK_EQ(7e-27, Strtod("7", -27));
    CHECK_EQ(37e-29, Strtod("37", -29));
    CHECK_EQ(743e-18, Strtod("743", -18));
    CHECK_EQ(7861e-33, Strtod("7861", -33));
    CHECK_EQ(46073e-30, Strtod("46073", -30));
    CHECK_EQ(774497e-34, Strtod("774497", -34));
    CHECK_EQ(8184513e-33, Strtod("8184513", -33));
    CHECK_EQ(89842219e-28, Strtod("89842219", -28));
    CHECK_EQ(449211095e-29, Strtod("449211095", -29));
    CHECK_EQ(8128913627e-40, Strtod("8128913627", -40));
    CHECK_EQ(87365670181e-18, Strtod("87365670181", -18));
    CHECK_EQ(436828350905e-19, Strtod("436828350905", -19));
    CHECK_EQ(5569902441849e-49, Strtod("5569902441849", -49));
    CHECK_EQ(60101945175297e-32, Strtod("60101945175297", -32));
    CHECK_EQ(754205928904091e-51, Strtod("754205928904091", -51));
    CHECK_EQ(5930988018823113e-37, Strtod("5930988018823113", -37));
    CHECK_EQ(51417459976130695e-27, Strtod("51417459976130695", -27));
    CheckStrtod(826224659167966417e-41);
    CheckStrtod(9612793100620708287e-57);
    CheckStrtod(93219542812847969081e-39);
    CheckStrtod(544579064588249633923e-48);
    CheckStrtod(4985301935905831716201e-48);

    //
    // Table 19:
    // Stress Inputs for Conversion to 56-bit Binary, > 1/2 ULP
    //

    CHECK_EQ(9e+26, Strtod("9", 26));
    CHECK_EQ(79e-8, Strtod("79", -8));
    CHECK_EQ(393e+26, Strtod("393", 26));
    CHECK_EQ(9171e-40, Strtod("9171", -40));
    CHECK_EQ(56257e-16, Strtod("56257", -16));
    CHECK_EQ(281285e-17, Strtod("281285", -17));
    CHECK_EQ(4691113e-43, Strtod("4691113", -43));
    CHECK_EQ(29994057e-15, Strtod("29994057", -15));
    CHECK_EQ(834548641e-46, Strtod("834548641", -46));
    CHECK_EQ(1058695771e-47, Strtod("1058695771", -47));
    CHECK_EQ(87365670181e-18, Strtod("87365670181", -18));
    CHECK_EQ(872580695561e-36, Strtod("872580695561", -36));
    CHECK_EQ(6638060417081e-51, Strtod("6638060417081", -51));
    CHECK_EQ(88473759402752e-52, Strtod("88473759402752", -52));
    CHECK_EQ(412413848938563e-27, Strtod("412413848938563", -27));
    CHECK_EQ(5592117679628511e-48, Strtod("5592117679628511", -48));
    CHECK_EQ(83881765194427665e-50, Strtod("83881765194427665", -50));
    CheckStrtod(638632866154697279e-35);
    CheckStrtod(3624461315401357483e-53);
    CheckStrtod(75831386216699428651e-30);
    CheckStrtod(356645068918103229683e-42);
    CheckStrtod(7022835002724438581513e-33);
}

TEST_CASE("Strtod - 2")
{
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

TEST_CASE("Strtof - all of them")
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

            const auto dec = ryu::ToDecimal(value);
            const auto m10 = dec.digits;
            const auto m10len = dtoa::impl::DecimalLength(m10);
            const auto e10 = dec.exponent;

            const auto value_out = RyuToBinary32(m10, m10len, e10);
            const auto bits_out = BitsFromFloat(value_out);

            if (bits != bits_out)
            {
                CHECK(false);
            }
        }
    }
}
#endif
