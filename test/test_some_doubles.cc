#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <random>

// #include "../src/grisu2.h"
// #include "../src/grisu3.h"
#include "../src/ryu.h"
// #include "../lib/floaxie.h"
// #include "../lib/fmt.h"

#include "double-conversion/double-conversion.h"

#include "scan_number.h"

static inline double FloatFromBits(uint64_t bits)
{
    double f;
    std::memcpy(&f, &bits, sizeof(uint64_t));
    return f;
}

static inline uint64_t BitsFromFloat(double f)
{
    uint64_t u;
    std::memcpy(&u, &f, sizeof(uint64_t));
    return u;
}

int main()
{
    constexpr int P = 53;
    constexpr uint64_t MaxF = (1ull << (P - 1)) - 1;
    constexpr int MinExp = 0;
    constexpr int MaxExp = 2047 - 1;
    constexpr int NumSignificandsPerExponent = 1 << 15;

    std::mt19937 random;

    uint64_t num_checked = 0;
    uint64_t num_optimal = 0;
    uint64_t num_short   = 0;

    for (int e = MinExp; e <= MaxExp; ++e)
    {
        uint32_t curr_num_checked = 0;
        uint32_t curr_num_optimal = 0;
        uint32_t curr_num_short   = 0;

        printf("e = %4d ... ", e);

        bool fail = false;
        for (int i = 0; i < NumSignificandsPerExponent; ++i)
        {
            ++curr_num_checked;

            std::uniform_int_distribution<uint64_t> gen(0, MaxF);
            const uint64_t f = gen(random);

            const uint64_t bits = (static_cast<uint64_t>(e) << (P - 1)) | f;
            const double value = FloatFromBits(bits);

            char buf[32];
            // char* end = grisu2::ToChars(buf, value);
            // char* end = grisu3::ToChars(buf, value);
            char* end = ryu::ToChars(buf, value);
            // char* end = floaxie_Dtoa(buf, 32, value);
            // char* end = fmt_Dtoa(buf, 32, value);
            end[0] = '\0';

            double_conversion::StringToDoubleConverter s2d(0, 0.0, 0.0, "inf", "nan");
            int unused;
            const double value_out = s2d.StringToDouble(buf, static_cast<int>(end - buf), &unused);

            const uint64_t bits_out = BitsFromFloat(value_out);
            if (bits != bits_out)
            {
                auto& conv = double_conversion::DoubleToStringConverter::EcmaScriptConverter();
                char tmp[32];
                double_conversion::StringBuilder builder(tmp, 32);
                conv.ToShortest(value, &builder);
                builder.Finalize();
                printf("\nFAIL: 0x%016llX [actual = %s] [expected = %s]\n", bits, buf, tmp);
                fail = true;
                break;
                // return -1;
            }
            {
                char tmp[32];
                double_conversion::StringBuilder builder(tmp, 32);
                auto& conv = double_conversion::DoubleToStringConverter::EcmaScriptConverter();
                conv.ToShortest(value, &builder);
                builder.Finalize();

                const auto num1 = ScanNumber(buf);
                const auto num2 = ScanNumber(tmp);
                assert(num1.digits.size() >= num2.digits.size());
                const auto len1 = static_cast<int>(num1.digits.size()) + num1.exponent;
                const auto len2 = static_cast<int>(num2.digits.size()) + num2.exponent;
                assert(len1 == len2);
                if (num1.digits.size() == num2.digits.size())
                {
                    ++curr_num_short;
                }
                else
                {
                    // printf("\nnum1 = %17s * 10^%d\n", num1.digits.c_str(), num1.exponent);
                    // printf("\nnum2 = %17s * 10^%d\n", num2.digits.c_str(), num2.exponent);
                }
                if (num1.digits == num2.digits)
                {
                    ++curr_num_optimal;
                }
                //else
                //{
                //    printf("\nNOT optimal: 0x%016llX [actual = %s] [expected = %s]\n", bits, buf, tmp);
                //    continue;
                //}
            }
        }

        if (!fail)
        {
            const uint32_t curr_not_short = curr_num_checked - curr_num_short;
            printf("optimal: %7.2f%%, not short: %7.2f%% (%u)\n",
                100.0 * (double)curr_num_optimal / (double)curr_num_checked,
                100.0 * (double)curr_not_short / (double)curr_num_checked,
                curr_not_short);
        }

        num_checked += curr_num_checked;
        num_optimal += curr_num_optimal;
        num_short   += curr_num_short;
    }

    printf("done.\n");
    printf("checked: %llu\n", num_checked);
    printf("optimal: %7.2f%% (%llu)\n", 100.0 * (double)num_optimal  / (double)num_checked, num_optimal);
    printf("short:   %7.2f%% (%llu)\n", 100.0 * (double)num_short / (double)num_checked, num_short);
    return 0;
}

#include "double-conversion/bignum-dtoa.cc"
#include "double-conversion/bignum.cc"
#include "double-conversion/cached-powers.cc"
#include "double-conversion/diy-fp.cc"
#include "double-conversion/double-conversion.cc"
#include "double-conversion/fast-dtoa.cc"
#include "double-conversion/fixed-dtoa.cc"
#include "double-conversion/strtod.cc"

// #include "../lib/floaxie.cc"
// #define FMT_HEADER_ONLY 1
// #include "../lib/fmt.cc"
