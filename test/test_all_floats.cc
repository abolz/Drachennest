#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <string>

#include "double-conversion/double-conversion.h"

#include "../src/grisu2.h"
// #include "../src/grisu3.h"
// #include "../src/ryu.h"
// #include "../lib/floaxie.h"

#include "scan_number.h"

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

int main()
{
    constexpr int P = 24;
    constexpr uint32_t MaxF = (1u << (P - 1)) - 1;
    constexpr int MinExp = 0;
    constexpr int MaxExp = 255 - 1;

    uint32_t num_checked = 0;
    uint32_t num_optimal = 0;
    uint32_t num_short   = 0;

    for (int e = MinExp; e <= MaxExp; ++e)
    {
        uint32_t curr_num_checked = 0;
        uint32_t curr_num_optimal = 0;
        uint32_t curr_num_short   = 0;

        printf("e = %3d ... ", e);

        bool fail = false;
        for (uint32_t f = 0; f <= MaxF; ++f)
        // uint32_t f = 0;
        {
            ++curr_num_checked;

            const uint32_t bits = ((uint32_t)e << (P-1)) | f;
            const float value = FloatFromBits(bits);

            char buf[256];
            char* end = grisu2::ToChars(buf, value);
            // char* end = grisu3::ToChars(buf, value);
            // char* end = ryu::ToChars(buf, value);
            // char* end = jkj_Ftoa(buf, 256, value);
            // char* end = floaxie_Ftoa(buf, 256, value);
            end[0] = '\0';

            double_conversion::StringToDoubleConverter s2f(0, 0.0, 0.0, "inf", "nan");
            int unused;
            const float value_out = s2f.StringToFloat(buf, (int)(end - buf), &unused);

            const uint32_t bits_out = BitsFromFloat(value_out);
            if (bits != bits_out)
            {
                auto& conv = double_conversion::DoubleToStringConverter::EcmaScriptConverter();
                char tmp[256];
                double_conversion::StringBuilder builder(tmp, 256);
                conv.ToShortestSingle(value, &builder);
                builder.Finalize();
                printf("\nFAIL: 0x%08X [actual = %s] [expected = %s]\n", bits, buf, tmp);
                fail = true;
                break;
                // return -1;
            }
            {
                char tmp[256];
                double_conversion::StringBuilder builder(tmp, 256);
                auto& conv = double_conversion::DoubleToStringConverter::EcmaScriptConverter();
                conv.ToShortestSingle(value, &builder);
                char* endtmp = tmp + builder.position();
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
                    // printf("\nnum1 = %9s * 10^%d\n", num1.digits.c_str(), num1.exponent);
                    // printf("\nnum2 = %9s * 10^%d\n", num2.digits.c_str(), num2.exponent);
                }
                if (num1.digits == num2.digits)
                {
                    ++curr_num_optimal;
                }
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
    printf("checked: %u\n", num_checked);
    printf("optimal: %7.2f%% (%u)\n", 100.0 * (double)num_optimal  / (double)num_checked, num_optimal);
    printf("short:   %7.2f%% (%u)\n", 100.0 * (double)num_short / (double)num_checked, num_short);
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
