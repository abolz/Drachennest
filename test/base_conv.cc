#include "base_conv.h"

#if 0
#ifdef _MSC_VER
#define DTOA_INLINE inline
#define STRTOD_INLINE inline
#endif
#endif
#if 1
#define DTOA_OPTIMIZE_SIZE 1
#define STRTOD_OPTIMIZE_SIZE 1
#include "../src/ryu.h"
#include "../src/strtod_0.h"
#else
#include "../src/dtoa.h"
#include "../src/strtod.h"
#endif

char* base_conv_Dtoa(char* next, char* last, double value, bool force_trailing_dot_zero)
{
    return base_conv::Dtoa(next, last, value, force_trailing_dot_zero);
}

double base_conv_DecimalToDouble(char const* digits, int num_digits, int exponent)
{
    return base_conv::DecimalToDouble(digits, num_digits, exponent);
}

bool base_conv_Strtod(double& value, char const* first, char const* last)
{
    auto const res = base_conv::Strtod(value, first, last);

    if (res == base_conv::StrtodStatus::success) {
        return true;
    }

    value = std::numeric_limits<double>::quiet_NaN();
    return false;
}
