#include "double-conversion.h"
#include "../ext/double-conversion/double-conversion.h"

char* double_conversion_Dtoa(char* buf, int buflen, double value)
{
    auto& converter = double_conversion::DoubleToStringConverter::EcmaScriptConverter();
    double_conversion::StringBuilder builder(buf, buflen);
    converter.ToShortest(value, &builder);
    return buf + builder.position();
}

char* double_conversion_Ftoa(char* buf, int buflen, float value)
{
    auto& converter = double_conversion::DoubleToStringConverter::EcmaScriptConverter();
    double_conversion::StringBuilder builder(buf, buflen);
    converter.ToShortestSingle(value, &builder);
    return buf + builder.position();
}

double double_conversion_Strtod(char* buf, int len)
{
    double_conversion::StringToDoubleConverter s2d(0, 0.0, 1.0, "inf", "nan");
    int processed_characters_count = 0;
    return s2d.StringToDouble(buf, len, &processed_characters_count);
}

float double_conversion_Strtof(char* buf, int len)
{
    double_conversion::StringToDoubleConverter s2d(0, 0.0, 1.0, "inf", "nan");
    int processed_characters_count = 0;
    return s2d.StringToFloat(buf, len, &processed_characters_count);
}

#include "../ext/double-conversion/bignum-dtoa.cc"
#include "../ext/double-conversion/bignum.cc"
#include "../ext/double-conversion/cached-powers.cc"
#include "../ext/double-conversion/diy-fp.cc"
#include "../ext/double-conversion/double-conversion.cc"
#include "../ext/double-conversion/fast-dtoa.cc"
#include "../ext/double-conversion/fixed-dtoa.cc"
#include "../ext/double-conversion/strtod.cc"
