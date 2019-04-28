#include "grisu2.h"

#include "../src/grisu2.h"

char* grisu2_Dtoa(char* buf, int /*buflen*/, double value)
{
    return grisu2::ToChars(buf, value);
}

char* grisu2_Ftoa(char* buf, int /*buflen*/, float value)
{
    return grisu2::ToChars(buf, value);
}

uint64_t grisu2_Dtoa(int& exponent, double value)
{
    const auto dec = grisu2::ToDecimal(value);
    exponent = dec.exponent;
    return dec.digits;
}

uint32_t grisu2_Ftoa(int& exponent, float value)
{
    const auto dec = grisu2::ToDecimal(value);
    exponent = dec.exponent;
    return dec.digits;
}
