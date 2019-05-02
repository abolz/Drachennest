#include "grisu3.h"

#include "../src/grisu3.h"

char* grisu3_Dtoa(char* buf, int /*buflen*/, double value)
{
    return grisu3::ToChars(buf, value);
}

char* grisu3_Ftoa(char* buf, int /*buflen*/, float value)
{
    return grisu3::ToChars(buf, value);
}

uint64_t grisu3_Dtoa(int& exponent, double value)
{
    const auto dec = grisu3::ToDecimal(value);
    exponent = dec.exponent;
    return dec.digits;
}

uint32_t grisu3_Ftoa(int& exponent, float value)
{
    const auto dec = grisu3::ToDecimal(value);
    exponent = dec.exponent;
    return dec.digits;
}
