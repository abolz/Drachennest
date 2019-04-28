#include "ryu.h"

#include "../src/ryu.h"

char* ryu_Dtoa(char* buf, int /*buflen*/, double value)
{
    return ryu::ToChars(buf, value);
}

char* ryu_Ftoa(char* buf, int /*buflen*/, float value)
{
    return ryu::ToChars(buf, value);
}

uint64_t ryu_Dtoa(int& exponent, double value)
{
    const auto dec = ryu::ToDecimal(value);
    exponent = dec.exponent;
    return dec.digits;
}

uint32_t ryu_Ftoa(int& exponent, float value)
{
    const auto dec = ryu::ToDecimal(value);
    exponent = dec.exponent;
    return dec.digits;
}
