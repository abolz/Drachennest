#include "ryu.h"

#include "../src/ryu.h"

#if 1
char* ryu_Dtoa(char* buf, int /*buflen*/, double value)
{
    return ryu::ToChars(buf, value);
}

char* ryu_Ftoa(char* buf, int /*buflen*/, float value)
{
    return ryu::ToChars(buf, value);
}
#else
char* ryu_Dtoa(char* buf, int /*buflen*/, double value)
{
    char tmp[64];
    const auto len = ryu::ToChars(tmp, value) - tmp;
    std::memcpy(buf, tmp, static_cast<size_t>(len));
    return buf + len;
}

char* ryu_Ftoa(char* buf, int /*buflen*/, float value)
{
    char tmp[64];
    const auto len = ryu::ToChars(tmp, value) - tmp;
    std::memcpy(buf, tmp, static_cast<size_t>(len));
    return buf + len;
}
#endif

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
