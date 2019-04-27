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
