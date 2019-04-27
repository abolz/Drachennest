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
