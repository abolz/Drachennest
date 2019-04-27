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
