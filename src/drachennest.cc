#include "drachennest.h"

#include "grisu2.h"
#include "grisu3.h"
#include "ryu.h"

char* drachennest::dtoa_grisu2(char* buf, double value)
{
    return grisu2::ToChars(buf, value);
}

char* drachennest::dtoa_grisu3(char* buf, double value)
{
    return grisu3::ToChars(buf, value);
}

char* drachennest::dtoa_ryu(char* buf, double value)
{
    return ryu::ToChars(buf, value);
}

char* drachennest::ftoa_grisu2(char* buf, float value)
{
    return grisu2::ToChars(buf, value);
}

char* drachennest::ftoa_grisu3(char* buf, float value)
{
    return grisu3::ToChars(buf, value);
}

char* drachennest::ftoa_ryu(char* buf, float value)
{
    return ryu::ToChars(buf, value);
}
