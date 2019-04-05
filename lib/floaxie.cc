#include "floaxie.h"

#include "../ext/floaxie/ftoa.h"

char* floaxie_Dtoa(char* next, char* /*last*/, double value)
{
    return floaxie::ftoa(value, next);
}

char* floaxie_Ftoa(char* next, char* /*last*/, float value)
{
    return floaxie::ftoa(value, next);
}
