#include "grisu3.h"
#include "../src/grisu3.h"

char* grisu3_Dtoa(char* next, char* last, double value, bool force_trailing_dot_zero)
{
    return grisu3::Dtoa(next, last, value, force_trailing_dot_zero);
}

char* grisu3_Ftoa(char* next, char* last, float value, bool force_trailing_dot_zero)
{
    return grisu3::Dtoa(next, last, value, force_trailing_dot_zero);
}
