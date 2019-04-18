#include "grisu2.h"
#include "../src/grisu2.h"

char* grisu2_Dtoa(char* next, char* last, double value, bool force_trailing_dot_zero)
{
    return grisu2::Dtoa(next, last, value, force_trailing_dot_zero);
}

char* grisu2_Ftoa(char* next, char* last, float value, bool force_trailing_dot_zero)
{
    return grisu2::Dtoa(next, last, value, force_trailing_dot_zero);
}
