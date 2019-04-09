#include "milo.h"

#include "../ext/milo/dtoa_milo.h"

char* milo_Dtoa(char* next, char* /*last*/, double value)
{
    return dtoa_milo(value, next);
}
