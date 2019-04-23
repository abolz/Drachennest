#include "grisu2.h"

#if _MSC_VER
#define GRISU_INLINE static __forceinline
#else
#define GRISU_INLINE static inline
#endif
#include "../src/grisu2.h"

char* grisu2_Dtoa(char* buf, int buflen, double value)
{
    return grisu2::Dtoa(buf, buf + buflen, value);
}

char* grisu2_Ftoa(char* buf, int buflen, float value)
{
    return grisu2::Dtoa(buf, buf + buflen, value);
}
