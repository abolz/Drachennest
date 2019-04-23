#include "grisu3.h"

#if _MSC_VER
#define GRISU_INLINE static __forceinline
#else
#define GRISU_INLINE static inline
#endif
#include "../src/grisu3.h"

char* grisu3_Dtoa(char* buf, int buflen, double value)
{
    return grisu3::Dtoa(buf, buf + buflen, value);
}

char* grisu3_Ftoa(char* buf, int buflen, float value)
{
    return grisu3::Dtoa(buf, buf + buflen, value);
}
