#include "ryu.h"

#if _MSC_VER
#define RYU_INLINE static __forceinline
#else
#define RYU_INLINE static inline
#endif
#include "../src/ryu.h"

char* ryu_Dtoa(char* buf, int buflen, double value)
{
    return ryu::Dtoa(buf, buf + buflen, value);
}

char* ryu_Ftoa(char* buf, int buflen, float value)
{
    return ryu::Dtoa(buf, buf + buflen, value);
}
