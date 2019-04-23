#include "charconv.h"

#if defined(_MSC_VER)
#include <charconv>

char* charconv_Dtoa(char* buf, int buflen, double value)
{
	return std::to_chars(buf, buf + buflen, value, std::chars_format::general).ptr;
}

char* charconv_Ftoa(char* buf, int buflen, float value)
{
	return std::to_chars(buf, buf + buflen, value, std::chars_format::general).ptr;
}
#else
char* charconv_Dtoa(char* buf, int /*buflen*/, double /*value*/)
{
	return buf;
}

char* charconv_Ftoa(char* buf, int /*buflen*/, float /*value*/)
{
	return buf;
}
#endif
