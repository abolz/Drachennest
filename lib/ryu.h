#pragma once

#include <cstdint>

char* ryu_Dtoa(char* buf, int buflen, double value);

char* ryu_Ftoa(char* buf, int buflen, float value);

uint64_t ryu_Dtoa(int& exponent, double value);

uint32_t ryu_Ftoa(int& exponent, float value);
