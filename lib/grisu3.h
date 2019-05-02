#pragma once

#include <cstdint>

char* grisu3_Dtoa(char* buf, int buflen, double value);

char* grisu3_Ftoa(char* buf, int buflen, float value);

uint64_t grisu3_Dtoa(int& exponent, double value);

uint32_t grisu3_Ftoa(int& exponent, float value);
