#pragma once

#include <cstdint>

char* grisu2_Dtoa(char* buf, int buflen, double value);

char* grisu2_Ftoa(char* buf, int buflen, float value);

uint64_t grisu2_Dtoa(int& exponent, double value);

uint32_t grisu2_Ftoa(int& exponent, float value);
