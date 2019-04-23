#pragma once

char* double_conversion_Dtoa(char* buf, int buflen, double value);

char* double_conversion_Ftoa(char* buf, int buflen, float value);

double double_conversion_Strtod(char* buf, int len);

float double_conversion_Strtof(char* buf, int len);
