#pragma once

char* base_conv_Dtoa(char* next, char* last, double value, bool force_trailing_dot_zero = false);

double base_conv_DecimalToDouble(char const* digits, int num_digits, int exponent);

bool base_conv_Strtod(double& value, char const* first, char const* last);
