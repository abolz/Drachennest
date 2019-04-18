#pragma once

char* grisu3_Dtoa(char* next, char* last, double value, bool force_trailing_dot_zero = false);

char* grisu3_Ftoa(char* next, char* last, float value, bool force_trailing_dot_zero = false);
