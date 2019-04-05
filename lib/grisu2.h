#pragma once

char* grisu2_Dtoa(char* next, char* last, double value, bool force_trailing_dot_zero = false);

char* grisu2_Ftoa(char* next, char* last, float value, bool force_trailing_dot_zero = false);
