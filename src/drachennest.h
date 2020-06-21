#pragma once

namespace drachennest {

char* dtoa_grisu2(char* buf, double value);
char* dtoa_grisu3(char* buf, double value);

char* ftoa_grisu2(char* buf, float value);
char* ftoa_grisu3(char* buf, float value);

} // namespace drachennest
