#pragma once

#include <cassert>
#include <cstdint>
#include <string>
#include <cstdio>

struct ScanNumberResult {
    std::string digits;
    int exponent;
};

inline bool IsDigit(char ch)
{
    return '0' <= ch && ch <= '9';
}

inline int DigitValue(char ch)
{
    assert(IsDigit(ch));
    return ch - '0';
}

// Assumes JSON number format ¯\_(ツ)_/¯
inline ScanNumberResult ScanNumber(char const* next, char const* last)
{
    // printf("ScanNumber(%.*s)\n", static_cast<int>(last - next), next);

    std::string digits;
    int exponent = 0;

    assert(next != last);
    assert(IsDigit(*next));

    if (*next == '0')
    {
        // Number is of the form 0[.xxx][e+nnn].
        // The leading zero here is not a significant digit.
        ++next;
    }
    else
    {
        for (;;)
        {
            digits += *next;
            ++next;
            if (next == last || !IsDigit(*next))
                break;
        }
    }

    if (next != last && *next == '.')
    {
        ++next;
        assert(next != last);
        assert(IsDigit(*next));

        if (digits.empty())
        {
            // Number is of the form 0.xxx[e+nnn].
            // Skip leading zeros in the fractional part and adjust the exponent.
            while (*next == '0')
            {
                --exponent;
                ++next;
                if (next == last)
                    return {"0", 0};
            }
        }

        while (IsDigit(*next))
        {
            digits += *next;
            --exponent;
            ++next;
            if (next == last)
                break;
        }
    }

    if (next != last && (*next == 'e' || *next == 'E'))
    {
        if (digits.empty())
        {
            // Number is of the form 0[.000]e+nnn.
            return {"0", 0};
        }

        ++next;
        assert(next != last);

        bool const exp_is_neg = (*next == '-');
        if (exp_is_neg || *next == '+')
        {
            ++next;
            assert(next != last);
        }

        // No overflow checks...
        int e = 0;
        while (next != last)
        {
            e = 10 * e + DigitValue(*next);
            ++next;
        }

        exponent += exp_is_neg ? -e : e;
    }

    // Move trailing zeros into the exponent
    while (!digits.empty() && digits.back() == '0')
    {
        digits.pop_back();
        exponent++;
    }

    // Normalize "0.0" and "0"
    if (digits.empty())
    {
        return {"0", 0};
    }

    return {digits, exponent};
}

inline ScanNumberResult ScanNumber(std::string const& str)
{
    char const* next = str.c_str();
    char const* last = str.c_str() + str.size();

    return ScanNumber(next, last);
}