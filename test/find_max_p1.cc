#include "../src/grisu2.h"
#include <stdio.h>
#include <limits.h>

using namespace grisu2::impl;

int main()
{
//  constexpr uint64_t MaxF = ((uint64_t{1} << 53) - 1) << 11;
    constexpr uint64_t MaxF = UINT64_MAX;
    constexpr int MinExp = -1137;
    constexpr int MaxExp =   960;

    uint32_t max_p1 = 0;
    for (int e = MinExp; e <= MaxExp; ++e)
    {
        auto const v = DiyFp(MaxF, e);
        auto const cached = GetCachedPowerForBinaryExponent(e);
        auto const c_minus_k = DiyFp(cached.f, cached.e);
        auto const w = Multiply(v, c_minus_k);
        auto const p1_64 = w.f >> -w.e;
        assert(p1_64 <= UINT32_MAX);
        auto const p1 = static_cast<uint32_t>(p1_64);
        if (max_p1 < p1)
            max_p1 = p1;
    }

    printf("max_p1 = %u [0x%08X]\n", max_p1, max_p1);
}
