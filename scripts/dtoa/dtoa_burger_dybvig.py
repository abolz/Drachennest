PRECISION = 53
EXPONENT_BITS = 11

# PRECISION = 24
# EXPONENT_BITS = 8

HIDDEN_BIT = 2**(PRECISION - 1)
BIAS = 2**(EXPONENT_BITS - 1) - 1 + (PRECISION - 1)
MIN_EXPONENT = 1 - BIAS
MAX_EXPONENT = 2**EXPONENT_BITS - 2 - BIAS

#===================================================================================================
# Dtoa
#===================================================================================================

# Burger, Dybvig, "Printing Floating-Point Numbers Quickly and Accurately",
# Proceedings of the ACM SIGPLAN 1996 Conference on Programming Language
# Design and Implementation, PLDI 1996

def EffectivePrecision(f):
    """Returns the effective precision of the significand, aka. f.bit_length()"""

    assert f > 0
    assert f < 2**PRECISION

    p = PRECISION
    while f < 2**(PRECISION - 1):
        f *= 2
        p -= 1

    return p

def CeilLog10Pow2(e):
    """Returns ceil(log_10(2^e))"""

    assert e >= -1650
    assert e <=  1650
    return (int(e) * 78913 + (2**18 - 1)) // 2**18 # floor-division (SAR)

def DtoaBurgerDybvig(f, e):
    assert f > 0
    assert f < 2**PRECISION
    assert e >= MIN_EXPONENT
    assert e <= MAX_EXPONENT

    #
    # Init
    #

    isEven = (f % 2 == 0)
    acceptBounds = isEven
    lowerBoundaryIsCloser = (f == HIDDEN_BIT and e != MIN_EXPONENT)

    if e >= 0:
        r, s, mp, mm = f * 2 * 2**e, 2, 2**e, 2**e
    else:
        r, s, mp, mm = f * 2, 2**(-e) * 2, 1, 1

    # Could do after scaling to keep the numbers a tiny bit smaller?!?!
    if lowerBoundaryIsCloser:
        r, s, mp = r * 2, s * 2, mp * 2

    #
    # Scale into the range [0.1, 1.0)
    #   aka: Find the smallest integer k, such that (r + m+) / s <= 10^k
    #   aka: k = ceil(log_10((r + m+) / s))
    #

    if False:
        k = 0
        while True:
            rp = r + mp
            if (rp >= s) if acceptBounds else (rp > s):
                s, k = s * 10, k + 1
            elif (rp * 10 < s) if acceptBounds else (rp * 10 <= s):
                r, mp, mm, k = r * 10, mp * 10, mm * 10, k - 1
            else:
                break
    else:
#       p = f.bit_length() # Effective precision
        p = EffectivePrecision(f)

        # Estimate:
        k = CeilLog10Pow2(e + (p - 1))
        if k >= 0:
            s *= 10**k
        else:
            r, mp, mm = r * 10**(-k), mp * 10**(-k), mm * 10**(-k)

        # Fixup:
        if (r + mp >= s) if acceptBounds else (r + mp > s):
            s, k = s * 10, k + 1

    assert (r + mp < s) if acceptBounds else (r + mp <= s)
    assert ((r + mp) * 10 >= s) if acceptBounds else ((r + mp) * 10 > s)

    #
    # Generate
    #

    d = []
    while True:
        assert r > 0

        r, mp, mm = r * 10, mp * 10, mm * 10
        q, r = divmod(r, s)
        assert q <= 9

        tc1 = (r <= mm) if acceptBounds else (r < mm)
        tc2 = (r + mp >= s) if acceptBounds else (r + mp > s)

        if tc1 and tc2:
            # Return the number closer to v. If the two are equidistant
            # from v, use **some** strategy to break the tie.
            if (r * 2 > s) or (r * 2 == s and q % 2 != 0):
                q += 1
        elif not tc1 and tc2:
            q += 1

        assert q >= 0
        assert q <= 9
        d.append(q) # d = 10 * d + q
        k -= 1

        if tc1 or tc2:
            break

    return long(''.join(map(str, d))), k
