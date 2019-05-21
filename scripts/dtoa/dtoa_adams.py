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

def DtoaAdams(m2, e2):
    isEven = (m2 % 2 == 0)
    acceptSmaller = isEven
    acceptLarger = isEven

    #
    # Step 2:
    # Determine the interval of information preserving outputs.
    #

    lowerBoundaryIsCloser = (m2 == HIDDEN_BIT and e2 > MIN_EXPONENT)

    u = 4 * m2 - 2 + (1 if lowerBoundaryIsCloser else 0)
    v = 4 * m2
    w = 4 * m2 + 2
    e2 -= 2
    # print 'u = {}'.format(u)
    # print 'v = {}'.format(v)
    # print 'w = {}'.format(w)
    # print 'e2 = {}'.format(e2)

    #
    # Step 3:
    # Convert (u,v,w) * 2^e2 to a decimal power base.
    #

    if e2 >= 0:
        e10 = 0
        a = u * 2**e2
        b = v * 2**e2
        c = w * 2**e2
    else:
        e10 = e2
        a = u * 5**(-e2)
        b = v * 5**(-e2)
        c = w * 5**(-e2)

    assert a < c - 1
    # print 'a = {}'.format(a)
    # print 'b = {}'.format(b)
    # print 'c = {}'.format(c)
    # print 'e10 = {}'.format(e10)

    #
    # Step 4:
    # Find a shortest, correctly-rounded decimal representations in the interval
    # of valid representations.
    #

    if not acceptLarger:
        c -= 1

    atz = True # a[0],...,a[i] are 0s
    btz = True # b[0],...,b[i-1] are 0s
    bi = 0
    while a / 10 < c / 10:
        atz = atz and a % 10 == 0
        btz = btz and bi == 0
        bi = b % 10
        a /= 10
        b /= 10
        c /= 10
        e10 += 1
        # print '--- a/10 < c/10 ---'
        # print 'atz = {}'.format(atz)
        # print 'btz = {}'.format(btz)
        # print 'bi = {}'.format(bi)
        # print 'a = {}'.format(a)
        # print 'b = {}'.format(b)
        # print 'c = {}'.format(c)
        # print 'e10 = {}'.format(e10)
    if acceptSmaller and atz:
        while a % 10 == 0:
            btz = btz and bi == 0
            bi = b % 10
            a /= 10
            b /= 10
            c /= 10
            e10 += 1
            # print '--- a%10 == 0 ---'
            # print 'atz = {}'.format(atz)
            # print 'btz = {}'.format(btz)
            # print 'bi = {}'.format(bi)
            # print 'a = {}'.format(a)
            # print 'b = {}'.format(b)
            # print 'c = {}'.format(c)
            # print 'e10 = {}'.format(e10)

    roundDown = bi < 5 or (bi == 5 and btz and b % 2 == 0)
    # print 'roundDown = {}'.format(roundDown)
    if (not roundDown and b < c) or (a == b and not atz):
        # print 'rounding up...'
        b += 1

    # print 'result = {} * 10^{}'.format(b, e10)
    return b, e10
