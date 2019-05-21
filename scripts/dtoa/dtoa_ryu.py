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

def FloorLog2Pow5(e):
    assert e >= -12654
    assert e <=  12654
    return (e * 38955489) >> 24

def Ceil(num, den):
    assert num >= 0
    assert den > 0
    q, r = divmod(num, den)
    if r > 0:
        q += 1
    return q

def ComputePow5(k, bits = 128):
    assert bits > 0
    e = FloorLog2Pow5(k) + 1 - bits
    if k >= 0:
        if e > 0:
            f = 5**k // 2**e
        else:
            f = 5**k * 2**(-e)
    else:
        f = Ceil(2**(-e), 5**(-k))
    assert f >= 2**(bits - 1)
    assert f < 2**bits
    return f

#---------------------------------------------------------------------------------------------------
#
#---------------------------------------------------------------------------------------------------

CachedPow5Bits = 128
CachedPow5MinDecExp = -290
CachedPow5MaxDecExp =  325
CachedPow5 = [ComputePow5(k, bits=CachedPow5Bits) for k in range(CachedPow5MinDecExp, CachedPow5MaxDecExp + 1)]

def FloorLog10Pow2(e):
    assert e >= -2620
    assert e <=  2620
    return (e * 315653) >> 20

def FloorLog10Pow5(e):
    assert e >= -2620
    assert e <=  2620
    return (e * 732923) >> 20

def MulShift(x, y, j):
    assert x >= 0
    assert y >= 0
    assert j >= 0
    return (x * y) >> j

def MulPow5DivPow2(u, v, w, e5, e2):
    assert e5 >= CachedPow5MinDecExp
    assert e5 <= CachedPow5MaxDecExp
    k = FloorLog2Pow5(e5) + 1 - CachedPow5Bits
    j = e2 - k
#   pow5 = ComputePow5(e5)
    pow5 = CachedPow5[e5 - CachedPow5MinDecExp]
    a = MulShift(u, pow5, j)
    b = MulShift(v, pow5, j)
    c = MulShift(w, pow5, j)
    return a, b, c

def MultipleOfPow5(value, e5):
    if e5 <= 0:
        return True
    return value % 5**e5 == 0
    # while e5 > 0:
    #     q, r = divmod(value, 5)
    #     if r != 0:
    #         return False
    #     value = q
    #     e5 -= 1
    # return True

def MultipleOfPow2(value, e2):
    if e2 <= 0:
        return True
    return value % 2**e2 == 0

def DtoaRyu(m2, e2):
    isEven = (m2 % 2 == 0)
    # acceptLower = isEven
    # acceptUpper = isEven
    acceptBounds = isEven

    #---------------------------------------------------------------------------
    # Step 2:
    #
    # Determine the interval of information preserving outputs.
    #

    lowerBoundaryIsCloser = (m2 == HIDDEN_BIT and e2 > MIN_EXPONENT)

    e2 -= 2
    u = 4 * m2 - 2 + (1 if lowerBoundaryIsCloser else 0)
    v = 4 * m2
    w = 4 * m2 + 2

    #---------------------------------------------------------------------------
    # Step 3:
    #
    # Convert (u,v,w) * 2^e2 to a decimal power base.
    #

    atz = False # a[0], ..., a[i-1], a[i] == 0
    btz = False # b[0], ..., b[i-1]       == 0
    ctz = False # c[0], ..., c[i-1], c[i] == 0
    # We only need to know atz iff acceptLower is true.
    # We only need to know ctz iff acceptUpper is false.
    # Since acceptLower == acceptUpper, we only need to know one or the other.

    if e2 >= 0:
        q_prime = FloorLog10Pow2(e2)
        q = max(0, q_prime - 1)
        e10 = q

        # Z(x,e2,q) = (x * 2^e2) % 10^q == 0
        #           = p10(x * 2^e2) >= q
        #           = min(p2(x) + p2(e2), p5(x)) >= q
        #           = p2(x) + e2 >= q and p5(x) >= q
        #           = p5(x) >= q
        #           = x % 5^q == 0

        #   atz = MultipleOfPow5(u, q)
        #   btz = MultipleOfPow5(v, q)
        #   ctz = MultipleOfPow5(w, q)

        #   btz = MultipleOfPow5(v, q)
        #   if acceptBounds:
        #       atz = MultipleOfPow5(u, q)
        #   if not acceptBounds:
        #       ctz = MultipleOfPow5(w, q)

        # Since w - u <= 4, only one of u,v,w can be divisible by 5, if any.
        if q <= 22:
            if v % 5 == 0:
                btz = MultipleOfPow5(v, q)
            elif acceptBounds:
                atz = MultipleOfPow5(u, q)
            else:
                ctz = MultipleOfPow5(w, q)
    else:
        q_prime = FloorLog10Pow5(-e2)
        q = max(0, q_prime - 1)
        e10 = q + e2

        # Z(x,e2,q) = (x * 5^-e2) % 10^q == 0
        #           = min(p2(x), p5(x) - e2) >= q
        #           = p2(x) >= q and p5(x) - e2 >= q
        #           = p2(x) >= q
        #           = u % 2^q == 0

        #   btz = MulPow5DivPow2(v, q - 1)
        #   if acceptBounds:
        #       atz = MulPow5DivPow2(u, q)
        #   else:
        #       ctz = MultipleOfPow2(w, q)

        if q <= 1:
            # v = 4 * f2 = 2^2 * f2 => p2(v) >= 2 >= q
            btz = True
            if acceptBounds:
                # if lowerBoundaryIsCloser:
                #   u = 4 * f2 - 1 = 2 * (2*f2) - 1 ==> p2(u) == 0
                # else:
                #   u = 4 * f2 - 2 = 2 * (2*f2 - 1) ==> p2(u) == 1
#               atz = (q == 0) if lowerBoundaryIsCloser else True
                atz = not lowerBoundaryIsCloser or q == 0
            else:
                # w = 4 * f2 + 2 = 2 * (2 * f2 + 1) ==> p2(w) == 1 >= q
                ctz = True
        elif q <= 53 + 2:
            btz = MultipleOfPow2(v, q - 1)
            # p2(u) <= 1 < q ==> atz = False
            # p2(w) == 1 < q ==> ctz = False

    # (a, b, c) = (u, v, w) * 2^e2 / 10^e10 = (u, v, w) * 5^(-e10) / 2^(e10 - e2)
    a, b, c = MulPow5DivPow2(u, v, w, -e10, e10 - e2)
    assert a < c - 1

    #---------------------------------------------------------------------------
    # Step 4:
    #
    # Find a shortest, correctly-rounded decimal representations in the interval
    # of valid representations.
    #

    # if not acceptBounds and ctz:
    #
    # Since we only assign to ctz iff acceptBounds == false, and ctz is
    # initialized to false, the test may be simplified here.
    if ctz:
        c -= 1

    bi = 0

    # If q_prime == 0, then bi = 0.
    # Otherwise we need to execute the loop at least once to assign a valid value bi.
    assert (q_prime == 0) or (a / 10 < c / 10)

    if atz or btz:
        while a / 10 < c / 10:
            atz = atz and a % 10 == 0
            btz = btz and bi == 0
            bi = b % 10
            a /= 10
            b /= 10
            c /= 10
            e10 += 1
        # if acceptBounds and atz:
        #
        # Since we only assign to atz iff acceptBounds == true, and atz is
        # initialized to false, the test may be simplified here.
        if atz:
            while a % 10 == 0:
                btz = btz and bi == 0
                bi = b % 10
                a /= 10
                b /= 10
                c /= 10
                e10 += 1

        roundDown = bi < 5 or (bi == 5 and btz and b % 2 == 0)
    else:
        roundDown = True
        while a / 10 < c / 10:
            roundDown = (b % 10 < 5)
            a /= 10
            b /= 10
            c /= 10
            e10 += 1

    # A return value of b is valid if and only if a != b or atz is true.
    # A return value of b + 1 is valid if and only if b + 1 <= c.
#   if (not roundDown and b < c) or (a == b and not (acceptBounds and atz)):
    if (not roundDown and b < c) or (a == b and not atz):
        b += 1

    return b, e10
