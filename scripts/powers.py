MIN_EXPONENT = -300
MAX_EXPONENT =  324
STEP = 8
BITS = 64

# MIN_EXPONENT = 0
# MAX_EXPONENT = 27
# STEP = 1
# BITS = 64

# MIN_EXPONENT = -37
# MAX_EXPONENT =  46
# STEP = 1
# BITS = 32

def FloorLog2Pow10(e):
    return (e * 1741647) // 2**19

def RoundUp(num, den):
    p, r = divmod(num, den)
    if 2 * r >= den:
        p += 1
    return p

for k in range(MIN_EXPONENT, MAX_EXPONENT + 1, STEP):
    e = FloorLog2Pow10(k) + 1 - BITS
    if k >= 0:
        if e > 0:
            p = RoundUp(10**k, 2**e)
        else:
            p = 10**k * 2**(-e)
    else:
        p = RoundUp(2**(-e), 10**(-k))
    assert p < 2**BITS
    assert p >= 2**(BITS - 1)
    print "0x{:X}, // e = {:5d}, k = {:4d}".format(p, e, k)
