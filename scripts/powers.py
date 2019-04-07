MIN_EXPONENT = -300
MAX_EXPONENT =  324
STEP = 8

def FloorLog2Pow10(e):
    return (e * 1741647) // 2**19

def RoundUp(num, den):
    p, r = divmod(num, den)
    if 2 * r >= den:
        p += 1
    return p

for k in range(MIN_EXPONENT, MAX_EXPONENT + 1, STEP):
    e = FloorLog2Pow10(k) + 1 - 64
    if k >= 0:
        if e > 0:
            p = RoundUp(10**k, 2**e)
        else:
            p = 10**k * 2**(-e)
    else:
        p = RoundUp(2**(-e), 10**(-k))
    assert p < 2**64
    assert p >= 2**(64 - 1)
    print "0x{:X}, // e = {:5d}, k = {:4d}".format(p, e, k)
