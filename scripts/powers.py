#===============================================================================
#
#===============================================================================

def FormatHexChunks(n, bits_per_chunk = 64):
    assert bits_per_chunk > 0
    assert bits_per_chunk % 4 == 0
    q, r = divmod(n, 2**bits_per_chunk)
    s = '0x{:0{}X}'.format(r, bits_per_chunk/4)
    chunks = 1
    while q != 0:
        q, r = divmod(q, 2**bits_per_chunk)
        s = '0x{:0{}X}, '.format(r, bits_per_chunk/4) + s
        chunks += 1
    if chunks > 1:
        s = '{' + s + '}'
    return s

#===============================================================================
# Grisu
#===============================================================================

def FloorLog2Pow10(e):
    return (e * 1741647) // 2**19

def RoundUp(num, den):
    p, r = divmod(num, den)
    if 2 * r >= den:
        p += 1
    return p

def ComputeGrisuPower(k, bits):
    e = FloorLog2Pow10(k) + 1 - bits
    if k >= 0:
        if e > 0:
            f = RoundUp(10**k, 2**e)
        else:
            f = 10**k * 2**(-e)
    else:
        f = RoundUp(2**(-e), 10**(-k))
    assert f >= 2**(bits - 1)
    assert f < 2**bits
    return f, e

def PrintGrisuPowers(bits, min_exponent, max_exponent, step = 1):
    print '// Let e = FloorLog2Pow10(k) + 1 - {}'.format(bits)
    print '// For k >= 0, stores 10^k in the form: round_up(10^k / 2^e )'
    print '// For k <= 0, stores 10^k in the form: round_up(2^-e / 10^-k)'
    for k in range(min_exponent, max_exponent + 1, step):
        f, e = ComputeGrisuPower(k, bits)
        print FormatHexChunks(f, bits_per_chunk=64) + ', // e = {:5d}, k = {:4d}'.format(e, k)

# For double-precision:
# PrintGrisuPowers(bits=64, min_exponent=-300, max_exponent=324, step=8)

# For single-precision:
# PrintGrisuPowers(bits=32, min_exponent=-37, max_exponent=46, step=1)

#===============================================================================
# Ryu
#===============================================================================

def FloorLog2Pow5(e):
    return (e * 1217359) // 2**19

def Ceil(num, den):
    q, r = divmod(num, den)
    if r > 0:
        q += 1
    return q

def ComputeRyuPower(k, bits):
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
    return f, e

def PrintRyuPowers(bits, min_exponent, max_exponent):
    print '// Let e = FloorLog2Pow5(k) + 1 - {}'.format(bits)
    print '// For k >= 0, stores 5^k in the form: floor( 5^k / 2^e )'
    print '// For k <= 0, stores 5^k in the form:  ceil(2^-e / 5^-k)'
    for k in range(min_exponent, max_exponent + 1):
        f, e = ComputeRyuPower(k, bits)
        print FormatHexChunks(f, bits_per_chunk=64) + ', // e = {:5d}, k = {:4d}'.format(e, k)

# For double-precision:
PrintRyuPowers(bits=128, min_exponent=-291, max_exponent=325)

# For single-precision:
# PrintRyuPowers(bits=64, min_exponent=-30, max_exponent=47)
