#===============================================================================
#
#===============================================================================

# def Sum(chunks, bits_per_chunk):
#     assert bits_per_chunk > 0
#     return sum( map(lambda (i, c): c * (2**bits_per_chunk)**i, enumerate(chunks)) )

def Split(n, bits_per_chunk):
    assert n >= 0
    assert bits_per_chunk > 0
    chunks = []
    while True:
        n, r = divmod(n, 2**bits_per_chunk)
        chunks.append(r)
        if n == 0:
            break
    return chunks

def ToHexString(n, bits):
    assert bits > 0
    p = (bits + (4 - 1)) // 4       # Round up to four bits per hexit
#   p = 2**((p - 1).bit_length())   # Round up to next power-of-2
    assert 4*p >= n.bit_length()
    return '0x{:0{}X}'.format(n, p)

def FormatHexChunks(n, bits_per_chunk = 64):
    chunks = Split(n, bits_per_chunk)
    s = ', '.join(map(lambda x: ToHexString(x, bits_per_chunk), reversed(chunks)))
    if len(chunks) > 1:
        s = '{' + s + '}'
    return s

#===============================================================================
# Grisu
#===============================================================================

def FloorLog2Pow10(e):
    assert e >= -1233
    assert e <= 1232
    return (e * 1741647) // 2**19

def RoundUp(num, den):
    assert num >= 0
    assert den > 0
    p, r = divmod(num, den)
    if 2 * r >= den:
        p += 1
    return p

def ComputeGrisuPower(k, bits):
    assert bits > 0
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
    assert e >= -1764
    assert e <= 1763
    return (e * 1217359) // 2**19

def Ceil(num, den):
    assert num >= 0
    assert den > 0
    q, r = divmod(num, den)
    if r > 0:
        q += 1
    return q

def ComputeRyuPower(k, bits):
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
    return f, e

def PrintRyuPowers(bits, min_exponent, max_exponent):
    print '// Let e = FloorLog2Pow5(k) + 1 - {}'.format(bits)
    print '// For k >= 0, stores 5^k in the form: floor( 5^k / 2^e )'
    print '// For k <= 0, stores 5^k in the form:  ceil(2^-e / 5^-k)'
    for k in range(min_exponent, max_exponent + 1):
        f, e = ComputeRyuPower(k, bits)
        print FormatHexChunks(f, bits_per_chunk=64) + ', // e = {:5d}, k = {:4d}'.format(e, k)

# For double-precision:
# PrintRyuPowers(bits=128, min_exponent=-291, max_exponent=325)

# For single-precision:
PrintRyuPowers(bits=64, min_exponent=-30, max_exponent=47)
