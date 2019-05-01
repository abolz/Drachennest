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
    return (e * 1741647) >> 19

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

def DivUp(num, den):
    return (num + (den - 1)) // den

def CeilLog10Pow2(e):
    assert e >= -2620
    assert e <=  2620
    return (e * 315653 + (2**20 - 1)) >> 20;

def FloorLog10Pow2(e):
    assert e >= -2620
    assert e <=  2620
    return (e * 315653) >> 20

def ComputeBinaryExponentRange(q, p, exponent_bits):
    assert 0 <= p and p + 3 <= q
    bias = 2**(exponent_bits - 1) - 1
    min_exp = (1                    - bias) - (p - 1) - (p - 1) - (q - p)
    max_exp = (2**exponent_bits - 2 - bias) - (p - 1)           - (q - p)
    return min_exp, max_exp

def PrintGrisuPowersForExponentRange(alpha, gamma):
    assert alpha + 3 <= gamma

    # DiyFp precision q = 64
    # For IEEE double-precision p = 53, exponent_bits = 11
    e_min, e_max = ComputeBinaryExponentRange(q=64, p=53, exponent_bits=11)

    k_del = max(1, FloorLog10Pow2(gamma - alpha))
    assert k_del >= 1
    k_min = CeilLog10Pow2(alpha - e_max - 1)
    k_max = CeilLog10Pow2(alpha - e_min - 1)

    num_cached = DivUp(k_max - k_min, k_del) + 1
    k_min_cached = k_min;
    k_max_cached = k_min + k_del * (num_cached - 1)

    print 'constexpr int kAlpha = {:3d};'.format(alpha)
    print 'constexpr int kGamma = {:3d};'.format(gamma)
    print '// k_min = {:4d}'.format(k_min)
    print '// k_max = {:4d}'.format(k_max)
#   print '// k_min (max) = {}'.format(k_min + (k_del - 1))
    print ''
    print 'constexpr int kCachedPowersSize       = {:>4d};'.format(num_cached)
    print 'constexpr int kCachedPowersMinDecExp  = {:>4d};'.format(k_min_cached)
    print 'constexpr int kCachedPowersMaxDecExp  = {:>4d};'.format(k_max_cached)
    print 'constexpr int kCachedPowersDecExpStep = {:>4d};'.format(k_del)
    print ''
    print 'inline CachedPower GetCachedPower(int index)'
    print '{'
    print '    static constexpr uint64_t kSignificands[] = {'
    for k in range(k_min_cached, k_max_cached + 1, k_del):
        f, e = ComputeGrisuPower(k, 64)
        print '        ' + FormatHexChunks(f, bits_per_chunk=64) + ', // e = {:5d}, k = {:4d}'.format(e, k)
    print '    };'
    print ''
    print '    GRISU_ASSERT(index >= 0);'
    print '    GRISU_ASSERT(index < kCachedPowersSize);'
    print ''
    print '    const int k = kCachedPowersMinDecExp + index * kCachedPowersDecExpStep;'
    print '    const int e = FloorLog2Pow10(k) + 1 - {};'.format(64)
    print ''
    print '    return {kSignificands[index], e, k};'
    print '}'

PrintGrisuPowersForExponentRange(-60, -32)
# PrintGrisuPowersForExponentRange(0, 3)
# PrintGrisuPowersForExponentRange(-35, -32)
# PrintGrisuPowersForExponentRange(-32, 0)
# PrintGrisuPowersForExponentRange(-31, 0)
# PrintGrisuPowersForExponentRange(-60, -45)
# PrintGrisuPowersForExponentRange(-60, -51)
# PrintGrisuPowersForExponentRange(-60, -33)
# PrintGrisuPowersForExponentRange(-60, -46)
# PrintGrisuPowersForExponentRange(-57, -43)
# PrintGrisuPowersForExponentRange(-3, 0)