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
    return('0x{:0{}X}'.format(n, p))

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
    print('// Let e = FloorLog2Pow10(k) + 1 - {}'.format(bits))
    print('// For k >= 0, stores 10^k in the form: round_up(10^k / 2^e )')
    print('// For k <= 0, stores 10^k in the form: round_up(2^-e / 10^-k)')
    for k in range(min_exponent, max_exponent + 1, step):
        f, e = ComputeGrisuPower(k, bits)
        print(FormatHexChunks(f, bits_per_chunk=64) + ', // e = {:5d}, k = {:4d}'.format(e, k))

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

def PrintGrisuPowersForExponentRange(alpha, gamma, q = 64, p = 53, exponent_bits = 11):
    assert alpha + 3 <= gamma

    # DiyFp precision q = 64
    # For IEEE double-precision p = 53, exponent_bits = 11
    # e_min, e_max = ComputeBinaryExponentRange(q=64, p=53, exponent_bits=11)
    # e_min, e_max = ComputeBinaryExponentRange(q=32, p=24, exponent_bits=8)
    e_min, e_max = ComputeBinaryExponentRange(q, p, exponent_bits)

    k_del = max(1, FloorLog10Pow2(gamma - alpha))
    # k_del = 1
    assert k_del >= 1
    k_min = CeilLog10Pow2(alpha - e_max - 1)
    # k_min += 7
    k_max = CeilLog10Pow2(alpha - e_min - 1)

    num_cached = DivUp(k_max - k_min, k_del) + 1
    k_min_cached = k_min;
    k_max_cached = k_min + k_del * (num_cached - 1)

    print('constexpr int kAlpha = {:3d};'.format(alpha))
    print('constexpr int kGamma = {:3d};'.format(gamma))
    print('// k_min = {:4d}'.format(k_min))
    print('// k_max = {:4d}'.format(k_max))
#   print('// k_del = {:4d}'.format(k_del))
#   print('// k_min (max) = {}'.format(k_min + (k_del - 1)))
    print('')
    print('constexpr int kCachedPowersSize       = {:>4d};'.format(num_cached))
    print('constexpr int kCachedPowersMinDecExp  = {:>4d};'.format(k_min_cached))
    print('constexpr int kCachedPowersMaxDecExp  = {:>4d};'.format(k_max_cached))
    print('constexpr int kCachedPowersDecExpStep = {:>4d};'.format(k_del))
    print('')

    # print('inline CachedPower GetCachedPower(int index)')
    # print('{')
    # print('    static constexpr uint{}_t kSignificands[] = {{'.format(q))
    # for k in range(k_min_cached, k_max_cached + 1, k_del):
    #     f, e = ComputeGrisuPower(k, q)
    #     print('        ' + FormatHexChunks(f, q) + ', // e = {:5d}, k = {:4d}'.format(e, k))
    # print('    };')
    # print('')
    # print('    GRISU_ASSERT(index >= 0);')
    # print('    GRISU_ASSERT(index < kCachedPowersSize);')
    # print('')
    # print('    const int k = kCachedPowersMinDecExp + index * kCachedPowersDecExpStep;')
    # print('    const int e = FloorLog2Pow10(k) + 1 - {};'.format(q))
    # print('')
    # print('    return {kSignificands[index], e, k};')
    # print('}')

    print('// For a normalized DiyFp w = f * 2^e, this function returns a (normalized)')
    print('// cached power-of-ten c = f_c * 2^e_c, such that the exponent of the product')
    print('// w * c satisfies')
    print('//')
    print('//      kAlpha <= e_c + e + q <= kGamma.')
    print('//')
    print('inline CachedPower GetCachedPowerForBinaryExponent(int e)')
    print('{')
    print('    static constexpr uint{}_t kSignificands[] = {{'.format(q))
    for k in range(k_min_cached, k_max_cached + 1, k_del):
        f, e = ComputeGrisuPower(k, q)
        print('        ' + FormatHexChunks(f, q) + ', // e = {:5d}, k = {:4d}'.format(e, k))
    print('    };')
    print('')
    print('    GRISU_ASSERT(e >= {:>5d});'.format(e_min))
    print('    GRISU_ASSERT(e <= {:>5d});'.format(e_max))
    print('')
    print('    const int k = CeilLog10Pow2(kAlpha - e - 1);')
    print('    GRISU_ASSERT(k >= kCachedPowersMinDecExp - (kCachedPowersDecExpStep - 1));')
    print('    GRISU_ASSERT(k <= kCachedPowersMaxDecExp);')
    print('')
    print('    const unsigned index = static_cast<unsigned>(k - (kCachedPowersMinDecExp - (kCachedPowersDecExpStep - 1))) / kCachedPowersDecExpStep;')
    print('    GRISU_ASSERT(index < kCachedPowersSize);')
    print('')
    print('    const int k_cached = kCachedPowersMinDecExp + static_cast<int>(index) * kCachedPowersDecExpStep;')
    print('    const int e_cached = FloorLog2Pow10(k_cached) + 1 - {};'.format(q))
    print('')
    print('    const CachedPower cached = {kSignificands[index], e_cached, k_cached};')
    print('    GRISU_ASSERT(kAlpha <= cached.e + e + {});'.format(q))
    print('    GRISU_ASSERT(kGamma >= cached.e + e + {});'.format(q))
    print('')
    print('    return cached;')
    print('}')

# PrintGrisuPowersForExponentRange(-60, -32, q=64, p=53, exponent_bits=11)
# PrintGrisuPowersForExponentRange(-59, -32, q=64, p=53, exponent_bits=11)
# PrintGrisuPowersForExponentRange(-56, -42, q=64, p=53, exponent_bits=11)
# PrintGrisuPowersForExponentRange(-3, 0, q=64, p=53, exponent_bits=11)
# PrintGrisuPowersForExponentRange(-28, 0, q=64, p=53, exponent_bits=11)
# PrintGrisuPowersForExponentRange(-53, -46, q=64, p=53, exponent_bits=11)
PrintGrisuPowersForExponentRange(-50, -36, q=64, p=53, exponent_bits=11)
# PrintGrisuPowersForExponentRange(-50, -41, q=64, p=53, exponent_bits=11)
# PrintGrisuPowersForExponentRange(-50, -44, q=64, p=53, exponent_bits=11)
# PrintGrisuPowersForExponentRange(-50, -47, q=64, p=53, exponent_bits=11)
# PrintGrisuPowersForExponentRange(-50, -36, q=64, p=53, exponent_bits=11)
# PrintGrisuPowersForExponentRange(-50, -41, q=64, p=53, exponent_bits=11)

# PrintGrisuPowersForExponentRange(0, 3, q=32, p=24, exponent_bits= 8)
# PrintGrisuPowersForExponentRange(0, 3, q=64, p=53, exponent_bits=11)
# PrintGrisuPowersForExponentRange(-25, -18, q=32, p=24, exponent_bits= 8)
