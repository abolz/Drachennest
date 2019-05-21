def ExtendedGCD(a, b):
    if a == 0:
        return b, 0, 1
    else:
        g, y, x = ExtendedGCD(b % a, a)
        return g, x - (b // a) * y, y

def ModularInverse(a, m):
    assert a >= 0
    assert m >= 1
    g, x, _ = ExtendedGCD(a, m)
    assert g == 1
    return x % m

def PrintTable(max_e5, bits):
    assert bits % 4 == 0
    print '    {{0x{:0{}X}u, 0x{:0{}X}u}},'.format(1, bits // 4, 2**bits - 1, bits // 4)
    for e in range(1, max_e5 + 1):
        a = 5**e
        m = ModularInverse(a, 2**bits)
        assert (a * m) % 2**bits == 1
        print '    {{0x{:0{}X}u, 0x{:0{}X}u}},'.format(m, bits // 4, 2**bits // a, bits // 4)
    print ''

PrintTable(max_e5=10, bits=32) # 10 = floor(log_5(2^26))
PrintTable(max_e5=22, bits=64) # 22 = floor(log_5(2^55))
