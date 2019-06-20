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

def PrintTable(base, max_exp, bits):
    assert bits % 4 == 0
    print('{{0x{:0{}X}u, 0x{:0{}X}u}}, // {}^{}'.format(1, bits // 4, 2**bits - 1, bits // 4, base, 0))
    e = 1
    while True:
        a = base**e
        if a >= 2**bits:
            break
        m = ModularInverse(a, 2**bits)
        assert (a * m) % 2**bits == 1
        print('{{0x{:0{}X}u, 0x{:0{}X}u}}, // {}^{}'.format(m, bits // 4, 2**bits // a, bits // 4, base, e))
        e += 1
    print('')

PrintTable(base=5, max_exp=10, bits=32) # 10 = floor(log_5(2^26))
PrintTable(base=5, max_exp=22, bits=64) # 22 = floor(log_5(2^55))
