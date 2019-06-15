#===================================================================================================
#
#===================================================================================================

class LogApprox:
    def __init__(self, b = 0, B = 0, mul = 0, shift = 0):
        self.b = b
        self.B = B
        self.mul = mul
        self.shift = shift

    def FloorOverflow(self, e, min_int, max_int):
        return e * self.mul > max_int or e * self.mul < min_int

    def Floor(self, e):
        return (e * self.mul) >> self.shift

    def IsFloor(self, k, e):
        # b^k <= B^E < b^(k+1)
        x = 1
        y = 1
        if k >= 0:
            x *= self.b**k
        else:
            y *= self.b**(-k)
        if e >= 0:
            y *= self.B**e
        else:
            x *= self.B**(-e)
        return x <= y and y < x * self.b

    def TestFloor(self, e, min_int, max_int):
        if self.FloorOverflow(e, min_int, max_int):
            # print "floor: Integer overflow at e = {}".format(e)
            return False
        k = self.Floor(e)
        if k < min_int or k > max_int:
            # print "floor: Result is out of range"
            return False
        if not self.IsFloor(k, e):
            # print "floor: test failed: k = {:6d}, b = {:3d}, B = {:3d}, e = {:6d}".format(k, self.b, self.B, e)
            return False
        return True

    def CeilOverflow(self, e, min_int, max_int):
        return e * self.mul + (2**self.shift - 1) > max_int or e * self.mul < min_int

    def Ceil(self, e):
        return (e * self.mul + (2**self.shift - 1)) >> self.shift

    def IsCeil(self, k, e):
        # b^(k-1) < B^E < b^k   <=>   b^k < B^E * b < b^k * b
        x = 1
        y = self.b
        if k >= 0:
            x *= self.b**k
        else:
            y *= self.b**(-k)
        if e >= 0:
            y *= self.B**e
        else:
            x *= self.B**(-e)
        return x < y and y <= x * self.b

    def TestCeil(self, e, min_int, max_int):
        if self.CeilOverflow(e, min_int, max_int):
            # print "ceil: Integer overflow at: e = {}".format(e)
            return False
        k = self.Ceil(e)
        if k < min_int or k > max_int:
            # print "ceil: Result is out of range"
            return False
        if not self.IsCeil(k, e):
            # print "ceil: test failed: k = {:6d}, b = {:3d}, B = {:3d}, e = {:6d}".format(k, self.b, self.B, e)
            return False
        return True

#===================================================================================================
#
#===================================================================================================

Log10_2 = 0.301029995663981198017467022509663365781307220458984375
Log10_5 = 0.69897000433601885749368420874816365540027618408203125
Log2_10 = 3.321928094887362181708567732130177319049835205078125
Log2_5  = 2.321928094887362181708567732130177319049835205078125
Log5_2  = 0.430676558073393056513822330089169554412364959716796875
Log5_10 = 1.4306765580733931120249735613469965755939483642578125

def FindFloorApprox(log_approx, b, B, int_width, min_exponent, max_exponent, signed = True):
    assert signed or min_exponent >= 0
    assert signed or max_exponent >= 0
    assert min_exponent <= max_exponent

    print('Find approximation for floor(log_{} {}^e)'.format(b, B))
    # print '  in [{}, {}]'.format(min_exponent, max_exponent)
    # print 'Integer bit width: {}'.format(int_width)

    # for bits in range(int_width // 2, int_width):
    for bits in range(0, int_width):
        # print '*** Bits = {}'.format(bits)

        L = LogApprox()
        L.b = b
        L.B = B
        L.mul = int( "{:.0f}".format(log_approx * 2**bits) )
        L.shift = bits

        if signed:
            min_int = -2**(int_width - 1)
            max_int =  2**(int_width - 1) - 1
        else:
            min_int = 0
            max_int = 2**int_width - 1

        max_e = 0
        min_e = 0
        if signed:
            for e in range(0, -min_exponent + 1):
                if not L.TestFloor(-e, min_int, max_int):
                    break
                min_e = -e
        for e in range(0, max_exponent + 1):
            if not L.TestFloor(e, min_int, max_int):
                break
            max_e = e

        if min_e <= min_exponent and max_exponent <= max_e:
            max_target_exponent = max_int # 262380
            if True:
                print('Found.')
                print('Computing exact valid range of target exponent...')
                real_min_e = -min_e
                if signed:
                    while True:
                        next_e = real_min_e + 1
                        if next_e > max_target_exponent:
                            break
                        if not L.TestFloor(-next_e, min_int, max_int):
                            break
                        real_min_e = next_e
                real_min_e = -real_min_e
                real_max_e = max_e
                while True:
                    next_e = real_max_e + 1
                    if next_e > max_target_exponent:
                        break
                    if not L.TestFloor(next_e, min_int, max_int):
                        break
                    real_max_e = next_e
            else:
                real_min_e = min_e
                real_max_e = max_e

            if signed:
                print('inline int{}_t FloorLog{}Pow{}(int{}_t e) {{'.format(int_width, b, B, int_width))
                print('    assert(e >= {:d});'.format(real_min_e))
                print('    assert(e <= {:d});'.format(real_max_e))
                print('    return (e * {}) >> {};'.format(L.mul, L.shift))
                print('}')
            else:
                print('inline uint{}_t FloorLog{}Pow{}(uint{}_t e) {{'.format(int_width, b, B, int_width))
                print('    assert(e <= {:d});'.format(real_max_e))
                print('    return (e * {}) >> {};'.format(L.mul, L.shift))
                print('}')
            return

    if int_width < 64:
        return FindFloorApprox(log_approx, b, B, int_width * 2, min_exponent, max_exponent)

    print('FAIL.')

INITIAL_INT_WIDTH = 32

def FindFloorApproxForPrec(float_width, exponent_bits, precision):
    print('Finding logs for float{}'.format(float_width))
    bias = 2**(exponent_bits - 1) - 1 + (precision - 1)
    min_exponent = 1 - bias
    max_exponent = 2**exponent_bits - 2 - bias
    min_exponent -= 2
    max_exponent -= 2
    print('min_exponent = {}'.format(min_exponent))
    print('max_exponent = {}'.format(max_exponent))
    FindFloorApprox(Log10_2, 10, 2, INITIAL_INT_WIDTH, 0, max_exponent)
    FindFloorApprox(Log10_5, 10, 5, INITIAL_INT_WIDTH, 0, -min_exponent)

# float16:
# FindFloorApproxForPrec(16, 5, 10)
# FindFloorApprox(Log2_5, 2, 5, INITIAL_INT_WIDTH, 0, 9, signed=True)

# float32:
# FindFloorApproxForPrec(32, 8, 24)
# FindFloorApprox(Log2_5, 2, 5, INITIAL_INT_WIDTH, -29, 47, signed=True)

# float64:
FindFloorApproxForPrec(64, 11, 53)
FindFloorApprox(Log2_5, 2, 5, INITIAL_INT_WIDTH, -290, 325, signed=True)

# float80:
# FindFloorApproxForPrec(80, 15, 64)
# FindFloorApprox(Log2_5, 2, 5, INITIAL_INT_WIDTH, -4911, 4953, signed=True)

#===================================================================================================
# Finding logs for float16
# min_exponent = -25
# max_exponent = 4
# Find approximation for floor(log_10 2^e)
# Found.
# Computing exact valid range of target exponent...
# inline int32_t FloorLog10Pow2(int32_t e) {
#     assert(e >= -3);
#     assert(e <= 6);
#     return (e * 1) >> 2;
# }
# Find approximation for floor(log_10 5^e)
# Found.
# Computing exact valid range of target exponent...
# inline int32_t FloorLog10Pow5(int32_t e) {
#     assert(e >= -102);
#     assert(e <= 102);
#     return (e * 179) >> 8;
# }
# Find approximation for floor(log_2 5^e)
# Found.
# Computing exact valid range of target exponent...
# inline int32_t FloorLog2Pow5(int32_t e) {
#     assert(e >= -15);
#     assert(e <= 18);
#     return (e * 37) >> 4;
# }

#===================================================================================================
# Finding logs for float32
# min_exponent = -151
# max_exponent = 102
# Find approximation for floor(log_10 2^e)
# Found.
# Computing exact valid range of target exponent...
# inline int32_t FloorLog10Pow2(int32_t e) {
#     assert(e >= -102);
#     assert(e <= 102);
#     return (e * 77) >> 8;
# }
# Find approximation for floor(log_10 5^e)
# Found.
# Computing exact valid range of target exponent...
# inline int32_t FloorLog10Pow5(int32_t e) {
#     assert(e >= -680);
#     assert(e <= 680);
#     return (e * 2863) >> 12;
# }
# Find approximation for floor(log_2 5^e)
# Found.
# Computing exact valid range of target exponent...
# inline int32_t FloorLog2Pow5(int32_t e) {
#     assert(e >= -58);
#     assert(e <= 58);
#     return (e * 1189) >> 9;
# }

#===================================================================================================
# Finding logs for float64
# min_exponent = -1076
# max_exponent = 969
# Find approximation for floor(log_10 2^e)
# Found.
# Computing exact valid range of target exponent...
# inline int32_t FloorLog10Pow2(int32_t e) {
#     assert(e >= -1650);
#     assert(e <= 1650);
#     return (e * 78913) >> 18;
# }
# Find approximation for floor(log_10 5^e)
# Found.
# Computing exact valid range of target exponent...
# inline int32_t FloorLog10Pow5(int32_t e) {
#     assert(e >= -1650);
#     assert(e <= 1650);
#     return (e * 183231) >> 18;
# }
# Find approximation for floor(log_2 5^e)
# Found.
# Computing exact valid range of target exponent...
# inline int32_t FloorLog2Pow5(int32_t e) {
#     assert(e >= -642);
#     assert(e <= 642);
#     return (e * 76085) >> 15;
# }

#===================================================================================================
# Finding logs for float80
# min_exponent = -16447
# max_exponent = 16318
# Find approximation for floor(log_10 2^e)
# Find approximation for floor(log_10 2^e)
# Found.
# Computing exact valid range of target exponent...
# inline int64_t FloorLog10Pow2(int64_t e) {
#     assert(e >= -28737);
#     assert(e <= 28737);
#     return (e * 20201781) >> 26;
# }
# Find approximation for floor(log_10 5^e)
# Find approximation for floor(log_10 5^e)
# Found.
# Computing exact valid range of target exponent...
# inline int64_t FloorLog10Pow5(int64_t e) {
#     assert(e >= -28737);
#     assert(e <= 28737);
#     return (e * 46907083) >> 26;
# }
# Find approximation for floor(log_2 5^e)
# Find approximation for floor(log_2 5^e)
# Found.
# Computing exact valid range of target exponent...
# inline int64_t FloorLog2Pow5(int64_t e) {
#     assert(e >= -12654);
#     assert(e <= 12654);
#     return (e * 38955489) >> 24;
# }
