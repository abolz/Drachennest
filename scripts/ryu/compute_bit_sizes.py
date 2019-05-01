# Translated from:
# https://github.com/ulfjack/ryu/blob/c9c3fb19791c44fbe35701ad3b8ca4dc0977eb08/src/main/java/info/adams/ryu/analysis/ComputeRequiredBitSizes.java

# Copyright 2018 Ulf Adams
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from euclid_min_max import EuclidMin, EuclidMax

def FloorLog10Pow2(e):
    assert e >= -28737
    assert e <=  28737
    return (e * 20201781) >> 26

def FloorLog10Pow5(e):
    assert e >= -28737
    assert e <=  28737
    return (e * 46907083) >> 26

def ComputeRequiredBitSizes(exponent_bits, explicit_mantissa_bits):
    """
    Computes appropriate values for B_0 and B_1 for a given floating point type.
    """

    bias = 2**(exponent_bits - 1) - 1
    mbits = explicit_mantissa_bits + 3
    # max(w) = 4 * ((1 << format.explicit_mantissa_bits()) * 2 - 1) + 2
    #        = (1 << (format.explicit_mantissa_bits() + 3)) - 2
    max_w = 2**mbits - 2

    #
    # For e2 >= 0:
    #
    min_e2 = 0
    max_e2 = (2**exponent_bits - 2) - bias - explicit_mantissa_bits - 2
    b0 = 0
    for e2 in range(min_e2, max_e2 + 1):
        q = max(0, FloorLog10Pow2(e2) - 1)
        pow5 = 5**q
        pow2 = 2**(e2 - q)
        euclid_max = EuclidMax(pow2, pow5, max_w - 1)
        bits = ((max_w * pow5 * pow2) // (pow5 - euclid_max)).bit_length()
        reqn = bits - pow5.bit_length()
        b0 = max(b0, reqn)

    #
    # For e2 < 0:
    #
    min_e2 = 0
    max_e2 = -(1 - bias - explicit_mantissa_bits - 2)
    b1 = 0
    for e2 in range(min_e2, max_e2 + 1):
        q = max(0, FloorLog10Pow5(e2) - 1)
        pow5 = 5**(e2 - q)
        pow2 = 2**q
        euclid_min = EuclidMin(pow5, pow2, max_w - 1)
        bits = (euclid_min // max_w).bit_length()
        reqn = pow5.bit_length() - bits
        b1 = max(b1, reqn)

    return b0, b1

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def PrintRequiredBitSizes(total_bits, exponent_bits, explicit_mantissa_bits):
    b = ComputeRequiredBitSizes(exponent_bits, explicit_mantissa_bits)
    print 'float{}: [{}, {}]'.format(total_bits, *b)

PrintRequiredBitSizes(16,  5, 10) # ==> [15, 21]
PrintRequiredBitSizes(32,  8, 23) # ==> [60, 63]
PrintRequiredBitSizes(64, 11, 52) # ==> [124, 124]
PrintRequiredBitSizes(80, 15, 63) # ==> [150, 152]
