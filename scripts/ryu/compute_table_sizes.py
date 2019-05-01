# Translated from:
# https://github.com/ulfjack/ryu/blob/c9c3fb19791c44fbe35701ad3b8ca4dc0977eb08/src/main/java/info/adams/ryu/analysis/ComputeTableSize.java

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

def FloorLog10Pow2(e):
    assert e >= -28737
    assert e <=  28737
    return (e * 20201781) >> 26

def FloorLog10Pow5(e):
    assert e >= -28737
    assert e <=  28737
    return (e * 46907083) >> 26

def ComputeTableSize(exponent_bits, explicit_mantissa_bits):
    bias = 2**(exponent_bits - 1) - 1

    min_e2 = 1 - bias - explicit_mantissa_bits - 2
    max_e2 = (2**exponent_bits - 2) - bias - explicit_mantissa_bits - 2
    # print 'min_e2 = {}'.format(min_e2)
    # print 'max_e2 = {}'.format(max_e2)

    q_for_pos_e2 = max(0, FloorLog10Pow2(max_e2) - 1)
    q_for_neg_e2 = max(0, FloorLog10Pow5(-min_e2) - 1)

    min_pow = -q_for_pos_e2
    max_pow = -min_e2 - q_for_neg_e2

    return min_pow, max_pow

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def PrintTableSize(total_bits, exponent_bits, explicit_mantissa_bits):
    b = ComputeTableSize(exponent_bits, explicit_mantissa_bits)
    print 'float{}: [{}, {}]'.format(total_bits, *b)

PrintTableSize(16,  5, 10) # ==> [0, 9]
PrintTableSize(32,  8, 23) # ==> [-29, 47]
PrintTableSize(64, 11, 52) # ==> [-290, 325]
PrintTableSize(80, 15, 63) # ==> [-4911, 4953]
