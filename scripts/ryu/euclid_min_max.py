# Translated from:
# https://github.com/ulfjack/ryu/blob/c9c3fb19791c44fbe35701ad3b8ca4dc0977eb08/src/main/java/info/adams/ryu/analysis/EuclidMinMax.java

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

from fractions import gcd

def EuclidMin(multiplier, modulo, maximum):
    """
    Computes the modular min using a modified version of Euclid's algorithm.
    """

    c = gcd(multiplier, modulo)
    b = modulo // c
    a = (multiplier // c) % b
    if maximum >= b:
        return 0
    s = 1
    t = 0
    u = 0
    v = 1
    while True:
        while b >= a:
            b = b - a
            u = u - s
            v = v - t
            if -u >= maximum:
                return a * c
        if b == 0:
            return 0
        while a >= b:
            old_a = a
            a = a - b
            s = s - u
            t = t - v
            if s >= maximum:
                if s > maximum:
                    return old_a * c
                else:
                    return a * c
        if a == 0:
            return 0

def EuclidMax(multiplier, modulo, maximum):
    """
    Computes the modular max using a modified version of Euclid's algorithm.
    """

    c = gcd(multiplier, modulo)
    b = modulo // c
    a = (multiplier // c) % b
    if maximum >= b:
        return modulo - c
    s = 1
    t = 0
    u = 0
    v = 1
    while True:
        while b >= a:
            q = b // a
            q = min(q, (maximum - (-u)) // s - 1)
            q = max(q, 1)
            old_b = b
            b = b - a * q
            u = u - s * q
            v = v - t * q
            if -u >= maximum:
                if -u > maximum:
                    return modulo - old_b * c
                else:
                    return modulo - b * c
        if b == 0:
            return modulo - c
        while a >= b:
            q = 1
            if u != 0:
                q = min(q, (maximum - s) // (-u + 1))
            q = max(q, 1)
            a = a - b * q
            s = s - u * q
            t = t - v * q
            if s >= maximum:
                return modulo - b * c
        if a == 0:
            return modulo - c
