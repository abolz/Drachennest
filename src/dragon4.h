// Copyright 2019 Alexander Bolz
//
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cassert>
#include <cstdint>
#include <cstring>
#include <limits>
#ifdef _MSC_VER
#include <intrin.h>
#endif

#ifndef DRAGON4_ASSERT
#define DRAGON4_ASSERT(X) assert(X)
#endif

namespace dragon4 {

//==================================================================================================
// Dragon4
//
// Implements the Dragon4 algorithm for (IEEE) binary to decimal floating-point conversion.
//
// References:
//
// [1]  Burger, Dybvig, "Printing Floating-Point Numbers Quickly and Accurately",
//      Proceedings of the ACM SIGPLAN 1996 Conference on Programming Language Design and Implementation, PLDI 1996
// [2]  Steele, White, "How to Print FloatingPoint Numbers Accurately",
//      Proceedings of the ACM SIGPLAN 1990 conference on Programming language design and implementation, PLDI 1990
//==================================================================================================
// Constant data: 56 bytes

namespace impl {

inline int Min(int x, int y) { return y < x ? y : x; }
inline int Max(int x, int y) { return y < x ? x : y; }

inline void Zeroize(uint32_t* dst, int n)
{
    DRAGON4_ASSERT(n >= 0);
    std::memset(dst, 0, static_cast<unsigned>(n) * sizeof(uint32_t));
}

inline void Memmove(uint32_t* dst, const uint32_t* src, int n)
{
    DRAGON4_ASSERT(n >= 0);
    std::memmove(dst, src, static_cast<unsigned>(n) * sizeof(uint32_t));
}

//template <int MaxBits>
struct DiyInt
{
    static constexpr int MaxBits = 1130;
    static constexpr int Capacity = (MaxBits + (32 - 1)) / 32;

    uint32_t bigits[Capacity]; // Significand stored in little-endian form.
    int      size = 0;

    DiyInt() = default;
    DiyInt(const DiyInt&) = delete;             // (not needed here)
    DiyInt& operator=(const DiyInt&) = delete;  // (not needed here)
};

inline void AssignU32(DiyInt& x, uint32_t value)
{
    x.bigits[0] = value;
    x.size = (value != 0) ? 1 : 0;
}

inline void AssignU64(DiyInt& x, uint64_t value)
{
    x.bigits[0] = static_cast<uint32_t>(value);
    x.bigits[1] = static_cast<uint32_t>(value >> 32);
    x.size = (x.bigits[1] != 0) ? 2 : ((x.bigits[0] != 0) ? 1 : 0);
}

// x := A * x + B
inline void MulAddU32(DiyInt& x, uint32_t A, uint32_t B = 0)
{
    DRAGON4_ASSERT(x.size >= 0);

    if (A == 1 && B == 0)
    {
        return;
    }
    if (A == 0 || x.size <= 0)
    {
        AssignU32(x, B);
        return;
    }

    uint32_t carry = B;
    for (int i = 0; i < x.size; ++i)
    {
        const uint64_t p = uint64_t{x.bigits[i]} * A + carry;
        x.bigits[i]      = static_cast<uint32_t>(p);
        carry            = static_cast<uint32_t>(p >> 32);
    }

    if (carry != 0)
    {
        DRAGON4_ASSERT(x.size < DiyInt::Capacity);
        x.bigits[x.size++] = carry;
    }
}

// x := x * 2^e2
inline void MulPow2(DiyInt& x, int e2) // aka left-shift
{
    DRAGON4_ASSERT(x.size >= 0);
    DRAGON4_ASSERT(e2 >= 0);

    if (x.size <= 0 || e2 == 0)
        return;

    const int bigit_shift = static_cast<int>(static_cast<unsigned>(e2) / 32);
    const int bit_shift   = static_cast<int>(static_cast<unsigned>(e2) % 32);

    if (bit_shift > 0)
    {
        uint32_t carry = 0;
        for (int i = 0; i < x.size; ++i)
        {
            const uint32_t h = x.bigits[i] >> (32 - bit_shift);
            x.bigits[i]      = x.bigits[i] << bit_shift | carry;
            carry            = h;
        }

        if (carry != 0)
        {
            DRAGON4_ASSERT(x.size < DiyInt::Capacity);
            x.bigits[x.size++] = carry;
        }
    }

    if (bigit_shift > 0)
    {
        DRAGON4_ASSERT(bigit_shift <= DiyInt::Capacity);
        DRAGON4_ASSERT(x.size <= DiyInt::Capacity - bigit_shift);

        Memmove(x.bigits + bigit_shift, x.bigits, x.size);
        Zeroize(x.bigits, bigit_shift);
        x.size += bigit_shift;
    }
}

// x := x * 5^e5
inline void MulPow5(DiyInt& x, int e5)
{
    // TODO:
    // Optimize for large powers???

    static constexpr uint32_t kPow5_32[] = {
        1, // (unused)
        5,
        25,
        125,
        625,
        3125,
        15625,
        78125,
        390625,
        1953125,
        9765625,
        48828125,
        244140625,
        1220703125, // 5^13
    };

    DRAGON4_ASSERT(x.size >= 0);
    DRAGON4_ASSERT(e5 >= 0);

    if (x.size <= 0)
        return;

    while (e5 > 0)
    {
        const int n = Min(e5, 13);
        MulAddU32(x, kPow5_32[n]);
        e5 -= n;
    }
}

inline void Mul2(DiyInt& x)
{
    uint32_t carry = 0;
    for (int i = 0; i < x.size; ++i)
    {
        const uint32_t h = x.bigits[i] >> 31;
        x.bigits[i]      = x.bigits[i] << 1 | carry;
        carry            = h;
    }
    if (carry != 0)
    {
        DRAGON4_ASSERT(x.size < DiyInt::Capacity);
        x.bigits[x.size++] = carry;
    }
}

inline void Mul10(DiyInt& x)
{
#if 0
    MulAddU32(x, 10, 0);
#else
    uint32_t carry = 0;
    for (int i = 0; i < x.size; ++i)
    {
        const uint64_t p = uint64_t{x.bigits[i]} * 10 + carry;
        x.bigits[i]      = static_cast<uint32_t>(p);
        carry            = static_cast<uint32_t>(p >> 32);
    }
    if (carry != 0)
    {
        DRAGON4_ASSERT(x.size < DiyInt::Capacity);
        x.bigits[x.size++] = carry;
    }
#endif
}

// x := 2^e2
inline void AssignPow2(DiyInt& x, int e2)
{
    DRAGON4_ASSERT(e2 >= 0);

    if (e2 == 0)
    {
        AssignU32(x, 1);
        return;
    }

    const int bigit_shift = static_cast<int>(static_cast<unsigned>(e2)) / 32;
    const int bit_shift   = static_cast<int>(static_cast<unsigned>(e2)) % 32;

    Zeroize(x.bigits, bigit_shift);
    x.bigits[bigit_shift] = 1u << bit_shift;
    x.size = bigit_shift + 1;
}

// x := 5^e5
inline void AssignPow5(DiyInt& x, int e5)
{
    // TODO:
    // Optimize?!

    AssignU32(x, 1);
    MulPow5(x, e5);
}

// x := 10^e10
inline void AssignPow10(DiyInt& x, int e10)
{
    // TODO:
    // Optimize?!

    AssignPow5(x, e10);
    MulPow2(x, e10);
}

// x := value * 2^e2
inline void AssignU64MulPow2(DiyInt& x, uint64_t value, int e2)
{
    DRAGON4_ASSERT(e2 >= 0);

    if (value == 0 || e2 == 0)
    {
        AssignU64(x, value);
        return;
    }

    const int bigit_shift = static_cast<int>(static_cast<unsigned>(e2)) / 32;
    const int bit_shift   = static_cast<int>(static_cast<unsigned>(e2)) % 32;

    Zeroize(x.bigits, bigit_shift);

    const uint32_t lo = static_cast<uint32_t>(value);
    const uint32_t hi = static_cast<uint32_t>(value >> 32);
    if (bit_shift == 0)
    {
        // Relax: only write to x.bigits if neccessary?!
        DRAGON4_ASSERT(DiyInt::Capacity >= bigit_shift + 2);

        x.bigits[bigit_shift + 0] = lo;
        x.bigits[bigit_shift + 1] = hi;
        x.size = bigit_shift + ((hi != 0) ? 2 : 1);
    }
    else
    {
        // Relax: only write to x.bigits if neccessary?!
        DRAGON4_ASSERT(DiyInt::Capacity >= bigit_shift + 3);

        const uint32_t v0 = lo << bit_shift;
        const uint32_t v1 = hi << bit_shift | lo >> (32 - bit_shift);
        const uint32_t v2 =                   hi >> (32 - bit_shift);
        x.bigits[bigit_shift + 0] = v0;
        x.bigits[bigit_shift + 1] = v1;
        x.bigits[bigit_shift + 2] = v2;
        x.size = bigit_shift + ((v2 != 0) ? 3 : ((v1 != 0) ? 2 : 1));
    }
}

// x := value * 10^e10
inline void AssignU64MulPow10(DiyInt& x, uint64_t value, int e10)
{
    AssignU64MulPow2(x, value, e10);
    MulPow5(x, e10);
}

// x := 2^e2 * 5^e5
inline void AssignPow2MulPow5(DiyInt& x, int e2, int e5)
{
    AssignPow2(x, e2);
    MulPow5(x, e5);

    //AssignPow5(x, e5);
    //MulPow2(x, e2);
}

// Returns the number of leading 0-bits in x, starting at the most significant bit position.
// If x is 0, the result is undefined.
inline int CountLeadingZeros32(uint32_t x)
{
    DRAGON4_ASSERT(x != 0);

#if defined(__GNUC__) || defined(__clang__)
    return __builtin_clz(x);
#elif defined(_MSC_VER) && (defined(_M_ARM) || defined(_M_ARM64))
    return static_cast<int>(_CountLeadingZeros(x));
#elif defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
    return static_cast<int>(__lzcnt(x));
#else
    int lz = 0;
    while ((x >> 31) == 0) {
        x <<= 1;
        ++lz;
    }
    return lz;
#endif
}

inline uint32_t DivModShort(DiyInt& u, uint32_t v)
{
    const int m = u.size;

    uint32_t q = 0;
    uint32_t r = 0;
    for (int i = m - 1; i >= 0; --i)
    {
        const uint64_t t = (uint64_t{r} << 32) | u.bigits[i];
        q = static_cast<uint32_t>(t / v);
        r = static_cast<uint32_t>(t % v);
    }
    AssignU32(u, r);
    return q;
}

// q, r = divmod(u, v)
// u := r
// return q
// PRE: 0 <= q <= 9
inline uint32_t DivMod(DiyInt& u, const DiyInt& v)
{
    DRAGON4_ASSERT(u.size > 0);
    DRAGON4_ASSERT(v.size > 0);
    DRAGON4_ASSERT(u.bigits[u.size - 1] != 0);
    DRAGON4_ASSERT(v.bigits[v.size - 1] != 0);

    const int m = u.size;
    const int n = v.size;
    if (m < n)
    {
        return 0;
    }

    DRAGON4_ASSERT(m >= n);
    DRAGON4_ASSERT(n >= 1);

    //--------------------------------------------------------------------------
    // D0.
    //
    // Handle the case of a single digit division first.
    // Note that this step is not only here for performance reasons. The
    // algorithm described below requires at least two digits in the
    // denominator.
    if (n == 1)
    {
        return DivModShort(u, v.bigits[0]);
    }

    DRAGON4_ASSERT(n >= 2);
    DRAGON4_ASSERT(DiyInt::Capacity >= m + 1);
    u.bigits[m] = 0;

    //
    // XXX:
    //
    // Is the normalization step only required once in Dragon4 before starting
    // the digit generation procedure?!?!?!
    // It might be more efficient to shift r and s and delta instead of
    // shifting the leading bits in each iteration...
    //

    //--------------------------------------------------------------------------
    // D1. [Normalize.]
    //
    // Set d := b / (v[n-1] + 1). Then set
    //    u' := (u[0] u[1] ... u[m-1] u[m])_b = d * (u[0] u[1] ... u[m-1])_b,
    //    v' := (v[0] v[1] ... v[n-1]     )_b = d * (v[0] v[1] ... v[n-1])_b.
    //
    // Note the introduction of a new digit position u[m] at the right of
    // u[m-1]; if d = 1, all we need to do in this step is set u[m] := 0.
    //
    // On a binary computer it may be preferable to choose d to be a power of 2
    // instead of using the value suggested here; any value of d that results in
    // v[n-1] >= b/2 will suffice.
    //

    // This normalization step is required only to efficiently estimate the
    // quotient q' (see below). Is is not necessary for the other steps of the
    // algorithm.
    // Instead of shifting both u and v into u' and v' resp., the required
    // digits of u' and v' are computed when they are needed.
    //
    // The variables vK here denote v'[n - K], where K = 1, 2, and v' denotes
    // the normalized value d * v.

    uint32_t v1 = v.bigits[n - 1];
    uint32_t v2 = v.bigits[n - 2];

    const int shift = CountLeadingZeros32(v1);
    if (shift > 0)
    {
        const uint32_t v3 = (n >= 3) ? v.bigits[n - 3] : 0;
        v1 = (v1 << shift) | (v2 >> (32 - shift));
        v2 = (v2 << shift) | (v3 >> (32 - shift));
    }
    // v1 and v2 now contain the leading digits of v'.

    //--------------------------------------------------------------------------
    // D2. [Initialize.]
    //
    // Set j := m - n.
    //
    // The loop on j, steps D2 through D7, will be essentially a division of
    // (u[j] u[j+1] ... u[j+n])_b by (v[0] v[1] ... v[n-1])_b to get a single
    // quotient digit.
    //

    //--------------------------------------------------------------------------
    // D3. [Calculate q'.]
    //
    // If u[j+n] = v[n-1], set
    //    q' := b - 1;
    // otherwise set
    //    q' := (u[j+n] * b + u[j+n-1]) / v[n-1].
    // Now test if
    //    q' * v[n-2] > (u[j+n] * b + u[j+n-1] - q' * v[n-1]) * b + u[j+n-2];
    // if so, decrease q' by 1 and repeat this test.
    //
    // The latter test determines at high speed most of the cases in which the
    // trial value q' is one too large, and it eliminates all cases where q' is
    // two too large.
    //

    // The variable uK here denotes u'[j + n - K], where K = 0, 1, 2, and u'
    // denotes the scaled value d * u.

    uint32_t u0 = u.bigits[n];
    uint32_t u1 = u.bigits[n - 1];
    uint32_t u2 = u.bigits[n - 2];

    if (shift > 0)
    {
        DRAGON4_ASSERT((u0 >> (32 - shift)) == 0);

        const uint32_t u3 = (n >= 3) ? u.bigits[n - 3] : 0;
        u0 = (u0 << shift) | (u1 >> (32 - shift));
        u1 = (u1 << shift) | (u2 >> (32 - shift));
        u2 = (u2 << shift) | (u3 >> (32 - shift));
    }
    // u0, u1 and u2 now contain the leading digits of u'.

    // NB: Use repeated subtraction for division to avoid a 64-bit div.
    uint64_t rp = uint64_t{u0} << 32 | u1;
    uint32_t qp = 0;
    while (rp >= v1)
    {
        rp -= v1;
        qp++;
    }
    DRAGON4_ASSERT(qp <= 10);

    if (uint64_t{qp} * v2 > (rp << 32 | u2))
    {
        DRAGON4_ASSERT(qp > 0);
        qp--;
    }
    DRAGON4_ASSERT(qp <= 9);

    //--------------------------------------------------------------------------
    // D4. [Multiply and subtract.]
    //
    // Replace
    //    (u[j] u[j+1] ... u[j+n])_b
    //      := (u[j] u[j+1] ... u[j+n])_b - q' * (v[0] v[1] ... v[n-1] 0)
    //
    // This step consists of a simple multiplication by a one-place number,
    // combined with subtraction. The digits (u[j] ... u[j+n])_b should be kept
    // positive; if the result of this step is actually negative,
    // (u[j] ... u[j+n])_b should be left as the true value plus b^(n+1), i.e.,
    // as the b's complement of the true value, and a "borrow" to the right
    // should be remembered.
    //

    if (qp == 0)
    {
        // No need to multiply.
        return 0;
    }

    uint32_t borrow = 0;
    for (int i = 0; i < n; ++i)
    {
        const uint32_t ui = u.bigits[i];
        const uint32_t vi = v.bigits[i];
        const uint64_t p  = uint64_t{qp} * vi + borrow;
        const uint32_t si = static_cast<uint32_t>(p);
        borrow            = static_cast<uint32_t>(p >> 32);
        const uint32_t di = ui - si;
        borrow           += di > ui;
        u.bigits[i]       = di;
    }
    // vn = 0:
    const uint32_t un = u.bigits[n];
    const uint32_t dn = un - borrow;
    u.bigits[n] = un;

    //--------------------------------------------------------------------------
    // D5. [Test remainder.]
    //
    // Set q[j] := q'. If the result of step D4 was negative, go to step D6;
    // otherwise go on to step D7.
    //

    const bool was_negative = (dn > un);
    if (was_negative)
    {
        //----------------------------------------------------------------------
        // D6. [Add back.]
        //
        // Decrease q[j] by 1, and add (v[0] v[1] ... v[n-1] 0)_b to
        // (u[j] u[j+1] ... u[j+n])_b. (A carry will occur to the right of
        // u[j+n], and it should be ignored since it cancels with the "borrow"
        // that occurred in D4.)
        //
        // The probability that this step is necessary is very small, on the
        // order of only 2/b; test data that activates this step should
        // therefore be specifically contrived when debugging.
        //

        qp--;

        uint32_t carry = 0;
        for (int i = 0; i < n; ++i)
        {
            const uint32_t ui = u.bigits[i];
            const uint32_t vi = v.bigits[i];
            const uint64_t s  = uint64_t{ui} + vi + carry;
            u.bigits[i]       = static_cast<uint32_t>(s);
            carry             = static_cast<uint32_t>(s >> 32);
        }
        // vn = 0:
        u.bigits[n] += carry;
    }

    //--------------------------------------------------------------------------
    // D7. [Loop on j.]
    //
    // Decrease j by one. Now if j >= 0, go back to D3.
    //

    //--------------------------------------------------------------------------
    // D8. [Unnormalize.]
    //
    // Now (q[0] q[1] ... q[m-n])_b is the desired quotient, and the desired
    // remainder may be obtained by dividing (u[0] u[1] ... u[n-1])_b by d.
    //

    // We didn't multiply in the first place, so we don't need to divide here.

    // Still need to clamp the remainder.
    int k = n;
    for ( ; k > 0 && u.bigits[k - 1] == 0; --k)
    {
    }
    u.size = k;

    return qp;
}

inline int Compare(const DiyInt& lhs, const DiyInt& rhs)
{
    const int n1 = lhs.size;
    const int n2 = rhs.size;

    if (n1 < n2) return -1;
    if (n1 > n2) return +1;

    for (int i = n1 - 1; i >= 0; --i)
    {
        const uint32_t b1 = lhs.bigits[i];
        const uint32_t b2 = rhs.bigits[i];

        if (b1 < b2) return -1;
        if (b1 > b2) return +1;
    }

    return 0;
}

// Returns Compare(a + b, c)
inline int CompareAdd(const DiyInt& a, const DiyInt& b, const DiyInt& c)
{
    // NB:
    // This function is only ever called with a >= c, which implies a.size >= c.size.
    DRAGON4_ASSERT(c.size >= a.size);

    const int na = a.size;
    const int nb = b.size;
    const int nc = c.size;

    const int m = Max(na, nb);
    if (m + 1 < nc)
        return -1; // s = (a + b) cannot be larger or equal to c
    if (m > nc)
        return +1; // max(a, b) > c

    // Perform a (partial) left-to-right subtraction, propagating a borrow digit
    // (base B = 2^32) along to the right, stopping as soon as s > c or s < c.

    uint32_t borrow = 0;
    for (int i = nc - 1; i >= 0; --i)
    {
        // Invariant:
        // The leading digits s[i+1],s[i+2],... of s and the leading digits
        // c[i+1],c[i+2],... (after possibly subtracting a borrow) are equal.

        DRAGON4_ASSERT(borrow == 0 || borrow == 1);
        const uint64_t ci = uint64_t{borrow} << 32 | c.bigits[i];
        const uint32_t ai = i < na ? a.bigits[i] : 0;
        const uint32_t bi = i < nb ? b.bigits[i] : 0;
        const uint64_t si = uint64_t{ai} + bi;
        const uint64_t di = ci - si;
//      if (ci < si)
        if (di > ci)
        {
            // Since all the leading digits are equal, this implies c < s,
            // or a + b > c.
            return +1;
        }
        if (di > 1)
        {
            // In this case, the trailing digits s[i-1],s[i-2],... cannot
            // possibly compensate the difference: therefore c > s, or a + b < c.
            return -1;
        }

        // di == 0 or di == 1.
        // If di == 1, borrow B = 2^32 from ci and add to c[i-1], which restores
        // the invariant.
        //  c:      1   2   9   9  ==>  1   1  19   9
        //  s:      1   1  12   3       1   1  12   3
        //              ^                   ^
        //              i                   i
        borrow = static_cast<uint32_t>(di);
    }

//  return borrow == 0 ? 0 : -1;
    return -static_cast<int>(borrow);
}

// Returns the number of leading 0-bits in x, starting at the most significant bit position.
// If x is 0, the result is undefined.
inline int CountLeadingZeros64(uint64_t x)
{
    DRAGON4_ASSERT(x != 0);

#if defined(__GNUC__) || defined(__clang__)
    return __builtin_clzll(x);
#elif defined(_MSC_VER) && (defined(_M_ARM) || defined(_M_ARM64))
    return static_cast<int>(_CountLeadingZeros64(x));
#elif defined(_MSC_VER) && defined(_M_X64)
    return static_cast<int>(__lzcnt64(x));
#elif defined(_MSC_VER) && defined(_M_IX86)
    int lz = static_cast<int>( __lzcnt(static_cast<uint32_t>(x >> 32)) );
    if (lz == 32) {
        lz += static_cast<int>( __lzcnt(static_cast<uint32_t>(x)) );
    }
    return lz;
#else
    int lz = 0;
    while ((x >> 63) == 0) {
        x <<= 1;
        ++lz;
    }
    return lz;
#endif
}

inline int EffectivePrecision(uint64_t f)
{
    DRAGON4_ASSERT(f != 0);
    return 64 - CountLeadingZeros64(f);
}

// Returns: floor(x / 2^n)
inline int SAR(int x, int n)
{
    // Technically, right-shift of negative integers is implementation defined...
    // Should easily get optimized into SAR (or equivalent) instruction.
#if 1
    return x < 0 ? ~(~x >> n) : (x >> n);
#else
    return x >> n;
#endif
}

// Returns: ceil(log_10(2^e))
inline int CeilLog10Pow2(int e)
{
    DRAGON4_ASSERT(e >= -2620);
    DRAGON4_ASSERT(e <=  2620);
    return SAR(e * 315653 + ((1 << 20) - 1), 20);
}

inline int ComputeInitialValuesAndEstimate(DiyInt& r, DiyInt& s, DiyInt& delta, uint64_t f, int e, bool lower_boundary_is_closer)
{
    const int boundary_shift = lower_boundary_is_closer ? 2 : 1;
    const int p = EffectivePrecision(f);
    DRAGON4_ASSERT(p >= 1);
    DRAGON4_ASSERT(p <= 53);
    const int k = CeilLog10Pow2(e + (p - 1));

//  const int cmpf = CompareEstimate((f << boundary_shift) + 1, boundary_shift - e, k);
    if (e >= 0)
    {
        DRAGON4_ASSERT(e >= 0);
        DRAGON4_ASSERT(e <= 971);
        DRAGON4_ASSERT(k >= 0);
        DRAGON4_ASSERT(k <= 308);

        // r = f * 2^(boundary_shift + e)
        AssignU64MulPow2(r, f << boundary_shift, e);
        // s = 2^boundary_shift * 10^k
        AssignPow2MulPow5(s, boundary_shift + k, k);
        // delta = 2^e
        AssignPow2(delta, e);
    }
    else if (k < 0)
    {
        DRAGON4_ASSERT(e >= -1074);
        DRAGON4_ASSERT(e <= -1);
        DRAGON4_ASSERT(k >= -323);
        DRAGON4_ASSERT(k <= -1);

        // r = f * 2^boundary_shift * 10^(-k)
        AssignU64MulPow10(r, f << boundary_shift, -k);
        // s = 2^(boundary_shift - e)
        AssignPow2(s, boundary_shift - e);
        // delta = 10^(-k)
        AssignPow10(delta, -k);
    }
    else
    {
        DRAGON4_ASSERT(e >= -55);
        DRAGON4_ASSERT(e <= -1);
        DRAGON4_ASSERT(k >= 0);
        DRAGON4_ASSERT(k <= 16);

        // r = f * 2^boundary_shift
        AssignU64(r, f << boundary_shift);
        // s = 2^(boundary_shift - e) * 10^k
        AssignPow2MulPow5(s, boundary_shift - e + k, k);
        // delta = 1
        AssignU32(delta, 1);
    }

    return k;
}

} // namespace impl

inline void Dragon4(uint64_t& digits, int& exponent, uint64_t f, int e, bool accept_bounds, bool lower_boundary_is_closer)
{
    using namespace dragon4::impl;

    DiyInt r;
    DiyInt s;
    DiyInt delta;

    //
    // Compute initial values.
    // Estimate k.
    //
    int k = ComputeInitialValuesAndEstimate(r, s, delta, f, e, lower_boundary_is_closer);

    //
    // Fixup, in case k is too low.
    //
    const int cmpf = CompareAdd(r, delta, s);
    if (accept_bounds ? (cmpf >= 0) : (cmpf > 0))
    {
        Mul10(s);
        k++;
    }

    //
    // Generate digits from left to right.
    //
    Mul10(r);       // (Move into init step above?)
    Mul10(delta);   // (Move into init step above?)

    uint64_t d10 = 0;
    for (;;)
    {
        DRAGON4_ASSERT(d10 <= 9999999999999999ull);
        DRAGON4_ASSERT(r.size > 0);

        // q = r / s
        // r = r % s
        uint32_t q = DivMod(r, s);
        DRAGON4_ASSERT(q <= 9);

        const int cmp1 = Compare(r, delta);
        if (lower_boundary_is_closer)
        {
            Mul2(delta);
        }
        const int cmp2 = CompareAdd(r, delta, s);

        const bool tc1 = accept_bounds ? (cmp1 <= 0) : (cmp1 < 0);
        const bool tc2 = accept_bounds ? (cmp2 >= 0) : (cmp2 > 0);
        if (tc1 && tc2)
        {
            // Return the number closer to v.
            // If the two are equidistant from v, use _some_ strategy to break
            // the tie.
            const int cmpr = CompareAdd(r, r, s);
            if (cmpr > 0 || (cmpr == 0 && q % 2 != 0))
            {
                q++;
            }
        }
        else if (!tc1 && tc2)
        {
            q++;
        }

        DRAGON4_ASSERT(q <= 9);
        d10 = d10 * 10 + q;
        k--;

        if (tc1 || tc2)
            break;

        Mul10(r);
        MulAddU32(delta, lower_boundary_is_closer ? 5 : 10);
    }

    digits = d10;
    exponent = k;
}

} // namespace dragon4
