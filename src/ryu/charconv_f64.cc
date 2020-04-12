// Copyright 2019 Ulf Adams
// Copyright 2019 Alexander Bolz
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "charconv_f64.h"

#include <cassert>
#include <climits>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#if _MSC_VER
#include <intrin.h>
#endif

#ifndef RYU_ASSERT
#define RYU_ASSERT(X) assert(X)
#endif

#ifndef RYU_NEVER_INLINE
#if _MSC_VER
#define RYU_NEVER_INLINE __declspec(noinline) inline
#elif __GNUC__
#define RYU_NEVER_INLINE __attribute__((noinline)) inline
#else
#define RYU_NEVER_INLINE inline
#endif
#endif

#ifndef RYU_FORCE_INLINE
#if _MSC_VER
#define RYU_FORCE_INLINE __forceinline
#elif __GNUC__
#define RYU_FORCE_INLINE __attribute__((always_inline)) inline
#else
#define RYU_FORCE_INLINE inline
#endif
#endif

//==================================================================================================
//
//==================================================================================================

template <typename Dest, typename Source>
static inline Dest ReinterpretBits(Source source)
{
    static_assert(sizeof(Dest) == sizeof(Source), "size mismatch");

    Dest dest;
    std::memcpy(&dest, &source, sizeof(Source));
    return dest;
}

namespace {
struct Double
{
    static_assert(std::numeric_limits<double>::is_iec559
               && std::numeric_limits<double>::digits == 53
               && std::numeric_limits<double>::max_exponent == 1024,
        "IEEE-754 double-precision implementation required");

    using value_type = double;
    using bits_type = uint64_t;

//  static constexpr int       MaxDigits10     = std::numeric_limits<value_type>::max_digits10;
    static constexpr int       SignificandSize = std::numeric_limits<value_type>::digits; // = p   (includes the hidden bit)
    static constexpr int       ExponentBias    = std::numeric_limits<value_type>::max_exponent - 1 + (SignificandSize - 1);
//  static constexpr int       MaxExponent     = std::numeric_limits<value_type>::max_exponent - 1 - (SignificandSize - 1);
//  static constexpr int       MinExponent     = std::numeric_limits<value_type>::min_exponent - 1 - (SignificandSize - 1);
    static constexpr bits_type HiddenBit       = bits_type{1} << (SignificandSize - 1);   // = 2^(p-1)
    static constexpr bits_type SignificandMask = HiddenBit - 1;                           // = 2^(p-1) - 1
    static constexpr bits_type ExponentMask    = (bits_type{2 * std::numeric_limits<value_type>::max_exponent - 1}) << (SignificandSize - 1);
    static constexpr bits_type SignMask        = ~(~bits_type{0} >> 1);

    bits_type bits;

    explicit Double(bits_type bits_) : bits(bits_) {}
    explicit Double(value_type value) : bits(ReinterpretBits<bits_type>(value)) {}

    bits_type PhysicalSignificand() const {
        return bits & SignificandMask;
    }

    bits_type PhysicalExponent() const {
        return (bits & ExponentMask) >> (SignificandSize - 1);
    }

    bool IsFinite() const {
        return (bits & ExponentMask) != ExponentMask;
    }

    bool IsInf() const {
        return (bits & ExponentMask) == ExponentMask && (bits & SignificandMask) == 0;
    }

    bool IsNaN() const {
        return (bits & ExponentMask) == ExponentMask && (bits & SignificandMask) != 0;
    }

    bool IsZero() const {
        return (bits & ~SignMask) == 0;
    }

    bool SignBit() const {
        return (bits & SignMask) != 0;
    }

    value_type Value() const {
        return ReinterpretBits<value_type>(bits);
    }

    value_type AbsValue() const {
        return ReinterpretBits<value_type>(bits & ~SignMask);
    }
};
} // namespace

//==================================================================================================
//
//==================================================================================================

static inline int Max(int x, int y)
{
    return y < x ? x : y;
}

// Returns floor(x / 2^n).
static inline int FloorDivPow2(int x, int n)
{
#if 1
    // Technically, right-shift of negative integers is implementation defined...
    // Should easily be optimized into SAR (or equivalent) instruction.
    return x < 0 ? ~(~x >> n) : (x >> n);
#else
    return x >> n;
#endif
}

static inline int FloorLog2Pow5(int e)
{
    RYU_ASSERT(e >= -1764);
    RYU_ASSERT(e <=  1763);
    return FloorDivPow2(e * 1217359, 19);
}

static inline int FloorLog10Pow2(int e)
{
    RYU_ASSERT(e >= -2620);
    RYU_ASSERT(e <=  2620);
    return FloorDivPow2(e * 315653, 20);
}

static inline int FloorLog10Pow5(int e)
{
    RYU_ASSERT(e >= -2620);
    RYU_ASSERT(e <=  2620);
    return FloorDivPow2(e * 732923, 20);
}

static inline uint32_t Lo32(uint64_t x)
{
    return static_cast<uint32_t>(x);
}

//[[maybe_unused]]
static inline uint32_t Hi32(uint64_t x)
{
    return static_cast<uint32_t>(x >> 32);
}

//==================================================================================================
// ToDecimal
//
// Double-precision implementation
//==================================================================================================
// Constant data = 10'656 (+ 400) bytes

static constexpr int BitsPerPow5_Double = 125;

namespace {
struct Uint64x2 {
    uint64_t hi;
    uint64_t lo;
};
}

static inline Uint64x2 ComputePow5_Double(int k)
{
    // Let e = FloorLog2Pow5(k) + 1 - 125
    // For k >= 0, stores 5^k in the form: ceil( 5^k / 2^e )
    // For k <= 0, stores 5^k in the form: ceil(2^-e / 5^-k)
    static constexpr int MinDecExp = -340;
    static constexpr int MaxDecExp =  325;
    static constexpr Uint64x2 Pow5[MaxDecExp - MinDecExp + 1] = {
        {0x1755DC2FF447D7EE, 0xCBAF379E01A5BECA}, // e =  -914, k = -340
        {0x1D2B533BF159CDEA, 0x7E9B0585820F2E7C}, // e =  -912, k = -339
        {0x123B140576D820B2, 0x8F20E37371497D0E}, // e =  -909, k = -338
        {0x16C9D906D48E28DF, 0x32E91C504D9BDC51}, // e =  -907, k = -337
        {0x1C7C4F4889B1B316, 0xFFA363646102D365}, // e =  -905, k = -336
        {0x11CDB18D560F0FEE, 0x5FC61E1EBCA1C41F}, // e =  -902, k = -335
        {0x16411DF0AB92D3E9, 0xF7B7A5A66BCA3527}, // e =  -900, k = -334
        {0x1BD1656CD67788E4, 0x75A58F1006BCC271}, // e =  -898, k = -333
        {0x1162DF64060AB58E, 0xC987796A0435F987}, // e =  -895, k = -332
        {0x15BB973D078D62F2, 0x7BE957C4854377E8}, // e =  -893, k = -331
        {0x1B2A7D0C4970BBAF, 0x1AE3ADB5A69455E2}, // e =  -891, k = -330
        {0x10FA8E27ADE6754D, 0x70CE4C91881CB5AE}, // e =  -888, k = -329
        {0x153931B1996012A0, 0xCD01DFB5EA23E319}, // e =  -886, k = -328
        {0x1A877E1DFFB81749, 0x004257A364ACDBDF}, // e =  -884, k = -327
        {0x1094AED2BFD30E8D, 0xA02976C61EEC096B}, // e =  -881, k = -326
        {0x14B9DA876FC7D231, 0x0833D477A6A70BC6}, // e =  -879, k = -325
        {0x19E851294BB9C6BD, 0x4A40C9959050CEB8}, // e =  -877, k = -324
        {0x103132B9CF541C36, 0x4E687DFD7A328133}, // e =  -874, k = -323
        {0x143D7F6843292343, 0xE2029D7CD8BF2180}, // e =  -872, k = -322
        {0x194CDF4253F36C14, 0xDA8344DC0EEEE9DF}, // e =  -870, k = -321
        {0x1FA01712E8F0471A, 0x1124161312AAA457}, // e =  -868, k = -320
        {0x13C40E6BD1962C70, 0x4AB68DCBEBAAA6B7}, // e =  -865, k = -319
        {0x18B51206C5FBB78C, 0x5D64313EE6955064}, // e =  -863, k = -318
        {0x1EE25688777AA56F, 0x74BD3D8EA03AA47D}, // e =  -861, k = -317
        {0x134D76154AACA765, 0xA8F646792424A6CE}, // e =  -858, k = -316
        {0x1820D39A9D57D13F, 0x1333D8176D2DD082}, // e =  -856, k = -315
        {0x1E29088144ADC58E, 0xD800CE1D487944A2}, // e =  -854, k = -314
        {0x12D9A550CAEC9B79, 0x470080D24D4BCAE6}, // e =  -851, k = -313
        {0x17900EA4FDA7C257, 0x98C0A106E09EBD9F}, // e =  -849, k = -312
        {0x1D74124E3D11B2ED, 0x7EF0C94898C66D06}, // e =  -847, k = -311
        {0x12688B70E62B0FD4, 0x6F567DCD5F7C0424}, // e =  -844, k = -310
        {0x1702AE4D1FB5D3C9, 0x8B2C1D40B75B052D}, // e =  -842, k = -309
        {0x1CC359E067A348BB, 0xEDF72490E531C678}, // e =  -840, k = -308
        {0x11FA182C40C60D75, 0x74BA76DA8F3F1C0B}, // e =  -837, k = -307
        {0x16789E3750F790D2, 0xD1E91491330EE30E}, // e =  -835, k = -306
        {0x1C16C5C525357507, 0x866359B57FD29BD1}, // e =  -833, k = -305
        {0x118E3B9B37416924, 0xB3FE18116FE3A163}, // e =  -830, k = -304
        {0x15F1CA820511C36D, 0xE0FD9E15CBDC89BC}, // e =  -828, k = -303
        {0x1B6E3D2286563449, 0x593D059B3ED3AC2B}, // e =  -826, k = -302
        {0x1124E63593F5E0AD, 0xD7C6238107444B9B}, // e =  -823, k = -301
        {0x156E1FC2F8F358D9, 0x4DB7AC6149155E81}, // e =  -821, k = -300
        {0x1AC9A7B3B7302F0F, 0xA12597799B5AB622}, // e =  -819, k = -299
        {0x10BE08D0527E1D69, 0xC4B77EAC0118B1D5}, // e =  -816, k = -298
        {0x14ED8B04671DA4C4, 0x35E55E57015EDE4A}, // e =  -814, k = -297
        {0x1A28EDC580E50DF5, 0x435EB5ECC1B695DD}, // e =  -812, k = -296
        {0x1059949B708F28B9, 0x4A1B31B3F9121DAA}, // e =  -809, k = -295
        {0x146FF9C24CB2F2E7, 0x9CA1FE20F756A515}, // e =  -807, k = -294
        {0x198BF832DFDFAFA1, 0x83CA7DA9352C4E5A}, // e =  -805, k = -293
        {0x1FEEF63F97D79B89, 0xE4BD1D13827761F0}, // e =  -803, k = -292
        {0x13F559E7BEE6C136, 0x2EF6322C318A9D36}, // e =  -800, k = -291
        {0x18F2B061AEA07183, 0xBAB3BEB73DED4483}, // e =  -798, k = -290
        {0x1F2F5C7A1A488DE4, 0xA960AE650D6895A4}, // e =  -796, k = -289
        {0x137D99CC506D58AE, 0xE9DC6CFF28615D87}, // e =  -793, k = -288
        {0x185D003F6488AEDA, 0xA453883EF279B4E8}, // e =  -791, k = -287
        {0x1E74404F3DAADA91, 0x4D686A4EAF182222}, // e =  -789, k = -286
        {0x1308A831868AC89A, 0xD06142712D6F1556}, // e =  -786, k = -285
        {0x17CAD23DE82D7AC1, 0x8479930D78CADAAB}, // e =  -784, k = -284
        {0x1DBD86CD6238D971, 0xE597F7D0D6FD9156}, // e =  -782, k = -283
        {0x129674405D6387E7, 0x2F7EFAE2865E7AD6}, // e =  -779, k = -282
        {0x173C115074BC69E0, 0xFB5EB99B27F6198B}, // e =  -777, k = -281
        {0x1D0B15A491EB8459, 0x3A366801F1F39FEE}, // e =  -775, k = -280
        {0x1226ED86DB3332B7, 0xC4620101373843F5}, // e =  -772, k = -279
        {0x16B0A8E891FFFF65, 0xB57A8141850654F2}, // e =  -770, k = -278
        {0x1C5CD322B67FFF3F, 0x22D92191E647EA2E}, // e =  -768, k = -277
        {0x11BA03F5B20FFF87, 0x75C7B4FB2FECF25D}, // e =  -765, k = -276
        {0x162884F31E93FF69, 0x5339A239FBE82EF4}, // e =  -763, k = -275
        {0x1BB2A62FE638FF43, 0xA8080AC87AE23AB1}, // e =  -761, k = -274
        {0x114FA7DDEFE39F8A, 0x490506BD4CCD64AF}, // e =  -758, k = -273
        {0x15A391D56BDC876C, 0xDB46486CA000BDDA}, // e =  -756, k = -272
        {0x1B0C764AC6D3A948, 0x1217DA87C800ED51}, // e =  -754, k = -271
        {0x10E7C9EEBC4449CD, 0x0B4EE894DD009453}, // e =  -751, k = -270
        {0x1521BC6A6B555C40, 0x4E22A2BA1440B967}, // e =  -749, k = -269
        {0x1A6A2B85062AB350, 0x61AB4B689950E7C1}, // e =  -747, k = -268
        {0x10825B3323DAB012, 0x3D0B0F215FD290D9}, // e =  -744, k = -267
        {0x14A2F1FFECD15C16, 0xCC4DD2E9B7C7350F}, // e =  -742, k = -266
        {0x19CBAE7FE805B31C, 0x7F6147A425B90252}, // e =  -740, k = -265
        {0x101F4D0FF1038FF1, 0xCF9CCCC69793A174}, // e =  -737, k = -264
        {0x14272053ED4473EE, 0x4383FFF83D7889D1}, // e =  -735, k = -263
        {0x1930E868E89590E9, 0xD464FFF64CD6AC45}, // e =  -733, k = -262
        {0x1F7D228322BAF524, 0x497E3FF3E00C5756}, // e =  -731, k = -261
        {0x13AE3591F5B4D936, 0xADEEE7F86C07B696}, // e =  -728, k = -260
        {0x1899C2F673220F84, 0x596AA1F68709A43B}, // e =  -726, k = -259
        {0x1EC033B40FEA9365, 0x6FC54A7428CC0D4A}, // e =  -724, k = -258
        {0x1338205089F29C1F, 0x65DB4E88997F884E}, // e =  -721, k = -257
        {0x18062864AC6F4327, 0x3F52222ABFDF6A62}, // e =  -719, k = -256
        {0x1E07B27DD78B13F1, 0x0F26AAB56FD744FA}, // e =  -717, k = -255
        {0x12C4CF8EA6B6EC76, 0xA9782AB165E68B1C}, // e =  -714, k = -254
        {0x177603725064A794, 0x53D6355DBF602DE3}, // e =  -712, k = -253
        {0x1D53844EE47DD179, 0x68CBC2B52F38395C}, // e =  -710, k = -252
        {0x125432B14ECEA2EB, 0xE17F59B13D8323DA}, // e =  -707, k = -251
        {0x16E93F5DA2824BA6, 0xD9DF301D8CE3ECD0}, // e =  -705, k = -250
        {0x1CA38F350B22DE90, 0x9056FC24F01CE804}, // e =  -703, k = -249
        {0x11E6398126F5CB1A, 0x5A365D9716121103}, // e =  -700, k = -248
        {0x165FC7E170B33DE0, 0xF0C3F4FCDB969543}, // e =  -698, k = -247
        {0x1BF7B9D9CCE00D59, 0x2CF4F23C127C3A94}, // e =  -696, k = -246
        {0x117AD428200C0857, 0xBC1917658B8DA49D}, // e =  -693, k = -245
        {0x15D98932280F0A6D, 0xAB1F5D3EEE710DC4}, // e =  -691, k = -244
        {0x1B4FEB7EB212CD09, 0x15E7348EAA0D5134}, // e =  -689, k = -243
        {0x1111F32F2F4BC025, 0xADB080D92A4852C1}, // e =  -686, k = -242
        {0x15566FFAFB1EB02F, 0x191CA10F74DA6771}, // e =  -684, k = -241
        {0x1AAC0BF9B9E65C3A, 0xDF63C9535211014D}, // e =  -682, k = -240
        {0x10AB877C142FF9A4, 0xCB9E5DD4134AA0D0}, // e =  -679, k = -239
        {0x14D6695B193BF80D, 0xFE85F549181D4904}, // e =  -677, k = -238
        {0x1A0C03B1DF8AF611, 0x7E27729B5E249B45}, // e =  -675, k = -237
        {0x1047824F2BB6D9CA, 0xEED8A7A11AD6E10C}, // e =  -672, k = -236
        {0x145962E2F6A4903D, 0xAA8ED189618C994E}, // e =  -670, k = -235
        {0x196FBB9BB44DB44D, 0x153285EBB9EFBFA2}, // e =  -668, k = -234
        {0x1FCBAA82A1612160, 0x5A7F2766A86BAF8A}, // e =  -666, k = -233
        {0x13DF4A91A4DCB4DC, 0x388F78A029434DB6}, // e =  -663, k = -232
        {0x18D71D360E13E213, 0x46B356C833942124}, // e =  -661, k = -231
        {0x1F0CE4839198DA98, 0x18602C7A4079296D}, // e =  -659, k = -230
        {0x13680ED23AFF889F, 0x0F3C1BCC684BB9E4}, // e =  -656, k = -229
        {0x18421286C9BF6AC6, 0xD30B22BF825EA85D}, // e =  -654, k = -228
        {0x1E5297287C2F4578, 0x87CDEB6F62F65274}, // e =  -652, k = -227
        {0x12F39E794D9D8B6B, 0x54E0B3259DD9F389}, // e =  -649, k = -226
        {0x17B08617A104EE46, 0x2A18DFEF0550706B}, // e =  -647, k = -225
        {0x1D9CA79D894629D7, 0xB49F17EAC6A48C86}, // e =  -645, k = -224
        {0x1281E8C275CBDA26, 0xD0E36EF2BC26D7D4}, // e =  -642, k = -223
        {0x172262F3133ED0B0, 0x851C4AAF6B308DC8}, // e =  -640, k = -222
        {0x1CEAFBAFD80E84DC, 0xA6635D5B45FCB13A}, // e =  -638, k = -221
        {0x1212DD4DE7091309, 0xE7FE1A590BBDEEC5}, // e =  -635, k = -220
        {0x169794A160CB57CC, 0x61FDA0EF4EAD6A76}, // e =  -633, k = -219
        {0x1C3D79C9B8FE2DBF, 0x7A7D092B2258C513}, // e =  -631, k = -218
        {0x11A66C1E139EDC97, 0xAC8E25BAF5777B2C}, // e =  -628, k = -217
        {0x16100725988693BD, 0x97B1AF29B2D559F7}, // e =  -626, k = -216
        {0x1B9408EEFEA838AC, 0xFD9E1AF41F8AB075}, // e =  -624, k = -215
        {0x113C85955F29236C, 0x1E82D0D893B6AE49}, // e =  -621, k = -214
        {0x158BA6FAB6F36C47, 0x2623850EB8A459DB}, // e =  -619, k = -213
        {0x1AEE90B964B04758, 0xEFAC665266CD7052}, // e =  -617, k = -212
        {0x10D51A73DEEE2C97, 0x95CBBFF380406633}, // e =  -614, k = -211
        {0x150A6110D6A9B7BD, 0x7B3EAFF060507FC0}, // e =  -612, k = -210
        {0x1A4CF9550C5425AC, 0xDA0E5BEC78649FB0}, // e =  -610, k = -209
        {0x10701BD527B4978C, 0x0848F973CB3EE3CE}, // e =  -607, k = -208
        {0x148C22CA71A1BD6F, 0x0A5B37D0BE0E9CC2}, // e =  -605, k = -207
        {0x19AF2B7D0E0A2CCA, 0xCCF205C4ED9243F2}, // e =  -603, k = -206
        {0x100D7B2E28C65BFE, 0xC017439B147B6A77}, // e =  -600, k = -205
        {0x1410D9F9B2F7F2FE, 0x701D1481D99A4515}, // e =  -598, k = -204
        {0x191510781FB5EFBE, 0x0C2459A25000D65A}, // e =  -596, k = -203
        {0x1F5A549627A36BAD, 0x8F2D700AE4010BF1}, // e =  -594, k = -202
        {0x139874DDD8C6234C, 0x797C6606CE80A777}, // e =  -591, k = -201
        {0x187E92154EF7AC1F, 0x97DB7F888220D154}, // e =  -589, k = -200
        {0x1E9E369AA2B59727, 0x7DD25F6AA2A905A9}, // e =  -587, k = -199
        {0x1322E220A5B17E78, 0xAEA37BA2A5A9A38A}, // e =  -584, k = -198
        {0x17EB9AA8CF1DDE16, 0xDA4C5A8B4F140C6C}, // e =  -582, k = -197
        {0x1DE6815302E5559C, 0x90DF712E22D90F87}, // e =  -580, k = -196
        {0x12B010D3E1CF5581, 0xDA8BA6BCD5C7A9B5}, // e =  -577, k = -195
        {0x175C1508DA432AE2, 0x512E906C0B399422}, // e =  -575, k = -194
        {0x1D331A4B10D3F59A, 0xE57A34870E07F92A}, // e =  -573, k = -193
        {0x123FF06EEA847980, 0xCF6C60D468C4FBBA}, // e =  -570, k = -192
        {0x16CFEC8AA52597E1, 0x0347790982F63AA9}, // e =  -568, k = -191
        {0x1C83E7AD4E6EFDD9, 0x4419574BE3B3C953}, // e =  -566, k = -190
        {0x11D270CC51055EA7, 0xCA8FD68F6E505DD4}, // e =  -563, k = -189
        {0x16470CFF6546B651, 0xBD33CC3349E47549}, // e =  -561, k = -188
        {0x1BD8D03F3E9863E6, 0x2C80BF401C5D929B}, // e =  -559, k = -187
        {0x11678227871F3E6F, 0xDBD0778811BA7BA1}, // e =  -556, k = -186
        {0x15C162B168E70E0B, 0xD2C4956A16291A89}, // e =  -554, k = -185
        {0x1B31BB5DC320D18E, 0xC775BAC49BB3612B}, // e =  -552, k = -184
        {0x10FF151A99F482F9, 0x3CA994BAE1501CBB}, // e =  -549, k = -183
        {0x153EDA614071A3B7, 0x8BD3F9E999A423EA}, // e =  -547, k = -182
        {0x1A8E90F9908E0CA5, 0x6EC8F864000D2CE4}, // e =  -545, k = -181
        {0x10991A9BFA58C7E7, 0x653D9B3E80083C0F}, // e =  -542, k = -180
        {0x14BF6142F8EEF9E1, 0x3E8D020E200A4B13}, // e =  -540, k = -179
        {0x19EF3993B72AB859, 0x8E304291A80CDDD7}, // e =  -538, k = -178
        {0x103583FC527AB337, 0xF8DE299B09080AA7}, // e =  -535, k = -177
        {0x1442E4FB67196005, 0xF715B401CB4A0D50}, // e =  -533, k = -176
        {0x19539E3A40DFB807, 0x74DB21023E1C90A4}, // e =  -531, k = -175
        {0x1FA885C8D117A609, 0x5211E942CDA3B4CD}, // e =  -529, k = -174
        {0x13C9539D82AEC7C5, 0xD34B31C9C0865100}, // e =  -526, k = -173
        {0x18BBA884E35A79B7, 0x481DFE3C30A7E540}, // e =  -524, k = -172
        {0x1EEA92A61C311825, 0x1A257DCB3CD1DE90}, // e =  -522, k = -171
        {0x13529BA7D19EAF17, 0x30576E9F06032B1A}, // e =  -519, k = -170
        {0x18274291C6065ADC, 0xFC6D4A46C783F5E1}, // e =  -517, k = -169
        {0x1E3113363787F194, 0x3B889CD87964F359}, // e =  -515, k = -168
        {0x12DEAC01E2B4F6FC, 0xA53562074BDF1818}, // e =  -512, k = -167
        {0x179657025B6234BB, 0xCE82BA891ED6DE1D}, // e =  -510, k = -166
        {0x1D7BECC2F23AC1EA, 0xC223692B668C95A5}, // e =  -508, k = -165
        {0x126D73F9D764B932, 0xB95621BB2017DD87}, // e =  -505, k = -164
        {0x1708D0F84D3DE77F, 0x67ABAA29E81DD4E9}, // e =  -503, k = -163
        {0x1CCB0536608D615F, 0x419694B462254A23}, // e =  -501, k = -162
        {0x11FEE341FC585CDB, 0x88FE1CF0BD574E56}, // e =  -498, k = -161
        {0x167E9C127B6E7412, 0x6B3DA42CECAD21EB}, // e =  -496, k = -160
        {0x1C1E43171A4A1117, 0x060D0D3827D86A66}, // e =  -494, k = -159
        {0x1192E9EE706E4AAE, 0x63C8284318E74280}, // e =  -491, k = -158
        {0x15F7A46A0C89DD59, 0xFCBA3253DF211320}, // e =  -489, k = -157
        {0x1B758D848FAC54B0, 0x7BE8BEE8D6E957E8}, // e =  -487, k = -156
        {0x11297872D9CBB4EE, 0x4D7177518651D6F1}, // e =  -484, k = -155
        {0x1573D68F903EA229, 0xE0CDD525E7E64CAD}, // e =  -482, k = -154
        {0x1AD0CC33744E4AB4, 0x59014A6F61DFDFD8}, // e =  -480, k = -153
        {0x10C27FA028B0EEB0, 0xB7A0CE859D2BEBE7}, // e =  -477, k = -152
        {0x14F31F8832DD2A5C, 0xE58902270476E6E1}, // e =  -475, k = -151
        {0x1A2FE76A3F9474F4, 0x1EEB42B0C594A099}, // e =  -473, k = -150
        {0x105DF0A267BCC918, 0x935309AE7B7CE460}, // e =  -470, k = -149
        {0x14756CCB01ABFB5E, 0xB827CC1A1A5C1D78}, // e =  -468, k = -148
        {0x1992C7FDC216FA36, 0x6631BF20A0F324D6}, // e =  -466, k = -147
        {0x1FF779FD329CB8C3, 0xFFBE2EE8C92FEE0B}, // e =  -464, k = -146
        {0x13FAAC3E3FA1F37A, 0x7FD6DD517DBDF4C7}, // e =  -461, k = -145
        {0x18F9574DCF8A7059, 0x1FCC94A5DD2D71F9}, // e =  -459, k = -144
        {0x1F37AD21436D0C6F, 0x67BFB9CF5478CE77}, // e =  -457, k = -143
        {0x1382CC34CA2427C5, 0xA0D7D42194CB810A}, // e =  -454, k = -142
        {0x18637F41FCAD31B7, 0x090DC929F9FE614D}, // e =  -452, k = -141
        {0x1E7C5F127BD87E24, 0xCB513B74787DF9A0}, // e =  -450, k = -140
        {0x130DBB6B8D674ED6, 0xFF12C528CB4EBC04}, // e =  -447, k = -139
        {0x17D12A4670C1228C, 0xBED77672FE226B05}, // e =  -445, k = -138
        {0x1DC574D80CF16B2F, 0xEE8D540FBDAB05C6}, // e =  -443, k = -137
        {0x129B69070816E2FD, 0xF5185489D68AE39C}, // e =  -440, k = -136
        {0x17424348CA1C9BBD, 0x725E69AC4C2D9C83}, // e =  -438, k = -135
        {0x1D12D41AFCA3C2AC, 0xCEF604175F3903A3}, // e =  -436, k = -134
        {0x122BC490DDE659AC, 0x0159C28E9B83A246}, // e =  -433, k = -133
        {0x16B6B5B5155FF017, 0x01B0333242648AD8}, // e =  -431, k = -132
        {0x1C6463225AB7EC1C, 0xC21C3FFED2FDAD8E}, // e =  -429, k = -131
        {0x11BEBDF578B2F391, 0xF951A7FF43DE8C79}, // e =  -426, k = -130
        {0x162E6D72D6DFB076, 0x77A611FF14D62F97}, // e =  -424, k = -129
        {0x1BBA08CF8C979C94, 0x158F967EDA0BBB7C}, // e =  -422, k = -128
        {0x11544581B7DEC1DC, 0x8D79BE0F4847552E}, // e =  -419, k = -127
        {0x15A956E225D67253, 0xB0D82D931A592A79}, // e =  -417, k = -126
        {0x1B13AC9AAF4C0EE8, 0x9D0E38F7E0EF7517}, // e =  -415, k = -125
        {0x10EC4BE0AD8F8951, 0x6228E39AEC95A92F}, // e =  -412, k = -124
        {0x15275ED8D8F36BA5, 0xBAB31C81A7BB137A}, // e =  -410, k = -123
        {0x1A71368F0F30468F, 0x295FE3A211A9D859}, // e =  -408, k = -122
        {0x1086C219697E2C19, 0x79DBEE454B0A2738}, // e =  -405, k = -121
        {0x14A8729FC3DDB71F, 0xD852E9D69DCCB106}, // e =  -403, k = -120
        {0x19D28F47B4D524E7, 0xCE67A44C453FDD47}, // e =  -401, k = -119
        {0x1023998CD1053710, 0xE100C6AFAB47EA4C}, // e =  -398, k = -118
        {0x142C7FF0054684D5, 0x1940F85B9619E4DF}, // e =  -396, k = -117
        {0x19379FEC0698260A, 0x5F9136727BA05E17}, // e =  -394, k = -116
        {0x1F8587E7083E2F8C, 0xF775840F1A88759D}, // e =  -392, k = -115
        {0x13B374F06526DDB8, 0x1AA9728970954982}, // e =  -389, k = -114
        {0x18A0522C7E709526, 0x2153CF2BCCBA9BE3}, // e =  -387, k = -113
        {0x1EC866B79E0CBA6F, 0xA9A8C2F6BFE942DB}, // e =  -385, k = -112
        {0x133D4032C2C7F485, 0xCA0979DA37F1C9C9}, // e =  -382, k = -111
        {0x180C903F7379F1A7, 0x3C8BD850C5EE3C3B}, // e =  -380, k = -110
        {0x1E0FB44F50586E11, 0x0BAECE64F769CB4A}, // e =  -378, k = -109
        {0x12C9D0B1923744CA, 0xA74D40FF1AA21F0E}, // e =  -375, k = -108
        {0x177C44DDF6C515FD, 0x5120913EE14AA6D2}, // e =  -373, k = -107
        {0x1D5B561574765B7C, 0xA568B58E999D5086}, // e =  -371, k = -106
        {0x125915CD68C9F92D, 0xE761717920025254}, // e =  -368, k = -105
        {0x16EF5B40C2FC7779, 0x6139CDD76802E6E9}, // e =  -366, k = -104
        {0x1CAB3210F3BB9557, 0xB988414D4203A0A3}, // e =  -364, k = -103
        {0x11EAFF4A98553D56, 0xD3F528D049424466}, // e =  -361, k = -102
        {0x1665BF1D3E6A8CAC, 0x88F273045B92D580}, // e =  -359, k = -101
        {0x1BFF2EE48E052FD7, 0xAB2F0FC572778ADF}, // e =  -357, k = -100
        {0x117F7D4ED8C33DE6, 0xCAFD69DB678AB6CC}, // e =  -354, k =  -99
        {0x15DF5CA28EF40D60, 0x7DBCC452416D647F}, // e =  -352, k =  -98
        {0x1B5733CB32B110B8, 0x9D2BF566D1C8BD9E}, // e =  -350, k =  -97
        {0x1116805EFFAEAA73, 0x623B7960431D7683}, // e =  -347, k =  -96
        {0x155C2076BF9A5510, 0x3ACA57B853E4D424}, // e =  -345, k =  -95
        {0x1AB328946F80EA54, 0x497CEDA668DE092C}, // e =  -343, k =  -94
        {0x10AFF95CC5B09274, 0xADEE1488018AC5BC}, // e =  -340, k =  -93
        {0x14DBF7B3F71CB711, 0xD96999AA01ED772B}, // e =  -338, k =  -92
        {0x1A12F5A0F4E3E4D6, 0x4FC400148268D4F5}, // e =  -336, k =  -91
        {0x104BD984990E6F05, 0xF1DA800CD181851A}, // e =  -333, k =  -90
        {0x145ECFE5BF520AC7, 0x6E51201005E1E660}, // e =  -331, k =  -89
        {0x197683DF2F268D79, 0x49E56814075A5FF8}, // e =  -329, k =  -88
        {0x1FD424D6FAF030D7, 0x9C5EC2190930F7F6}, // e =  -327, k =  -87
        {0x13E497065CD61E86, 0xC1BB394FA5BE9AFA}, // e =  -324, k =  -86
        {0x18DDBCC7F40BA628, 0x722A07A38F2E41B8}, // e =  -322, k =  -85
        {0x1F152BF9F10E8FB2, 0x8EB4898C72F9D226}, // e =  -320, k =  -84
        {0x136D3B7C36A919CF, 0x9930D5F7C7DC2358}, // e =  -317, k =  -83
        {0x18488A5B44536043, 0x7F7D0B75B9D32C2E}, // e =  -315, k =  -82
        {0x1E5AACF215683854, 0x5F5C4E532847F739}, // e =  -313, k =  -81
        {0x12F8AC174D612334, 0xBB99B0F3F92CFA84}, // e =  -310, k =  -80
        {0x17B6D71D20B96C01, 0xEA801D30F7783925}, // e =  -308, k =  -79
        {0x1DA48CE468E7C702, 0x6520247D3556476E}, // e =  -306, k =  -78
        {0x1286D80EC190DC61, 0x7F3416CE4155ECA5}, // e =  -303, k =  -77
        {0x17288E1271F51379, 0xDF011C81D1AB67CE}, // e =  -301, k =  -76
        {0x1CF2B1970E725858, 0x56C163A2461641C1}, // e =  -299, k =  -75
        {0x1217AEFE69077737, 0x3638DE456BCDE919}, // e =  -296, k =  -74
        {0x169D9ABE03495505, 0x03C715D6C6C1635F}, // e =  -294, k =  -73
        {0x1C45016D841BAA46, 0x44B8DB4C7871BC37}, // e =  -292, k =  -72
        {0x11AB20E472914A6B, 0xEAF3890FCB4715A2}, // e =  -289, k =  -71
        {0x1615E91D8F359D06, 0xE5B06B53BE18DB0B}, // e =  -287, k =  -70
        {0x1B9B6364F3030448, 0x9F1C8628AD9F11CD}, // e =  -285, k =  -69
        {0x11411E1F17E1E2AD, 0x6371D3D96C836B20}, // e =  -282, k =  -68
        {0x159165A6DDDA5B58, 0xBC4E48CFC7A445E8}, // e =  -280, k =  -67
        {0x1AF5BF109550F22E, 0xEB61DB03B98D5762}, // e =  -278, k =  -66
        {0x10D9976A5D52975D, 0x531D28E253F8569E}, // e =  -275, k =  -65
        {0x150FFD44F4A73D34, 0xA7E4731AE8F66C45}, // e =  -273, k =  -64
        {0x1A53FC9631D10C81, 0xD1DD8FE1A3340756}, // e =  -271, k =  -63
        {0x10747DDDDF22A7D1, 0x232A79ED06008496}, // e =  -268, k =  -62
        {0x14919D5556EB51C5, 0x6BF518684780A5BB}, // e =  -266, k =  -61
        {0x19B604AAACA62636, 0xC6F25E825960CF2A}, // e =  -264, k =  -60
        {0x1011C2EAABE7D7E2, 0x3C577B1177DC817B}, // e =  -261, k =  -59
        {0x141633A556E1CDDA, 0xCB6D59D5D5D3A1D9}, // e =  -259, k =  -58
        {0x191BC08EAC9A4151, 0x7E48B04B4B488A4F}, // e =  -257, k =  -57
        {0x1F62B0B257C0D1A5, 0xDDDADC5E1E1AACE3}, // e =  -255, k =  -56
        {0x139DAE6F76D88307, 0xAAA8C9BAD2D0AC0E}, // e =  -252, k =  -55
        {0x18851A0B548EA3C9, 0x9552FC298784D711}, // e =  -250, k =  -54
        {0x1EA6608E29B24CBB, 0xFAA7BB33E9660CD6}, // e =  -248, k =  -53
        {0x1327FC58DA0F6FF5, 0x7CA8D50071DFC806}, // e =  -245, k =  -52
        {0x17F1FB6F10934BF2, 0xDBD30A408E57BA07}, // e =  -243, k =  -51
        {0x1DEE7A4AD4B81EEF, 0x92C7CCD0B1EDA889}, // e =  -241, k =  -50
        {0x12B50C6EC4F31355, 0xBBBCE0026F348956}, // e =  -238, k =  -49
        {0x17624F8A762FD82B, 0x2AAC18030B01ABAB}, // e =  -236, k =  -48
        {0x1D3AE36D13BBCE35, 0xF5571E03CDC21695}, // e =  -234, k =  -47
        {0x1244CE242C5560E1, 0xB95672C260994E1E}, // e =  -231, k =  -46
        {0x16D601AD376AB91A, 0x27AC0F72F8BFA1A5}, // e =  -229, k =  -45
        {0x1C8B821885456760, 0xB197134FB6EF8A0E}, // e =  -227, k =  -44
        {0x11D7314F534B609C, 0x6EFE6C11D255B649}, // e =  -224, k =  -43
        {0x164CFDA3281E38C3, 0x8ABE071646EB23DB}, // e =  -222, k =  -42
        {0x1BE03D0BF225C6F4, 0x6D6D88DBD8A5ECD2}, // e =  -220, k =  -41
        {0x116C262777579C58, 0xC46475896767B403}, // e =  -217, k =  -40
        {0x15C72FB1552D836E, 0xF57D92EBC141A104}, // e =  -215, k =  -39
        {0x1B38FB9DAA78E44A, 0xB2DCF7A6B1920945}, // e =  -213, k =  -38
        {0x11039D428A8B8EAE, 0xAFCA1AC82EFB45CB}, // e =  -210, k =  -37
        {0x154484932D2E725A, 0x5BBCA17A3ABA173E}, // e =  -208, k =  -36
        {0x1A95A5B7F87A0EF0, 0xF2ABC9D8C9689D0D}, // e =  -206, k =  -35
        {0x109D8792FB4C4956, 0x97AB5E277DE16228}, // e =  -203, k =  -34
        {0x14C4E977BA1F5BAC, 0x3D9635B15D59BAB2}, // e =  -201, k =  -33
        {0x19F623D5A8A73297, 0x4CFBC31DB4B0295F}, // e =  -199, k =  -32
        {0x1039D66589687F9E, 0x901D59F290EE19DB}, // e =  -196, k =  -31
        {0x14484BFEEBC29F86, 0x3424B06F3529A052}, // e =  -194, k =  -30
        {0x195A5EFEA6B34767, 0xC12DDC8B02740867}, // e =  -192, k =  -29
        {0x1FB0F6BE50601941, 0xB17953ADC3110A80}, // e =  -190, k =  -28
        {0x13CE9A36F23C0FC9, 0x0EEBD44C99EAA690}, // e =  -187, k =  -27
        {0x18C240C4AECB13BB, 0x52A6C95FC0655034}, // e =  -185, k =  -26
        {0x1EF2D0F5DA7DD8AA, 0x27507BB7B07EA441}, // e =  -183, k =  -25
        {0x1357C299A88EA76A, 0x58924D52CE4F26A9}, // e =  -180, k =  -24
        {0x182DB34012B25144, 0xEEB6E0A781E2F053}, // e =  -178, k =  -23
        {0x1E392010175EE596, 0x2A6498D1625BAC68}, // e =  -176, k =  -22
        {0x12E3B40A0E9B4F7D, 0xDA7EDF82DD794BC1}, // e =  -173, k =  -21
        {0x179CA10C9242235D, 0x511E976394D79EB1}, // e =  -171, k =  -20
        {0x1D83C94FB6D2AC34, 0xA5663D3C7A0D865D}, // e =  -169, k =  -19
        {0x12725DD1D243ABA0, 0xE75FE645CC4873FA}, // e =  -166, k =  -18
        {0x170EF54646D49689, 0x2137DFD73F5A90F9}, // e =  -164, k =  -17
        {0x1CD2B297D889BC2B, 0x6985D7CD0F313537}, // e =  -162, k =  -16
        {0x1203AF9EE756159B, 0x21F3A6E0297EC143}, // e =  -159, k =  -15
        {0x16849B86A12B9B01, 0xEA70909833DE7193}, // e =  -157, k =  -14
        {0x1C25C268497681C2, 0x650CB4BE40D60DF8}, // e =  -155, k =  -13
        {0x119799812DEA1119, 0x7F27F0F6E885C8BB}, // e =  -152, k =  -12
        {0x15FD7FE17964955F, 0xDEF1ED34A2A73AEA}, // e =  -150, k =  -11
        {0x1B7CDFD9D7BDBAB7, 0xD6AE6881CB5109A4}, // e =  -148, k =  -10
        {0x112E0BE826D694B2, 0xE62D01511F12A607}, // e =  -145, k =   -9
        {0x15798EE2308C39DF, 0x9FB841A566D74F88}, // e =  -143, k =   -8
        {0x1AD7F29ABCAF4857, 0x87A6520EC08D236A}, // e =  -141, k =   -7
        {0x10C6F7A0B5ED8D36, 0xB4C7F34938583622}, // e =  -138, k =   -6
        {0x14F8B588E368F084, 0x61F9F01B866E43AB}, // e =  -136, k =   -5
        {0x1A36E2EB1C432CA5, 0x7A786C226809D496}, // e =  -134, k =   -4
        {0x10624DD2F1A9FBE7, 0x6C8B4395810624DE}, // e =  -131, k =   -3
        {0x147AE147AE147AE1, 0x47AE147AE147AE15}, // e =  -129, k =   -2
        {0x1999999999999999, 0x999999999999999A}, // e =  -127, k =   -1
        {0x1000000000000000, 0x0000000000000000}, // e =  -124, k =    0
        {0x1400000000000000, 0x0000000000000000}, // e =  -122, k =    1
        {0x1900000000000000, 0x0000000000000000}, // e =  -120, k =    2
        {0x1F40000000000000, 0x0000000000000000}, // e =  -118, k =    3
        {0x1388000000000000, 0x0000000000000000}, // e =  -115, k =    4
        {0x186A000000000000, 0x0000000000000000}, // e =  -113, k =    5
        {0x1E84800000000000, 0x0000000000000000}, // e =  -111, k =    6
        {0x1312D00000000000, 0x0000000000000000}, // e =  -108, k =    7
        {0x17D7840000000000, 0x0000000000000000}, // e =  -106, k =    8
        {0x1DCD650000000000, 0x0000000000000000}, // e =  -104, k =    9
        {0x12A05F2000000000, 0x0000000000000000}, // e =  -101, k =   10
        {0x174876E800000000, 0x0000000000000000}, // e =   -99, k =   11
        {0x1D1A94A200000000, 0x0000000000000000}, // e =   -97, k =   12
        {0x12309CE540000000, 0x0000000000000000}, // e =   -94, k =   13
        {0x16BCC41E90000000, 0x0000000000000000}, // e =   -92, k =   14
        {0x1C6BF52634000000, 0x0000000000000000}, // e =   -90, k =   15
        {0x11C37937E0800000, 0x0000000000000000}, // e =   -87, k =   16
        {0x16345785D8A00000, 0x0000000000000000}, // e =   -85, k =   17
        {0x1BC16D674EC80000, 0x0000000000000000}, // e =   -83, k =   18
        {0x1158E460913D0000, 0x0000000000000000}, // e =   -80, k =   19
        {0x15AF1D78B58C4000, 0x0000000000000000}, // e =   -78, k =   20
        {0x1B1AE4D6E2EF5000, 0x0000000000000000}, // e =   -76, k =   21
        {0x10F0CF064DD59200, 0x0000000000000000}, // e =   -73, k =   22
        {0x152D02C7E14AF680, 0x0000000000000000}, // e =   -71, k =   23
        {0x1A784379D99DB420, 0x0000000000000000}, // e =   -69, k =   24
        {0x108B2A2C28029094, 0x0000000000000000}, // e =   -66, k =   25
        {0x14ADF4B7320334B9, 0x0000000000000000}, // e =   -64, k =   26
        {0x19D971E4FE8401E7, 0x4000000000000000}, // e =   -62, k =   27
        {0x1027E72F1F128130, 0x8800000000000000}, // e =   -59, k =   28
        {0x1431E0FAE6D7217C, 0xAA00000000000000}, // e =   -57, k =   29
        {0x193E5939A08CE9DB, 0xD480000000000000}, // e =   -55, k =   30
        {0x1F8DEF8808B02452, 0xC9A0000000000000}, // e =   -53, k =   31
        {0x13B8B5B5056E16B3, 0xBE04000000000000}, // e =   -50, k =   32
        {0x18A6E32246C99C60, 0xAD85000000000000}, // e =   -48, k =   33
        {0x1ED09BEAD87C0378, 0xD8E6400000000000}, // e =   -46, k =   34
        {0x13426172C74D822B, 0x878FE80000000000}, // e =   -43, k =   35
        {0x1812F9CF7920E2B6, 0x6973E20000000000}, // e =   -41, k =   36
        {0x1E17B84357691B64, 0x03D0DA8000000000}, // e =   -39, k =   37
        {0x12CED32A16A1B11E, 0x8262889000000000}, // e =   -36, k =   38
        {0x178287F49C4A1D66, 0x22FB2AB400000000}, // e =   -34, k =   39
        {0x1D6329F1C35CA4BF, 0xABB9F56100000000}, // e =   -32, k =   40
        {0x125DFA371A19E6F7, 0xCB54395CA0000000}, // e =   -29, k =   41
        {0x16F578C4E0A060B5, 0xBE2947B3C8000000}, // e =   -27, k =   42
        {0x1CB2D6F618C878E3, 0x2DB399A0BA000000}, // e =   -25, k =   43
        {0x11EFC659CF7D4B8D, 0xFC90400474400000}, // e =   -22, k =   44
        {0x166BB7F0435C9E71, 0x7BB4500591500000}, // e =   -20, k =   45
        {0x1C06A5EC5433C60D, 0xDAA16406F5A40000}, // e =   -18, k =   46
        {0x118427B3B4A05BC8, 0xA8A4DE8459868000}, // e =   -15, k =   47
        {0x15E531A0A1C872BA, 0xD2CE16256FE82000}, // e =   -13, k =   48
        {0x1B5E7E08CA3A8F69, 0x87819BAECBE22800}, // e =   -11, k =   49
        {0x111B0EC57E6499A1, 0xF4B1014D3F6D5900}, // e =    -8, k =   50
        {0x1561D276DDFDC00A, 0x71DD41A08F48AF40}, // e =    -6, k =   51
        {0x1ABA4714957D300D, 0x0E549208B31ADB10}, // e =    -4, k =   52
        {0x10B46C6CDD6E3E08, 0x28F4DB456FF0C8EA}, // e =    -1, k =   53
        {0x14E1878814C9CD8A, 0x33321216CBECFB25}, // e =     1, k =   54
        {0x1A19E96A19FC40EC, 0xBFFE969C7EE839EE}, // e =     3, k =   55
        {0x105031E2503DA893, 0xF7FF1E21CF512435}, // e =     6, k =   56
        {0x14643E5AE44D12B8, 0xF5FEE5AA43256D42}, // e =     8, k =   57
        {0x197D4DF19D605767, 0x337E9F14D3EEC893}, // e =    10, k =   58
        {0x1FDCA16E04B86D41, 0x005E46DA08EA7AB7}, // e =    12, k =   59
        {0x13E9E4E4C2F34448, 0xA03AEC4845928CB3}, // e =    15, k =   60
        {0x18E45E1DF3B0155A, 0xC849A75A56F72FDF}, // e =    17, k =   61
        {0x1F1D75A5709C1AB1, 0x7A5C1130ECB4FBD7}, // e =    19, k =   62
        {0x13726987666190AE, 0xEC798ABE93F11D66}, // e =    22, k =   63
        {0x184F03E93FF9F4DA, 0xA797ED6E38ED64C0}, // e =    24, k =   64
        {0x1E62C4E38FF87211, 0x517DE8C9C728BDF0}, // e =    26, k =   65
        {0x12FDBB0E39FB474A, 0xD2EEB17E1C7976B6}, // e =    29, k =   66
        {0x17BD29D1C87A191D, 0x87AA5DDDA397D463}, // e =    31, k =   67
        {0x1DAC74463A989F64, 0xE994F5550C7DC97C}, // e =    33, k =   68
        {0x128BC8ABE49F639F, 0x11FD195527CE9DEE}, // e =    36, k =   69
        {0x172EBAD6DDC73C86, 0xD67C5FAA71C24569}, // e =    38, k =   70
        {0x1CFA698C95390BA8, 0x8C1B77950E32D6C3}, // e =    40, k =   71
        {0x121C81F7DD43A749, 0x57912ABD28DFC63A}, // e =    43, k =   72
        {0x16A3A275D494911B, 0xAD75756C7317B7C9}, // e =    45, k =   73
        {0x1C4C8B1349B9B562, 0x98D2D2C78FDDA5BB}, // e =    47, k =   74
        {0x11AFD6EC0E14115D, 0x9F83C3BCB9EA8795}, // e =    50, k =   75
        {0x161BCCA7119915B5, 0x0764B4ABE865297A}, // e =    52, k =   76
        {0x1BA2BFD0D5FF5B22, 0x493DE1D6E27E73D8}, // e =    54, k =   77
        {0x1145B7E285BF98F5, 0x6DC6AD264D8F0867}, // e =    57, k =   78
        {0x159725DB272F7F32, 0xC938586FE0F2CA81}, // e =    59, k =   79
        {0x1AFCEF51F0FB5EFF, 0x7B866E8BD92F7D21}, // e =    61, k =   80
        {0x10DE1593369D1B5F, 0xAD34051767BDAE35}, // e =    64, k =   81
        {0x15159AF804446237, 0x9881065D41AD19C2}, // e =    66, k =   82
        {0x1A5B01B605557AC5, 0x7EA147F492186033}, // e =    68, k =   83
        {0x1078E111C3556CBB, 0x6F24CCF8DB4F3C20}, // e =    71, k =   84
        {0x14971956342AC7EA, 0x4AEE003712230B28}, // e =    73, k =   85
        {0x19BCDFABC13579E4, 0xDDA98044D6ABCDF1}, // e =    75, k =   86
        {0x10160BCB58C16C2F, 0x0A89F02B062B60B7}, // e =    78, k =   87
        {0x141B8EBE2EF1C73A, 0xCD2C6C35C7B638E5}, // e =    80, k =   88
        {0x1922726DBAAE3909, 0x8077874339A3C71E}, // e =    82, k =   89
        {0x1F6B0F092959C74B, 0xE0956914080CB8E5}, // e =    84, k =   90
        {0x13A2E965B9D81C8F, 0x6C5D61AC8507F38F}, // e =    87, k =   91
        {0x188BA3BF284E23B3, 0x4774BA17A649F073}, // e =    89, k =   92
        {0x1EAE8CAEF261ACA0, 0x1951E89D8FDC6C90}, // e =    91, k =   93
        {0x132D17ED577D0BE4, 0x0FD3316279E9C3DA}, // e =    94, k =   94
        {0x17F85DE8AD5C4EDD, 0x13C7FDBB186434D0}, // e =    96, k =   95
        {0x1DF67562D8B36294, 0x58B9FD29DE7D4204}, // e =    98, k =   96
        {0x12BA095DC7701D9C, 0xB7743E3A2B0E4943}, // e =   101, k =   97
        {0x17688BB5394C2503, 0xE5514DC8B5D1DB93}, // e =   103, k =   98
        {0x1D42AEA2879F2E44, 0xDEA5A13AE3465278}, // e =   105, k =   99
        {0x1249AD2594C37CEB, 0x0B2784C4CE0BF38B}, // e =   108, k =  100
        {0x16DC186EF9F45C25, 0xCDF165F6018EF06E}, // e =   110, k =  101
        {0x1C931E8AB871732F, 0x416DBF7381F2AC89}, // e =   112, k =  102
        {0x11DBF316B346E7FD, 0x88E497A83137ABD6}, // e =   115, k =  103
        {0x1652EFDC6018A1FC, 0xEB1DBD923D8596CB}, // e =   117, k =  104
        {0x1BE7ABD3781ECA7C, 0x25E52CF6CCE6FC7E}, // e =   119, k =  105
        {0x1170CB642B133E8D, 0x97AF3C1A40105DCF}, // e =   122, k =  106
        {0x15CCFE3D35D80E30, 0xFD9B0B20D0147543}, // e =   124, k =  107
        {0x1B403DCC834E11BD, 0x3D01CDE904199293}, // e =   126, k =  108
        {0x1108269FD210CB16, 0x462120B1A28FFB9C}, // e =   129, k =  109
        {0x154A3047C694FDDB, 0xD7A968DE0B33FA83}, // e =   131, k =  110
        {0x1A9CBC59B83A3D52, 0xCD93C3158E00F924}, // e =   133, k =  111
        {0x10A1F5B813246653, 0xC07C59ED78C09BB7}, // e =   136, k =  112
        {0x14CA732617ED7FE8, 0xB09B7068D6F0C2A4}, // e =   138, k =  113
        {0x19FD0FEF9DE8DFE2, 0xDCC24C830CACF34D}, // e =   140, k =  114
        {0x103E29F5C2B18BED, 0xC9F96FD1E7EC1810}, // e =   143, k =  115
        {0x144DB473335DEEE9, 0x3C77CBC661E71E14}, // e =   145, k =  116
        {0x1961219000356AA3, 0x8B95BEB7FA60E599}, // e =   147, k =  117
        {0x1FB969F40042C54C, 0x6E7B2E65F8F91EFF}, // e =   149, k =  118
        {0x13D3E2388029BB4F, 0xC50CFCFFBB9BB360}, // e =   152, k =  119
        {0x18C8DAC6A0342A23, 0xB6503C3FAA82A038}, // e =   154, k =  120
        {0x1EFB1178484134AC, 0xA3E44B4F95234845}, // e =   156, k =  121
        {0x135CEAEB2D28C0EB, 0xE66EAF11BD360D2C}, // e =   159, k =  122
        {0x183425A5F872F126, 0xE00A5AD62C839076}, // e =   161, k =  123
        {0x1E412F0F768FAD70, 0x980CF18BB7A47494}, // e =   163, k =  124
        {0x12E8BD69AA19CC66, 0x5F0816F752C6C8DD}, // e =   166, k =  125
        {0x17A2ECC414A03F7F, 0xF6CA1CB527787B14}, // e =   168, k =  126
        {0x1D8BA7F519C84F5F, 0xF47CA3E2715699D8}, // e =   170, k =  127
        {0x127748F9301D319B, 0xF8CDE66D86D62027}, // e =   173, k =  128
        {0x17151B377C247E02, 0xF7016008E88BA831}, // e =   175, k =  129
        {0x1CDA62055B2D9D83, 0xB4C1B80B22AE923D}, // e =   177, k =  130
        {0x12087D4358FC8272, 0x50F91306F5AD1B66}, // e =   180, k =  131
        {0x168A9C942F3BA30E, 0xE53757C8B3186240}, // e =   182, k =  132
        {0x1C2D43B93B0A8BD2, 0x9E852DBADFDE7AD0}, // e =   184, k =  133
        {0x119C4A53C4E69763, 0xA3133C94CBEB0CC2}, // e =   187, k =  134
        {0x16035CE8B6203D3C, 0x8BD80BB9FEE5CFF2}, // e =   189, k =  135
        {0x1B843422E3A84C8B, 0xAECE0EA87E9F43EF}, // e =   191, k =  136
        {0x1132A095CE492FD7, 0x4D40C9294F238A76}, // e =   194, k =  137
        {0x157F48BB41DB7BCD, 0x2090FB73A2EC6D13}, // e =   196, k =  138
        {0x1ADF1AEA12525AC0, 0x68B53A508BA78857}, // e =   198, k =  139
        {0x10CB70D24B7378B8, 0x417144725748B537}, // e =   201, k =  140
        {0x14FE4D06DE5056E6, 0x51CD958EED1AE284}, // e =   203, k =  141
        {0x1A3DE04895E46C9F, 0xE640FAF2A8619B25}, // e =   205, k =  142
        {0x1066AC2D5DAEC3E3, 0xEFE89CD7A93D00F8}, // e =   208, k =  143
        {0x14805738B51A74DC, 0xEBE2C40D938C4135}, // e =   210, k =  144
        {0x19A06D06E2611214, 0x26DB7510F86F5182}, // e =   212, k =  145
        {0x100444244D7CAB4C, 0x9849292A9B4592F2}, // e =   215, k =  146
        {0x1405552D60DBD61F, 0xBE5B73754216F7AE}, // e =   217, k =  147
        {0x1906AA78B912CBA7, 0xADF25052929CB599}, // e =   219, k =  148
        {0x1F485516E7577E91, 0x996EE4673743E300}, // e =   221, k =  149
        {0x138D352E5096AF1A, 0xFFE54EC0828A6DE0}, // e =   224, k =  150
        {0x18708279E4BC5AE1, 0xBFDEA270A32D0958}, // e =   226, k =  151
        {0x1E8CA3185DEB719A, 0x2FD64B0CCBF84BAE}, // e =   228, k =  152
        {0x1317E5EF3AB32700, 0x5DE5EEE7FF7B2F4D}, // e =   231, k =  153
        {0x17DDDF6B095FF0C0, 0x755F6AA1FF59FB20}, // e =   233, k =  154
        {0x1DD55745CBB7ECF0, 0x92B7454A7F3079E8}, // e =   235, k =  155
        {0x12A5568B9F52F416, 0x5BB28B4E8F7E4C31}, // e =   238, k =  156
        {0x174EAC2E8727B11B, 0xF29F2E22335DDF3D}, // e =   240, k =  157
        {0x1D22573A28F19D62, 0xEF46F9AAC035570C}, // e =   242, k =  158
        {0x123576845997025D, 0xD58C5C0AB8215668}, // e =   245, k =  159
        {0x16C2D4256FFCC2F5, 0x4AEF730D6629AC02}, // e =   247, k =  160
        {0x1C73892ECBFBF3B2, 0x9DAB4FD0BFB41702}, // e =   249, k =  161
        {0x11C835BD3F7D784F, 0xA28B11E277D08E61}, // e =   252, k =  162
        {0x163A432C8F5CD663, 0x8B2DD65B15C4B1FA}, // e =   254, k =  163
        {0x1BC8D3F7B3340BFC, 0x6DF94BF1DB35DE78}, // e =   256, k =  164
        {0x115D847AD000877D, 0xC4BBCF772901AB0B}, // e =   259, k =  165
        {0x15B4E5998400A95D, 0x35EAC354F34215CE}, // e =   261, k =  166
        {0x1B221EFFE500D3B4, 0x8365742A30129B41}, // e =   263, k =  167
        {0x10F5535FEF208450, 0xD21F689A5E0BA109}, // e =   266, k =  168
        {0x1532A837EAE8A565, 0x06A742C0F58E894B}, // e =   268, k =  169
        {0x1A7F5245E5A2CEBE, 0x4851137132F22B9E}, // e =   270, k =  170
        {0x108F936BAF85C136, 0xED32AC26BFD75B43}, // e =   273, k =  171
        {0x14B378469B673184, 0xA87F57306FCD3213}, // e =   275, k =  172
        {0x19E056584240FDE5, 0xD29F2CFC8BC07E98}, // e =   277, k =  173
        {0x102C35F729689EAF, 0xA3A37C1DD7584F1F}, // e =   280, k =  174
        {0x14374374F3C2C65B, 0x8C8C5B254D2E62E7}, // e =   282, k =  175
        {0x1945145230B377F2, 0x6FAF71EEA079FBA0}, // e =   284, k =  176
        {0x1F965966BCE055EF, 0x0B9B4E6A48987A88}, // e =   286, k =  177
        {0x13BDF7E0360C35B5, 0x674111026D5F4C95}, // e =   289, k =  178
        {0x18AD75D8438F4322, 0xC111554308B71FBB}, // e =   291, k =  179
        {0x1ED8D34E547313EB, 0x7155AA93CAE4E7A9}, // e =   293, k =  180
        {0x13478410F4C7EC73, 0x26D58A9C5ECF10CA}, // e =   296, k =  181
        {0x1819651531F9E78F, 0xF08AED437682D4FC}, // e =   298, k =  182
        {0x1E1FBE5A7E786173, 0xECADA89454238A3B}, // e =   300, k =  183
        {0x12D3D6F88F0B3CE8, 0x73EC895CB4963665}, // e =   303, k =  184
        {0x1788CCB6B2CE0C22, 0x90E7ABB3E1BBC3FE}, // e =   305, k =  185
        {0x1D6AFFE45F818F2B, 0x352196A0DA2AB4FE}, // e =   307, k =  186
        {0x1262DFEEBBB0F97B, 0x0134FE24885AB11F}, // e =   310, k =  187
        {0x16FB97EA6A9D37D9, 0xC1823DADAA715D66}, // e =   312, k =  188
        {0x1CBA7DE5054485D0, 0x31E2CD19150DB4C0}, // e =   314, k =  189
        {0x11F48EAF234AD3A2, 0x1F2DC02FAD2890F8}, // e =   317, k =  190
        {0x1671B25AEC1D888A, 0xA6F9303B9872B536}, // e =   319, k =  191
        {0x1C0E1EF1A724EAAD, 0x50B77C4A7E8F6283}, // e =   321, k =  192
        {0x1188D357087712AC, 0x5272ADAE8F199D92}, // e =   324, k =  193
        {0x15EB082CCA94D757, 0x670F591A32E004F7}, // e =   326, k =  194
        {0x1B65CA37FD3A0D2D, 0x40D32F60BF980634}, // e =   328, k =  195
        {0x111F9E62FE44483C, 0x4883FD9C77BF03E1}, // e =   331, k =  196
        {0x156785FBBDD55A4B, 0x5AA4FD0395AEC4D9}, // e =   333, k =  197
        {0x1AC1677AAD4AB0DE, 0x314E3C447B1A760F}, // e =   335, k =  198
        {0x10B8E0ACAC4EAE8A, 0xDED0E5AACCF089CA}, // e =   338, k =  199
        {0x14E718D7D7625A2D, 0x96851F15802CAC3C}, // e =   340, k =  200
        {0x1A20DF0DCD3AF0B8, 0xFC2666DAE037D74B}, // e =   342, k =  201
        {0x10548B68A044D673, 0x9D980048CC22E68F}, // e =   345, k =  202
        {0x1469AE42C8560C10, 0x84FE005AFF2BA033}, // e =   347, k =  203
        {0x198419D37A6B8F14, 0xA63D8071BEF6883F}, // e =   349, k =  204
        {0x1FE52048590672D9, 0xCFCCE08E2EB42A4F}, // e =   351, k =  205
        {0x13EF342D37A407C8, 0x21E00C58DD309A71}, // e =   354, k =  206
        {0x18EB0138858D09BA, 0x2A580F6F147CC10E}, // e =   356, k =  207
        {0x1F25C186A6F04C28, 0xB4EE134AD99BF151}, // e =   358, k =  208
        {0x137798F428562F99, 0x7114CC0EC80176D3}, // e =   361, k =  209
        {0x18557F31326BBB7F, 0xCD59FF127A01D487}, // e =   363, k =  210
        {0x1E6ADEFD7F06AA5F, 0xC0B07ED7188249A9}, // e =   365, k =  211
        {0x1302CB5E6F642A7B, 0xD86E4F466F516E0A}, // e =   368, k =  212
        {0x17C37E360B3D351A, 0xCE89E3180B25C98C}, // e =   370, k =  213
        {0x1DB45DC38E0C8261, 0x822C5BDE0DEF3BEF}, // e =   372, k =  214
        {0x1290BA9A38C7D17C, 0xF15BB96AC8B58576}, // e =   375, k =  215
        {0x1734E940C6F9C5DC, 0x2DB2A7C57AE2E6D3}, // e =   377, k =  216
        {0x1D022390F8B83753, 0x391F51B6D99BA087}, // e =   379, k =  217
        {0x1221563A9B732294, 0x03B3931248014455}, // e =   382, k =  218
        {0x16A9ABC9424FEB39, 0x04A077D6DA01956A}, // e =   384, k =  219
        {0x1C5416BB92E3E607, 0x45C895CC9081FAC4}, // e =   386, k =  220
        {0x11B48E353BCE6FC4, 0x8B9D5D9FDA513CBB}, // e =   389, k =  221
        {0x1621B1C28AC20BB5, 0xAE84B507D0E58BE9}, // e =   391, k =  222
        {0x1BAA1E332D728EA3, 0x1A25E249C51EEEE4}, // e =   393, k =  223
        {0x114A52DFFC679925, 0xF057AD6E1B33554E}, // e =   396, k =  224
        {0x159CE797FB817F6F, 0x6C6D98C9A2002AA2}, // e =   398, k =  225
        {0x1B04217DFA61DF4B, 0x4788FEFC0A80354A}, // e =   400, k =  226
        {0x10E294EEBC7D2B8F, 0x0CB59F5D8690214F}, // e =   403, k =  227
        {0x151B3A2A6B9C7672, 0xCFE30734E83429A2}, // e =   405, k =  228
        {0x1A6208B50683940F, 0x83DBC9022241340B}, // e =   407, k =  229
        {0x107D457124123C89, 0xB2695DA15568C087}, // e =   410, k =  230
        {0x149C96CD6D16CBAC, 0x1F03B509AAC2F0A8}, // e =   412, k =  231
        {0x19C3BC80C85C7E97, 0x26C4A24C1573ACD2}, // e =   414, k =  232
        {0x101A55D07D39CF1E, 0x783AE56F8D684C04}, // e =   417, k =  233
        {0x1420EB449C8842E6, 0x16499ECB70C25F04}, // e =   419, k =  234
        {0x19292615C3AA539F, 0x9BDC067E4CF2F6C5}, // e =   421, k =  235
        {0x1F736F9B3494E887, 0x82D3081DE02FB477}, // e =   423, k =  236
        {0x13A825C100DD1154, 0xB1C3E512AC1DD0CA}, // e =   426, k =  237
        {0x18922F31411455A9, 0xDE34DE57572544FD}, // e =   428, k =  238
        {0x1EB6BAFD91596B14, 0x55C215ED2CEE963C}, // e =   430, k =  239
        {0x133234DE7AD7E2EC, 0xB5994DB43C151DE6}, // e =   433, k =  240
        {0x17FEC216198DDBA7, 0xE2FFA1214B1A655F}, // e =   435, k =  241
        {0x1DFE729B9FF15291, 0xDBBF89699DE0FEB7}, // e =   437, k =  242
        {0x12BF07A143F6D39B, 0x2957B5E202AC9F32}, // e =   440, k =  243
        {0x176EC98994F48881, 0xF3ADA35A8357C6FF}, // e =   442, k =  244
        {0x1D4A7BEBFA31AAA2, 0x70990C31242DB8BE}, // e =   444, k =  245
        {0x124E8D737C5F0AA5, 0x865FA79EB69C9377}, // e =   447, k =  246
        {0x16E230D05B76CD4E, 0xE7F791866443B855}, // e =   449, k =  247
        {0x1C9ABD04725480A2, 0xA1F575E7FD54A66A}, // e =   451, k =  248
        {0x11E0B622C774D065, 0xA53969B0FE54E802}, // e =   454, k =  249
        {0x1658E3AB7952047F, 0x0E87C41D3DEA2203}, // e =   456, k =  250
        {0x1BEF1C9657A6859E, 0xD229B5248D64AA83}, // e =   458, k =  251
        {0x117571DDF6C81383, 0x435A1136D85EEA92}, // e =   461, k =  252
        {0x15D2CE55747A1864, 0x143095848E76A537}, // e =   463, k =  253
        {0x1B4781EAD1989E7D, 0x193CBAE5B2144E84}, // e =   465, k =  254
        {0x110CB132C2FF630E, 0x2FC5F4CF8F4CB113}, // e =   468, k =  255
        {0x154FDD7F73BF3BD1, 0xBBB77203731FDD57}, // e =   470, k =  256
        {0x1AA3D4DF50AF0AC6, 0x2AA54E844FE7D4AD}, // e =   472, k =  257
        {0x10A6650B926D66BB, 0xDAA75112B1F0E4EC}, // e =   475, k =  258
        {0x14CFFE4E7708C06A, 0xD15125575E6D1E27}, // e =   477, k =  259
        {0x1A03FDE214CAF085, 0x85A56EAD360865B1}, // e =   479, k =  260
        {0x10427EAD4CFED653, 0x7387652C41C53F8F}, // e =   482, k =  261
        {0x14531E58A03E8BE8, 0x50693E7752368F72}, // e =   484, k =  262
        {0x1967E5EEC84E2EE2, 0x64838E1526C4334F}, // e =   486, k =  263
        {0x1FC1DF6A7A61BA9A, 0xFDA4719A70754023}, // e =   488, k =  264
        {0x13D92BA28C7D14A0, 0xDE86C70086494816}, // e =   491, k =  265
        {0x18CF768B2F9C59C9, 0x162878C0A7DB9A1B}, // e =   493, k =  266
        {0x1F03542DFB83703B, 0x5BB296F0D1D280A2}, // e =   495, k =  267
        {0x1362149CBD322625, 0x194F9E5683239065}, // e =   498, k =  268
        {0x183A99C3EC7EAFAE, 0x5FA385EC23EC747F}, // e =   500, k =  269
        {0x1E494034E79E5B99, 0xF78C67672CE7919E}, // e =   502, k =  270
        {0x12EDC82110C2F940, 0x3AB7C0A07C10BB03}, // e =   505, k =  271
        {0x17A93A2954F3B790, 0x4965B0C89B14E9C4}, // e =   507, k =  272
        {0x1D9388B3AA30A574, 0x5BBF1CFAC1DA2434}, // e =   509, k =  273
        {0x127C35704A5E6768, 0xB957721CB92856A1}, // e =   512, k =  274
        {0x171B42CC5CF60142, 0xE7AD4EA3E7726C49}, // e =   514, k =  275
        {0x1CE2137F74338193, 0xA198A24CE14F075B}, // e =   516, k =  276
        {0x120D4C2FA8A030FC, 0x44FF65700CD16499}, // e =   519, k =  277
        {0x16909F3B92C83D3B, 0x563F3ECC1005BDBF}, // e =   521, k =  278
        {0x1C34C70A777A4C8A, 0x2BCF0E7F14072D2F}, // e =   523, k =  279
        {0x11A0FC668AAC6FD6, 0x5B61690F6C847C3E}, // e =   526, k =  280
        {0x16093B802D578BCB, 0xF239C35347A59B4D}, // e =   528, k =  281
        {0x1B8B8A6038AD6EBE, 0xEEC83428198F0220}, // e =   530, k =  282
        {0x1137367C236C6537, 0x553D20990FF96154}, // e =   533, k =  283
        {0x1585041B2C477E85, 0x2A8C68BF53F7B9A9}, // e =   535, k =  284
        {0x1AE64521F7595E26, 0x752F82EF28F5A813}, // e =   537, k =  285
        {0x10CFEB353A97DAD8, 0x093DB1D57999890C}, // e =   540, k =  286
        {0x1503E602893DD18E, 0x0B8D1E4AD7FFEB4F}, // e =   542, k =  287
        {0x1A44DF832B8D45F1, 0x8E7065DD8DFFE623}, // e =   544, k =  288
        {0x106B0BB1FB384BB6, 0xF9063FAA78BFEFD6}, // e =   547, k =  289
        {0x1485CE9E7A065EA4, 0xB747CF9516EFEBCB}, // e =   549, k =  290
        {0x19A742461887F64D, 0xE519C37A5CABE6BE}, // e =   551, k =  291
        {0x1008896BCF54F9F0, 0xAF301A2C79EB7037}, // e =   554, k =  292
        {0x140AABC6C32A386C, 0xDAFC20B798664C44}, // e =   556, k =  293
        {0x190D56B873F4C688, 0x11BB28E57E7FDF55}, // e =   558, k =  294
        {0x1F50AC6690F1F82A, 0x1629F31EDE1FD72B}, // e =   560, k =  295
        {0x13926BC01A973B1A, 0x4DDA37F34AD3E67B}, // e =   563, k =  296
        {0x187706B0213D09E0, 0xE150C5F01D88E01A}, // e =   565, k =  297
        {0x1E94C85C298C4C59, 0x19A4F76C24EB1820}, // e =   567, k =  298
        {0x131CFD3999F7AFB7, 0xB0071AA39712EF14}, // e =   570, k =  299
        {0x17E43C8800759BA5, 0x9C08E14C7CD7AAD9}, // e =   572, k =  300
        {0x1DDD4BAA0093028F, 0x030B199F9C0D958F}, // e =   574, k =  301
        {0x12AA4F4A405BE199, 0x61E6F003C1887D7A}, // e =   577, k =  302
        {0x1754E31CD072D9FF, 0xBA60AC04B1EA9CD8}, // e =   579, k =  303
        {0x1D2A1BE4048F907F, 0xA8F8D705DE65440E}, // e =   581, k =  304
        {0x123A516E82D9BA4F, 0xC99B8663AAFF4A89}, // e =   584, k =  305
        {0x16C8E5CA239028E3, 0xBC0267FC95BF1D2B}, // e =   586, k =  306
        {0x1C7B1F3CAC74331C, 0xAB0301FBBB2EE475}, // e =   588, k =  307
        {0x11CCF385EBC89FF1, 0xEAE1E13D54FD4ECA}, // e =   591, k =  308
        {0x1640306766BAC7EE, 0x659A598CAA3CA27C}, // e =   593, k =  309
        {0x1BD03C81406979E9, 0xFF00EFEFD4CBCB1B}, // e =   595, k =  310
        {0x116225D0C841EC32, 0x3F6095F5E4FF5EF1}, // e =   598, k =  311
        {0x15BAAF44FA52673E, 0xCF38BB735E3F36AD}, // e =   600, k =  312
        {0x1B295B1638E7010E, 0x8306EA5035CF0458}, // e =   602, k =  313
        {0x10F9D8EDE39060A9, 0x11E4527221A162B7}, // e =   605, k =  314
        {0x15384F295C7478D3, 0x565D670EAA09BB65}, // e =   607, k =  315
        {0x1A8662F3B3919708, 0x2BF4C0D2548C2A3E}, // e =   609, k =  316
        {0x1093FDD8503AFE65, 0x1B78F88374D79A67}, // e =   612, k =  317
        {0x14B8FD4E6449BDFE, 0x625736A4520D8101}, // e =   614, k =  318
        {0x19E73CA1FD5C2D7D, 0xFAED044D6690E141}, // e =   616, k =  319
        {0x103085E53E599C6E, 0xBCD422B0601A8CC9}, // e =   619, k =  320
        {0x143CA75E8DF0038A, 0x6C092B5C78212FFB}, // e =   621, k =  321
        {0x194BD136316C046D, 0x070B763396297BF9}, // e =   623, k =  322
        {0x1F9EC583BDC70588, 0x48CE53C07BB3DAF7}, // e =   625, k =  323
        {0x13C33B72569C6375, 0x2D80F4584D5068DB}, // e =   628, k =  324
        {0x18B40A4EEC437C52, 0x78E1316E60A48311}, // e =   630, k =  325
    };

    RYU_ASSERT(k >= MinDecExp);
    RYU_ASSERT(k <= MaxDecExp);
    return Pow5[static_cast<unsigned>(k - MinDecExp)];
}

#if defined(__SIZEOF_INT128__)

static inline uint64_t MulShift(uint64_t m, const Uint64x2* mul, int j)
{
    __extension__ using uint128_t = unsigned __int128;

    RYU_ASSERT(j >= 64);
    RYU_ASSERT(j <= 127);

    const uint128_t b0 = uint128_t{m} * mul->lo;
    const uint128_t b2 = uint128_t{m} * mul->hi;

    // We need shift = j - 64 here.
    // Since 64 < j < 128, this is equivalent to shift = (j - 64) % 64 = j % 64.
    // When written as shift = j & 63, clang can optimize the 128-bit shift into
    // a simple funnel shift.
    const int shift = j & 63;
    return static_cast<uint64_t>((b2 + static_cast<uint64_t>(b0 >> 64)) >> shift);
}

#elif defined(_MSC_VER) && defined(_M_X64)

static inline uint64_t MulShift(uint64_t m, const Uint64x2* mul, int j)
{
    RYU_ASSERT(j >= 64);
    RYU_ASSERT(j <= 127);

    uint64_t b0_hi;
    uint64_t b0_lo = _umul128(m, mul->lo, &b0_hi);
    uint64_t b2_hi;
    uint64_t b2_lo = _umul128(m, mul->hi, &b2_hi);
    static_cast<void>(b0_lo);

    // b2 + (b0 >> 64)
    // b2_lo += b0_hi;
    // b2_hi += b2_lo < b0_hi;
    _addcarry_u64(_addcarry_u64(0, b2_lo, b0_hi, &b2_lo), b2_hi, 0, &b2_hi);

    // We need shift = j - 64 here.
    // For the __shiftright128 intrinsic, the shift value is always modulo 64.
    // Since (j - 64) % 64 = j, we can simply use j here.
    return __shiftright128(b2_lo, b2_hi, static_cast<unsigned char>(j));
}

#else

static inline Uint64x2 Mul128(uint64_t a, uint64_t b)
{
    const uint64_t b00 = uint64_t{Lo32(a)} * Lo32(b);
    const uint64_t b01 = uint64_t{Lo32(a)} * Hi32(b);
    const uint64_t b10 = uint64_t{Hi32(a)} * Lo32(b);
    const uint64_t b11 = uint64_t{Hi32(a)} * Hi32(b);

    const uint64_t mid1 = b10 + Hi32(b00);
    const uint64_t mid2 = b01 + Lo32(mid1);

    const uint64_t hi = b11 + Hi32(mid1) + Hi32(mid2);
    const uint64_t lo = Lo32(b00) | uint64_t{Lo32(mid2)} << 32;
    return {hi, lo};
}

static inline uint64_t ShiftRight128(uint64_t lo, uint64_t hi, int n)
{
    // For the __shiftright128 intrinsic, the shift value is always modulo 64.
    // In the current implementation of the double-precision version of Ryu, the
    // shift value is always < 64.
    // Check this here in case a future change requires larger shift values. In
    // this case this function needs to be adjusted.
    RYU_ASSERT(n >= 1);
    RYU_ASSERT(n <= 63);

    const int lshift = -n & 63;
    const int rshift =  n & 63;
    return (hi << lshift) | (lo >> rshift);
}

static inline uint64_t MulShift(uint64_t m, const Uint64x2* mul, int j)
{
    auto b0 = Mul128(m, mul->lo);
    auto b2 = Mul128(m, mul->hi);

    // b2 + (b0 >> 64)
    b2.lo += b0.hi;
    b2.hi += b2.lo < b0.hi;

    const int shift = j & 63;
    return ShiftRight128(b2.lo, b2.hi, shift);
}

#endif

static inline void MulPow5DivPow2_Double(uint64_t u, uint64_t v, uint64_t w, int e5, int e2, uint64_t& a, uint64_t& b, uint64_t& c)
{
    // j >= 118 and m has at most 53 + 2 = 55 bits.
    // The product along with the subsequent shift therefore requires
    // 55 + 125 - 118 = 62 bits.

    const auto k = FloorLog2Pow5(e5) + 1 - BitsPerPow5_Double;
    const auto j = e2 - k;
    RYU_ASSERT(j >= BitsPerPow5_Double - 7); // 118 - 64 = 54
    RYU_ASSERT(j <= BitsPerPow5_Double - 1); // 124 - 64 = 60

    const auto pow5 = ComputePow5_Double(e5);

    a = MulShift(u, &pow5, j);
    b = MulShift(v, &pow5, j);
    c = MulShift(w, &pow5, j);
}

// Returns whether value is divisible by 5^e5
static inline bool MultipleOfPow5(uint64_t value, int e5)
{
    struct MulCmp {
        uint64_t mul;
        uint64_t cmp;
    };

    static constexpr MulCmp Mod5[] = {
        {0x0000000000000001u, 0xFFFFFFFFFFFFFFFFu}, // 5^0
        {0xCCCCCCCCCCCCCCCDu, 0x3333333333333333u}, // 5^1
        {0x8F5C28F5C28F5C29u, 0x0A3D70A3D70A3D70u}, // 5^2
        {0x1CAC083126E978D5u, 0x020C49BA5E353F7Cu}, // 5^3
        {0xD288CE703AFB7E91u, 0x0068DB8BAC710CB2u}, // 5^4
        {0x5D4E8FB00BCBE61Du, 0x0014F8B588E368F0u}, // 5^5
        {0x790FB65668C26139u, 0x000431BDE82D7B63u}, // 5^6
        {0xE5032477AE8D46A5u, 0x0000D6BF94D5E57Au}, // 5^7
        {0xC767074B22E90E21u, 0x00002AF31DC46118u}, // 5^8
        {0x8E47CE423A2E9C6Du, 0x0000089705F4136Bu}, // 5^9
        {0x4FA7F60D3ED61F49u, 0x000001B7CDFD9D7Bu}, // 5^10
        {0x0FEE64690C913975u, 0x00000057F5FF85E5u}, // 5^11
        {0x3662E0E1CF503EB1u, 0x000000119799812Du}, // 5^12
        {0xA47A2CF9F6433FBDu, 0x0000000384B84D09u}, // 5^13
        {0x54186F653140A659u, 0x00000000B424DC35u}, // 5^14
        {0x7738164770402145u, 0x0000000024075F3Du}, // 5^15
        {0xE4A4D1417CD9A041u, 0x000000000734ACA5u}, // 5^16
        {0xC75429D9E5C5200Du, 0x000000000170EF54u}, // 5^17
        {0xC1773B91FAC10669u, 0x000000000049C977u}, // 5^18
        {0x26B172506559CE15u, 0x00000000000EC1E4u}, // 5^19
        {0xD489E3A9ADDEC2D1u, 0x000000000002F394u}, // 5^20
        {0x90E860BB892C8D5Du, 0x000000000000971Du}, // 5^21
        {0x502E79BF1B6F4F79u, 0x0000000000001E39u}, // 5^22
        {0xDCD618596BE30FE5u, 0x000000000000060Bu}, // 5^23
        {0x2C2AD1AB7BFA3661u, 0x0000000000000135u}, // 5^24
    };

    RYU_ASSERT(e5 >= 0);
    RYU_ASSERT(e5 <= 24);
    const auto m5 = Mod5[static_cast<unsigned>(e5)];

    return value * m5.mul <= m5.cmp;
}

// Returns whether value is divisible by 2^e2
static inline bool MultipleOfPow2(uint64_t value, int e2)
{
    RYU_ASSERT(e2 >= 0);
    RYU_ASSERT(e2 <= 63);

    return (value & ((uint64_t{1} << e2) - 1)) == 0;
}

namespace {
struct ToDecimalResultDouble {
    uint64_t digits; // num_digits <= 17
    int exponent;
};
}

static inline ToDecimalResultDouble ToDecimal(double value)
{
    RYU_ASSERT(Double(value).IsFinite());
    RYU_ASSERT(value > 0);

    //
    // Step 1:
    // Decode the floating point number, and unify normalized and subnormal cases.
    //

    const Double ieee_value(value);

    // Decode bits into mantissa, and exponent.
    const uint64_t ieee_mantissa = ieee_value.PhysicalSignificand();
    const uint64_t ieee_exponent = ieee_value.PhysicalExponent();

    uint64_t m2;
    int e2;
    if (ieee_exponent == 0)
    {
        m2 = ieee_mantissa;
        e2 = 1 - Double::ExponentBias;
    }
    else
    {
        m2 = Double::HiddenBit | ieee_mantissa;
        e2 = static_cast<int>(ieee_exponent) - Double::ExponentBias;

        if /*unlikely*/ ((0 <= -e2 && -e2 < Double::SignificandSize) && MultipleOfPow2(m2, -e2))
        {
            // Since 2^52 <= m2 < 2^53 and 0 <= -e2 <= 52:
            //  1 <= value = m2 / 2^-e2 < 2^53.
            // Since m2 is divisible by 2^-e2, value is an integer.
            return {m2 >> -e2, 0};
        }
    }

    const bool is_even = (m2 % 2) == 0;
    const bool accept_lower = is_even;
    const bool accept_upper = is_even;

    //
    // Step 2:
    // Determine the interval of valid decimal representations.
    //

    const uint32_t lower_boundary_is_closer = (ieee_mantissa == 0 && ieee_exponent > 1);

    e2 -= 2;
    const uint64_t u = 4 * m2 - 2 + lower_boundary_is_closer;
    const uint64_t v = 4 * m2;
    const uint64_t w = 4 * m2 + 2;

    //
    // Step 3:
    // Convert to a decimal power base.
    //

    int e10;

    bool za = false; // a[0, ..., i-1] == 0
    bool zb = false; // b[0, ..., i-1] == 0
    bool zc = false; // c[0, ..., i-1] == 0

    if (e2 >= 0)
    {
        // We need
        //  (a,b,c) = (u,v,w) * 2^e2
        // and we need to remove at least q' = log_10(2^e2) digits from the
        // scaled values a,b,c, i.e. we want to compute
        //  (a,b,c) = (u,v,w) * 2^e2 / 10^(q')
        //          = (u,v,w) * 2^e2 / 10^(e10)
        //          = (u,v,w) * 5^(-e10) / 2^(e10 - e2)
        //
        // However, to correctly round the result we need to know the value of
        // the last removed digit. We therefore remove only q = q' - 1 digits in
        // the first step and make sure that we execute the loop below at least
        // once and determine the correct value of the last removed digit.

        const int q = FloorLog10Pow2(e2) - (e2 > 3); // == max(0, q' - 1)
        RYU_ASSERT(q >= 0);

        e10 = q;
        RYU_ASSERT(e10 >= 0);
        RYU_ASSERT(e10 - e2 <= 0);

        // Determine whether all the removed digits are 0.
        //
        // Z(x,e2,q) = (x * 2^e2) % 10^q == 0
        //           = p10(x * 2^e2) >= q
        //           = min(p2(x) + p2(e2), p5(x)) >= q
        //           = p2(x) + e2 >= q and p5(x) >= q
        //           = p5(x) >= q
        //           = x % 5^q == 0

        if (q <= 22) // 22 = floor(log_5(2^53))
        {
            za = MultipleOfPow5(u, q);
            zb = MultipleOfPow5(v, q);
            zc = MultipleOfPow5(w, q);
        }
    }
    else
    {
        // We need
        //  (a,b,c) = (u,v,w) * 2^e2 / 10^e2
        // and we need to remove at least q' = log_10(5^-e2) digits from the
        // scaled values a,b,c, i.e. we want to compute
        //  (a,b,c) = (u,v,w) * 2^e2 / 10^(e2 + q')
        //          = (u,v,w) * 2^e2 / 10^(e10),
        //          = (u,v,w) * 5^(-e10) / 2^(e10 - e2)

        const int q = FloorLog10Pow5(-e2) - (-e2 > 1); // == max(0, q' - 1)
        RYU_ASSERT(q >= 0);

        e10 = q + e2;
        RYU_ASSERT(e10 < 0);
        RYU_ASSERT(e10 - e2 >= 0);

        // Determine whether all the removed digits are 0.
        //
        // Z(x,e2,q) = (x * 5^-e2) % 10^q == 0
        //           = min(p2(x), p5(x) - e2) >= q
        //           = p2(x) >= q and p5(x) - e2 >= q
        //           = p2(x) >= q
        //           = x % 2^q == 0

        if (q <= Double::SignificandSize + 2)
        {
            za = MultipleOfPow2(u, q);
            zb = MultipleOfPow2(v, q);
            zc = MultipleOfPow2(w, q);
        }
    }

    uint64_t aq;
    uint64_t bq;
    uint64_t cq;
    MulPow5DivPow2_Double(u, v, w, -e10, e10 - e2, aq, bq, cq);

    //
    // Step 4:
    // Find the shortest decimal representation in the interval of valid representations.
    //

    cq -= !accept_upper && zc;

    // mask = 10^(number of digits removed),
    // i.e., (bq % mask) contains the actual digits removed from bq.
    // Since c < 2^62, which has 19 decimal digits, we remove at most 18 decimal digits.
    uint64_t mask = 1;

    uint64_t a = aq;
    uint64_t b = bq;
    uint64_t c = cq;

#if 0
    while (a / 10000 < c / 10000)
    {
        mask *= 10000;
        a /= 10000;
        b /= 10000;
        c /= 10000;
        e10 += 4;
    }
#else
    if (a / 10000 < c / 10000) // 4
    {
        mask = 10000;
        a /= 10000;
        b /= 10000;
        c /= 10000;
        e10 += 4;
        if (a / 10000 < c / 10000) // 8
        {
            mask = 100000000;
            a /= 10000;
            b /= 10000;
            c /= 10000;
            e10 += 4;
            if (a / 10000 < c / 10000) // 12
            {
                mask = 1000000000000;
                a /= 10000;
                b /= 10000;
                c /= 10000;
                e10 += 4;
                if (a / 10000 < c / 10000) // 16
                {
                    mask = 10000000000000000;
                    a /= 10000;
                    b /= 10000;
                    c /= 10000;
                    e10 += 4;
                }
            }
        }
    }
#endif

    if (a / 100 < c / 100)
    {
        mask *= 100;
        a /= 100;
        b /= 100;
        c /= 100;
        e10 += 2;
    }

    if (a / 10 < c / 10)
    {
        mask *= 10;
        a /= 10;
        b /= 10;
//      c /= 10;
        ++e10;
    }

    if /*likely*/ (!za && !zb)
    {
        const uint64_t br = bq - b * mask; // Digits removed from bq
        const uint64_t half = mask / 2;

        b += (a == b || br >= half);
    }
    else
    {
        // za currently determines whether the first q removed digits were all
        // 0's. Still need to check whether the digits removed in the loop above
        // are all 0's.
        const bool can_use_lower = accept_lower && za && (aq - a * mask == 0);
        if (can_use_lower)
        {
            // If the loop is executed at least once, we have a == b == c when
            // the loop terminates.
            // We only remove 0's from a, so ar and za don't change.
            RYU_ASSERT(a != 0);
            for (;;)
            {
                const uint64_t q = a / 10;
                const uint32_t r = Lo32(a) - 10 * Lo32(q); // = a % 10
                if (r != 0)
                    break;
                mask *= 10;
                a = q;
                b = q;
//              c = q;
                ++e10;
            }
        }

        const uint64_t br = bq - b * mask; // Digits removed from bq
        const uint64_t half = mask / 2;

        // A return value of b is valid if and only if a != b or za == true.
        // A return value of b + 1 is valid if and only if b + 1 <= c.
        const bool round_up = (a == b && !can_use_lower) // out of range
            || (br > half)
            || (br == half && (!zb || b % 2 != 0));

//      RYU_ASSERT(!round_up || b < c);
        b += round_up;
    }

    return {b, e10};
}

//==================================================================================================
// ToChars
//==================================================================================================

static inline char* Utoa_2Digits(char* buf, uint32_t digits)
{
    static constexpr char Digits100[200] = {
        '0','0','0','1','0','2','0','3','0','4','0','5','0','6','0','7','0','8','0','9',
        '1','0','1','1','1','2','1','3','1','4','1','5','1','6','1','7','1','8','1','9',
        '2','0','2','1','2','2','2','3','2','4','2','5','2','6','2','7','2','8','2','9',
        '3','0','3','1','3','2','3','3','3','4','3','5','3','6','3','7','3','8','3','9',
        '4','0','4','1','4','2','4','3','4','4','4','5','4','6','4','7','4','8','4','9',
        '5','0','5','1','5','2','5','3','5','4','5','5','5','6','5','7','5','8','5','9',
        '6','0','6','1','6','2','6','3','6','4','6','5','6','6','6','7','6','8','6','9',
        '7','0','7','1','7','2','7','3','7','4','7','5','7','6','7','7','7','8','7','9',
        '8','0','8','1','8','2','8','3','8','4','8','5','8','6','8','7','8','8','8','9',
        '9','0','9','1','9','2','9','3','9','4','9','5','9','6','9','7','9','8','9','9',
    };

    RYU_ASSERT(digits <= 99);
    std::memcpy(buf, &Digits100[2 * digits], 2 * sizeof(char));
    return buf + 2;
}

static inline char* Utoa_4Digits(char* buf, uint32_t digits)
{
    RYU_ASSERT(digits <= 9999);
    const uint32_t q = digits / 100;
    const uint32_t r = digits % 100;
    Utoa_2Digits(buf + 0, q);
    Utoa_2Digits(buf + 2, r);
    return buf + 4;
}

static inline char* Utoa_8Digits(char* buf, uint32_t digits)
{
    RYU_ASSERT(digits <= 99999999);
    const uint32_t q = digits / 10000;
    const uint32_t r = digits % 10000;
    Utoa_4Digits(buf + 0, q);
    Utoa_4Digits(buf + 4, r);
    return buf + 8;
}

static inline int DecimalLength(uint64_t v)
{
    RYU_ASSERT(v >= 1);
    RYU_ASSERT(v <= 99999999999999999ull);

    if (v >= 10000000000000000ull) { return 17; }
    if (v >= 1000000000000000ull) { return 16; }
    if (v >= 100000000000000ull) { return 15; }
    if (v >= 10000000000000ull) { return 14; }
    if (v >= 1000000000000ull) { return 13; }
    if (v >= 100000000000ull) { return 12; }
    if (v >= 10000000000ull) { return 11; }
    if (v >= 1000000000ull) { return 10; }
    if (v >= 100000000ull) { return 9; }
    if (v >= 10000000ull) { return 8; }
    if (v >= 1000000ull) { return 7; }
    if (v >= 100000ull) { return 6; }
    if (v >= 10000ull) { return 5; }
    if (v >= 1000ull) { return 4; }
    if (v >= 100ull) { return 3; }
    if (v >= 10ull) { return 2; }
    return 1;
}

static inline void PrintDecimalDigits(char* buf, uint64_t output, int output_length)
{
    // We prefer 32-bit operations, even on 64-bit platforms.
    // We have at most 17 digits, and uint32_t can store 9 digits.
    // If output doesn't fit into uint32_t, we cut off 8 digits,
    // so the rest will fit into uint32_t.
    if (static_cast<uint32_t>(output >> 32) != 0)
    {
        RYU_ASSERT(output_length > 8);
        const uint64_t q = output / 100000000;
        const uint32_t r = static_cast<uint32_t>(output % 100000000);
        output = q;
        output_length -= 8;
        Utoa_8Digits(buf + output_length, r);
    }

    RYU_ASSERT(output <= UINT32_MAX);
    uint32_t output2 = static_cast<uint32_t>(output);

    while (output2 >= 10000)
    {
        RYU_ASSERT(output_length > 4);
        const uint32_t q = output2 / 10000;
        const uint32_t r = output2 % 10000;
        output2 = q;
        output_length -= 4;
        Utoa_4Digits(buf + output_length, r);
    }

    if (output2 >= 100)
    {
        RYU_ASSERT(output_length > 2);
        const uint32_t q = output2 / 100;
        const uint32_t r = output2 % 100;
        output2 = q;
        output_length -= 2;
        Utoa_2Digits(buf + output_length, r);
    }

    if (output2 >= 10)
    {
        RYU_ASSERT(output_length == 2);
        Utoa_2Digits(buf, output2);
    }
    else
    {
        RYU_ASSERT(output_length == 1);
        buf[0] = static_cast<char>('0' + output2);
    }
}

static inline char* FormatDigits(char* buffer, uint64_t digits, int decimal_exponent, bool force_trailing_dot_zero = false)
{
    RYU_ASSERT(digits >= 1);
    RYU_ASSERT(digits <= 99999999999999999ull);
    RYU_ASSERT(decimal_exponent >= -999);
    RYU_ASSERT(decimal_exponent <=  999);

    const int num_digits = DecimalLength(digits);
    const int decimal_point = num_digits + decimal_exponent;

    // In order to successfully parse all numbers output by Dtoa using the Strtod implementation
    // below, we have to make sure to never emit more than 17 (significant) digits.
    static constexpr int MaxFixedDecimalPoint =  17;
    static constexpr int MinFixedDecimalPoint = -6;

    const bool use_fixed = MinFixedDecimalPoint <= decimal_point && decimal_point <= MaxFixedDecimalPoint;

    // Prepare the buffer.
    // Avoid calling memset/memcpy with variable arguments below...

    int decimal_digits_position;
    if (use_fixed)
    {
        if (decimal_point <= 0)
        {
            // 0.[000]digits
            // -6 <= decimal_point <= 0
            //  ==> 2 <= 2 + -decimal_point <= 8
            // Pre-filling the buffer with 8 '0's is therefore sufficient.
            std::memset(buffer, '0', 8);
            decimal_digits_position = 2 + (-decimal_point);
        }
        else if (decimal_point < num_digits)
        {
            // dig.its
            // 0 < decimal_point <= Min(17 - 1, MaxExp)
            // We need to move at most 16 bytes to the right.
            decimal_digits_position = 0;
        }
        else
        {
            // digits[000]
            // 1 <= num_digits <= 17 decimal_point <= 21.
            // Pre-filling buffer with 21 '0's is therefore sufficient.
            static_assert(MaxFixedDecimalPoint <= 24, "invalid parameter");
            std::memset(buffer, '0', 24);
            decimal_digits_position = 0;
        }
    }
    else
    {
        // dE+123 or d.igitsE+123
        // We only need to copy the first digit one position to the left.
        decimal_digits_position = 1;
    }

    PrintDecimalDigits(buffer + decimal_digits_position, digits, num_digits);

    if (use_fixed)
    {
        if (decimal_point <= 0)
        {
            // 0.[000]digits
            buffer[1] = '.';
            buffer += 2 + (-decimal_point) + num_digits;
        }
        else if (decimal_point < num_digits)
        {
            // dig.its
            // We need to move at most 16 bytes one place to the right.
            std::memmove(buffer + (decimal_point + 1), buffer + decimal_point, 16);
            buffer[decimal_point] = '.';
            buffer += num_digits + 1;
        }
        else // 0 < num_digits <= decimal_point
        {
            // digits[000]
            buffer += decimal_point;
            if (force_trailing_dot_zero)
            {
                std::memcpy(buffer, ".0", 2);
                buffer += 2;
            }
        }
    }
    else
    {
        // Copy the first digit one place to the left.
        buffer[0] = buffer[1];
        if (num_digits == 1)
        {
            // dE+123
            buffer += 1;
        }
        else
        {
            // d.igitsE+123
            buffer[1] = '.';
            buffer += 1 + num_digits;
        }

        const auto scientific_exponent = decimal_point - 1;
//      RYU_ASSERT(scientific_exponent != 0);

        std::memcpy(buffer, scientific_exponent < 0 ? "e-" : "e+", 2);
        buffer += 2;
        const uint32_t k = static_cast<uint32_t>(scientific_exponent < 0 ? -scientific_exponent : scientific_exponent);
        if (k < 10)
        {
            *buffer++ = static_cast<char>('0' + k);
        }
        else if (k < 100)
        {
            buffer = Utoa_2Digits(buffer, k);
        }
        else
        {
            const uint32_t r = k % 10;
            const uint32_t q = k / 10;
            buffer = Utoa_2Digits(buffer, q);
            *buffer++ = static_cast<char>('0' + r);
        }
    }

    return buffer;
}

static inline char* ToChars(char* buffer, double value, bool force_trailing_dot_zero = false)
{
    const Double v(value);

    if (!v.IsFinite())
    {
        if (v.IsNaN())
        {
            std::memcpy(buffer, "nan ", 4);
            return buffer + 3;
        }
        if (v.SignBit())
        {
            *buffer++ = '-';
        }
        std::memcpy(buffer, "inf ", 4);
        return buffer + 3;
    }


    if (v.SignBit())
    {
        value = v.AbsValue();
        *buffer++ = '-';
    }

    if (v.IsZero())
    {
        std::memcpy(buffer, "0.0 ", 4);
        buffer += 1 + (force_trailing_dot_zero ? 2 : 0);
        return buffer;
    }

    const auto dec = ToDecimal(value);
    return FormatDigits(buffer, dec.digits, dec.exponent, force_trailing_dot_zero);
}

//==================================================================================================
//
//==================================================================================================

char* charconv::Dtoa(char* buffer, double value)
{
    return ToChars(buffer, value);
}

//==================================================================================================
// ToBinary64
//==================================================================================================

// Maximum number of decimal digits in the significand the fast ToBinary method can handle.
// Inputs with more significant digits must be processed using another algorithm.
static constexpr int ToBinaryMaxDecimalDigits = 17;

// Any input <= 10^MinDecimalExponent is interpreted as 0.
// Any input >  10^MaxDecimalExponent is interpreted as +Infinity.
static constexpr int MinDecimalExponent = -324; // denorm_min = 4.9406564584124654e-324 >=  1 * 10^-324
static constexpr int MaxDecimalExponent =  309; //        max = 1.7976931348623158e+308 <= 10 * 10^+308

static inline int FloorLog2(uint64_t x)
{
    RYU_ASSERT(x != 0);

#if defined(__GNUC__) || defined(__clang__)
    return 63 - __builtin_clzll(x);
#elif defined(_MSC_VER) && defined(_M_X64)
    unsigned long index;
    _BitScanReverse64(&index, x);
    return static_cast<int>(index);
#else
    int l2 = 0;
    for (;;)
    {
        x >>= 1;
        if (x == 0)
            break;
        ++l2;
    }
    return l2;
#endif
}

static inline int FloorLog2Pow10(int e)
{
    RYU_ASSERT(e >= -1233);
    RYU_ASSERT(e <= 1233);
    return FloorDivPow2(e * 1741647, 19);
}

static inline int ExtractBit(uint64_t x, int n)
{
    RYU_ASSERT(n >= 0);
    RYU_ASSERT(n <= 63);
    return (x & (uint64_t{1} << n)) != 0;
}

static inline double ToBinary64(uint64_t m10, int m10_digits, int e10)
{
    static constexpr int MantissaBits = Double::SignificandSize - 1;
    static constexpr int ExponentBias = Double::ExponentBias - (Double::SignificandSize - 1);

    RYU_ASSERT(m10 > 0);
    RYU_ASSERT(m10_digits == DecimalLength(m10));
    RYU_ASSERT(m10_digits <= ToBinaryMaxDecimalDigits);
    RYU_ASSERT(e10 >  MinDecimalExponent - m10_digits);
    RYU_ASSERT(e10 <= MaxDecimalExponent - m10_digits);
    static_cast<void>(m10_digits);

#if defined(_M_X64) || defined(__x86_64__)
    // If the significand fits into a double (m10 <= 2^53) and the exponent 10^e10 (or 10^-e10)
    // fits into a double too (-22 <= e10 <= 22), we can compute the result simply by multiplying,
    // resp. dividing, the two numbers.
    // This is possible since IEEE guarantees that the result is correctly rounded in this case.
    //
    // If there are less than 15 decimal digits in the significand, we can extend the upper bound
    // for e10 slightly to 22 + (15 - m10_digits) by moving some zeros from the exponent into the
    // significand before the multiplication. E.g., we can write
    //      123 * 10^25 = (123 * 10^3) * 10^22,
    // since the product inside the parantheses will be < 2^53 and therefore will be exactly
    // representable as double.

    static constexpr double ExactPowersOfTen[23] = {
        1e+00,
        1e+01,
        1e+02,
        1e+03,
        1e+04,
        1e+05,
        1e+06,
        1e+07,
        1e+08,
        1e+09,
        1e+10,
        1e+11,
        1e+12,
        1e+13,
        1e+14,
        1e+15, // 10^15 < 9007199254740992 = 2^53
        1e+16, // 10^16 = 5000000000000000 * 2^1  = (10^15 * 5^1 ) * 2^1
        1e+17, // 10^17 = 6250000000000000 * 2^4  = (10^13 * 5^4 ) * 2^4
        1e+18, // 10^18 = 7812500000000000 * 2^7  = (10^11 * 5^7 ) * 2^7
        1e+19, // 10^19 = 4882812500000000 * 2^11 = (10^8  * 5^11) * 2^11
        1e+20, // 10^20 = 6103515625000000 * 2^14 = (10^6  * 5^14) * 2^14
        1e+21, // 10^21 = 7629394531250000 * 2^17 = (10^4  * 5^17) * 2^17
        1e+22, // 10^22 = 4768371582031250 * 2^21 = (10^1  * 5^21) * 2^21
    };

#if 0
    if (m10 <= (uint64_t{1} << 53) && -22 <= e10 && e10 <= 22 + Max(0, 15 - m10_digits))
    {
        double flt = static_cast<double>(static_cast<int64_t>(m10));
        if (e10 < 0)
        {
            flt /= ExactPowersOfTen[static_cast<uint32_t>(-e10)];
        }
        else if (e10 <= 22)
        {
            flt *= ExactPowersOfTen[static_cast<uint32_t>( e10)];
        }
        else
        {
            RYU_ASSERT(m10_digits < 15);

            const int zeros = 15 - m10_digits;
            flt *= ExactPowersOfTen[static_cast<uint32_t>(zeros)];
            flt *= ExactPowersOfTen[static_cast<uint32_t>(e10 - zeros)];
        }

        return flt;
    }
#else
    // NB:
    // num_digits is unused...

    if (m10 <= (uint64_t{1} << 53) && -22 <= e10 && e10 <= 22)
    {
        double flt = static_cast<double>(static_cast<int64_t>(m10));
        if (e10 < 0)
            flt /= ExactPowersOfTen[static_cast<uint32_t>(-e10)];
        else
            flt *= ExactPowersOfTen[static_cast<uint32_t>( e10)];

        return flt;
    }
#endif
#endif

    // Convert to binary float m2 * 2^e2, while retaining information about whether the conversion
    // was exact.

    const auto log2_m10 = FloorLog2(m10);
    RYU_ASSERT(log2_m10 >= 0);
    RYU_ASSERT(log2_m10 <= 56);

    // Let b = floor(log_2(m10))
    // Let n = floor(log_2(5^e10))
    // Then
    //  j = ( e2 - e10 ) - ( n + 1 - BitsPerPow5 )
    //    = ( ( b + e10 + n - (MantissaBits + 1) ) - e10 ) - ( n + 1 - BitsPerPow5 )
    //    = b + BitsPerPow5 - MantissaBits - 2
    //    = b + 125 - 52 - 2
    //    = b + 71
    // Since 0 <= b <= 56, we have
    //    71 <= j <= 127
    // The product along with the subsequent shift therefore has (at most)
    //  b + 125 - (125 - 54 + b) = 54
    // bits.

    const auto log2_10_e10 = FloorLog2Pow10(e10);
    const auto e2 = log2_m10 + log2_10_e10 - (MantissaBits + 1);

    const auto pow5 = ComputePow5_Double(e10);
    const auto j = log2_m10 + (BitsPerPow5_Double - MantissaBits - 2);
    const auto m2 = MulShift(m10, &pow5, j);

    const auto log2_m2 = FloorLog2(m2);
    RYU_ASSERT(log2_m2 >= 53);
    RYU_ASSERT(log2_m2 <= 54);

    bool is_exact;
    if (e10 >= 0)
    {
        // 56 = floor(log_2(10^17))
        is_exact = (e2 < e10) || (e2 - e10 < 64 && MultipleOfPow2(m10, e2 - e10));
    }
    else
    {
        // 57 = ceil(log_2(10^17))
        // 24 = floor(log_5(2^57))
        is_exact = -e10 <= 24 && MultipleOfPow5(m10, -e10);
    }

    // Compute the final IEEE exponent.
    int ieee_e2 = Max(0, log2_m2 + e2 + ExponentBias);
    if (ieee_e2 >= 2 * std::numeric_limits<double>::max_exponent - 1)
    {
        // Overflow:
        // Final IEEE exponent is larger than the maximum representable.
        return std::numeric_limits<double>::infinity();
    }

    // We need to figure out how much we need to shift m2.
    // The tricky part is that we need to take the final IEEE exponent into account, so we need to
    // reverse the bias and also special-case the value 0.
    const auto shift = (ieee_e2 == 0 ? 1 : ieee_e2) - e2 - ExponentBias - MantissaBits;
    RYU_ASSERT(shift > 0);

    // We need to round up if the exact value is more than 0.5 above the value we computed. That's
    // equivalent to checking if the last removed bit was 1 and either the value was not just
    // trailing zeros or the result would otherwise be odd.
    const auto trailing_zeros
        = is_exact && MultipleOfPow2(m2, shift - 1);
    const auto last_removed_bit
        = ExtractBit(m2, shift - 1);
    const auto round_up
        = last_removed_bit != 0 && (!trailing_zeros || ExtractBit(m2, shift) != 0);

    auto significand = (m2 >> shift) + round_up;
    RYU_ASSERT(significand <= 2 * Double::HiddenBit);

    if (significand == 2 * Double::HiddenBit)
    {
        // Due to how the IEEE represents +/-Infinity, we don't need to check for overflow here.
        significand >>= 1;
        ++ieee_e2;
    }
    else if (significand >= 1 * Double::HiddenBit && ieee_e2 == 0)
    {
        RYU_ASSERT((significand & 1) == 0);
        ++ieee_e2;
    }

    RYU_ASSERT(ieee_e2 <= 2 * std::numeric_limits<double>::max_exponent - 1);
    const auto ieee = static_cast<uint64_t>(ieee_e2) << MantissaBits | (significand & Double::SignificandMask);
    return ReinterpretBits<double>(ieee);
}

//==================================================================================================
// Strtod
//==================================================================================================

#define RYU_ASSUME_NULL_TERMINATED_INPUT() 0

using charconv::StrtodStatus;
using charconv::StrtodResult;

static inline bool IsDigit(char ch)
{
    return static_cast<unsigned>(ch - '0') <= 9u;
}

static inline int DigitValue(char ch)
{
    RYU_ASSERT(IsDigit(ch));
    return ch - '0';
}

static inline bool IsLowerASCII(char ch)
{
    return 'a' <= ch && ch <= 'z';
}

static inline bool IsUpperASCII(char ch)
{
    return 'A' <= ch && ch <= 'Z';
}

static inline char ToLowerASCII(char ch)
{
//  return IsUpperASCII(ch) ? static_cast<char>(ch - 'A' + 'a') : ch;
    return static_cast<char>(static_cast<unsigned char>(ch) | 0x20);
}

static inline bool StartsWith(const char* next, const char* last, const char* lower_case_prefix)
{
    for ( ; next != last && *lower_case_prefix != '\0'; ++next, ++lower_case_prefix)
    {
        RYU_ASSERT(IsLowerASCII(*lower_case_prefix));
        if (ToLowerASCII(*next) != *lower_case_prefix)
            return false;
    }

    return *lower_case_prefix == '\0';
}

static inline StrtodResult ParseInfinity(const char* next, const char* last)
{
    RYU_ASSERT((*next == 'i' || *next == 'I'));

    if (!StartsWith(next + 1, last, "nf"))
        return {next, StrtodStatus::invalid};

    next += 3;
    if (StartsWith(next, last, "inity"))
        next += 5;

    return {next, StrtodStatus::ok};
}

// FIXME:
// Don't ignore the nan-sequence!!!
static inline StrtodResult ParseNaN(const char* next, const char* last)
{
    RYU_ASSERT((*next == 'n' || *next == 'N'));

    if (!StartsWith(next + 1, last, "an"))
        return {next, StrtodStatus::invalid};

    next += 3;
    if (next != last && *next == '(')
    {
        const char* const first = next;
        for (const char* p = next + 1; p != last; ++p)
        {
            if (*p == ')')
                return {p + 1, StrtodStatus::ok};

            if (*p == '_' || IsDigit(*p) || IsUpperASCII(*p) || IsLowerASCII(*p))
                continue;

            return {first, StrtodStatus::invalid}; // invalid/incomplete nan-sequence
        }
    }

    return {next, StrtodStatus::ok};
}

static RYU_NEVER_INLINE StrtodResult ParseSpecial(bool is_negative, const char* next, const char* last, double& value)
{
    if (*next == 'i' || *next == 'I')
    {
        const auto res = ParseInfinity(next, last);
        if (res.status != StrtodStatus::invalid)
        {
            value = is_negative ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
        }
        return res;
    }

    if (*next == 'n' || *next == 'N')
    {
        const auto res = ParseNaN(next, last);
        if (res.status != StrtodStatus::invalid)
        {
            value = std::numeric_limits<double>::quiet_NaN();
        }
        return res;
    }

    return {next, StrtodStatus::invalid};
}

static RYU_NEVER_INLINE StrtodResult StdStrtod(const char* next, const char* last, double& value)
{
    //
    // FIXME:
    // _strtod_l( ..., C_LOCALE )
    //

#if RYU_ASSUME_NULL_TERMINATED_INPUT()
    const char* const ptr = next;
#else
    // std::strtod expects null-terminated inputs. So we need to make a copy and null-terminate the input.
    // This function is actually almost never going to be called, so that should be ok.
    const std::string inp(next, last);

    const char* const ptr = inp.c_str();
#endif
    char* end;
    const auto flt = ::strtod(ptr, &end);

    if (ptr == end)
        return {next, StrtodStatus::invalid}; // invalid argument (this shouldn't happen here...)
#if 0
    if (errno == ERANGE)
        return {next, StrtodStatus::invalid};
#endif

    // std::strtod should have consumed all of the input.
    RYU_ASSERT(last - next == end - ptr);

    value = flt;
    return {last, StrtodStatus::ok};
}

StrtodResult charconv::Strtod(const char* next, const char* last, double& value)
{
    const char* const start = next;

    if (next == last)
        return {next, StrtodStatus::invalid};

    // Decompose the input into the form significand * 10^exponent,
    // where significand has num_digits decimal digits.

    uint64_t significand = 0; // only valid iff num_digits <= 19
    int64_t  num_digits  = 0; // 64-bit to avoid overflow...
    int64_t  exponent    = 0; // 64-bit to avoid overflow...

// [-]

    const bool is_negative = (*next == '-');
    if (is_negative || *next == '+')
    {
        ++next;
        if (next == last)
            return {next, StrtodStatus::invalid};
    }

// int

    const bool has_leading_zero = (*next == '0');
    const bool has_leading_dot  = (*next == '.');

    if (has_leading_zero)
    {
        for (;;)
        {
            ++next;
            if (next == last || *next != '0')
                break;
        }
    }

    if (next != last && IsDigit(*next)) // non-0
    {
        const char* const p = next;

        significand = DigitValue(*next);
        ++next;
        while (next != last && IsDigit(*next))
        {
            significand = 10 * significand + DigitValue(*next);
            ++next;
        }

        num_digits = next - p;
    }
    else if (!has_leading_zero && !has_leading_dot)
    {
        return ParseSpecial(is_negative, next, last, value);
    }

// frac

    if (has_leading_dot || (next != last && *next == '.'))
    {
        ++next; // skip '.'
        if (next != last && IsDigit(*next))
        {
            // TODO:
            // Ignore trailing zeros...

            const char* const p = next;

            significand = 10 * significand + DigitValue(*next);
            ++next;
            while (next != last && IsDigit(*next))
            {
                significand = 10 * significand + DigitValue(*next);
                ++next;
            }

            const char* nz = p;
            if (num_digits == 0)
            {
                // Number is of the form "0.xxx...".
                // Move the leading zeros in the fractional part into the exponent.
                while (nz != next && *nz == '0')
                    ++nz;
            }

            num_digits += next - nz;
            exponent = -(next - p);
        }
        else
        {
            // No digits in the fractional part.
            // But at least one digit must appear in either the integral or the fractional part.
            if (has_leading_dot)
                return {next, StrtodStatus::invalid};
        }
    }

// exp

    // Exponents larger than this limit will be treated as +Infinity.
    // But we must still scan all the digits if this happens to be the case.
    static constexpr int MaxExp = 999999;
    static_assert(MaxExp >= 999, "invalid parameter");
    static_assert(MaxExp <= (INT_MAX - 9) / 10, "invalid parameter");

    int parsed_exponent = 0;
    if (next != last && (*next == 'e' || *next == 'E'))
    {
        // Possibly the start of an exponent...
        // We accept (and ignore!) invalid or incomplete exponents.
        // The 'next' pointer is updated if and only if a valid exponent has been found.
        const char* p = next;

        ++p; // skip 'e' or 'E'
        if (p != last)
        {
            const bool parsed_exponent_is_negative = (*p == '-');
            if (parsed_exponent_is_negative || *p == '+')
                ++p;

            if (p != last && IsDigit(*p))
            {
                next = p; // Found a valid exponent.

                parsed_exponent = DigitValue(*next);
                ++next;
                if (next != last && IsDigit(*next))
                {
                    parsed_exponent = 10 * parsed_exponent + DigitValue(*next);
                    ++next;
                }
                if (next != last && IsDigit(*next))
                {
                    parsed_exponent = 10 * parsed_exponent + DigitValue(*next);
                    ++next;
                }
                while (next != last && IsDigit(*next))
                {
                    if (parsed_exponent <= MaxExp)
                        parsed_exponent = 10 * parsed_exponent + DigitValue(*next);
                    ++next;
                }

                parsed_exponent = parsed_exponent_is_negative ? -parsed_exponent : parsed_exponent;

                // (Assume overflow does not happen here...)
                exponent += parsed_exponent;
            }
        }
    }

    RYU_ASSERT(num_digits >= 0);

    double flt;
    if (num_digits == 0)
    {
        flt = 0;
    }
    else if (parsed_exponent < -MaxExp || exponent + num_digits <= MinDecimalExponent)
    {
        // input = x * 10^-inf = 0
        // or
        // input < 10^MinDecimalExponent, which rounds to +-0.
        flt = 0;
    }
    else if (parsed_exponent > +MaxExp || exponent + num_digits > MaxDecimalExponent)
    {
        // input = x * 10^+inf = +inf
        // or
        // input >= 10^MaxDecimalExponent, which rounds to +-infinity.
        flt = std::numeric_limits<double>::infinity();
    }
    else if (num_digits <= ToBinaryMaxDecimalDigits)
    {
        RYU_ASSERT(exponent >= INT_MIN);
        RYU_ASSERT(exponent <= INT_MAX);
        flt = ToBinary64(significand, static_cast<int>(num_digits), static_cast<int>(exponent));
    }
    else
    {
        // We need to fall back to another algorithm if the input is too long.
        return StdStrtod(start, next, value);
    }

    value = is_negative ? -flt : flt;
    return {next, StrtodStatus::ok};
}
