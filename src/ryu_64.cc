// Copyright 2020 Ulf Adams
// Copyright 2020 Alexander Bolz
//
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#include "ryu_64.h"

#if defined(__has_include) && __has_include(<version>)
#include <version>
#else
#include <cassert>
#endif

#if defined(__cpp_lib_to_chars)
#define HAS_CHARCONV() 1
#else
#define HAS_CHARCONV() 0
#endif

//#undef NDEBUG
#include <cassert>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#if HAS_CHARCONV()
#include <charconv>
#else
#include <string>
#endif
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

//  static constexpr int32_t   MaxDigits10     = std::numeric_limits<value_type>::max_digits10;
    static constexpr int32_t   SignificandSize = std::numeric_limits<value_type>::digits; // = p   (includes the hidden bit)
    static constexpr int32_t   ExponentBias    = std::numeric_limits<value_type>::max_exponent - 1 + (SignificandSize - 1);
//  static constexpr int32_t   MaxExponent     = std::numeric_limits<value_type>::max_exponent - 1 - (SignificandSize - 1);
//  static constexpr int32_t   MinExponent     = std::numeric_limits<value_type>::min_exponent - 1 - (SignificandSize - 1);
    static constexpr int32_t   MaxIeeeExponent = 2 * std::numeric_limits<value_type>::max_exponent - 1;
    static constexpr bits_type HiddenBit       = bits_type{1} << (SignificandSize - 1);   // = 2^(p-1)
    static constexpr bits_type SignificandMask = HiddenBit - 1;                           // = 2^(p-1) - 1
    static constexpr bits_type ExponentMask    = bits_type{MaxIeeeExponent} << (SignificandSize - 1);
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
};
} // namespace

//==================================================================================================
//
//==================================================================================================

static inline int32_t Max(int32_t x, int32_t y)
{
    return y < x ? x : y;
}

// Returns floor(x / 2^n).
static inline int32_t FloorDivPow2(int32_t x, int32_t n)
{
#if 0
    // Technically, right-shift of negative integers is implementation defined...
    // Should easily be optimized into SAR (or equivalent) instruction.
    return x < 0 ? ~(~x >> n) : (x >> n);
#else
    return x >> n;
#endif
}

static inline int32_t FloorLog2Pow5(int32_t e)
{
    RYU_ASSERT(e >= -1764);
    RYU_ASSERT(e <=  1763);
    return FloorDivPow2(e * 1217359, 19);
}

static inline int32_t FloorLog10Pow2(int32_t e)
{
    RYU_ASSERT(e >= -2620);
    RYU_ASSERT(e <=  2620);
    return FloorDivPow2(e * 315653, 20);
}

static inline int32_t FloorLog10Pow5(int32_t e)
{
    RYU_ASSERT(e >= -2620);
    RYU_ASSERT(e <=  2620);
    return FloorDivPow2(e * 732923, 20);
}

static inline uint32_t Lo32(uint64_t x)
{
    return static_cast<uint32_t>(x & 0xFFFFFFFFu);
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

static constexpr int32_t BitsPerPow5_Double = 128;

namespace {
struct uint64x2 {
    uint64_t hi;
    uint64_t lo;
};
}

static inline uint64x2 ComputePow5_Double(int32_t k)
{
    // Let e = FloorLog2Pow5(k) + 1 - BitsPerPow5_Double
    // For k <  0, stores 5^k in the form:  ceil(2^-e / 5^-k)
    // For k >= 0, stores 5^k in the form: floor( 5^k / 2^e )
    static constexpr int32_t MinDecExp = -340;
    static constexpr int32_t MaxDecExp =  325;
    static constexpr uint64x2 Pow5[MaxDecExp - MinDecExp + 1] = {
        {0xBAAEE17FA23EBF76, 0x5D79BCF00D2DF64A}, // e =  -917, k = -340
        {0xE95A99DF8ACE6F53, 0xF4D82C2C107973DD}, // e =  -915, k = -339
        {0x91D8A02BB6C10594, 0x79071B9B8A4BE86A}, // e =  -912, k = -338
        {0xB64EC836A47146F9, 0x9748E2826CDEE285}, // e =  -910, k = -337
        {0xE3E27A444D8D98B7, 0xFD1B1B2308169B26}, // e =  -908, k = -336
        {0x8E6D8C6AB0787F72, 0xFE30F0F5E50E20F8}, // e =  -905, k = -335
        {0xB208EF855C969F4F, 0xBDBD2D335E51A936}, // e =  -903, k = -334
        {0xDE8B2B66B3BC4723, 0xAD2C788035E61383}, // e =  -901, k = -333
        {0x8B16FB203055AC76, 0x4C3BCB5021AFCC32}, // e =  -898, k = -332
        {0xADDCB9E83C6B1793, 0xDF4ABE242A1BBF3E}, // e =  -896, k = -331
        {0xD953E8624B85DD78, 0xD71D6DAD34A2AF0E}, // e =  -894, k = -330
        {0x87D4713D6F33AA6B, 0x8672648C40E5AD69}, // e =  -891, k = -329
        {0xA9C98D8CCB009506, 0x680EFDAF511F18C3}, // e =  -889, k = -328
        {0xD43BF0EFFDC0BA48, 0x0212BD1B2566DEF3}, // e =  -887, k = -327
        {0x84A57695FE98746D, 0x014BB630F7604B58}, // e =  -884, k = -326
        {0xA5CED43B7E3E9188, 0x419EA3BD35385E2E}, // e =  -882, k = -325
        {0xCF42894A5DCE35EA, 0x52064CAC828675BA}, // e =  -880, k = -324
        {0x818995CE7AA0E1B2, 0x7343EFEBD1940994}, // e =  -877, k = -323
        {0xA1EBFB4219491A1F, 0x1014EBE6C5F90BF9}, // e =  -875, k = -322
        {0xCA66FA129F9B60A6, 0xD41A26E077774EF7}, // e =  -873, k = -321
        {0xFD00B897478238D0, 0x8920B098955522B5}, // e =  -871, k = -320
        {0x9E20735E8CB16382, 0x55B46E5F5D5535B1}, // e =  -868, k = -319
        {0xC5A890362FDDBC62, 0xEB2189F734AA831E}, // e =  -866, k = -318
        {0xF712B443BBD52B7B, 0xA5E9EC7501D523E5}, // e =  -864, k = -317
        {0x9A6BB0AA55653B2D, 0x47B233C92125366F}, // e =  -861, k = -316
        {0xC1069CD4EABE89F8, 0x999EC0BB696E840B}, // e =  -859, k = -315
        {0xF148440A256E2C76, 0xC00670EA43CA250E}, // e =  -857, k = -314
        {0x96CD2A865764DBCA, 0x380406926A5E5729}, // e =  -854, k = -313
        {0xBC807527ED3E12BC, 0xC605083704F5ECF3}, // e =  -852, k = -312
        {0xEBA09271E88D976B, 0xF7864A44C633682F}, // e =  -850, k = -311
        {0x93445B8731587EA3, 0x7AB3EE6AFBE0211E}, // e =  -847, k = -310
        {0xB8157268FDAE9E4C, 0x5960EA05BAD82965}, // e =  -845, k = -309
        {0xE61ACF033D1A45DF, 0x6FB92487298E33BE}, // e =  -843, k = -308
        {0x8FD0C16206306BAB, 0xA5D3B6D479F8E057}, // e =  -840, k = -307
        {0xB3C4F1BA87BC8696, 0x8F48A4899877186D}, // e =  -838, k = -306
        {0xE0B62E2929ABA83C, 0x331ACDABFE94DE88}, // e =  -836, k = -305
        {0x8C71DCD9BA0B4925, 0x9FF0C08B7F1D0B15}, // e =  -833, k = -304
        {0xAF8E5410288E1B6F, 0x07ECF0AE5EE44DDA}, // e =  -831, k = -303
        {0xDB71E91432B1A24A, 0xC9E82CD9F69D6151}, // e =  -829, k = -302
        {0x892731AC9FAF056E, 0xBE311C083A225CD3}, // e =  -826, k = -301
        {0xAB70FE17C79AC6CA, 0x6DBD630A48AAF407}, // e =  -824, k = -300
        {0xD64D3D9DB981787D, 0x092CBBCCDAD5B109}, // e =  -822, k = -299
        {0x85F0468293F0EB4E, 0x25BBF56008C58EA6}, // e =  -819, k = -298
        {0xA76C582338ED2621, 0xAF2AF2B80AF6F24F}, // e =  -817, k = -297
        {0xD1476E2C07286FAA, 0x1AF5AF660DB4AEE2}, // e =  -815, k = -296
        {0x82CCA4DB847945CA, 0x50D98D9FC890ED4E}, // e =  -812, k = -295
        {0xA37FCE126597973C, 0xE50FF107BAB528A1}, // e =  -810, k = -294
        {0xCC5FC196FEFD7D0C, 0x1E53ED49A96272C9}, // e =  -808, k = -293
        {0xFF77B1FCBEBCDC4F, 0x25E8E89C13BB0F7B}, // e =  -806, k = -292
        {0x9FAACF3DF73609B1, 0x77B191618C54E9AD}, // e =  -803, k = -291
        {0xC795830D75038C1D, 0xD59DF5B9EF6A2418}, // e =  -801, k = -290
        {0xF97AE3D0D2446F25, 0x4B0573286B44AD1E}, // e =  -799, k = -289
        {0x9BECCE62836AC577, 0x4EE367F9430AEC33}, // e =  -796, k = -288
        {0xC2E801FB244576D5, 0x229C41F793CDA740}, // e =  -794, k = -287
        {0xF3A20279ED56D48A, 0x6B43527578C11110}, // e =  -792, k = -286
        {0x9845418C345644D6, 0x830A13896B78AAAA}, // e =  -789, k = -285
        {0xBE5691EF416BD60C, 0x23CC986BC656D554}, // e =  -787, k = -284
        {0xEDEC366B11C6CB8F, 0x2CBFBE86B7EC8AA9}, // e =  -785, k = -283
        {0x94B3A202EB1C3F39, 0x7BF7D71432F3D6AA}, // e =  -782, k = -282
        {0xB9E08A83A5E34F07, 0xDAF5CCD93FB0CC54}, // e =  -780, k = -281
        {0xE858AD248F5C22C9, 0xD1B3400F8F9CFF69}, // e =  -778, k = -280
        {0x91376C36D99995BE, 0x23100809B9C21FA2}, // e =  -775, k = -279
        {0xB58547448FFFFB2D, 0xABD40A0C2832A78B}, // e =  -773, k = -278
        {0xE2E69915B3FFF9F9, 0x16C90C8F323F516D}, // e =  -771, k = -277
        {0x8DD01FAD907FFC3B, 0xAE3DA7D97F6792E4}, // e =  -768, k = -276
        {0xB1442798F49FFB4A, 0x99CD11CFDF41779D}, // e =  -766, k = -275
        {0xDD95317F31C7FA1D, 0x40405643D711D584}, // e =  -764, k = -274
        {0x8A7D3EEF7F1CFC52, 0x482835EA666B2573}, // e =  -761, k = -273
        {0xAD1C8EAB5EE43B66, 0xDA3243650005EED0}, // e =  -759, k = -272
        {0xD863B256369D4A40, 0x90BED43E40076A83}, // e =  -757, k = -271
        {0x873E4F75E2224E68, 0x5A7744A6E804A292}, // e =  -754, k = -270
        {0xA90DE3535AAAE202, 0x711515D0A205CB37}, // e =  -752, k = -269
        {0xD3515C2831559A83, 0x0D5A5B44CA873E04}, // e =  -750, k = -268
        {0x8412D9991ED58091, 0xE858790AFE9486C3}, // e =  -747, k = -267
        {0xA5178FFF668AE0B6, 0x626E974DBE39A873}, // e =  -745, k = -266
        {0xCE5D73FF402D98E3, 0xFB0A3D212DC81290}, // e =  -743, k = -265
        {0x80FA687F881C7F8E, 0x7CE66634BC9D0B9A}, // e =  -740, k = -264
        {0xA139029F6A239F72, 0x1C1FFFC1EBC44E81}, // e =  -738, k = -263
        {0xC987434744AC874E, 0xA327FFB266B56221}, // e =  -736, k = -262
        {0xFBE9141915D7A922, 0x4BF1FF9F0062BAA9}, // e =  -734, k = -261
        {0x9D71AC8FADA6C9B5, 0x6F773FC3603DB4AA}, // e =  -731, k = -260
        {0xC4CE17B399107C22, 0xCB550FB4384D21D4}, // e =  -729, k = -259
        {0xF6019DA07F549B2B, 0x7E2A53A146606A49}, // e =  -727, k = -258
        {0x99C102844F94E0FB, 0x2EDA7444CBFC426E}, // e =  -724, k = -257
        {0xC0314325637A1939, 0xFA911155FEFB5309}, // e =  -722, k = -256
        {0xF03D93EEBC589F88, 0x793555AB7EBA27CB}, // e =  -720, k = -255
        {0x96267C7535B763B5, 0x4BC1558B2F3458DF}, // e =  -717, k = -254
        {0xBBB01B9283253CA2, 0x9EB1AAEDFB016F17}, // e =  -715, k = -253
        {0xEA9C227723EE8BCB, 0x465E15A979C1CADD}, // e =  -713, k = -252
        {0x92A1958A7675175F, 0x0BFACD89EC191ECA}, // e =  -710, k = -251
        {0xB749FAED14125D36, 0xCEF980EC671F667C}, // e =  -708, k = -250
        {0xE51C79A85916F484, 0x82B7E12780E7401B}, // e =  -706, k = -249
        {0x8F31CC0937AE58D2, 0xD1B2ECB8B0908811}, // e =  -703, k = -248
        {0xB2FE3F0B8599EF07, 0x861FA7E6DCB4AA16}, // e =  -701, k = -247
        {0xDFBDCECE67006AC9, 0x67A791E093E1D49B}, // e =  -699, k = -246
        {0x8BD6A141006042BD, 0xE0C8BB2C5C6D24E1}, // e =  -696, k = -245
        {0xAECC49914078536D, 0x58FAE9F773886E19}, // e =  -694, k = -244
        {0xDA7F5BF590966848, 0xAF39A475506A899F}, // e =  -692, k = -243
        {0x888F99797A5E012D, 0x6D8406C952429604}, // e =  -689, k = -242
        {0xAAB37FD7D8F58178, 0xC8E5087BA6D33B84}, // e =  -687, k = -241
        {0xD5605FCDCF32E1D6, 0xFB1E4A9A90880A65}, // e =  -685, k = -240
        {0x855C3BE0A17FCD26, 0x5CF2EEA09A550680}, // e =  -682, k = -239
        {0xA6B34AD8C9DFC06F, 0xF42FAA48C0EA481F}, // e =  -680, k = -238
        {0xD0601D8EFC57B08B, 0xF13B94DAF124DA27}, // e =  -678, k = -237
        {0x823C12795DB6CE57, 0x76C53D08D6B70859}, // e =  -675, k = -236
        {0xA2CB1717B52481ED, 0x54768C4B0C64CA6F}, // e =  -673, k = -235
        {0xCB7DDCDDA26DA268, 0xA9942F5DCF7DFD0A}, // e =  -671, k = -234
        {0xFE5D54150B090B02, 0xD3F93B35435D7C4D}, // e =  -669, k = -233
        {0x9EFA548D26E5A6E1, 0xC47BC5014A1A6DB0}, // e =  -666, k = -232
        {0xC6B8E9B0709F109A, 0x359AB6419CA1091C}, // e =  -664, k = -231
        {0xF867241C8CC6D4C0, 0xC30163D203C94B63}, // e =  -662, k = -230
        {0x9B407691D7FC44F8, 0x79E0DE63425DCF1E}, // e =  -659, k = -229
        {0xC21094364DFB5636, 0x985915FC12F542E5}, // e =  -657, k = -228
        {0xF294B943E17A2BC4, 0x3E6F5B7B17B2939E}, // e =  -655, k = -227
        {0x979CF3CA6CEC5B5A, 0xA705992CEECF9C43}, // e =  -652, k = -226
        {0xBD8430BD08277231, 0x50C6FF782A838354}, // e =  -650, k = -225
        {0xECE53CEC4A314EBD, 0xA4F8BF5635246429}, // e =  -648, k = -224
        {0x940F4613AE5ED136, 0x871B7795E136BE9A}, // e =  -645, k = -223
        {0xB913179899F68584, 0x28E2557B59846E40}, // e =  -643, k = -222
        {0xE757DD7EC07426E5, 0x331AEADA2FE589D0}, // e =  -641, k = -221
        {0x9096EA6F3848984F, 0x3FF0D2C85DEF7622}, // e =  -638, k = -220
        {0xB4BCA50B065ABE63, 0x0FED077A756B53AA}, // e =  -636, k = -219
        {0xE1EBCE4DC7F16DFB, 0xD3E8495912C62895}, // e =  -634, k = -218
        {0x8D3360F09CF6E4BD, 0x64712DD7ABBBD95D}, // e =  -631, k = -217
        {0xB080392CC4349DEC, 0xBD8D794D96AACFB4}, // e =  -629, k = -216
        {0xDCA04777F541C567, 0xECF0D7A0FC5583A1}, // e =  -627, k = -215
        {0x89E42CAAF9491B60, 0xF41686C49DB57245}, // e =  -624, k = -214
        {0xAC5D37D5B79B6239, 0x311C2875C522CED6}, // e =  -622, k = -213
        {0xD77485CB25823AC7, 0x7D633293366B828C}, // e =  -620, k = -212
        {0x86A8D39EF77164BC, 0xAE5DFF9C02033198}, // e =  -617, k = -211
        {0xA8530886B54DBDEB, 0xD9F57F830283FDFD}, // e =  -615, k = -210
        {0xD267CAA862A12D66, 0xD072DF63C324FD7C}, // e =  -613, k = -209
        {0x8380DEA93DA4BC60, 0x4247CB9E59F71E6E}, // e =  -610, k = -208
        {0xA46116538D0DEB78, 0x52D9BE85F074E609}, // e =  -608, k = -207
        {0xCD795BE870516656, 0x67902E276C921F8C}, // e =  -606, k = -206
        {0x806BD9714632DFF6, 0x00BA1CD8A3DB53B7}, // e =  -603, k = -205
        {0xA086CFCD97BF97F3, 0x80E8A40ECCD228A5}, // e =  -601, k = -204
        {0xC8A883C0FDAF7DF0, 0x6122CD128006B2CE}, // e =  -599, k = -203
        {0xFAD2A4B13D1B5D6C, 0x796B805720085F82}, // e =  -597, k = -202
        {0x9CC3A6EEC6311A63, 0xCBE3303674053BB1}, // e =  -594, k = -201
        {0xC3F490AA77BD60FC, 0xBEDBFC4411068A9D}, // e =  -592, k = -200
        {0xF4F1B4D515ACB93B, 0xEE92FB5515482D45}, // e =  -590, k = -199
        {0x991711052D8BF3C5, 0x751BDD152D4D1C4B}, // e =  -587, k = -198
        {0xBF5CD54678EEF0B6, 0xD262D45A78A0635E}, // e =  -585, k = -197
        {0xEF340A98172AACE4, 0x86FB897116C87C35}, // e =  -583, k = -196
        {0x9580869F0E7AAC0E, 0xD45D35E6AE3D4DA1}, // e =  -580, k = -195
        {0xBAE0A846D2195712, 0x8974836059CCA10A}, // e =  -578, k = -194
        {0xE998D258869FACD7, 0x2BD1A438703FC94C}, // e =  -576, k = -193
        {0x91FF83775423CC06, 0x7B6306A34627DDD0}, // e =  -573, k = -192
        {0xB67F6455292CBF08, 0x1A3BC84C17B1D543}, // e =  -571, k = -191
        {0xE41F3D6A7377EECA, 0x20CABA5F1D9E4A94}, // e =  -569, k = -190
        {0x8E938662882AF53E, 0x547EB47B7282EE9D}, // e =  -566, k = -189
        {0xB23867FB2A35B28D, 0xE99E619A4F23AA44}, // e =  -564, k = -188
        {0xDEC681F9F4C31F31, 0x6405FA00E2EC94D5}, // e =  -562, k = -187
        {0x8B3C113C38F9F37E, 0xDE83BC408DD3DD05}, // e =  -559, k = -186
        {0xAE0B158B4738705E, 0x9624AB50B148D446}, // e =  -557, k = -185
        {0xD98DDAEE19068C76, 0x3BADD624DD9B0958}, // e =  -555, k = -184
        {0x87F8A8D4CFA417C9, 0xE54CA5D70A80E5D7}, // e =  -552, k = -183
        {0xA9F6D30A038D1DBC, 0x5E9FCF4CCD211F4D}, // e =  -550, k = -182
        {0xD47487CC8470652B, 0x7647C32000696720}, // e =  -548, k = -181
        {0x84C8D4DFD2C63F3B, 0x29ECD9F40041E074}, // e =  -545, k = -180
        {0xA5FB0A17C777CF09, 0xF468107100525891}, // e =  -543, k = -179
        {0xCF79CC9DB955C2CC, 0x7182148D4066EEB5}, // e =  -541, k = -178
        {0x81AC1FE293D599BF, 0xC6F14CD848405531}, // e =  -538, k = -177
        {0xA21727DB38CB002F, 0xB8ADA00E5A506A7D}, // e =  -536, k = -176
        {0xCA9CF1D206FDC03B, 0xA6D90811F0E4851D}, // e =  -534, k = -175
        {0xFD442E4688BD304A, 0x908F4A166D1DA664}, // e =  -532, k = -174
        {0x9E4A9CEC15763E2E, 0x9A598E4E043287FF}, // e =  -529, k = -173
        {0xC5DD44271AD3CDBA, 0x40EFF1E1853F29FE}, // e =  -527, k = -172
        {0xF7549530E188C128, 0xD12BEE59E68EF47D}, // e =  -525, k = -171
        {0x9A94DD3E8CF578B9, 0x82BB74F8301958CF}, // e =  -522, k = -170
        {0xC13A148E3032D6E7, 0xE36A52363C1FAF02}, // e =  -520, k = -169
        {0xF18899B1BC3F8CA1, 0xDC44E6C3CB279AC2}, // e =  -518, k = -168
        {0x96F5600F15A7B7E5, 0x29AB103A5EF8C0BA}, // e =  -515, k = -167
        {0xBCB2B812DB11A5DE, 0x7415D448F6B6F0E8}, // e =  -513, k = -166
        {0xEBDF661791D60F56, 0x111B495B3464AD22}, // e =  -511, k = -165
        {0x936B9FCEBB25C995, 0xCAB10DD900BEEC35}, // e =  -508, k = -164
        {0xB84687C269EF3BFB, 0x3D5D514F40EEA743}, // e =  -506, k = -163
        {0xE65829B3046B0AFA, 0x0CB4A5A3112A5113}, // e =  -504, k = -162
        {0x8FF71A0FE2C2E6DC, 0x47F0E785EABA72AC}, // e =  -501, k = -161
        {0xB3F4E093DB73A093, 0x59ED216765690F57}, // e =  -499, k = -160
        {0xE0F218B8D25088B8, 0x306869C13EC3532D}, // e =  -497, k = -159
        {0x8C974F7383725573, 0x1E414218C73A13FC}, // e =  -494, k = -158
        {0xAFBD2350644EEACF, 0xE5D1929EF90898FB}, // e =  -492, k = -157
        {0xDBAC6C247D62A583, 0xDF45F746B74ABF3A}, // e =  -490, k = -156
        {0x894BC396CE5DA772, 0x6B8BBA8C328EB784}, // e =  -487, k = -155
        {0xAB9EB47C81F5114F, 0x066EA92F3F326565}, // e =  -485, k = -154
        {0xD686619BA27255A2, 0xC80A537B0EFEFEBE}, // e =  -483, k = -153
        {0x8613FD0145877585, 0xBD06742CE95F5F37}, // e =  -480, k = -152
        {0xA798FC4196E952E7, 0x2C48113823B73705}, // e =  -478, k = -151
        {0xD17F3B51FCA3A7A0, 0xF75A15862CA504C6}, // e =  -476, k = -150
        {0x82EF85133DE648C4, 0x9A984D73DBE722FC}, // e =  -473, k = -149
        {0xA3AB66580D5FDAF5, 0xC13E60D0D2E0EBBB}, // e =  -471, k = -148
        {0xCC963FEE10B7D1B3, 0x318DF905079926A9}, // e =  -469, k = -147
        {0xFFBBCFE994E5C61F, 0xFDF17746497F7053}, // e =  -467, k = -146
        {0x9FD561F1FD0F9BD3, 0xFEB6EA8BEDEFA634}, // e =  -464, k = -145
        {0xC7CABA6E7C5382C8, 0xFE64A52EE96B8FC1}, // e =  -462, k = -144
        {0xF9BD690A1B68637B, 0x3DFDCE7AA3C673B1}, // e =  -460, k = -143
        {0x9C1661A651213E2D, 0x06BEA10CA65C084F}, // e =  -457, k = -142
        {0xC31BFA0FE5698DB8, 0x486E494FCFF30A63}, // e =  -455, k = -141
        {0xF3E2F893DEC3F126, 0x5A89DBA3C3EFCCFB}, // e =  -453, k = -140
        {0x986DDB5C6B3A76B7, 0xF89629465A75E01D}, // e =  -450, k = -139
        {0xBE89523386091465, 0xF6BBB397F1135824}, // e =  -448, k = -138
        {0xEE2BA6C0678B597F, 0x746AA07DED582E2D}, // e =  -446, k = -137
        {0x94DB483840B717EF, 0xA8C2A44EB4571CDD}, // e =  -443, k = -136
        {0xBA121A4650E4DDEB, 0x92F34D62616CE414}, // e =  -441, k = -135
        {0xE896A0D7E51E1566, 0x77B020BAF9C81D18}, // e =  -439, k = -134
        {0x915E2486EF32CD60, 0x0ACE1474DC1D122F}, // e =  -436, k = -133
        {0xB5B5ADA8AAFF80B8, 0x0D819992132456BB}, // e =  -434, k = -132
        {0xE3231912D5BF60E6, 0x10E1FFF697ED6C6A}, // e =  -432, k = -131
        {0x8DF5EFABC5979C8F, 0xCA8D3FFA1EF463C2}, // e =  -429, k = -130
        {0xB1736B96B6FD83B3, 0xBD308FF8A6B17CB3}, // e =  -427, k = -129
        {0xDDD0467C64BCE4A0, 0xAC7CB3F6D05DDBDF}, // e =  -425, k = -128
        {0x8AA22C0DBEF60EE4, 0x6BCDF07A423AA96C}, // e =  -422, k = -127
        {0xAD4AB7112EB3929D, 0x86C16C98D2C953C7}, // e =  -420, k = -126
        {0xD89D64D57A607744, 0xE871C7BF077BA8B8}, // e =  -418, k = -125
        {0x87625F056C7C4A8B, 0x11471CD764AD4973}, // e =  -415, k = -124
        {0xA93AF6C6C79B5D2D, 0xD598E40D3DD89BD0}, // e =  -413, k = -123
        {0xD389B47879823479, 0x4AFF1D108D4EC2C4}, // e =  -411, k = -122
        {0x843610CB4BF160CB, 0xCEDF722A585139BB}, // e =  -408, k = -121
        {0xA54394FE1EEDB8FE, 0xC2974EB4EE658829}, // e =  -406, k = -120
        {0xCE947A3DA6A9273E, 0x733D226229FEEA33}, // e =  -404, k = -119
        {0x811CCC668829B887, 0x0806357D5A3F5260}, // e =  -401, k = -118
        {0xA163FF802A3426A8, 0xCA07C2DCB0CF26F8}, // e =  -399, k = -117
        {0xC9BCFF6034C13052, 0xFC89B393DD02F0B6}, // e =  -397, k = -116
        {0xFC2C3F3841F17C67, 0xBBAC2078D443ACE3}, // e =  -395, k = -115
        {0x9D9BA7832936EDC0, 0xD54B944B84AA4C0E}, // e =  -392, k = -114
        {0xC5029163F384A931, 0x0A9E795E65D4DF12}, // e =  -390, k = -113
        {0xF64335BCF065D37D, 0x4D4617B5FF4A16D6}, // e =  -388, k = -112
        {0x99EA0196163FA42E, 0x504BCED1BF8E4E46}, // e =  -385, k = -111
        {0xC06481FB9BCF8D39, 0xE45EC2862F71E1D7}, // e =  -383, k = -110
        {0xF07DA27A82C37088, 0x5D767327BB4E5A4D}, // e =  -381, k = -109
        {0x964E858C91BA2655, 0x3A6A07F8D510F870}, // e =  -378, k = -108
        {0xBBE226EFB628AFEA, 0x890489F70A55368C}, // e =  -376, k = -107
        {0xEADAB0ABA3B2DBE5, 0x2B45AC74CCEA842F}, // e =  -374, k = -106
        {0x92C8AE6B464FC96F, 0x3B0B8BC90012929E}, // e =  -371, k = -105
        {0xB77ADA0617E3BBCB, 0x09CE6EBB40173745}, // e =  -369, k = -104
        {0xE55990879DDCAABD, 0xCC420A6A101D0516}, // e =  -367, k = -103
        {0x8F57FA54C2A9EAB6, 0x9FA946824A12232E}, // e =  -364, k = -102
        {0xB32DF8E9F3546564, 0x47939822DC96ABFA}, // e =  -362, k = -101
        {0xDFF9772470297EBD, 0x59787E2B93BC56F8}, // e =  -360, k = -100
        {0x8BFBEA76C619EF36, 0x57EB4EDB3C55B65B}, // e =  -357, k =  -99
        {0xAEFAE51477A06B03, 0xEDE622920B6B23F2}, // e =  -355, k =  -98
        {0xDAB99E59958885C4, 0xE95FAB368E45ECEE}, // e =  -353, k =  -97
        {0x88B402F7FD75539B, 0x11DBCB0218EBB415}, // e =  -350, k =  -96
        {0xAAE103B5FCD2A881, 0xD652BDC29F26A11A}, // e =  -348, k =  -95
        {0xD59944A37C0752A2, 0x4BE76D3346F04960}, // e =  -346, k =  -94
        {0x857FCAE62D8493A5, 0x6F70A4400C562DDC}, // e =  -343, k =  -93
        {0xA6DFBD9FB8E5B88E, 0xCB4CCD500F6BB953}, // e =  -341, k =  -92
        {0xD097AD07A71F26B2, 0x7E2000A41346A7A8}, // e =  -339, k =  -91
        {0x825ECC24C873782F, 0x8ED400668C0C28C9}, // e =  -336, k =  -90
        {0xA2F67F2DFA90563B, 0x728900802F0F32FB}, // e =  -334, k =  -89
        {0xCBB41EF979346BCA, 0x4F2B40A03AD2FFBA}, // e =  -332, k =  -88
        {0xFEA126B7D78186BC, 0xE2F610C84987BFA9}, // e =  -330, k =  -87
        {0x9F24B832E6B0F436, 0x0DD9CA7D2DF4D7CA}, // e =  -327, k =  -86
        {0xC6EDE63FA05D3143, 0x91503D1C79720DBC}, // e =  -325, k =  -85
        {0xF8A95FCF88747D94, 0x75A44C6397CE912B}, // e =  -323, k =  -84
        {0x9B69DBE1B548CE7C, 0xC986AFBE3EE11ABB}, // e =  -320, k =  -83
        {0xC24452DA229B021B, 0xFBE85BADCE996169}, // e =  -318, k =  -82
        {0xF2D56790AB41C2A2, 0xFAE27299423FB9C4}, // e =  -316, k =  -81
        {0x97C560BA6B0919A5, 0xDCCD879FC967D41B}, // e =  -313, k =  -80
        {0xBDB6B8E905CB600F, 0x5400E987BBC1C921}, // e =  -311, k =  -79
        {0xED246723473E3813, 0x290123E9AAB23B69}, // e =  -309, k =  -78
        {0x9436C0760C86E30B, 0xF9A0B6720AAF6522}, // e =  -306, k =  -77
        {0xB94470938FA89BCE, 0xF808E40E8D5B3E6A}, // e =  -304, k =  -76
        {0xE7958CB87392C2C2, 0xB60B1D1230B20E05}, // e =  -302, k =  -75
        {0x90BD77F3483BB9B9, 0xB1C6F22B5E6F48C3}, // e =  -299, k =  -74
        {0xB4ECD5F01A4AA828, 0x1E38AEB6360B1AF4}, // e =  -297, k =  -73
        {0xE2280B6C20DD5232, 0x25C6DA63C38DE1B1}, // e =  -295, k =  -72
        {0x8D590723948A535F, 0x579C487E5A38AD0F}, // e =  -292, k =  -71
        {0xB0AF48EC79ACE837, 0x2D835A9DF0C6D852}, // e =  -290, k =  -70
        {0xDCDB1B2798182244, 0xF8E431456CF88E66}, // e =  -288, k =  -69
        {0x8A08F0F8BF0F156B, 0x1B8E9ECB641B5900}, // e =  -285, k =  -68
        {0xAC8B2D36EED2DAC5, 0xE272467E3D222F40}, // e =  -283, k =  -67
        {0xD7ADF884AA879177, 0x5B0ED81DCC6ABB10}, // e =  -281, k =  -66
        {0x86CCBB52EA94BAEA, 0x98E947129FC2B4EA}, // e =  -278, k =  -65
        {0xA87FEA27A539E9A5, 0x3F2398D747B36225}, // e =  -276, k =  -64
        {0xD29FE4B18E88640E, 0x8EEC7F0D19A03AAE}, // e =  -274, k =  -63
        {0x83A3EEEEF9153E89, 0x1953CF68300424AD}, // e =  -271, k =  -62
        {0xA48CEAAAB75A8E2B, 0x5FA8C3423C052DD8}, // e =  -269, k =  -61
        {0xCDB02555653131B6, 0x3792F412CB06794E}, // e =  -267, k =  -60
        {0x808E17555F3EBF11, 0xE2BBD88BBEE40BD1}, // e =  -264, k =  -59
        {0xA0B19D2AB70E6ED6, 0x5B6ACEAEAE9D0EC5}, // e =  -262, k =  -58
        {0xC8DE047564D20A8B, 0xF245825A5A445276}, // e =  -260, k =  -57
        {0xFB158592BE068D2E, 0xEED6E2F0F0D56713}, // e =  -258, k =  -56
        {0x9CED737BB6C4183D, 0x55464DD69685606C}, // e =  -255, k =  -55
        {0xC428D05AA4751E4C, 0xAA97E14C3C26B887}, // e =  -253, k =  -54
        {0xF53304714D9265DF, 0xD53DD99F4B3066A9}, // e =  -251, k =  -53
        {0x993FE2C6D07B7FAB, 0xE546A8038EFE402A}, // e =  -248, k =  -52
        {0xBF8FDB78849A5F96, 0xDE98520472BDD034}, // e =  -246, k =  -51
        {0xEF73D256A5C0F77C, 0x963E66858F6D4441}, // e =  -244, k =  -50
        {0x95A8637627989AAD, 0xDDE7001379A44AA9}, // e =  -241, k =  -49
        {0xBB127C53B17EC159, 0x5560C018580D5D53}, // e =  -239, k =  -48
        {0xE9D71B689DDE71AF, 0xAAB8F01E6E10B4A7}, // e =  -237, k =  -47
        {0x9226712162AB070D, 0xCAB3961304CA70E9}, // e =  -234, k =  -46
        {0xB6B00D69BB55C8D1, 0x3D607B97C5FD0D23}, // e =  -232, k =  -45
        {0xE45C10C42A2B3B05, 0x8CB89A7DB77C506B}, // e =  -230, k =  -44
        {0x8EB98A7A9A5B04E3, 0x77F3608E92ADB243}, // e =  -227, k =  -43
        {0xB267ED1940F1C61C, 0x55F038B237591ED4}, // e =  -225, k =  -42
        {0xDF01E85F912E37A3, 0x6B6C46DEC52F6689}, // e =  -223, k =  -41
        {0x8B61313BBABCE2C6, 0x2323AC4B3B3DA016}, // e =  -220, k =  -40
        {0xAE397D8AA96C1B77, 0xABEC975E0A0D081B}, // e =  -218, k =  -39
        {0xD9C7DCED53C72255, 0x96E7BD358C904A22}, // e =  -216, k =  -38
        {0x881CEA14545C7575, 0x7E50D64177DA2E55}, // e =  -213, k =  -37
        {0xAA242499697392D2, 0xDDE50BD1D5D0B9EA}, // e =  -211, k =  -36
        {0xD4AD2DBFC3D07787, 0x955E4EC64B44E865}, // e =  -209, k =  -35
        {0x84EC3C97DA624AB4, 0xBD5AF13BEF0B113F}, // e =  -206, k =  -34
        {0xA6274BBDD0FADD61, 0xECB1AD8AEACDD58F}, // e =  -204, k =  -33
        {0xCFB11EAD453994BA, 0x67DE18EDA5814AF3}, // e =  -202, k =  -32
        {0x81CEB32C4B43FCF4, 0x80EACF948770CED8}, // e =  -199, k =  -31
        {0xA2425FF75E14FC31, 0xA1258379A94D028E}, // e =  -197, k =  -30
        {0xCAD2F7F5359A3B3E, 0x096EE45813A04331}, // e =  -195, k =  -29
        {0xFD87B5F28300CA0D, 0x8BCA9D6E188853FD}, // e =  -193, k =  -28
        {0x9E74D1B791E07E48, 0x775EA264CF55347E}, // e =  -190, k =  -27
        {0xC612062576589DDA, 0x95364AFE032A819E}, // e =  -188, k =  -26
        {0xF79687AED3EEC551, 0x3A83DDBD83F52205}, // e =  -186, k =  -25
        {0x9ABE14CD44753B52, 0xC4926A9672793543}, // e =  -183, k =  -24
        {0xC16D9A0095928A27, 0x75B7053C0F178294}, // e =  -181, k =  -23
        {0xF1C90080BAF72CB1, 0x5324C68B12DD6339}, // e =  -179, k =  -22
        {0x971DA05074DA7BEE, 0xD3F6FC16EBCA5E04}, // e =  -176, k =  -21
        {0xBCE5086492111AEA, 0x88F4BB1CA6BCF585}, // e =  -174, k =  -20
        {0xEC1E4A7DB69561A5, 0x2B31E9E3D06C32E6}, // e =  -172, k =  -19
        {0x9392EE8E921D5D07, 0x3AFF322E62439FD0}, // e =  -169, k =  -18
        {0xB877AA3236A4B449, 0x09BEFEB9FAD487C3}, // e =  -167, k =  -17
        {0xE69594BEC44DE15B, 0x4C2EBE687989A9B4}, // e =  -165, k =  -16
        {0x901D7CF73AB0ACD9, 0x0F9D37014BF60A11}, // e =  -162, k =  -15
        {0xB424DC35095CD80F, 0x538484C19EF38C95}, // e =  -160, k =  -14
        {0xE12E13424BB40E13, 0x2865A5F206B06FBA}, // e =  -158, k =  -13
        {0x8CBCCC096F5088CB, 0xF93F87B7442E45D4}, // e =  -155, k =  -12
        {0xAFEBFF0BCB24AAFE, 0xF78F69A51539D749}, // e =  -153, k =  -11
        {0xDBE6FECEBDEDD5BE, 0xB573440E5A884D1C}, // e =  -151, k =  -10
        {0x89705F4136B4A597, 0x31680A88F8953031}, // e =  -148, k =   -9
        {0xABCC77118461CEFC, 0xFDC20D2B36BA7C3E}, // e =  -146, k =   -8
        {0xD6BF94D5E57A42BC, 0x3D32907604691B4D}, // e =  -144, k =   -7
        {0x8637BD05AF6C69B5, 0xA63F9A49C2C1B110}, // e =  -141, k =   -6
        {0xA7C5AC471B478423, 0x0FCF80DC33721D54}, // e =  -139, k =   -5
        {0xD1B71758E219652B, 0xD3C36113404EA4A9}, // e =  -137, k =   -4
        {0x83126E978D4FDF3B, 0x645A1CAC083126EA}, // e =  -134, k =   -3
        {0xA3D70A3D70A3D70A, 0x3D70A3D70A3D70A4}, // e =  -132, k =   -2
        {0xCCCCCCCCCCCCCCCC, 0xCCCCCCCCCCCCCCCD}, // e =  -130, k =   -1
        {0x8000000000000000, 0x0000000000000000}, // e =  -127, k =    0
        {0xA000000000000000, 0x0000000000000000}, // e =  -125, k =    1
        {0xC800000000000000, 0x0000000000000000}, // e =  -123, k =    2
        {0xFA00000000000000, 0x0000000000000000}, // e =  -121, k =    3
        {0x9C40000000000000, 0x0000000000000000}, // e =  -118, k =    4
        {0xC350000000000000, 0x0000000000000000}, // e =  -116, k =    5
        {0xF424000000000000, 0x0000000000000000}, // e =  -114, k =    6
        {0x9896800000000000, 0x0000000000000000}, // e =  -111, k =    7
        {0xBEBC200000000000, 0x0000000000000000}, // e =  -109, k =    8
        {0xEE6B280000000000, 0x0000000000000000}, // e =  -107, k =    9
        {0x9502F90000000000, 0x0000000000000000}, // e =  -104, k =   10
        {0xBA43B74000000000, 0x0000000000000000}, // e =  -102, k =   11
        {0xE8D4A51000000000, 0x0000000000000000}, // e =  -100, k =   12
        {0x9184E72A00000000, 0x0000000000000000}, // e =   -97, k =   13
        {0xB5E620F480000000, 0x0000000000000000}, // e =   -95, k =   14
        {0xE35FA931A0000000, 0x0000000000000000}, // e =   -93, k =   15
        {0x8E1BC9BF04000000, 0x0000000000000000}, // e =   -90, k =   16
        {0xB1A2BC2EC5000000, 0x0000000000000000}, // e =   -88, k =   17
        {0xDE0B6B3A76400000, 0x0000000000000000}, // e =   -86, k =   18
        {0x8AC7230489E80000, 0x0000000000000000}, // e =   -83, k =   19
        {0xAD78EBC5AC620000, 0x0000000000000000}, // e =   -81, k =   20
        {0xD8D726B7177A8000, 0x0000000000000000}, // e =   -79, k =   21
        {0x878678326EAC9000, 0x0000000000000000}, // e =   -76, k =   22
        {0xA968163F0A57B400, 0x0000000000000000}, // e =   -74, k =   23
        {0xD3C21BCECCEDA100, 0x0000000000000000}, // e =   -72, k =   24
        {0x84595161401484A0, 0x0000000000000000}, // e =   -69, k =   25
        {0xA56FA5B99019A5C8, 0x0000000000000000}, // e =   -67, k =   26
        {0xCECB8F27F4200F3A, 0x0000000000000000}, // e =   -65, k =   27
        {0x813F3978F8940984, 0x4000000000000000}, // e =   -62, k =   28
        {0xA18F07D736B90BE5, 0x5000000000000000}, // e =   -60, k =   29
        {0xC9F2C9CD04674EDE, 0xA400000000000000}, // e =   -58, k =   30
        {0xFC6F7C4045812296, 0x4D00000000000000}, // e =   -56, k =   31
        {0x9DC5ADA82B70B59D, 0xF020000000000000}, // e =   -53, k =   32
        {0xC5371912364CE305, 0x6C28000000000000}, // e =   -51, k =   33
        {0xF684DF56C3E01BC6, 0xC732000000000000}, // e =   -49, k =   34
        {0x9A130B963A6C115C, 0x3C7F400000000000}, // e =   -46, k =   35
        {0xC097CE7BC90715B3, 0x4B9F100000000000}, // e =   -44, k =   36
        {0xF0BDC21ABB48DB20, 0x1E86D40000000000}, // e =   -42, k =   37
        {0x96769950B50D88F4, 0x1314448000000000}, // e =   -39, k =   38
        {0xBC143FA4E250EB31, 0x17D955A000000000}, // e =   -37, k =   39
        {0xEB194F8E1AE525FD, 0x5DCFAB0800000000}, // e =   -35, k =   40
        {0x92EFD1B8D0CF37BE, 0x5AA1CAE500000000}, // e =   -32, k =   41
        {0xB7ABC627050305AD, 0xF14A3D9E40000000}, // e =   -30, k =   42
        {0xE596B7B0C643C719, 0x6D9CCD05D0000000}, // e =   -28, k =   43
        {0x8F7E32CE7BEA5C6F, 0xE4820023A2000000}, // e =   -25, k =   44
        {0xB35DBF821AE4F38B, 0xDDA2802C8A800000}, // e =   -23, k =   45
        {0xE0352F62A19E306E, 0xD50B2037AD200000}, // e =   -21, k =   46
        {0x8C213D9DA502DE45, 0x4526F422CC340000}, // e =   -18, k =   47
        {0xAF298D050E4395D6, 0x9670B12B7F410000}, // e =   -16, k =   48
        {0xDAF3F04651D47B4C, 0x3C0CDD765F114000}, // e =   -14, k =   49
        {0x88D8762BF324CD0F, 0xA5880A69FB6AC800}, // e =   -11, k =   50
        {0xAB0E93B6EFEE0053, 0x8EEA0D047A457A00}, // e =    -9, k =   51
        {0xD5D238A4ABE98068, 0x72A4904598D6D880}, // e =    -7, k =   52
        {0x85A36366EB71F041, 0x47A6DA2B7F864750}, // e =    -4, k =   53
        {0xA70C3C40A64E6C51, 0x999090B65F67D924}, // e =    -2, k =   54
        {0xD0CF4B50CFE20765, 0xFFF4B4E3F741CF6D}, // e =     0, k =   55
        {0x82818F1281ED449F, 0xBFF8F10E7A8921A4}, // e =     3, k =   56
        {0xA321F2D7226895C7, 0xAFF72D52192B6A0D}, // e =     5, k =   57
        {0xCBEA6F8CEB02BB39, 0x9BF4F8A69F764490}, // e =     7, k =   58
        {0xFEE50B7025C36A08, 0x02F236D04753D5B4}, // e =     9, k =   59
        {0x9F4F2726179A2245, 0x01D762422C946590}, // e =    12, k =   60
        {0xC722F0EF9D80AAD6, 0x424D3AD2B7B97EF5}, // e =    14, k =   61
        {0xF8EBAD2B84E0D58B, 0xD2E0898765A7DEB2}, // e =    16, k =   62
        {0x9B934C3B330C8577, 0x63CC55F49F88EB2F}, // e =    19, k =   63
        {0xC2781F49FFCFA6D5, 0x3CBF6B71C76B25FB}, // e =    21, k =   64
        {0xF316271C7FC3908A, 0x8BEF464E3945EF7A}, // e =    23, k =   65
        {0x97EDD871CFDA3A56, 0x97758BF0E3CBB5AC}, // e =    26, k =   66
        {0xBDE94E8E43D0C8EC, 0x3D52EEED1CBEA317}, // e =    28, k =   67
        {0xED63A231D4C4FB27, 0x4CA7AAA863EE4BDD}, // e =    30, k =   68
        {0x945E455F24FB1CF8, 0x8FE8CAA93E74EF6A}, // e =    33, k =   69
        {0xB975D6B6EE39E436, 0xB3E2FD538E122B44}, // e =    35, k =   70
        {0xE7D34C64A9C85D44, 0x60DBBCA87196B616}, // e =    37, k =   71
        {0x90E40FBEEA1D3A4A, 0xBC8955E946FE31CD}, // e =    40, k =   72
        {0xB51D13AEA4A488DD, 0x6BABAB6398BDBE41}, // e =    42, k =   73
        {0xE264589A4DCDAB14, 0xC696963C7EED2DD1}, // e =    44, k =   74
        {0x8D7EB76070A08AEC, 0xFC1E1DE5CF543CA2}, // e =    47, k =   75
        {0xB0DE65388CC8ADA8, 0x3B25A55F43294BCB}, // e =    49, k =   76
        {0xDD15FE86AFFAD912, 0x49EF0EB713F39EBE}, // e =    51, k =   77
        {0x8A2DBF142DFCC7AB, 0x6E3569326C784337}, // e =    54, k =   78
        {0xACB92ED9397BF996, 0x49C2C37F07965404}, // e =    56, k =   79
        {0xD7E77A8F87DAF7FB, 0xDC33745EC97BE906}, // e =    58, k =   80
        {0x86F0AC99B4E8DAFD, 0x69A028BB3DED71A3}, // e =    61, k =   81
        {0xA8ACD7C0222311BC, 0xC40832EA0D68CE0C}, // e =    63, k =   82
        {0xD2D80DB02AABD62B, 0xF50A3FA490C30190}, // e =    65, k =   83
        {0x83C7088E1AAB65DB, 0x792667C6DA79E0FA}, // e =    68, k =   84
        {0xA4B8CAB1A1563F52, 0x577001B891185938}, // e =    70, k =   85
        {0xCDE6FD5E09ABCF26, 0xED4C0226B55E6F86}, // e =    72, k =   86
        {0x80B05E5AC60B6178, 0x544F8158315B05B4}, // e =    75, k =   87
        {0xA0DC75F1778E39D6, 0x696361AE3DB1C721}, // e =    77, k =   88
        {0xC913936DD571C84C, 0x03BC3A19CD1E38E9}, // e =    79, k =   89
        {0xFB5878494ACE3A5F, 0x04AB48A04065C723}, // e =    81, k =   90
        {0x9D174B2DCEC0E47B, 0x62EB0D64283F9C76}, // e =    84, k =   91
        {0xC45D1DF942711D9A, 0x3BA5D0BD324F8394}, // e =    86, k =   92
        {0xF5746577930D6500, 0xCA8F44EC7EE36479}, // e =    88, k =   93
        {0x9968BF6ABBE85F20, 0x7E998B13CF4E1ECB}, // e =    91, k =   94
        {0xBFC2EF456AE276E8, 0x9E3FEDD8C321A67E}, // e =    93, k =   95
        {0xEFB3AB16C59B14A2, 0xC5CFE94EF3EA101E}, // e =    95, k =   96
        {0x95D04AEE3B80ECE5, 0xBBA1F1D158724A12}, // e =    98, k =   97
        {0xBB445DA9CA61281F, 0x2A8A6E45AE8EDC97}, // e =   100, k =   98
        {0xEA1575143CF97226, 0xF52D09D71A3293BD}, // e =   102, k =   99
        {0x924D692CA61BE758, 0x593C2626705F9C56}, // e =   105, k =  100
        {0xB6E0C377CFA2E12E, 0x6F8B2FB00C77836C}, // e =   107, k =  101
        {0xE498F455C38B997A, 0x0B6DFB9C0F956447}, // e =   109, k =  102
        {0x8EDF98B59A373FEC, 0x4724BD4189BD5EAC}, // e =   112, k =  103
        {0xB2977EE300C50FE7, 0x58EDEC91EC2CB657}, // e =   114, k =  104
        {0xDF3D5E9BC0F653E1, 0x2F2967B66737E3ED}, // e =   116, k =  105
        {0x8B865B215899F46C, 0xBD79E0D20082EE74}, // e =   119, k =  106
        {0xAE67F1E9AEC07187, 0xECD8590680A3AA11}, // e =   121, k =  107
        {0xDA01EE641A708DE9, 0xE80E6F4820CC9495}, // e =   123, k =  108
        {0x884134FE908658B2, 0x3109058D147FDCDD}, // e =   126, k =  109
        {0xAA51823E34A7EEDE, 0xBD4B46F0599FD415}, // e =   128, k =  110
        {0xD4E5E2CDC1D1EA96, 0x6C9E18AC7007C91A}, // e =   130, k =  111
        {0x850FADC09923329E, 0x03E2CF6BC604DDB0}, // e =   133, k =  112
        {0xA6539930BF6BFF45, 0x84DB8346B786151C}, // e =   135, k =  113
        {0xCFE87F7CEF46FF16, 0xE612641865679A63}, // e =   137, k =  114
        {0x81F14FAE158C5F6E, 0x4FCB7E8F3F60C07E}, // e =   140, k =  115
        {0xA26DA3999AEF7749, 0xE3BE5E330F38F09D}, // e =   142, k =  116
        {0xCB090C8001AB551C, 0x5CADF5BFD3072CC5}, // e =   144, k =  117
        {0xFDCB4FA002162A63, 0x73D9732FC7C8F7F6}, // e =   146, k =  118
        {0x9E9F11C4014DDA7E, 0x2867E7FDDCDD9AFA}, // e =   149, k =  119
        {0xC646D63501A1511D, 0xB281E1FD541501B8}, // e =   151, k =  120
        {0xF7D88BC24209A565, 0x1F225A7CA91A4226}, // e =   153, k =  121
        {0x9AE757596946075F, 0x3375788DE9B06958}, // e =   156, k =  122
        {0xC1A12D2FC3978937, 0x0052D6B1641C83AE}, // e =   158, k =  123
        {0xF209787BB47D6B84, 0xC0678C5DBD23A49A}, // e =   160, k =  124
        {0x9745EB4D50CE6332, 0xF840B7BA963646E0}, // e =   163, k =  125
        {0xBD176620A501FBFF, 0xB650E5A93BC3D898}, // e =   165, k =  126
        {0xEC5D3FA8CE427AFF, 0xA3E51F138AB4CEBE}, // e =   167, k =  127
        {0x93BA47C980E98CDF, 0xC66F336C36B10137}, // e =   170, k =  128
        {0xB8A8D9BBE123F017, 0xB80B0047445D4184}, // e =   172, k =  129
        {0xE6D3102AD96CEC1D, 0xA60DC059157491E5}, // e =   174, k =  130
        {0x9043EA1AC7E41392, 0x87C89837AD68DB2F}, // e =   177, k =  131
        {0xB454E4A179DD1877, 0x29BABE4598C311FB}, // e =   179, k =  132
        {0xE16A1DC9D8545E94, 0xF4296DD6FEF3D67A}, // e =   181, k =  133
        {0x8CE2529E2734BB1D, 0x1899E4A65F58660C}, // e =   184, k =  134
        {0xB01AE745B101E9E4, 0x5EC05DCFF72E7F8F}, // e =   186, k =  135
        {0xDC21A1171D42645D, 0x76707543F4FA1F73}, // e =   188, k =  136
        {0x899504AE72497EBA, 0x6A06494A791C53A8}, // e =   191, k =  137
        {0xABFA45DA0EDBDE69, 0x0487DB9D17636892}, // e =   193, k =  138
        {0xD6F8D7509292D603, 0x45A9D2845D3C42B6}, // e =   195, k =  139
        {0x865B86925B9BC5C2, 0x0B8A2392BA45A9B2}, // e =   198, k =  140
        {0xA7F26836F282B732, 0x8E6CAC7768D7141E}, // e =   200, k =  141
        {0xD1EF0244AF2364FF, 0x3207D795430CD926}, // e =   202, k =  142
        {0x8335616AED761F1F, 0x7F44E6BD49E807B8}, // e =   205, k =  143
        {0xA402B9C5A8D3A6E7, 0x5F16206C9C6209A6}, // e =   207, k =  144
        {0xCD036837130890A1, 0x36DBA887C37A8C0F}, // e =   209, k =  145
        {0x802221226BE55A64, 0xC2494954DA2C9789}, // e =   212, k =  146
        {0xA02AA96B06DEB0FD, 0xF2DB9BAA10B7BD6C}, // e =   214, k =  147
        {0xC83553C5C8965D3D, 0x6F92829494E5ACC7}, // e =   216, k =  148
        {0xFA42A8B73ABBF48C, 0xCB772339BA1F17F9}, // e =   218, k =  149
        {0x9C69A97284B578D7, 0xFF2A760414536EFB}, // e =   221, k =  150
        {0xC38413CF25E2D70D, 0xFEF5138519684ABA}, // e =   223, k =  151
        {0xF46518C2EF5B8CD1, 0x7EB258665FC25D69}, // e =   225, k =  152
        {0x98BF2F79D5993802, 0xEF2F773FFBD97A61}, // e =   228, k =  153
        {0xBEEEFB584AFF8603, 0xAAFB550FFACFD8FA}, // e =   230, k =  154
        {0xEEAABA2E5DBF6784, 0x95BA2A53F983CF38}, // e =   232, k =  155
        {0x952AB45CFA97A0B2, 0xDD945A747BF26183}, // e =   235, k =  156
        {0xBA756174393D88DF, 0x94F971119AEEF9E4}, // e =   237, k =  157
        {0xE912B9D1478CEB17, 0x7A37CD5601AAB85D}, // e =   239, k =  158
        {0x91ABB422CCB812EE, 0xAC62E055C10AB33A}, // e =   242, k =  159
        {0xB616A12B7FE617AA, 0x577B986B314D6009}, // e =   244, k =  160
        {0xE39C49765FDF9D94, 0xED5A7E85FDA0B80B}, // e =   246, k =  161
        {0x8E41ADE9FBEBC27D, 0x14588F13BE847307}, // e =   249, k =  162
        {0xB1D219647AE6B31C, 0x596EB2D8AE258FC8}, // e =   251, k =  163
        {0xDE469FBD99A05FE3, 0x6FCA5F8ED9AEF3BB}, // e =   253, k =  164
        {0x8AEC23D680043BEE, 0x25DE7BB9480D5854}, // e =   256, k =  165
        {0xADA72CCC20054AE9, 0xAF561AA79A10AE6A}, // e =   258, k =  166
        {0xD910F7FF28069DA4, 0x1B2BA1518094DA04}, // e =   260, k =  167
        {0x87AA9AFF79042286, 0x90FB44D2F05D0842}, // e =   263, k =  168
        {0xA99541BF57452B28, 0x353A1607AC744A53}, // e =   265, k =  169
        {0xD3FA922F2D1675F2, 0x42889B8997915CE8}, // e =   267, k =  170
        {0x847C9B5D7C2E09B7, 0x69956135FEBADA11}, // e =   270, k =  171
        {0xA59BC234DB398C25, 0x43FAB9837E699095}, // e =   272, k =  172
        {0xCF02B2C21207EF2E, 0x94F967E45E03F4BB}, // e =   274, k =  173
        {0x8161AFB94B44F57D, 0x1D1BE0EEBAC278F5}, // e =   277, k =  174
        {0xA1BA1BA79E1632DC, 0x6462D92A69731732}, // e =   279, k =  175
        {0xCA28A291859BBF93, 0x7D7B8F7503CFDCFE}, // e =   281, k =  176
        {0xFCB2CB35E702AF78, 0x5CDA735244C3D43E}, // e =   283, k =  177
        {0x9DEFBF01B061ADAB, 0x3A0888136AFA64A7}, // e =   286, k =  178
        {0xC56BAEC21C7A1916, 0x088AAA1845B8FDD0}, // e =   288, k =  179
        {0xF6C69A72A3989F5B, 0x8AAD549E57273D45}, // e =   290, k =  180
        {0x9A3C2087A63F6399, 0x36AC54E2F678864B}, // e =   293, k =  181
        {0xC0CB28A98FCF3C7F, 0x84576A1BB416A7DD}, // e =   295, k =  182
        {0xF0FDF2D3F3C30B9F, 0x656D44A2A11C51D5}, // e =   297, k =  183
        {0x969EB7C47859E743, 0x9F644AE5A4B1B325}, // e =   300, k =  184
        {0xBC4665B596706114, 0x873D5D9F0DDE1FEE}, // e =   302, k =  185
        {0xEB57FF22FC0C7959, 0xA90CB506D155A7EA}, // e =   304, k =  186
        {0x9316FF75DD87CBD8, 0x09A7F12442D588F2}, // e =   307, k =  187
        {0xB7DCBF5354E9BECE, 0x0C11ED6D538AEB2F}, // e =   309, k =  188
        {0xE5D3EF282A242E81, 0x8F1668C8A86DA5FA}, // e =   311, k =  189
        {0x8FA475791A569D10, 0xF96E017D694487BC}, // e =   314, k =  190
        {0xB38D92D760EC4455, 0x37C981DCC395A9AC}, // e =   316, k =  191
        {0xE070F78D3927556A, 0x85BBE253F47B1417}, // e =   318, k =  192
        {0x8C469AB843B89562, 0x93956D7478CCEC8E}, // e =   321, k =  193
        {0xAF58416654A6BABB, 0x387AC8D1970027B2}, // e =   323, k =  194
        {0xDB2E51BFE9D0696A, 0x06997B05FCC0319E}, // e =   325, k =  195
        {0x88FCF317F22241E2, 0x441FECE3BDF81F03}, // e =   328, k =  196
        {0xAB3C2FDDEEAAD25A, 0xD527E81CAD7626C3}, // e =   330, k =  197
        {0xD60B3BD56A5586F1, 0x8A71E223D8D3B074}, // e =   332, k =  198
        {0x85C7056562757456, 0xF6872D5667844E49}, // e =   335, k =  199
        {0xA738C6BEBB12D16C, 0xB428F8AC016561DB}, // e =   337, k =  200
        {0xD106F86E69D785C7, 0xE13336D701BEBA52}, // e =   339, k =  201
        {0x82A45B450226B39C, 0xECC0024661173473}, // e =   342, k =  202
        {0xA34D721642B06084, 0x27F002D7F95D0190}, // e =   344, k =  203
        {0xCC20CE9BD35C78A5, 0x31EC038DF7B441F4}, // e =   346, k =  204
        {0xFF290242C83396CE, 0x7E67047175A15271}, // e =   348, k =  205
        {0x9F79A169BD203E41, 0x0F0062C6E984D386}, // e =   351, k =  206
        {0xC75809C42C684DD1, 0x52C07B78A3E60868}, // e =   353, k =  207
        {0xF92E0C3537826145, 0xA7709A56CCDF8A82}, // e =   355, k =  208
        {0x9BBCC7A142B17CCB, 0x88A66076400BB691}, // e =   358, k =  209
        {0xC2ABF989935DDBFE, 0x6ACFF893D00EA435}, // e =   360, k =  210
        {0xF356F7EBF83552FE, 0x0583F6B8C4124D43}, // e =   362, k =  211
        {0x98165AF37B2153DE, 0xC3727A337A8B704A}, // e =   365, k =  212
        {0xBE1BF1B059E9A8D6, 0x744F18C0592E4C5C}, // e =   367, k =  213
        {0xEDA2EE1C7064130C, 0x1162DEF06F79DF73}, // e =   369, k =  214
        {0x9485D4D1C63E8BE7, 0x8ADDCB5645AC2BA8}, // e =   372, k =  215
        {0xB9A74A0637CE2EE1, 0x6D953E2BD7173692}, // e =   374, k =  216
        {0xE8111C87C5C1BA99, 0xC8FA8DB6CCDD0437}, // e =   376, k =  217
        {0x910AB1D4DB9914A0, 0x1D9C9892400A22A2}, // e =   379, k =  218
        {0xB54D5E4A127F59C8, 0x2503BEB6D00CAB4B}, // e =   381, k =  219
        {0xE2A0B5DC971F303A, 0x2E44AE64840FD61D}, // e =   383, k =  220
        {0x8DA471A9DE737E24, 0x5CEAECFED289E5D2}, // e =   386, k =  221
        {0xB10D8E1456105DAD, 0x7425A83E872C5F47}, // e =   388, k =  222
        {0xDD50F1996B947518, 0xD12F124E28F77719}, // e =   390, k =  223
        {0x8A5296FFE33CC92F, 0x82BD6B70D99AAA6F}, // e =   393, k =  224
        {0xACE73CBFDC0BFB7B, 0x636CC64D1001550B}, // e =   395, k =  225
        {0xD8210BEFD30EFA5A, 0x3C47F7E05401AA4E}, // e =   397, k =  226
        {0x8714A775E3E95C78, 0x65ACFAEC34810A71}, // e =   400, k =  227
        {0xA8D9D1535CE3B396, 0x7F1839A741A14D0D}, // e =   402, k =  228
        {0xD31045A8341CA07C, 0x1EDE48111209A050}, // e =   404, k =  229
        {0x83EA2B892091E44D, 0x934AED0AAB460432}, // e =   407, k =  230
        {0xA4E4B66B68B65D60, 0xF81DA84D5617853F}, // e =   409, k =  231
        {0xCE1DE40642E3F4B9, 0x36251260AB9D668E}, // e =   411, k =  232
        {0x80D2AE83E9CE78F3, 0xC1D72B7C6B426019}, // e =   414, k =  233
        {0xA1075A24E4421730, 0xB24CF65B8612F81F}, // e =   416, k =  234
        {0xC94930AE1D529CFC, 0xDEE033F26797B627}, // e =   418, k =  235
        {0xFB9B7CD9A4A7443C, 0x169840EF017DA3B1}, // e =   420, k =  236
        {0x9D412E0806E88AA5, 0x8E1F289560EE864E}, // e =   423, k =  237
        {0xC491798A08A2AD4E, 0xF1A6F2BAB92A27E2}, // e =   425, k =  238
        {0xF5B5D7EC8ACB58A2, 0xAE10AF696774B1DB}, // e =   427, k =  239
        {0x9991A6F3D6BF1765, 0xACCA6DA1E0A8EF29}, // e =   430, k =  240
        {0xBFF610B0CC6EDD3F, 0x17FD090A58D32AF3}, // e =   432, k =  241
        {0xEFF394DCFF8A948E, 0xDDFC4B4CEF07F5B0}, // e =   434, k =  242
        {0x95F83D0A1FB69CD9, 0x4ABDAF101564F98E}, // e =   437, k =  243
        {0xBB764C4CA7A4440F, 0x9D6D1AD41ABE37F1}, // e =   439, k =  244
        {0xEA53DF5FD18D5513, 0x84C86189216DC5ED}, // e =   441, k =  245
        {0x92746B9BE2F8552C, 0x32FD3CF5B4E49BB4}, // e =   444, k =  246
        {0xB7118682DBB66A77, 0x3FBC8C33221DC2A1}, // e =   446, k =  247
        {0xE4D5E82392A40515, 0x0FABAF3FEAA5334A}, // e =   448, k =  248
        {0x8F05B1163BA6832D, 0x29CB4D87F2A7400E}, // e =   451, k =  249
        {0xB2C71D5BCA9023F8, 0x743E20E9EF511012}, // e =   453, k =  250
        {0xDF78E4B2BD342CF6, 0x914DA9246B255416}, // e =   455, k =  251
        {0x8BAB8EEFB6409C1A, 0x1AD089B6C2F7548E}, // e =   458, k =  252
        {0xAE9672ABA3D0C320, 0xA184AC2473B529B1}, // e =   460, k =  253
        {0xDA3C0F568CC4F3E8, 0xC9E5D72D90A2741E}, // e =   462, k =  254
        {0x8865899617FB1871, 0x7E2FA67C7A658892}, // e =   465, k =  255
        {0xAA7EEBFB9DF9DE8D, 0xDDBB901B98FEEAB7}, // e =   467, k =  256
        {0xD51EA6FA85785631, 0x552A74227F3EA565}, // e =   469, k =  257
        {0x8533285C936B35DE, 0xD53A88958F87275F}, // e =   472, k =  258
        {0xA67FF273B8460356, 0x8A892ABAF368F137}, // e =   474, k =  259
        {0xD01FEF10A657842C, 0x2D2B7569B0432D85}, // e =   476, k =  260
        {0x8213F56A67F6B29B, 0x9C3B29620E29FC73}, // e =   479, k =  261
        {0xA298F2C501F45F42, 0x8349F3BA91B47B8F}, // e =   481, k =  262
        {0xCB3F2F7642717713, 0x241C70A936219A73}, // e =   483, k =  263
        {0xFE0EFB53D30DD4D7, 0xED238CD383AA0110}, // e =   485, k =  264
        {0x9EC95D1463E8A506, 0xF4363804324A40AA}, // e =   488, k =  265
        {0xC67BB4597CE2CE48, 0xB143C6053EDCD0D5}, // e =   490, k =  266
        {0xF81AA16FDC1B81DA, 0xDD94B7868E94050A}, // e =   492, k =  267
        {0x9B10A4E5E9913128, 0xCA7CF2B4191C8326}, // e =   495, k =  268
        {0xC1D4CE1F63F57D72, 0xFD1C2F611F63A3F0}, // e =   497, k =  269
        {0xF24A01A73CF2DCCF, 0xBC633B39673C8CEC}, // e =   499, k =  270
        {0x976E41088617CA01, 0xD5BE0503E085D813}, // e =   502, k =  271
        {0xBD49D14AA79DBC82, 0x4B2D8644D8A74E18}, // e =   504, k =  272
        {0xEC9C459D51852BA2, 0xDDF8E7D60ED1219E}, // e =   506, k =  273
        {0x93E1AB8252F33B45, 0xCABB90E5C942B503}, // e =   509, k =  274
        {0xB8DA1662E7B00A17, 0x3D6A751F3B936243}, // e =   511, k =  275
        {0xE7109BFBA19C0C9D, 0x0CC512670A783AD4}, // e =   513, k =  276
        {0x906A617D450187E2, 0x27FB2B80668B24C5}, // e =   516, k =  277
        {0xB484F9DC9641E9DA, 0xB1F9F660802DEDF6}, // e =   518, k =  278
        {0xE1A63853BBD26451, 0x5E7873F8A0396973}, // e =   520, k =  279
        {0x8D07E33455637EB2, 0xDB0B487B6423E1E8}, // e =   523, k =  280
        {0xB049DC016ABC5E5F, 0x91CE1A9A3D2CDA62}, // e =   525, k =  281
        {0xDC5C5301C56B75F7, 0x7641A140CC7810FB}, // e =   527, k =  282
        {0x89B9B3E11B6329BA, 0xA9E904C87FCB0A9D}, // e =   530, k =  283
        {0xAC2820D9623BF429, 0x546345FA9FBDCD44}, // e =   532, k =  284
        {0xD732290FBACAF133, 0xA97C177947AD4095}, // e =   534, k =  285
        {0x867F59A9D4BED6C0, 0x49ED8EABCCCC485D}, // e =   537, k =  286
        {0xA81F301449EE8C70, 0x5C68F256BFFF5A74}, // e =   539, k =  287
        {0xD226FC195C6A2F8C, 0x73832EEC6FFF3111}, // e =   541, k =  288
        {0x83585D8FD9C25DB7, 0xC831FD53C5FF7EAB}, // e =   544, k =  289
        {0xA42E74F3D032F525, 0xBA3E7CA8B77F5E55}, // e =   546, k =  290
        {0xCD3A1230C43FB26F, 0x28CE1BD2E55F35EB}, // e =   548, k =  291
        {0x80444B5E7AA7CF85, 0x7980D163CF5B81B3}, // e =   551, k =  292
        {0xA0555E361951C366, 0xD7E105BCC332621F}, // e =   553, k =  293
        {0xC86AB5C39FA63440, 0x8DD9472BF3FEFAA7}, // e =   555, k =  294
        {0xFA856334878FC150, 0xB14F98F6F0FEB951}, // e =   557, k =  295
        {0x9C935E00D4B9D8D2, 0x6ED1BF9A569F33D3}, // e =   560, k =  296
        {0xC3B8358109E84F07, 0x0A862F80EC4700C8}, // e =   562, k =  297
        {0xF4A642E14C6262C8, 0xCD27BB612758C0FA}, // e =   564, k =  298
        {0x98E7E9CCCFBD7DBD, 0x8038D51CB897789C}, // e =   567, k =  299
        {0xBF21E44003ACDD2C, 0xE0470A63E6BD56C3}, // e =   569, k =  300
        {0xEEEA5D5004981478, 0x1858CCFCE06CAC74}, // e =   571, k =  301
        {0x95527A5202DF0CCB, 0x0F37801E0C43EBC8}, // e =   574, k =  302
        {0xBAA718E68396CFFD, 0xD30560258F54E6BA}, // e =   576, k =  303
        {0xE950DF20247C83FD, 0x47C6B82EF32A2069}, // e =   578, k =  304
        {0x91D28B7416CDD27E, 0x4CDC331D57FA5441}, // e =   581, k =  305
        {0xB6472E511C81471D, 0xE0133FE4ADF8E952}, // e =   583, k =  306
        {0xE3D8F9E563A198E5, 0x58180FDDD97723A6}, // e =   585, k =  307
        {0x8E679C2F5E44FF8F, 0x570F09EAA7EA7648}, // e =   588, k =  308
        {0xB201833B35D63F73, 0x2CD2CC6551E513DA}, // e =   590, k =  309
        {0xDE81E40A034BCF4F, 0xF8077F7EA65E58D1}, // e =   592, k =  310
        {0x8B112E86420F6191, 0xFB04AFAF27FAF782}, // e =   595, k =  311
        {0xADD57A27D29339F6, 0x79C5DB9AF1F9B563}, // e =   597, k =  312
        {0xD94AD8B1C7380874, 0x18375281AE7822BC}, // e =   599, k =  313
        {0x87CEC76F1C830548, 0x8F2293910D0B15B5}, // e =   602, k =  314
        {0xA9C2794AE3A3C69A, 0xB2EB3875504DDB22}, // e =   604, k =  315
        {0xD433179D9C8CB841, 0x5FA60692A46151EB}, // e =   606, k =  316
        {0x849FEEC281D7F328, 0xDBC7C41BA6BCD333}, // e =   609, k =  317
        {0xA5C7EA73224DEFF3, 0x12B9B522906C0800}, // e =   611, k =  318
        {0xCF39E50FEAE16BEF, 0xD768226B34870A00}, // e =   613, k =  319
        {0x81842F29F2CCE375, 0xE6A1158300D46640}, // e =   616, k =  320
        {0xA1E53AF46F801C53, 0x60495AE3C1097FD0}, // e =   618, k =  321
        {0xCA5E89B18B602368, 0x385BB19CB14BDFC4}, // e =   620, k =  322
        {0xFCF62C1DEE382C42, 0x46729E03DD9ED7B5}, // e =   622, k =  323
        {0x9E19DB92B4E31BA9, 0x6C07A2C26A8346D1}, // e =   625, k =  324
        {0xC5A05277621BE293, 0xC7098B7305241885}, // e =   627, k =  325
    };

    RYU_ASSERT(k >= MinDecExp);
    RYU_ASSERT(k <= MaxDecExp);
    return Pow5[static_cast<unsigned>(k - MinDecExp)];
}

#if defined(__SIZEOF_INT128__)

static inline uint64_t MulShift(uint64_t m, const uint64x2* mul, int32_t j)
{
    __extension__ using uint128_t = unsigned __int128;

    RYU_ASSERT(j >= 64 + 1);
    RYU_ASSERT(j <= 64 + 127);

    const uint128_t b0 = uint128_t{m} * mul->lo;
    const uint128_t b2 = uint128_t{m} * mul->hi;

#if 1
    const int32_t shift = j - 64;
#else
    RYU_ASSERT(j <= 64 + 63);
    // We need shift = j - 64 here.
    // Since 64 < j < 128, this is equivalent to shift = (j - 64) % 64 = j % 64.
    // When written as shift = j & 63, clang can optimize the 128-bit shift into
    // a simple funnel shift.
    const int32_t shift = j & 63;
#endif
    return static_cast<uint64_t>((b2 + static_cast<uint64_t>(b0 >> 64)) >> shift);
}

#elif defined(_MSC_VER) && defined(_M_X64)

static inline uint64_t MulShift(uint64_t m, const uint64x2* mul, int32_t j)
{
    RYU_ASSERT(j >= 64 + 1);
    RYU_ASSERT(j <= 64 + 127);

    uint64_t b0_hi;
    uint64_t b0_lo = _umul128(m, mul->lo, &b0_hi);
    uint64_t b2_hi;
    uint64_t b2_lo = _umul128(m, mul->hi, &b2_hi);
    static_cast<void>(b0_lo);

    // b2 + (b0 >> 64)
    // b2_lo += b0_hi;
    // b2_hi += b2_lo < b0_hi;
    _addcarry_u64(_addcarry_u64(0, b2_lo, b0_hi, &b2_lo), b2_hi, 0, &b2_hi);

#if 1
    const int32_t shift = j - 64;
    const uint64_t l = __shiftright128(b2_lo, b2_hi, static_cast<unsigned char>(shift));
    const uint64_t h = __ull_rshift(b2_hi, shift);
    return (shift & 64) ? h : l;
#else
    RYU_ASSERT(j <= 64 + 63);
    // We need shift = j - 64 here.
    // For the __shiftright128 intrinsic, the shift value is always modulo 64.
    // Since (j - 64) % 64 = j, we can simply use j here.
    return __shiftright128(b2_lo, b2_hi, static_cast<unsigned char>(j));
#endif
}

#else

static inline uint64x2 Mul128(uint64_t a, uint64_t b)
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

static inline uint64_t ShiftRight128(uint64_t lo, uint64_t hi, int32_t n)
{
    // For the __shiftright128 intrinsic, the shift value is always modulo 64.
    // In the current implementation of the double-precision version of Ryu, the
    // shift value is always < 64.
    // Check this here in case a future change requires larger shift values. In
    // this case this function needs to be adjusted.
    RYU_ASSERT(n >= 1);
    RYU_ASSERT(n <= 63);

    const int32_t lshift = -n & 63;
    const int32_t rshift =  n;
    return (hi << lshift) | (lo >> rshift);
}

static inline uint64_t MulShift(uint64_t m, const uint64x2* mul, int32_t j)
{
    RYU_ASSERT(j >= 65);
    RYU_ASSERT(j <= 64 + 96 - 1);

    auto b0 = Mul128(m, mul->lo);
    auto b2 = Mul128(m, mul->hi);

    // b2 + (b0 >> 64)
    b2.lo += b0.hi;
    b2.hi += b2.lo < b0.hi;

#if 1
    const int32_t shift = j - 64; // [0, 128)
    // NB:
    // shift >= 64 ==> 0 <= shift - 64 <= 31 (use 32-bit intrinsics?)
    return shift <= 63 ? ShiftRight128(b2.lo, b2.hi, shift) : b2.hi >> (shift - 64);
#else
    RYU_ASSERT(j <= 64 + 63);
    const int32_t shift = j - 64;
    return ShiftRight128(b2.lo, b2.hi, shift);
#endif
}

#endif

static inline void MulPow5DivPow2_Double(uint64_t u, uint64_t v, uint64_t w, int32_t e5, int32_t e2, uint64_t& a, uint64_t& b, uint64_t& c)
{
    // j >= 121 and m has at most 53 + 2 = 55 bits.
    // The product along with the subsequent shift therefore requires
    // 55 + 128 - 121 = 62 bits.

    const auto k = FloorLog2Pow5(e5) + 1 - BitsPerPow5_Double;
    const auto j = e2 - k;
    RYU_ASSERT(j >= BitsPerPow5_Double - 7); // 121 - 64 = 57
    RYU_ASSERT(j <= BitsPerPow5_Double - 1); // 127 - 64 = 63

    const auto pow5 = ComputePow5_Double(e5);

    a = MulShift(u, &pow5, j);
    b = MulShift(v, &pow5, j);
    c = MulShift(w, &pow5, j);
}

// Returns whether value is divisible by 5^e5
static inline bool MultipleOfPow5(uint64_t value, int32_t e5)
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
static inline bool MultipleOfPow2(uint64_t value, int32_t e2)
{
    RYU_ASSERT(e2 >= 0);
    RYU_ASSERT(e2 <= 63);

    return (value & ((uint64_t{1} << e2) - 1)) == 0;
}

namespace {
struct FloatingDecimal64 {
    uint64_t digits; // num_digits <= 17
    int32_t exponent;
};
}

// TODO:
// Export?!
static inline FloatingDecimal64 ToDecimal64(uint64_t ieee_significand, uint64_t ieee_exponent)
{
    //
    // Step 1:
    // Decode the floating point number, and unify normalized and subnormal cases.
    //

    uint64_t m2;
    int32_t e2;
    if (ieee_exponent == 0)
    {
        m2 = ieee_significand;
        e2 = 1 - Double::ExponentBias;
    }
    else
    {
        m2 = Double::HiddenBit | ieee_significand;
        e2 = static_cast<int32_t>(ieee_exponent) - Double::ExponentBias;

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

    const uint32_t lower_boundary_is_closer = (ieee_significand == 0 && ieee_exponent > 1);

    e2 -= 2;
    const uint64_t u = 4 * m2 - 2 + lower_boundary_is_closer;
    const uint64_t v = 4 * m2;
    const uint64_t w = 4 * m2 + 2;

    //
    // Step 3:
    // Convert to a decimal power base.
    //

    int32_t e10;

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

        const int32_t q = FloorLog10Pow2(e2) - (e2 > 3); // == max(0, q' - 1)
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

        const int32_t q = FloorLog10Pow5(-e2) - (-e2 > 1); // == max(0, q' - 1)
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

static inline char* PrintDecimalDigitsBackwards(char* buf, uint64_t output)
{
    // We prefer 32-bit operations, even on 64-bit platforms.
    // We have at most 17 digits, and uint32_t can store 9 digits.
    // If output doesn't fit into uint32_t, we cut off 8 digits,
    // so the rest will fit into uint32_t.
    if (static_cast<uint32_t>(output >> 32) != 0)
    {
        const uint64_t q = output / 100000000;
        const uint32_t r = static_cast<uint32_t>(output % 100000000);
        output = q;
        buf -= 8;
        Utoa_8Digits(buf, r);
    }

    RYU_ASSERT(output <= UINT32_MAX);
    uint32_t output2 = static_cast<uint32_t>(output);

    // (Runs up to 4 times...)
    while (output2 >= 100)
    {
        const uint32_t q = output2 / 100;
        const uint32_t r = output2 % 100;
        output2 = q;
        buf -= 2;
        Utoa_2Digits(buf, r);
    }

    if (output2 >= 10)
    {
        buf -= 2;
        Utoa_2Digits(buf, output2);
    }
    else
    {
        *--buf = static_cast<char>('0' + output2);
    }

    return buf;
}

static inline int32_t DecimalLength(uint64_t v)
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

static inline char* FormatDigits(char* buffer, uint64_t digits, int32_t decimal_exponent, bool force_trailing_dot_zero = false)
{
    static constexpr int32_t MinFixedDecimalPoint = -6;
    static constexpr int32_t MaxFixedDecimalPoint =  17;
    static_assert(MinFixedDecimalPoint <= -1, "internal error");
    static_assert(MaxFixedDecimalPoint >= 17, "internal error");

    RYU_ASSERT(digits >= 1);
    RYU_ASSERT(digits <= 99999999999999999ull);
    RYU_ASSERT(decimal_exponent >= -999);
    RYU_ASSERT(decimal_exponent <=  999);

    const int32_t num_digits = DecimalLength(digits);
    const int32_t decimal_point = num_digits + decimal_exponent;

    const bool use_fixed = MinFixedDecimalPoint <= decimal_point && decimal_point <= MaxFixedDecimalPoint;

    // Prepare the buffer.
    // Avoid calling memset/memcpy with variable arguments below...

    int32_t decimal_digits_position;
    if (use_fixed)
    {
        if (decimal_point <= 0)
        {
            // 0.[000]digits
            decimal_digits_position = 2 - decimal_point;
            static_assert(MinFixedDecimalPoint >= -14, "internal error");
            std::memcpy(buffer, "0.00000000000000", 16);
        }
        else if (decimal_point < num_digits)
        {
            // dig.its
            decimal_digits_position = 0;
        }
        else
        {
            // digits[000]
            decimal_digits_position = 0;
            static_assert(MaxFixedDecimalPoint <= 32, "internal error");
            std::memset(buffer +  0, '0', 16);
            std::memset(buffer + 16, '0', 16);
        }
    }
    else
    {
        // dE+123 or d.igitsE+123
        decimal_digits_position = 1;
    }

    char* const digits_end = buffer + decimal_digits_position + num_digits;
    PrintDecimalDigitsBackwards(digits_end, digits);

    if (use_fixed)
    {
        if (decimal_point <= 0)
        {
            // 0.[000]digits
            buffer = digits_end;
        }
        else if (decimal_point < num_digits)
        {
            // dig.its
#if defined(_MSC_VER) && !defined(__clang__)
            // VC does not inline the memmove call below. (Even if compiled with /arch:AVX2.)
            // However, memcpy will be inlined.
            uint8_t tmp[16];
            char* const src = buffer + decimal_point;
            char* const dst = src + 1;
            std::memcpy(tmp, src, 16);
            std::memcpy(dst, tmp, 16);
#else
            std::memmove(buffer + decimal_point + 1, buffer + decimal_point, 16);
#endif
            buffer[decimal_point] = '.';
            buffer = digits_end + 1;
        }
        else
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
            ++buffer;
        }
        else
        {
            // d.igitsE+123
            buffer[1] = '.';
            buffer = digits_end;
        }

        const int32_t scientific_exponent = decimal_point - 1;
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
            const uint32_t q = k / 100;
            const uint32_t r = k % 100;
            *buffer++ = static_cast<char>('0' + q);
            buffer = Utoa_2Digits(buffer, r);
        }
    }

    return buffer;
}

static inline char* ToChars(char* buffer, double value, bool force_trailing_dot_zero = false)
{
    const Double v(value);

    const uint64_t significand = v.PhysicalSignificand();
    const uint64_t exponent = v.PhysicalExponent();

    if (exponent != Double::MaxIeeeExponent) // [[likely]]
    {
        // Finite

        buffer[0] = '-';
        buffer += v.SignBit();

        if (exponent != 0 || significand != 0) // [[likely]]
        {
            // != 0

            const auto dec = ToDecimal64(significand, exponent);
            return FormatDigits(buffer, dec.digits, dec.exponent, force_trailing_dot_zero);
        }
        else
        {
            std::memcpy(buffer, "0.0 ", 4);
            buffer += force_trailing_dot_zero ? 3 : 1;
            return buffer;
        }
    }

    if (significand == 0)
    {
        buffer[0] = '-';
        buffer += v.SignBit();

        std::memcpy(buffer, "inf ", 4);
        return buffer + 3;
    }
    else
    {
        std::memcpy(buffer, "nan ", 4);
        return buffer + 3;
    }
}

//==================================================================================================
//
//==================================================================================================

char* ryu::Dtoa(char* buffer, double value)
{
    return ToChars(buffer, value);
}

//==================================================================================================
// ToBinary64
//==================================================================================================

// Maximum number of decimal digits in the significand the fast ToBinary method can handle.
// Inputs with more significant digits must be processed using another algorithm.
static constexpr int32_t ToBinaryMaxDecimalDigits = 17;

// Any input <= 10^MinDecimalExponent is interpreted as 0.
// Any input >  10^MaxDecimalExponent is interpreted as +Infinity.
static constexpr int32_t MinDecimalExponent = -324; // denorm_min / 2 = 2.4703282292062327e-324 >= 10^-324
static constexpr int32_t MaxDecimalExponent =  309; //            max = 1.7976931348623158e+308 <= 10^+309

static inline int32_t FloorLog2(uint64_t x)
{
    RYU_ASSERT(x != 0);

#if defined(__GNUC__) || defined(__clang__)
    return 63 - __builtin_clzll(x);
#elif defined(_MSC_VER) && defined(_M_X64)
    unsigned long index;
    _BitScanReverse64(&index, x);
    return static_cast<int32_t>(index);
#elif defined(_MSC_VER) && defined(_M_IX86)
    unsigned long index;
    if (_BitScanReverse(&index, static_cast<uint32_t>(x >> 32)))
        return 32 + static_cast<int32_t>(index);
    _BitScanReverse(&index, static_cast<uint32_t>(x));
    return static_cast<int32_t>(index);
#else
    int32_t l2 = 0;
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

static inline int32_t FloorLog2Pow10(int32_t e)
{
    RYU_ASSERT(e >= -1233);
    RYU_ASSERT(e <=  1233);
    return FloorDivPow2(e * 1741647, 19);
}

static inline int32_t ExtractBit(uint64_t x, int32_t n)
{
    RYU_ASSERT(n >= 0);
    RYU_ASSERT(n <= 63);
    return (x & (uint64_t{1} << n)) != 0;
}

static inline double ToBinary64(uint64_t m10, int32_t m10_digits, int32_t e10)
{
    static constexpr int32_t MantissaBits = Double::SignificandSize - 1;
    static constexpr int32_t ExponentBias = Double::ExponentBias - (Double::SignificandSize - 1);

    RYU_ASSERT(m10 > 0);
    RYU_ASSERT(m10_digits == DecimalLength(m10));
    RYU_ASSERT(m10_digits <= ToBinaryMaxDecimalDigits);
    RYU_ASSERT(e10 >  MinDecimalExponent - m10_digits);
    RYU_ASSERT(e10 <= MaxDecimalExponent - m10_digits);
    static_cast<void>(m10_digits);

    // e10 >= MinDecimalExponent - m10_digits + 1 >= -324 - 17 + 1 = -340
    // e10 <= MaxDecimalExponent - m10_digits     <=  309 -  1     =  308

#if defined(__x86_64__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
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

            const int32_t zeros = 15 - m10_digits;
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
    RYU_ASSERT(log2_m10 <= 56); // 56 = floor(log_2(10^17))

    // The length of m10 * 10^e10 in bits is: log2(m10 * 10^e10) = log2(m10) + log2(10^e10).
    // We want to compute the (MantissaBits + 1) top-most bits (+1 for the implicit leading
    // one in IEEE format). We therefore choose a binary output exponent of
    //   e2 = log2(m10 * 10^e10) - (MantissaBits + 1).
    //
    // We use floor(log2(5^e10)) so that we get at least this many bits; better to have an
    // additional bit than to not have enough bits.

    // We compute [m10 * 10^e10 / 2^e2] == [m10 * 5^e10 / 2^(e2 - e10)]
    //
    // Let b = floor(log_2(m10))
    // Let n = floor(log_2(5^e10))
    // Then
    //  j = ( e2 - e10 ) - ( n + 1 - BitsPerPow5 )
    //    = ( ( b + e10 + n - (MantissaBits + 1) ) - e10 ) - ( n + 1 - BitsPerPow5 )
    //    = b + BitsPerPow5 - MantissaBits - 2
    //    = b + 128 - 52 - 2
    //    = b + 74
    // Since 0 <= b <= 56, we have
    //    74 <= j <= 130
    // The product along with the subsequent shift therefore has (at most)
    //  b + 130 - (130 - 54 + b) = 54
    // bits.

    const auto log2_10_e10 = FloorLog2Pow10(e10);
    const auto e2 = log2_m10 + log2_10_e10 - (MantissaBits + 1);

    const auto pow5 = ComputePow5_Double(e10);
    const auto j = log2_m10 + (BitsPerPow5_Double - MantissaBits - 2);
    RYU_ASSERT(j >= 74);
    RYU_ASSERT(j <= 130);
    const auto m2 = MulShift(m10, &pow5, j);

    const auto log2_m2 = FloorLog2(m2);
    RYU_ASSERT(log2_m2 >= 53);
    RYU_ASSERT(log2_m2 <= 54);

    // We also compute if the result is exact, i.e., [m10 * 10^e10 / 2^e2] == m10 * 10^e10 / 2^e2.
    //  (See: Ryu Revisited, Section 4.3)

    bool is_exact = (e2 <= e10) || (e2 - e10 < 64 && MultipleOfPow2(m10, e2 - e10));
    if (e10 >= 0)
    {
        // 2^(e2 - e10) | m10 5^e10
        //  <==> p2(m10 5^e10)       >= e2 - e10
        //  <==> p2(m10) + e10 p2(5) >= e2 - e10
        //  <==> p2(m10)             >= e2 - e10

        // is_exact
        //  <==>   (e2 <= e10   OR   p2(m10) >= e2 - e10)

        // is_exact = (e2 <= e10) || (e2 - e10 < 64 && MultipleOfPow2(m10, e2 - e10));
    }
    else
    {
        // e2 <= e10:
        //
        // m10 10^e10 / 2^e2
        //  == m10 2^e10 5^e10 / 2^e2
        //  == m10 2^(e10 - e2) / 5^(-e10)
        //
        // 5^(-e10) | m10 2^(e10 - e2)
        //  <==> p5(m10 2^(e10 - e2))       >= -e10
        //  <==> p5(m10) + (e10 - e2) p5(2) >= -e10   (p5(2) == 0)
        //  <==> p5(m10)                    >= -e10

        // e2 > e10:
        //
        // m10 10^e10 / 2^e2
        //  == m10 (2^e10 5^e10) / 2^e2
        //  == m10 / (5^(-e10) 2^(e2 - e10))
        //  == m10 / (10^(-e10) 2^e2)
        //
        // 5^(-e10) 2^(e2 - e10) | m10
        //  <==> 5^(-e10) | m10   AND   2^(e2 - e10) | m10
        //  <==> p5(m10) >= -e10   AND   p2(m10) >= e2 - e10

        // is_exact
        //  <==>   (e2 <= e10   OR   p2(m10) >= e2 - e10)   AND   p5(m10) >= -e10

        // e2 <= e10 ==> is_exact = true
        // In this case we need to check p5(m10) >= -e10.
        // Check that the test below works.
        RYU_ASSERT(e2 > e10 || is_exact);

        // 57 = ceil(log_2(10^17))
        // 24 = floor(log_5(2^57))
        is_exact = is_exact && (-e10 <= 24 && MultipleOfPow5(m10, -e10));
    }

    // Compute the final IEEE exponent.
    int32_t ieee_e2 = Max(0, log2_m2 + e2 + ExponentBias);
    if (ieee_e2 >= Double::MaxIeeeExponent)
    {
        // Overflow:
        // Final IEEE exponent is larger than the maximum representable.
        return std::numeric_limits<double>::infinity();
    }

    // We need to figure out how much we need to shift m2.
    // The tricky part is that we need to take the final IEEE exponent into account, so we need to
    // reverse the bias and also special-case the value 0.
    const int32_t shift = (ieee_e2 == 0 ? 1 : ieee_e2) - e2 - (ExponentBias + MantissaBits);
    RYU_ASSERT(shift > 0);
    RYU_ASSERT(shift < 64);

    // We need to round up if the exact value is more than 0.5 above the value we computed. That's
    // equivalent to checking if the last removed bit was 1 and either the value was not just
    // trailing zeros or the result would otherwise be odd.
    const auto trailing_zeros
        = is_exact && MultipleOfPow2(m2, shift - 1);
    const auto last_removed_bit
        = ExtractBit(m2, shift - 1);
    const auto round_up
        = last_removed_bit != 0 && (!trailing_zeros || ExtractBit(m2, shift) != 0);

    uint64_t significand = (m2 >> shift) + round_up;
    RYU_ASSERT(significand <= 2 * Double::HiddenBit); // significand <= 2^p = 2^53

    significand &= Double::SignificandMask;

    // Rounding up may cause overflow...
    if (significand == 0 && round_up)
    {
        // Rounding up did overflow the p-bit significand.
        // Move a trailing zero of the significand into the exponent.
        // Due to how the IEEE represents +/-Infinity, we don't need to check for overflow here.
        ++ieee_e2;
    }

    RYU_ASSERT(ieee_e2 <= Double::MaxIeeeExponent);
    const uint64_t ieee = static_cast<uint64_t>(ieee_e2) << MantissaBits | significand;
    return ReinterpretBits<double>(ieee);
}

//==================================================================================================
// Strtod
//==================================================================================================

using ryu::StrtodStatus;
using ryu::StrtodResult;

static inline bool IsDigit(char ch)
{
    return static_cast<unsigned>(ch - '0') <= 9u;
}

static inline int32_t DigitValue(char ch)
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
    RYU_ASSERT(*next == 'i' || *next == 'I');

    if (!StartsWith(next + 1, last, "nf"))
        return {next, StrtodStatus::invalid};

    next += 3;
    if (StartsWith(next, last, "inity"))
        next += 5;

    return {next, StrtodStatus::inf};
}

static inline bool IsNaNSequenceChar(char ch)
{
    return ch == '_' || IsDigit(ch) || IsUpperASCII(ch) || IsLowerASCII(ch);
}

// FIXME:
// Don't ignore the nan-sequence!!!
static inline StrtodResult ParseNaN(const char* next, const char* last)
{
    RYU_ASSERT(*next == 'n' || *next == 'N');

    if (!StartsWith(next + 1, last, "an"))
        return {next, StrtodStatus::invalid};

    next += 3;
    if (next != last && *next == '(')
    {
        for (const char* p = next + 1; p != last; ++p)
        {
            if (*p == ')')
                return {p + 1, StrtodStatus::nan};

            if (!IsNaNSequenceChar(*p))
                break; // invalid/incomplete nan-sequence
        }
    }

    return {next, StrtodStatus::nan};
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

#if RYU_STRTOD_FALLBACK()
static RYU_NEVER_INLINE double ToBinary64Slow(const char* next, const char* last)
{
#if HAS_CHARCONV()
    double flt = 0;
    std::from_chars(next, last, flt);
    return flt;
#else
    //
    // FIXME:
    // _strtod_l( ..., C_LOCALE )
    //

    // std::strtod expects null-terminated inputs. So we need to make a copy and null-terminate the input.
    // This function is actually almost never going to be called, so that should be ok.
    const std::string inp(next, last);
    const char* const ptr = inp.c_str();

    char* end;
    const auto flt = ::strtod(ptr, &end);

    // std::strtod should have consumed all of the input.
    RYU_ASSERT(ptr != end);
    RYU_ASSERT(last - next == end - ptr);

    return flt;
#endif
}
#endif

StrtodResult ryu::Strtod(const char* next, const char* last, double& value)
{
    if (next == last)
        return {next, StrtodStatus::invalid};

    // Decompose the input into the form significand * 10^exponent,
    // where significand has num_digits decimal digits.

    uint64_t significand = 0; // only valid iff num_digits <= 19
    int64_t  num_digits  = 0; // 64-bit to avoid overflow...
    int64_t  exponent    = 0; // 64-bit to avoid overflow...
    StrtodStatus status = StrtodStatus::integer;

// [-]

    const bool is_negative = (*next == '-');
    if (is_negative || *next == '+')
    {
        ++next;
        if (next == last)
            return {next, StrtodStatus::invalid};
    }

// int32_t

#if RYU_STRTOD_FALLBACK()
    const char* const start = next;
#endif

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
        status = StrtodStatus::fixed;

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
    static constexpr int32_t MaxExp = 999999;
    static_assert(MaxExp >= 999, "invalid parameter");
    static_assert(MaxExp <= (INT_MAX - 9) / 10, "invalid parameter");

    int32_t parsed_exponent = 0;
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
                // Found a valid exponent.
                status = StrtodStatus::scientific;
                next = p;

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
        flt = ToBinary64(significand, static_cast<int32_t>(num_digits), static_cast<int32_t>(exponent));
    }
    else
    {
        // We need to fall back to another algorithm if the input is too long.
#if RYU_STRTOD_FALLBACK()
        flt = ToBinary64Slow(start, next);
#else
        return {next, StrtodStatus::input_too_long};
#endif
    }

    value = is_negative ? -flt : flt;
    return {next, status};
}

//==================================================================================================
// Round10
//==================================================================================================

static inline uint64_t SmallPow10(int32_t e10)
{
    static constexpr uint64_t Pow10Table[] = {
        1,
        10,
        100,
        1000,
        10000,
        100000,
        1000000,
        10000000,
        100000000,
        1000000000,
        10000000000,
        100000000000,
        1000000000000,
        10000000000000,
        100000000000000,
        1000000000000000,
        10000000000000000,
        100000000000000000,
    };

    RYU_ASSERT(e10 >= 0);
    RYU_ASSERT(e10 <= 17);
    return Pow10Table[static_cast<uint32_t>(e10)];
}

static double MulRoundDiv(const double value, const int32_t mul_e10, const int32_t div_e10)
{
    const Double v(value);

    const uint64_t F = v.PhysicalSignificand();
    const uint64_t E = v.PhysicalExponent();

    if (E == Double::MaxIeeeExponent || (E == 0 && F == 0))
    {
        // +-0, or Infinity, or NaN
        // Multiplying by 10^n does not change the value.
        return value;
    }

    // Convert to decimal
    const FloatingDecimal64 dec = ToDecimal64(F, E);

    uint64_t digits     = dec.digits;
    int32_t  num_digits = DecimalLength(dec.digits);
    int32_t  exponent   = dec.exponent; // (Using 64 bits to avoid checking for overflow.)

    // Multiply by 10^mul_e10
    exponent += mul_e10;

    // Round x = digits * 10^exponent to the nearest integer.

    // We have
    // x = digits * 10^exponent
    //   = digits / 10^e10
    const int32_t e10 = -exponent;
    if (e10 <= 0)
    {
        // x = digits * 10^exponent, where exponent >= 0.
        // Nothing to do.
    }
    else if (e10 < num_digits)
    {
        // 1 <= x < D

        const uint64_t pow10 = SmallPow10(e10);

        RYU_ASSERT(digits >= pow10);
        const uint64_t i = digits / pow10;
        const uint64_t f = digits % pow10;

        // Round to int (towards +inf)
        digits      = i + (f >= pow10 / 2);
        num_digits -= e10;
        exponent    = 0;
    }
    else if (e10 == num_digits)
    {
        // 1/10 <= x < 1

        // x < 1/2 <==> 10x < 5
        //         <==> 10 (digits / 10^e10) < 5
        //         <==> digits < 5 * 10^(e10 - 1)

        digits      = (digits >= 5 * SmallPow10(e10 - 1)) ? 1 : 0;
        num_digits  = 1;
        exponent    = 0;
    }
    else
    {
        // x < 1/10
        // This definitely rounds to 0.
        digits      = 0;
        num_digits  = 1;
        exponent    = 0;
    }

    // Divide by 10^div_e10
    exponent -= div_e10;

    // And convert back to binary.
    double flt;
    if (digits == 0)
    {
        flt = 0;
    }
    else if (exponent + num_digits <= MinDecimalExponent)
    {
        // x * 10^-inf = 0
        flt = 0;
    }
    else if (exponent + num_digits > MaxDecimalExponent)
    {
        // x * 10^+inf = +inf
        flt = std::numeric_limits<double>::infinity();
    }
    else
    {
        flt = ToBinary64(digits, num_digits, exponent);
    }

    return value < 0 ? -flt : flt;
}

double ryu::Round10(const double value, const int n)
{
    if (n < -1000 || n > +1000) // (Not supported yet)
        return value;

    return MulRoundDiv(value, -n, -n);
}
