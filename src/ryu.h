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

#pragma once

#include <cassert>
#include <climits>
#include <cstdint>
#include <cstring>
#include <limits>

#ifndef RYU_ASSERT
#define RYU_ASSERT(X) assert(X)
#endif

#ifndef RYU_INLINE
#define RYU_INLINE inline
#endif

#if defined(_M_IX86) || defined(_M_ARM) || defined(__i386__) || defined(__arm__)
#define RYU_32_BIT_PLATFORM 1
#endif

#if defined(__SIZEOF_INT128__)
#define RYU_HAS_UINT128 1
#elif defined(_M_X64) || defined(__x86_64__)
#define RYU_HAS_X64_INTRINSICS 1
#define RYU_HAS_X86_INTRINSICS 1
#elif defined(_M_IX86) || defined(__i386__)
#define RYU_HAS_X86_INTRINSICS 1
#endif

#if RYU_HAS_X64_INTRINSICS || RYU_HAS_X86_INTRINSICS
#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
#endif

namespace ryu {
namespace impl {

//==================================================================================================
//
//==================================================================================================

template <typename Dest, typename Source>
RYU_INLINE Dest ReinterpretBits(Source source)
{
    static_assert(sizeof(Dest) == sizeof(Source), "size mismatch");

    Dest dest;
    std::memcpy(&dest, &source, sizeof(Source));
    return dest;
}

template <int Precision> struct BitsType;
template <> struct BitsType<24> { using type = uint32_t; };
template <> struct BitsType<53> { using type = uint64_t; };

template <typename Float>
struct IEEE
{
    // NB:
    // Works for double == long double.
    static_assert(std::numeric_limits<Float>::is_iec559 &&
                  ((std::numeric_limits<Float>::digits == 24 && std::numeric_limits<Float>::max_exponent == 128) ||
                   (std::numeric_limits<Float>::digits == 53 && std::numeric_limits<Float>::max_exponent == 1024)),
        "IEEE-754 single- or double-precision implementation required");

    using value_type = Float;
    using bits_type = typename BitsType<std::numeric_limits<Float>::digits>::type;

    static constexpr int       SignificandSize = std::numeric_limits<value_type>::digits; // = p   (includes the hidden bit)
    static constexpr int       ExponentBias    = std::numeric_limits<value_type>::max_exponent - 1 + (SignificandSize - 1);
    static constexpr int       MaxExponent     = std::numeric_limits<value_type>::max_exponent - 1 - (SignificandSize - 1);
    static constexpr int       MinExponent     = std::numeric_limits<value_type>::min_exponent - 1 - (SignificandSize - 1);
    static constexpr bits_type HiddenBit       = bits_type{1} << (SignificandSize - 1);   // = 2^(p-1)
    static constexpr bits_type SignificandMask = HiddenBit - 1;                           // = 2^(p-1) - 1
    static constexpr bits_type ExponentMask    = bits_type{2 * std::numeric_limits<value_type>::max_exponent - 1} << (SignificandSize - 1);
    static constexpr bits_type SignMask        = ~(~bits_type{0} >> 1);

    bits_type bits;

    explicit IEEE(bits_type bits_) : bits(bits_) {}
    explicit IEEE(value_type value) : bits(ReinterpretBits<bits_type>(value)) {}

    bits_type PhysicalSignificand() const {
        return bits & SignificandMask;
    }

    bits_type PhysicalExponent() const {
        return (bits & ExponentMask) >> (SignificandSize - 1);
    }

    // Returns the significand for a normalized double.
    bits_type NormalizedSignificand() const {
        return HiddenBit | PhysicalSignificand();
    }

    // Returns the exponent for a normalized double.
    int NormalizedExponent() const {
        return static_cast<int>(PhysicalExponent()) - ExponentBias;
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

using Double = IEEE<double>;
using Single = IEEE<float>;

//==================================================================================================
// ToDecimal
//
// Double-precision implementation
//==================================================================================================
// Constant data = 9872 bytes

struct Uint64x2 {
    uint64_t hi;
    uint64_t lo;
};

RYU_INLINE Uint64x2 ComputePow5Double(int k)
{
    // Let e = FloorLog2Pow5(k) + 1 - 128
    // For k >= 0, stores 5^k in the form: floor( 5^k / 2^e )
    // For k <= 0, stores 5^k in the form:  ceil(2^-e / 5^-k)
    static constexpr int MinDecExp = -291;
    static constexpr int MaxDecExp =  325;
    static constexpr Uint64x2 Pow5[MaxDecExp - MinDecExp + 1] = {
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
    return Pow5[k - MinDecExp];
}

// Returns floor(x / 2^n).
RYU_INLINE int FloorDivPow2(int x, int n)
{
    // Technically, right-shift of negative integers is implementation defined...
    // Should easily be optimized into SAR (or equivalent) instruction.
#if 1
    return x < 0 ? ~(~x >> n) : (x >> n);
#else
    return x >> n;
#endif
}

RYU_INLINE int FloorLog2Pow5(int e)
{
    RYU_ASSERT(e >= -1764);
    RYU_ASSERT(e <=  1763);
    return FloorDivPow2(e * 1217359, 19);
}

RYU_INLINE int FloorLog10Pow2(int e)
{
    RYU_ASSERT(e >= -2620);
    RYU_ASSERT(e <=  2620);
    return FloorDivPow2(e * 315653, 20);
}

RYU_INLINE int FloorLog10Pow5(int e)
{
    RYU_ASSERT(e >= -2620);
    RYU_ASSERT(e <=  2620);
    return FloorDivPow2(e * 732923, 20);
}

RYU_INLINE uint32_t Lo32(uint64_t x) {
    return static_cast<uint32_t>(x);
}

RYU_INLINE uint32_t Hi32(uint64_t x) {
    return static_cast<uint32_t>(x >> 32);
}

RYU_INLINE uint64_t Load64(uint32_t lo, uint32_t hi) {
    return lo | uint64_t{hi} << 32;
}

// Add unsigned 64-bit integers X and Y with unsigned 8-bit carry-in C, and stores the unsigned 64-bit result in OUT.
// Returns the carry-out.
RYU_INLINE unsigned char Addc64(unsigned char c, uint64_t x, uint64_t y, uint64_t* out)
{
#if RYU_HAS_X64_INTRINSICS
    return _addcarry_u64(c, x, y, reinterpret_cast<unsigned long long*>(out));
#elif RYU_HAS_X86_INTRINSICS
    unsigned int lo;
    unsigned int hi;
    c = _addcarry_u32(c, Lo32(x), Lo32(y), &lo);
    c = _addcarry_u32(c, Hi32(x), Hi32(y), &hi);
    out[0] = Load64(lo, hi);
    return c;
#else
    const uint64_t s = x + (y + c);
    out[0] = s;
    return c ? (s <= x) : (s < x);
#endif
}

// Add unsigned 8-bit borrow C (carry flag) to unsigned 64-bit integer Y, and subtract the result from unsigned 64-bit integer X.
// Stores the unsigned 64-bit result in OUT.
// Returns the carry-out.
RYU_INLINE unsigned char Subb64(unsigned char c, uint64_t x, uint64_t y, uint64_t* out)
{
#if RYU_HAS_X64_INTRINSICS
    return _subborrow_u64(c, x, y, reinterpret_cast<unsigned long long*>(out));
#elif RYU_HAS_X86_INTRINSICS
    unsigned int lo;
    unsigned int hi;
    c = _subborrow_u32(c, Lo32(x), Lo32(y), &lo);
    c = _subborrow_u32(c, Hi32(x), Hi32(y), &hi);
    out[0] = Load64(lo, hi);
    return c;
#else
    const uint64_t s = x - (y + c);
    out[0] = s;
    return c ? (s >= x) : (s > x);
#endif
}

RYU_INLINE Uint64x2 Mul128(uint64_t a, uint64_t b)
{
#if RYU_HAS_UINT128
    __extension__ using uint128_t = unsigned __int128;

    const uint128_t product = uint128_t{a} * b;

    const uint64_t lo = static_cast<uint64_t>(product);
    const uint64_t hi = static_cast<uint64_t>(product >> 64);
    return {hi, lo};
#elif RYU_HAS_X64_INTRINSICS
    uint64_t hi;
    uint64_t lo = _umul128(a, b, &hi);
    return {hi, lo};
#else
    const uint64_t b00 = uint64_t{Lo32(a)} * Lo32(b);
    const uint64_t b01 = uint64_t{Lo32(a)} * Hi32(b);
    const uint64_t b10 = uint64_t{Hi32(a)} * Lo32(b);
    const uint64_t b11 = uint64_t{Hi32(a)} * Hi32(b);

    const uint64_t mid1 = b10 + Hi32(b00);
    const uint64_t mid2 = b01 + Lo32(mid1);

    const uint64_t hi = b11 + Hi32(mid1) + Hi32(mid2);
    const uint64_t lo = Load64(Lo32(b00), Lo32(mid2));
    return {hi, lo};
#endif
}

RYU_INLINE uint64_t ShiftRight128(uint64_t lo, uint64_t hi, int dist)
{
    // For the __shiftright128 intrinsic, the shift value is always modulo 64.
    // In the current implementation of the double-precision version of Ryu, the
    // shift value is always < 64.
    // Check this here in case a future change requires larger shift values. In
    // this case this function needs to be adjusted.
    RYU_ASSERT(dist >= 56); // 56: MulShiftAll fallback, 57: otherwise.
    RYU_ASSERT(dist <= 63);

#if RYU_HAS_UINT128
    __extension__ using uint128_t = unsigned __int128;

    return static_cast<uint64_t>(((uint128_t{hi} << 64) | lo) >> (dist & 63));
#elif RYU_HAS_X64_INTRINSICS
    return __shiftright128(lo, hi, dist);
#else
#if RYU_32_BIT_PLATFORM
    // Avoid a 64-bit shift by taking advantage of the range of shift values.
    // We know that 0 < 64 - dist < 32 and 0 < dist - 32 < 32.
    // Use (64 - dist) = (64 - dist) % 32 = -dist & 31
    // and (dist - 32) = (dist - 32) % 32 =  dist & 31.
#if defined(_MSC_VER) && defined(_M_IX86) && !defined(__clang__)
    const int lshift = -dist; // % 32 is implicit in __ll_lshift
    const int rshift =  dist & 31;
    return __ll_lshift(hi, lshift) | (Hi32(lo) >> rshift);
#else
    const int lshift = -dist & 31;
    const int rshift =  dist & 31;
    return (hi << lshift) | (Hi32(lo) >> rshift);
#endif
#else
    const int lshift = -dist & 31;
    const int rshift =  dist & 63;
    return (hi << lshift) | (lo >> rshift);
#endif
#endif
}

#if RYU_32_BIT_PLATFORM

RYU_INLINE void MulShiftAll(uint64_t mv, uint64_t mp, uint64_t mm, const Uint64x2* mul, int j, uint64_t& vr, uint64_t& vp, uint64_t& vm)
{
    // m is maximum 55 bits
    unsigned char c;

    // m = 2*m2
    const uint64_t m = mv >> 1;
    const uint32_t mmShift = static_cast<uint32_t>(mv - 1 - mm);

    const Uint64x2 b0 = Mul128(m, mul->lo);
    const Uint64x2 b2 = Mul128(m, mul->hi);

    // vr = 2*m2 * mul
    uint64_t vr0 = b0.lo;
    uint64_t vr1;
    uint64_t vr2;
    c = Addc64(0, b2.lo, b0.hi, &vr1);
    c = Addc64(c, b2.hi, 0,     &vr2);
    vr = ShiftRight128(vr1, vr2, j - 64 - 1);

    // vp = 2*m2 * mul + mul
    uint64_t vp0;
    uint64_t vp1;
    uint64_t vp2;
    c = Addc64(0, vr0, mul->lo, &vp0);
    c = Addc64(c, vr1, mul->hi, &vp1);
    c = Addc64(c, vr2, 0,       &vp2);
    vp = ShiftRight128(vp1, vp2, j - 64 - 1);

    // vm = 2*m2 * mul - mul
    //   or 4*m2 * mul - mul (not common)
    uint64_t vm0;
    uint64_t vm1;
    uint64_t vm2;
    if (mmShift == 0)
    {
        // vr = 2 * vr
        c = Addc64(0, vr0, vr0, &vr0);
        c = Addc64(c, vr1, vr1, &vr1);
        c = Addc64(c, vr2, vr2, &vr2);
    }
    c = Subb64(0, vr0, mul->lo, &vm0);
    c = Subb64(c, vr1, mul->hi, &vm1);
    c = Subb64(c, vr2, 0,       &vm2);
    vm = ShiftRight128(vm1, vm2, (j - 64 - 1) + (mmShift == 0));
}

#else

RYU_INLINE uint64_t MulShift(uint64_t m, const Uint64x2* mul, int j)
{
    RYU_ASSERT((m >> 55) == 0); // m is maximum 55 bits

#if RYU_HAS_UINT128
    __extension__ using uint128_t = unsigned __int128;

    const uint128_t b0 = uint128_t{m} * mul->lo;
    const uint128_t b2 = uint128_t{m} * mul->hi;

    // We need shift = j - 64 here.
    // Since 64 < j < 128, this is equivalent to shift = (j - 64) % 64 = j % 64.
    // When written as shift = j & 63, clang and gcc (9+) can optimize the
    // 128-bit shift into a simple funnel shift.
    const int shift = j & 63;
    return static_cast<uint64_t>((b2 + (b0 >> 64)) >> shift);
#elif RYU_HAS_X64_INTRINSICS
    uint64_t b0Hi;
    uint64_t b0Lo = _umul128(m, mul->lo, &b0Hi);
    uint64_t b2Hi;
    uint64_t b2Lo = _umul128(m, mul->hi, &b2Hi);
    static_cast<void>(b0Lo);

    // b2 + (b0 >> 64)
    // b2Lo += b0Hi;
    // b2Hi += b2Lo < b0Hi;
    _addcarry_u64(_addcarry_u64(0, b2Lo, b0Hi, &b2Lo), b2Hi, 0, &b2Hi);

    // We need shift = j - 64 here.
    // For the __shiftright128 intrinsic, the shift value is always modulo 64.
    // Since (j - 64) % 64 = j, we can simply use j here.
    return __shiftright128(b2Lo, b2Hi, static_cast<unsigned char>(j));
#else
    auto b0 = Mul128(m, mul->lo);
    auto b2 = Mul128(m, mul->hi);

    // b2 + (b0 >> 64)
    b2.lo += b0.hi;
    b2.hi += b2.lo < b0.hi;

    return ShiftRight128(b2.lo, b2.hi, j & 63);
#endif
}

RYU_INLINE void MulShiftAll(uint64_t mv, uint64_t mp, uint64_t mm, const Uint64x2* mul, int j, uint64_t& vr, uint64_t& vp, uint64_t& vm)
{
    vr = MulShift(mv, mul, j);
    vp = MulShift(mp, mul, j);
    vm = MulShift(mm, mul, j);
}

#endif

#if RYU_32_BIT_PLATFORM

// On x86 platforms, compilers typically generate calls to library
// functions for 64-bit divisions, even if the divisor is a constant.
//
// E.g.:
// https://bugs.llvm.org/show_bug.cgi?id=37932
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=17958
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=37443
//
// The functions here perform division-by-constant using multiplications
// in the same way as 64-bit compilers would do.
//
// NB:
// The multipliers and shift values are the ones generated by clang x64
// for expressions like x/5, x/10, etc.

RYU_INLINE uint64_t Div5(uint64_t x) {
    return Mul128(x, 0xCCCCCCCCCCCCCCCDu).hi >> 2;
}

RYU_INLINE uint64_t Div10(uint64_t x) {
    return Mul128(x, 0xCCCCCCCCCCCCCCCDu).hi >> 3;
}

RYU_INLINE uint64_t Div100(uint64_t x) {
    return Mul128(x >> 2, 0x28F5C28F5C28F5C3u).hi >> 2;
}

RYU_INLINE uint64_t Div1e4(uint64_t x) {
    return Mul128(x, 0x346DC5D63886594Bu).hi >> 11;
}

RYU_INLINE uint64_t Div1e8(uint64_t x) {
    return Mul128(x, 0xABCC77118461CEFDu).hi >> 26;
}

// Avoid 64-bit math as much as possible.
// For D <= 10^9, returning (uint32_t) (x - D * Div<D>(x)) would perform
// 32x64-bit multiplication and 64-bit subtraction. x and D * Div<D>(x)
// are guaranteed to differ by less than D, so their highest 32 bits
// must be identical, so we can truncate both sides to uint32_t before
// subtracting.
// We can also simplify (uint32_t) (D * Div<D>(x)).
// We can truncate before multiplying instead of after, as multiplying
// the highest 32 bits of Div<D>(x) can't affect the lowest 32 bits.

RYU_INLINE uint32_t Mod5(uint64_t x, uint64_t q) {
    return static_cast<uint32_t>(x) - 5 * static_cast<uint32_t>(q);
}

RYU_INLINE uint32_t Mod10(uint64_t x, uint64_t q) {
    return static_cast<uint32_t>(x) - 10 * static_cast<uint32_t>(q);
}

RYU_INLINE uint32_t Mod100(uint64_t x, uint64_t q) {
    return static_cast<uint32_t>(x) - 100 * static_cast<uint32_t>(q);
}

RYU_INLINE uint32_t Mod1e4(uint64_t x, uint64_t q) {
    return static_cast<uint32_t>(x) - 10000 * static_cast<uint32_t>(q);
}

RYU_INLINE uint32_t Mod1e8(uint64_t x, uint64_t q) {
    return static_cast<uint32_t>(x) - 100000000 * static_cast<uint32_t>(q);
}

#else // RYU_32_BIT_PLATFORM

RYU_INLINE uint64_t Div5(uint64_t x) {
    return x / 5;
}

RYU_INLINE uint64_t Div10(uint64_t x) {
    return x / 10;
}

RYU_INLINE uint64_t Div100(uint64_t x) {
    return x / 100;
}

RYU_INLINE uint64_t Div1e4(uint64_t x) {
    return x / 10000;
}

RYU_INLINE uint64_t Div1e8(uint64_t x) {
    return x / 100000000;
}

RYU_INLINE uint32_t Mod5(uint64_t x, uint64_t q) {
    return static_cast<uint32_t>(x - 5 * q);
}

RYU_INLINE uint32_t Mod10(uint64_t x, uint64_t q) {
    return static_cast<uint32_t>(x - 10 * q);
}

RYU_INLINE uint32_t Mod100(uint64_t x, uint64_t q) {
    return static_cast<uint32_t>(x - 100 * q);
}

RYU_INLINE uint32_t Mod1e4(uint64_t x, uint64_t q) {
    return static_cast<uint32_t>(x - 10000 * q);
}

RYU_INLINE uint32_t Mod1e8(uint64_t x, uint64_t q) {
    return static_cast<uint32_t>(x - 100000000 * q);
}

#endif // RYU_32_BIT_PLATFORM

RYU_INLINE int Pow5Factor(uint64_t value)
{
    // For 64-bit integers: result <= 27
    // Since value here has at most 55-bits: result <= 23

    int factor = 0;
    for (;;) {
        RYU_ASSERT(value != 0);
        RYU_ASSERT(factor <= 23);

        const uint64_t q = Div5(value);
        const uint32_t r = Mod5(value, q);
        if (r != 0)
            return factor;
        value = q;
        ++factor;
    }
}

RYU_INLINE bool MultipleOfPow5(uint64_t value, int p)
{
    return Pow5Factor(value) >= p;
}

RYU_INLINE bool MultipleOfPow2(uint64_t value, int p)
{
    RYU_ASSERT(p >= 0);
    RYU_ASSERT(p <= 63);

    //return (value << (64 - p)) == 0;
    return (value & ((uint64_t{1} << p) - 1)) == 0;
}

} // namespace impl

struct F64ToDecimalResult {
    uint64_t digits;
    int exponent;
};

RYU_INLINE F64ToDecimalResult ToDecimal(double value)
{
    using namespace ryu::impl;

    RYU_ASSERT(Double(value).IsFinite());
    RYU_ASSERT(value > 0);

    //
    // Step 1:
    // Decode the floating point number, and unify normalized and subnormal cases.
    //

    const Double ieeeValue(value);

    // Decode bits into mantissa, and exponent.
    const uint64_t ieeeMantissa = ieeeValue.PhysicalSignificand();
    const uint64_t ieeeExponent = ieeeValue.PhysicalExponent();

    uint64_t m2;
    int e2;
    if (ieeeExponent == 0) {
        m2 = ieeeMantissa;
        e2 = 1;
    } else {
        m2 = Double::HiddenBit | ieeeMantissa;
        e2 = static_cast<int>(ieeeExponent);
    }

    const bool even = (m2 & 1) == 0;
    const bool acceptBounds = even;

    //
    // Step 2:
    // Determine the interval of legal decimal representations.
    //

    // We subtract 2 so that the bounds computation has 2 additional bits.
    e2 -= Double::ExponentBias + 2;

    const uint64_t mv = 4 * m2;
    const uint64_t mp = mv + 2;
    const uint32_t mmShift = (ieeeMantissa != 0 || ieeeExponent <= 1) ? 1 : 0;
    const uint64_t mm = mv - 1 - mmShift;

    //
    // Step 3:
    // Convert to a decimal power base using 128-bit arithmetic.
    //

    int e10;

    uint64_t vm;
    uint64_t vr;
    uint64_t vp;

    bool vmIsTrailingZeros = false;
    bool vrIsTrailingZeros = false;
//  bool vpIsTrailingZeros = false;

    if (e2 >= 0)
    {
        // TODO:
        // Do the math and simplify!?

        // I tried special-casing q == 0, but there was no effect on performance.
        // q = max(0, log_10(2^e2) - 1)
        const int q = FloorLog10Pow2(e2) - (e2 > 3); // exponent <= 0
        RYU_ASSERT(q >= 0);
        const int k = FloorLog2Pow5(-q) + 1 - 128;
        const int j = -e2 + q - k; // shift
        RYU_ASSERT(j >= 115);

        e10 = q;

        // mul = ceil(2^-k / 5^q)
        const auto mul = ComputePow5Double(-q);
        MulShiftAll(mv, mp, mm, &mul, j, vr, vp, vm);

        // 22 = floor(log_5(2^53))
        // 23 = floor(log_5(2^(53+2)))
        if (q <= 22)
        {
            // This should use q <= 22, but I think 21 is also safe. Smaller values
            // may still be safe, but it's more difficult to reason about them.
            // Only one of mp, mv, and mm can be a multiple of 5, if any.
            if (Mod5(mv, Div5(mv)) == 0)
            {
                vrIsTrailingZeros = MultipleOfPow5(mv, q);
            }
            else if (acceptBounds)
            {
                // Same as min(e2 + (~mm & 1), Pow5Factor(mm)) >= q
                // <=> e2 + (~mm & 1) >= q && Pow5Factor(mm) >= q
                // <=> true && Pow5Factor(mm) >= q, since e2 >= q.
                vmIsTrailingZeros = MultipleOfPow5(mm, q);
            }
            else
            {
                // Same as min(e2 + 1, Pow5Factor(mp)) >= q.
//              vpIsTrailingZeros = MultipleOfPow5(mp, q);
                vp -= MultipleOfPow5(mp, q);
            }
        }
    }
    else
    {
        // TODO:
        // Do the math and simplify!?

        // q = max(0, log_10(5^-e2) - 1)
        const int q = FloorLog10Pow5(-e2) - (-e2 > 1);
        RYU_ASSERT(q >= 0);
        const int i = -e2 - q; // -exponent > 0
        RYU_ASSERT(i > 0);
        const int k = FloorLog2Pow5(i) + 1 - 128;
        const int j = q - k; // shift
        RYU_ASSERT(j >= 114);

        e10 = -i;

        // mul = floor(5^i / 2^-k)
        const auto mul = ComputePow5Double(i);
        MulShiftAll(mv, mp, mm, &mul, j, vr, vp, vm);

        if (q <= 1)
        {
            // {vr,vp,vm} is trailing zeros if {mv,mp,mm} has at least q trailing 0 bits.
            // mv = 4 * m2, so it always has at least two trailing 0 bits.
            vrIsTrailingZeros = true;

            if (acceptBounds)
            {
                // mm = mv - 1 - mmShift, so it has 1 trailing 0 bit iff mmShift == 1.
                vmIsTrailingZeros = (mmShift == 1);
            }
            else
            {
                // mp = mv + 2, so it always has at least one trailing 0 bit.
//              vpIsTrailingZeros = true;
                --vp;
            }
        }
        else if (q <= Double::SignificandSize + 2)
        {
            // TODO(ulfjack): Use a tighter bound here.

            // We need to compute min(ntz(mv), Pow5Factor(mv) - e2) >= q-1
            // <=> ntz(mv) >= q-1  &&  Pow5Factor(mv) - e2 >= q-1
            // <=> ntz(mv) >= q-1
            // <=> mv & ((1 << (q-1)) - 1) == 0
            // We also need to make sure that the left shift does not overflow.
            vrIsTrailingZeros = MultipleOfPow2(mv, q - 1);
        }
    }

    //
    // Step 4:
    // Find the shortest decimal representation in the interval of legal representations.
    //

//  vp -= vpIsTrailingZeros;

    uint64_t output;
    if (vmIsTrailingZeros || vrIsTrailingZeros)
    {
        // General case, which happens rarely (<1%).

        uint32_t lastRemovedDigit = 0;

        bool vrPrevIsTrailingZeros = vrIsTrailingZeros;

        for (;;)
        {
            const uint64_t vmDiv10 = Div10(vm);
            const uint64_t vpDiv10 = Div10(vp);
            if (vmDiv10 >= vpDiv10)
                break;

            const uint32_t vmMod10 = Mod10(vm, vmDiv10);
            vmIsTrailingZeros &= (vmMod10 == 0);
            vrPrevIsTrailingZeros &= (lastRemovedDigit == 0);

            const uint64_t vrDiv10 = Div10(vr);
            const uint32_t vrMod10 = Mod10(vr, vrDiv10);
            lastRemovedDigit = vrMod10;

            vm = vmDiv10;
            vr = vrDiv10;
            vp = vpDiv10;
            ++e10;
        }

        if (vmIsTrailingZeros)
        {
            for (;;)
            {
                const uint64_t vmDiv10 = Div10(vm);
                const uint32_t vmMod10 = Mod10(vm, vmDiv10);
                if (vmMod10 != 0)
                    break;

                vrPrevIsTrailingZeros &= (lastRemovedDigit == 0);

                const uint64_t vrDiv10 = Div10(vr);
                const uint32_t vrMod10 = Mod10(vr, vrDiv10);
                lastRemovedDigit = vrMod10;

                vm = vmDiv10;
                vr = vrDiv10;
                //vp = Div10(vp);
                ++e10;
            }
        }

        bool roundUp = (lastRemovedDigit >= 5);
        if (lastRemovedDigit == 5 && vrPrevIsTrailingZeros)
        {
            // Halfway case: The number ends in ...500...00.
            roundUp = (static_cast<uint32_t>(vr) % 2 != 0);
        }

        // We need to take vr+1 if vr is outside bounds...
        // or we need to round up.
        const bool inc = (vr == vm && !(acceptBounds && vmIsTrailingZeros)) || roundUp;

        output = vr + (inc ? 1 : 0);
    }
    else
    {
        // Specialized for the common case (>99%).

        bool roundUp = false;

        // Remove 4 digits in each iteration.
        // This loop runs at most 20/4 = 5 times.
        for (;;)
        {
            const uint64_t vmDiv1e4 = Div1e4(vm);
            const uint64_t vpDiv1e4 = Div1e4(vp);
            if (vmDiv1e4 >= vpDiv1e4)
                break;

            const uint64_t vrDiv1e4 = Div1e4(vr);
            const uint32_t vrMod1e4 = Mod1e4(vr, vrDiv1e4);
            roundUp = (vrMod1e4 >= 10000 / 2);

            vm = vmDiv1e4;
            vr = vrDiv1e4;
            vp = vpDiv1e4;
            e10 += 4;
        }

        for (;;)
        {
            const uint64_t vmDiv10 = Div10(vm);
            const uint64_t vpDiv10 = Div10(vp);
            if (vmDiv10 >= vpDiv10)
                break;

            const uint64_t vrDiv10 = Div10(vr);
            const uint32_t vrMod10 = Mod10(vr, vrDiv10);
            roundUp = (vrMod10 >= 10 / 2);

            vm = vmDiv10;
            vr = vrDiv10;
            vp = vpDiv10;
            ++e10;
        }

        // We need to take vr+1 if vr is outside bounds...
        // or we need to round up.
        const bool inc = vr == vm || roundUp;

        output = vr + (inc ? 1 : 0);
    }

    return {output, e10};
}

//==================================================================================================
// ToDecimal
//
// Single-precision implementation
//==================================================================================================
// Constant data: 624 bytes

namespace impl {

RYU_INLINE uint64_t ComputePow5Single(int k)
{
#if 0
    // May use this if the double-precision table is available.
    // Works for all (required) decimal exponents.
    const auto p = ComputePow5Double(k);
    return k < 0
        ? p.hi + (p.lo != 0)
        : p.hi;
#else
    // Let e = FloorLog2Pow5(k) + 1 - 64
    // For k >= 0, stores 5^k in the form: floor( 5^k / 2^e )
    // For k <= 0, stores 5^k in the form:  ceil(2^-e / 5^-k)
    static constexpr int MinDecExp = -30;
    static constexpr int MaxDecExp =  47;
    static constexpr uint64_t Pow5[MaxDecExp - MinDecExp + 1] = {
        0xA2425FF75E14FC32, // e =  -133, k =  -30
        0xCAD2F7F5359A3B3F, // e =  -131, k =  -29
        0xFD87B5F28300CA0E, // e =  -129, k =  -28
        0x9E74D1B791E07E49, // e =  -126, k =  -27
        0xC612062576589DDB, // e =  -124, k =  -26
        0xF79687AED3EEC552, // e =  -122, k =  -25
        0x9ABE14CD44753B53, // e =  -119, k =  -24
        0xC16D9A0095928A28, // e =  -117, k =  -23
        0xF1C90080BAF72CB2, // e =  -115, k =  -22
        0x971DA05074DA7BEF, // e =  -112, k =  -21
        0xBCE5086492111AEB, // e =  -110, k =  -20
        0xEC1E4A7DB69561A6, // e =  -108, k =  -19
        0x9392EE8E921D5D08, // e =  -105, k =  -18
        0xB877AA3236A4B44A, // e =  -103, k =  -17
        0xE69594BEC44DE15C, // e =  -101, k =  -16
        0x901D7CF73AB0ACDA, // e =   -98, k =  -15
        0xB424DC35095CD810, // e =   -96, k =  -14
        0xE12E13424BB40E14, // e =   -94, k =  -13
        0x8CBCCC096F5088CC, // e =   -91, k =  -12
        0xAFEBFF0BCB24AAFF, // e =   -89, k =  -11
        0xDBE6FECEBDEDD5BF, // e =   -87, k =  -10
        0x89705F4136B4A598, // e =   -84, k =   -9
        0xABCC77118461CEFD, // e =   -82, k =   -8
        0xD6BF94D5E57A42BD, // e =   -80, k =   -7
        0x8637BD05AF6C69B6, // e =   -77, k =   -6
        0xA7C5AC471B478424, // e =   -75, k =   -5
        0xD1B71758E219652C, // e =   -73, k =   -4
        0x83126E978D4FDF3C, // e =   -70, k =   -3
        0xA3D70A3D70A3D70B, // e =   -68, k =   -2
        0xCCCCCCCCCCCCCCCD, // e =   -66, k =   -1
        0x8000000000000000, // e =   -63, k =    0
        0xA000000000000000, // e =   -61, k =    1
        0xC800000000000000, // e =   -59, k =    2
        0xFA00000000000000, // e =   -57, k =    3
        0x9C40000000000000, // e =   -54, k =    4
        0xC350000000000000, // e =   -52, k =    5
        0xF424000000000000, // e =   -50, k =    6
        0x9896800000000000, // e =   -47, k =    7
        0xBEBC200000000000, // e =   -45, k =    8
        0xEE6B280000000000, // e =   -43, k =    9
        0x9502F90000000000, // e =   -40, k =   10
        0xBA43B74000000000, // e =   -38, k =   11
        0xE8D4A51000000000, // e =   -36, k =   12
        0x9184E72A00000000, // e =   -33, k =   13
        0xB5E620F480000000, // e =   -31, k =   14
        0xE35FA931A0000000, // e =   -29, k =   15
        0x8E1BC9BF04000000, // e =   -26, k =   16
        0xB1A2BC2EC5000000, // e =   -24, k =   17
        0xDE0B6B3A76400000, // e =   -22, k =   18
        0x8AC7230489E80000, // e =   -19, k =   19
        0xAD78EBC5AC620000, // e =   -17, k =   20
        0xD8D726B7177A8000, // e =   -15, k =   21
        0x878678326EAC9000, // e =   -12, k =   22
        0xA968163F0A57B400, // e =   -10, k =   23
        0xD3C21BCECCEDA100, // e =    -8, k =   24
        0x84595161401484A0, // e =    -5, k =   25
        0xA56FA5B99019A5C8, // e =    -3, k =   26
        0xCECB8F27F4200F3A, // e =    -1, k =   27
        0x813F3978F8940984, // e =     2, k =   28
        0xA18F07D736B90BE5, // e =     4, k =   29
        0xC9F2C9CD04674EDE, // e =     6, k =   30
        0xFC6F7C4045812296, // e =     8, k =   31
        0x9DC5ADA82B70B59D, // e =    11, k =   32
        0xC5371912364CE305, // e =    13, k =   33
        0xF684DF56C3E01BC6, // e =    15, k =   34
        0x9A130B963A6C115C, // e =    18, k =   35
        0xC097CE7BC90715B3, // e =    20, k =   36
        0xF0BDC21ABB48DB20, // e =    22, k =   37
        0x96769950B50D88F4, // e =    25, k =   38
        0xBC143FA4E250EB31, // e =    27, k =   39
        0xEB194F8E1AE525FD, // e =    29, k =   40
        0x92EFD1B8D0CF37BE, // e =    32, k =   41
        0xB7ABC627050305AD, // e =    34, k =   42
        0xE596B7B0C643C719, // e =    36, k =   43
        0x8F7E32CE7BEA5C6F, // e =    39, k =   44
        0xB35DBF821AE4F38B, // e =    41, k =   45
        0xE0352F62A19E306E, // e =    43, k =   46
        0x8C213D9DA502DE45, // e =    46, k =   47
    };

    RYU_ASSERT(k >= MinDecExp);
    RYU_ASSERT(k <= MaxDecExp);
    return Pow5[k - MinDecExp];
#endif
}

RYU_INLINE uint32_t MulShift(uint32_t m, uint64_t mul, int j)
{
    RYU_ASSERT(j >= 59);
    RYU_ASSERT(j <= 63);

    const uint64_t bits0 = uint64_t{m} * Lo32(mul);
    const uint64_t bits1 = uint64_t{m} * Hi32(mul);

#if RYU_32_BIT_PLATFORM
    // On 32-bit platforms we can avoid a 64-bit shift-right since we only
    // need the upper 32 bits of the result and the shift value is > 32.
    const uint32_t bits0Hi = Hi32(bits0);
    uint32_t bits1Lo = Lo32(bits1);
    uint32_t bits1Hi = Hi32(bits1);
    bits1Lo += bits0Hi;
    bits1Hi += bits1Lo < bits0Hi;
    const int lshift = -j & 31; // == (32 - (j - 32)) % 32
    const int rshift =  j & 31; // == (     (j - 32)) % 32
    return (bits1Hi << lshift) | (bits1Lo >> rshift);
#else
    const uint64_t sum = bits1 + Hi32(bits0);
    const uint64_t shiftedSum = sum >> (j - 32);
    RYU_ASSERT(shiftedSum <= UINT32_MAX);
    return static_cast<uint32_t>(shiftedSum);
#endif
}

RYU_INLINE int Pow5Factor(uint32_t value)
{
    int factor = 0;
    for (;;) {
        RYU_ASSERT(value != 0);
        RYU_ASSERT(factor <= 13);

        if (value % 5 != 0) {
            return factor;
        }
        value /= 5;
        ++factor;
    }
}

RYU_INLINE bool MultipleOfPow5(uint32_t value, int p)
{
    return Pow5Factor(value) >= p;
}

RYU_INLINE bool MultipleOfPow2(uint32_t value, int p)
{
    RYU_ASSERT(p >= 0);
    RYU_ASSERT(p <= 31);

//  return (value << (32 - p)) == 0;
    return (value & ((uint32_t{1} << p) - 1)) == 0;
}

} // namespace impl

struct F32ToDecimalResult {
    uint32_t digits;
    int exponent;
};

RYU_INLINE F32ToDecimalResult ToDecimal(float value)
{
    using namespace ryu::impl;

    RYU_ASSERT(Single(value).IsFinite());
    RYU_ASSERT(value > 0);

    //
    // Step 1:
    // Decode the floating point number, and unify normalized and subnormal cases.
    //

    const Single ieee_value(value);

    // Decode bits into mantissa, and exponent.
    const uint32_t ieeeMantissa = ieee_value.PhysicalSignificand();
    const uint32_t ieeeExponent = ieee_value.PhysicalExponent();

    uint32_t m2;
    int e2;
    if (ieeeExponent == 0) {
        m2 = ieeeMantissa;
        e2 = 1;
    } else {
        m2 = Single::HiddenBit | ieeeMantissa;
        e2 = static_cast<int>(ieeeExponent);
    }

    const bool even = (m2 & 1) == 0;
    const bool acceptBounds = even;

    //
    // Step 2:
    // Determine the interval of legal decimal representations.
    //

    // We subtract 2 so that the bounds computation has 2 additional bits.
    e2 -= Single::ExponentBias + 2;

    const uint32_t mv = 4 * m2;
    const uint32_t mp = mv + 2;
    const uint32_t mmShift = (ieeeMantissa != 0 || ieeeExponent <= 1) ? 1 : 0;
    const uint32_t mm = mv - 1 - mmShift;

    //
    // Step 3:
    // Convert to a decimal power base using 128-bit arithmetic.
    //

    int e10;

    uint32_t vm;
    uint32_t vr;
    uint32_t vp;

    bool vmIsTrailingZeros = false;
    bool vrIsTrailingZeros = false;
//  bool vpIsTrailingZeros = false;

    uint32_t lastRemovedDigit = 0;

    if (e2 >= 0)
    {
        // TODO:
        // Do the math and simplify!?

        const int q = FloorLog10Pow2(e2);
        RYU_ASSERT(q >= 0);
        const int k = FloorLog2Pow5(-q) + 1 - 64;
        const int j = -e2 + q - k; // shift

        e10 = q;

        const uint64_t mul = ComputePow5Single(-q);
        vr = MulShift(mv, mul, j);
        vp = MulShift(mp, mul, j);
        vm = MulShift(mm, mul, j);

        if (q != 0 && (vp - 1) / 10 <= vm / 10)
        {
            // We need to know one removed digit even if we are not going to loop below. We could use
            // q = X - 1 above, except that would require 33 bits for the result, and we've found that
            // 32-bit arithmetic is faster even on 64-bit machines.

            const int q1 = q - 1;
            RYU_ASSERT(q1 >= 0);
            const int k1 = FloorLog2Pow5(-q1) + 1 - 64;
            const int j1 = -e2 + q1 - k1; // shift

            const uint64_t mul1 = ComputePow5Single(-q1);
            lastRemovedDigit = MulShift(mv, mul1, j1) % 10;
        }

        // 10 = floor(log_5(2^24))
        // 11 = floor(log_5(2^(24+2)))
        if (q <= 10)
        {
            // The largest power of 5 that fits in 24 bits is 5^10, but q <= 9 seems to be safe as well.
            // Only one of mp, mv, and mm can be a multiple of 5, if any.
            if (mv % 5 == 0)
            {
                vrIsTrailingZeros = MultipleOfPow5(mv, q);
            }
            else if (acceptBounds)
            {
                // Same as min(e2 + (~mm & 1), Pow5Factor(mm)) >= q
                // <=> e2 + (~mm & 1) >= q && Pow5Factor(mm) >= q
                // <=> true && Pow5Factor(mm) >= q, since e2 >= q.
                vmIsTrailingZeros = MultipleOfPow5(mm, q);
            }
            else
            {
                // Same as min(e2 + 1, Pow5Factor(mp)) >= q.
//              vpIsTrailingZeros = MultipleOfPow5(mp, q);
                vp -= MultipleOfPow5(mp, q);
            }
        }
    }
    else
    {
        // TODO:
        // Do the math and simplify!?

        const int q = FloorLog10Pow5(-e2);
        RYU_ASSERT(q >= 0);
        const int i = -e2 - q;
        RYU_ASSERT(i >= 0);
        const int k = FloorLog2Pow5(i) + 1 - 64;
        const int j = q - k; // shift

        e10 = q + e2;

        const uint64_t mul = ComputePow5Single(i);
        vr = MulShift(mv, mul, j);
        vp = MulShift(mp, mul, j);
        vm = MulShift(mm, mul, j);

        if (q != 0 && (vp - 1) / 10 <= vm / 10)
        {
            // We need to know one removed digit even if we are not going to loop below. We could use
            // q = X - 1 above, except that would require 33 bits for the result, and we've found that
            // 32-bit arithmetic is faster even on 64-bit machines.

            const int q1 = q - 1;
            RYU_ASSERT(q1 >= 0);
            const int i1 = i + 1; // = -e2 - q1
            RYU_ASSERT(i1 >= 0);
            const int k1 = FloorLog2Pow5(i1) + 1 - 64;
            const int j1 = q1 - k1; // shift

            const uint64_t mul1 = ComputePow5Single(i1);
            lastRemovedDigit = MulShift(mv, mul1, j1) % 10;
        }

        if (q <= 1)
        {
            // {vr,vp,vm} is trailing zeros if {mv,mp,mm} has at least q trailing 0 bits.
            // mv = 4 * m2, so it always has at least two trailing 0 bits.
            vrIsTrailingZeros = true;

            if (acceptBounds)
            {
                // mm = mv - 1 - mmShift, so it has 1 trailing 0 bit iff mmShift == 1.
                vmIsTrailingZeros = mmShift == 1;
            }
            else
            {
                // mp = mv + 2, so it always has at least one trailing 0 bit.
//              vpIsTrailingZeros = true;
                vp--;
            }
        }
        else if (q <= Single::SignificandSize + 2)
        {
            // We need to compute min(ntz(mv), Pow5Factor(mv) - e2) >= q-1
            // <=> ntz(mv) >= q-1  &&  Pow5Factor(mv) - e2 >= q-1
            // <=> ntz(mv) >= q-1
            // <=> mv & ((1 << (q-1)) - 1) == 0
            // We also need to make sure that the left shift does not overflow.
            vrIsTrailingZeros = MultipleOfPow2(mv, q - 1);
        }
    }

    //
    // Step 4:
    // Find the shortest decimal representation in the interval of legal representations.
    //

//  vp -= vpIsTrailingZeros;

    uint32_t output;

    if (vmIsTrailingZeros || vrIsTrailingZeros)
    {
        // General case, which happens rarely (~4.0%).

        bool vrPrevIsTrailingZeros = vrIsTrailingZeros;

        while (vm / 10 < vp / 10)
        {
            vmIsTrailingZeros &= (vm % 10 == 0);
            vrPrevIsTrailingZeros &= (lastRemovedDigit == 0);

            lastRemovedDigit = vr % 10;

            vm /= 10;
            vr /= 10;
            vp /= 10;
            ++e10;
        }

        if (vmIsTrailingZeros)
        {
            while (vm % 10 == 0)
            {
                vrPrevIsTrailingZeros &= (lastRemovedDigit == 0);

                lastRemovedDigit = vr % 10;

                vm /= 10;
                vr /= 10;
                //vp /= 10;
                ++e10;
            }
        }

        bool roundUp = (lastRemovedDigit >= 5);
        if (lastRemovedDigit == 5 && vrPrevIsTrailingZeros)
        {
            // Halfway case: The number ends in ...500...00.
            roundUp = vr % 2 != 0;
        }

        // We need to take vr+1 if vr is outside bounds...
        // or we need to round up.
        const bool inc = (vr == vm && !(acceptBounds && vmIsTrailingZeros)) || roundUp;

        output = vr + (inc ? 1 : 0);
    }
    else
    {
        // Specialized for the common case (~96.0%).

        while (vm / 10 < vp / 10)
        {
            lastRemovedDigit = vr % 10;
            vm /= 10;
            vr /= 10;
            vp /= 10;
            ++e10;
        }

        // We need to take vr+1 if vr is outside bounds...
        // or we need to round up.
        const bool inc = (vr == vm || lastRemovedDigit >= 5);

        output = vr + (inc ? 1 : 0);
    }

    return {output, e10};
}

//==================================================================================================
// ToDigits
//==================================================================================================
// Constant data: 200 bytes

namespace impl {

RYU_INLINE char* Utoa_2Digits(char* buf, uint32_t digits)
{
    static constexpr char kDigits100[200] = {
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

    RYU_ASSERT(digits < 100);
    std::memcpy(buf, &kDigits100[2*digits], 2*sizeof(char));
    return buf + 2;
}

RYU_INLINE char* Utoa_4Digits(char* buf, uint32_t digits)
{
    RYU_ASSERT(digits < 10000);
    const uint32_t q = digits / 100;
    const uint32_t r = digits % 100;
    Utoa_2Digits(buf + 0, q);
    Utoa_2Digits(buf + 2, r);
    return buf + 4;
}

RYU_INLINE char* Utoa_8Digits(char* buf, uint32_t digits)
{
    RYU_ASSERT(digits < 100000000);
    const uint32_t q = digits / 10000;
    const uint32_t r = digits % 10000;
    Utoa_4Digits(buf + 0, q);
    Utoa_4Digits(buf + 4, r);
    return buf + 8;
}

RYU_INLINE int DecimalLength(uint64_t v)
{
    RYU_ASSERT(v >= 1);
    RYU_ASSERT(v < 100000000000000000ull); // 10^17 = 0x0163'4578'5D8A'0000

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

RYU_INLINE int PrintDecimalDigits(char* buf, uint64_t output)
{
    const int output_length = DecimalLength(output);
    int i = output_length;

    // We prefer 32-bit operations, even on 64-bit platforms.
    // We have at most 17 digits, and uint32_t can store 9 digits.
    // If output doesn't fit into uint32_t, we cut off 8 digits,
    // so the rest will fit into uint32_t.
    if (static_cast<uint32_t>(output >> 32) != 0)
    {
        RYU_ASSERT(i > 8);
        const uint64_t q = Div1e8(output);
        const uint32_t r = Mod1e8(output, q);
        output = q;
        i -= 8;
        Utoa_8Digits(buf + i, r);
    }

    RYU_ASSERT(output <= UINT32_MAX);
    uint32_t output2 = static_cast<uint32_t>(output);

    while (output2 >= 10000)
    {
        RYU_ASSERT(i > 4);
        const uint32_t q = output2 / 10000;
        const uint32_t r = output2 % 10000;
        output2 = q;
        i -= 4;
        Utoa_4Digits(buf + i, r);
    }

    if (output2 >= 100)
    {
        RYU_ASSERT(i > 2);
        const uint32_t q = output2 / 100;
        const uint32_t r = output2 % 100;
        output2 = q;
        i -= 2;
        Utoa_2Digits(buf + i, r);
    }

    if (output2 >= 10)
    {
        RYU_ASSERT(i == 2);
        Utoa_2Digits(buf, output2);
    }
    else
    {
        RYU_ASSERT(i == 1);
        buf[0] = static_cast<char>('0' + output2);
    }

    return output_length;
}

RYU_INLINE int DecimalLength(uint32_t v)
{
    RYU_ASSERT(v >= 1);
    RYU_ASSERT(v < 1000000000);

    if (v >= 100000000) { return 9; }
    if (v >= 10000000) { return 8; }
    if (v >= 1000000) { return 7; }
    if (v >= 100000) { return 6; }
    if (v >= 10000) { return 5; }
    if (v >= 1000) { return 4; }
    if (v >= 100) { return 3; }
    if (v >= 10) { return 2; }
    return 1;
}

RYU_INLINE int PrintDecimalDigits(char* buf, uint32_t output)
{
    const int output_length = DecimalLength(output);
    int i = output_length;

    while (output >= 10000)
    {
        RYU_ASSERT(i > 4);
        const uint32_t q = output / 10000;
        const uint32_t r = output % 10000;
        output = q;
        i -= 4;
        Utoa_4Digits(buf + i, r);
    }

    if (output >= 100)
    {
        RYU_ASSERT(i > 2);
        const uint32_t q = output / 100;
        const uint32_t r = output % 100;
        output = q;
        i -= 2;
        Utoa_2Digits(buf + i, r);
    }

    if (output >= 10)
    {
        RYU_ASSERT(i == 2);
        Utoa_2Digits(buf, output);
    }
    else
    {
        RYU_ASSERT(i == 1);
        buf[0] = static_cast<char>('0' + output);
    }

    return output_length;
}

} // namespace impl

template <typename Float>
RYU_INLINE void ToDigits(char* buffer, int& num_digits, int& exponent, Float value)
{
    const auto dec = ryu::ToDecimal(value);

    num_digits = ryu::impl::PrintDecimalDigits(buffer, dec.digits);
    exponent = dec.exponent;
}

//==================================================================================================
// ToChars
//==================================================================================================

namespace impl {

// Appends a decimal representation of 'value' to buffer.
// Returns a pointer to the element following the digits.
//
// PRE: -1000 < value < 1000
RYU_INLINE char* ExponentToString(char* buffer, int value)
{
    RYU_ASSERT(value > -1000);
    RYU_ASSERT(value <  1000);

    int n = 0;

    if (value < 0)
    {
        buffer[n++] = '-';
        value = -value;
    }
    else
    {
        buffer[n++] = '+';
    }

    const uint32_t k = static_cast<uint32_t>(value);
    if (k < 10)
    {
        buffer[n++] = static_cast<char>('0' + k);
    }
    else if (k < 100)
    {
        Utoa_2Digits(buffer + n, k);
        n += 2;
    }
    else
    {
        const uint32_t r = k % 10;
        const uint32_t q = k / 10;
        Utoa_2Digits(buffer + n, q);
        n += 2;
        buffer[n++] = static_cast<char>('0' + r);
    }

    return buffer + n;
}

RYU_INLINE char* FormatFixed(char* buffer, intptr_t num_digits, intptr_t decimal_point, bool force_trailing_dot_zero)
{
    RYU_ASSERT(buffer != nullptr);
    RYU_ASSERT(num_digits >= 1);

    if (num_digits <= decimal_point)
    {
        // digits[000]
        // GRISU_ASSERT(buffer_capacity >= decimal_point + (force_trailing_dot_zero ? 2 : 0));

        std::memset(buffer + num_digits, '0', static_cast<size_t>(decimal_point - num_digits));
        buffer += decimal_point;
        if (force_trailing_dot_zero)
        {
            *buffer++ = '.';
            *buffer++ = '0';
        }
        return buffer;
    }
    else if (0 < decimal_point)
    {
        // dig.its
        // GRISU_ASSERT(buffer_capacity >= length + 1);

        std::memmove(buffer + (decimal_point + 1), buffer + decimal_point, static_cast<size_t>(num_digits - decimal_point));
        buffer[decimal_point] = '.';
        return buffer + (num_digits + 1);
    }
    else // decimal_point <= 0
    {
        // 0.[000]digits
        // GRISU_ASSERT(buffer_capacity >= 2 + (-decimal_point) + length);

        std::memmove(buffer + (2 + -decimal_point), buffer, static_cast<size_t>(num_digits));
        buffer[0] = '0';
        buffer[1] = '.';
        std::memset(buffer + 2, '0', static_cast<size_t>(-decimal_point));
        return buffer + (2 + (-decimal_point) + num_digits);
    }
}

RYU_INLINE char* FormatScientific(char* buffer, intptr_t num_digits, int exponent, bool force_trailing_dot_zero)
{
    RYU_ASSERT(buffer != nullptr);
    RYU_ASSERT(num_digits >= 1);

    if (num_digits == 1)
    {
        // dE+123
        // GRISU_ASSERT(buffer_capacity >= num_digits + 5);

        buffer += 1;
        if (force_trailing_dot_zero)
        {
            *buffer++ = '.';
            *buffer++ = '0';
        }
    }
    else
    {
        // d.igitsE+123
        // GRISU_ASSERT(buffer_capacity >= num_digits + 1 + 5);

        std::memmove(buffer + 2, buffer + 1, static_cast<size_t>(num_digits - 1));
        buffer[1] = '.';
        buffer += 1 + num_digits;
    }

    buffer[0] = 'e';
    buffer = ExponentToString(buffer + 1, exponent);

    return buffer;
}

// Format the digits similar to printf's %g style.
template <typename UnsignedInteger>
RYU_INLINE char* FormatDigits(char* buffer, UnsignedInteger digits, int exponent, bool force_trailing_dot_zero)
{
    const int num_digits = PrintDecimalDigits(buffer, digits);
    const int decimal_point = num_digits + exponent;

    // NB:
    // These are the values used by JavaScript's ToString applied to Number
    // type. Printf uses the values -4 and max_digits10 resp. (sort of).
    constexpr int kMinExp = -6;
    constexpr int kMaxExp = 21;

    const bool use_fixed = kMinExp < decimal_point && decimal_point <= kMaxExp;

    return use_fixed
        ? FormatFixed(buffer, num_digits, decimal_point, force_trailing_dot_zero)
        : FormatScientific(buffer, num_digits, decimal_point - 1, force_trailing_dot_zero);
}

} // namespace impl

// Generates a decimal representation of the floating-point number `value` in 'buffer'.
// Note: The result is _not_ null-terminated.
//
// PRE: The buffer must be large enough (32 bytes is sufficient).
template <typename Float>
RYU_INLINE char* ToChars(char* buffer, Float value, bool force_trailing_dot_zero = false)
{
    using Fp = ryu::impl::IEEE<Float>;
    const Fp v(value);

    const bool is_neg = v.SignBit();

    if (!v.IsFinite())
    {
        if (v.IsNaN())
        {
            std::memcpy(buffer, "NaN", 3);
            return buffer + 3;
        }
        if (is_neg)
        {
            *buffer++ = '-';
        }
        std::memcpy(buffer, "Infinity", 8);
        return buffer + 8;
    }

    if (v.SignBit())
    {
        value = v.AbsValue();
        *buffer++ = '-';
    }

    if (v.IsZero())
    {
        *buffer++ = '0';
        if (force_trailing_dot_zero)
        {
            *buffer++ = '.';
            *buffer++ = '0';
        }
        return buffer;
    }

    const auto dec = ryu::ToDecimal(value);
    return ryu::impl::FormatDigits(buffer, dec.digits, dec.exponent, force_trailing_dot_zero);
}

} // namespace ryu

//char* Dtoa(char* buffer, double value)
//{
//    return ryu::ToChars(buffer, value);
//}

//char* Ftoa(char* buffer, float value)
//{
//    return ryu::ToChars(buffer, value);
//}
