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

#define RYU_KEEP_TRAILING_ZEROS_IN_SMALL_INT()  1
#define RYU_SCIENTIFIC_NOTATION_ONLY()          0
#define RYU_USE_INTRINSICS()                    1

#include <cassert>
#include <climits>
#include <cstdint>
#include <cstring>
#include <limits>
#if RYU_USE_INTRINSICS() && defined(_MSC_VER)
#include <intrin.h>
#endif

#ifndef RYU_ASSERT
#define RYU_ASSERT(X) assert(X)
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

    bool IsFinite() const {
        return (bits & ExponentMask) != ExponentMask;
    }

//  bool IsInf() const {
//      return (bits & ExponentMask) == ExponentMask && (bits & SignificandMask) == 0;
//  }

    bool IsNaN() const {
        return (bits & ExponentMask) == ExponentMask && (bits & SignificandMask) != 0;
    }

    bool IsZero() const {
        return (bits & ~SignMask) == 0;
    }

    bool SignBit() const {
        return (bits & SignMask) != 0;
    }

//  value_type Value() const {
//      return ReinterpretBits<value_type>(bits);
//  }

    value_type AbsValue() const {
        return ReinterpretBits<value_type>(bits & ~SignMask);
    }
};

//==================================================================================================
//
//==================================================================================================

// Returns floor(x / 2^n).
static inline int FloorDivPow2(int x, int n)
{
    // Technically, right-shift of negative integers is implementation defined...
    // Should easily be optimized into SAR (or equivalent) instruction.
#if 1
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

static inline uint32_t Hi32(uint64_t x)
{
    return static_cast<uint32_t>(x >> 32);
}

//==================================================================================================
// ToDecimal
//
// Double-precision implementation
//==================================================================================================
// Constant data = 9872 (+ 368) bytes

static constexpr int BitsPerPow5_Double = 124;

struct Uint64x2 {
    uint64_t hi;
    uint64_t lo;
};

static inline Uint64x2 ComputePow5_Double(int k)
{
    // Let e = FloorLog2Pow5(k) + 1 - 124
    // For k >= 0, stores 5^k in the form: ceil( 5^k / 2^e )
    // For k <= 0, stores 5^k in the form: ceil(2^-e / 5^-k)
    static constexpr int MinDecExp = -290;
    static constexpr int MaxDecExp =  325;
    static constexpr Uint64x2 Pow5[MaxDecExp - MinDecExp + 1] = {
        {0x0C795830D75038C1, 0xDD59DF5B9EF6A242}, // e =  -797, k = -290
        {0x0F97AE3D0D2446F2, 0x54B0573286B44AD2}, // e =  -795, k = -289
        {0x09BECCE62836AC57, 0x74EE367F9430AEC4}, // e =  -792, k = -288
        {0x0C2E801FB244576D, 0x5229C41F793CDA74}, // e =  -790, k = -287
        {0x0F3A20279ED56D48, 0xA6B43527578C1111}, // e =  -788, k = -286
        {0x09845418C345644D, 0x6830A13896B78AAB}, // e =  -785, k = -285
        {0x0BE5691EF416BD60, 0xC23CC986BC656D56}, // e =  -783, k = -284
        {0x0EDEC366B11C6CB8, 0xF2CBFBE86B7EC8AB}, // e =  -781, k = -283
        {0x094B3A202EB1C3F3, 0x97BF7D71432F3D6B}, // e =  -778, k = -282
        {0x0B9E08A83A5E34F0, 0x7DAF5CCD93FB0CC6}, // e =  -776, k = -281
        {0x0E858AD248F5C22C, 0x9D1B3400F8F9CFF7}, // e =  -774, k = -280
        {0x091376C36D99995B, 0xE23100809B9C21FB}, // e =  -771, k = -279
        {0x0B58547448FFFFB2, 0xDABD40A0C2832A79}, // e =  -769, k = -278
        {0x0E2E69915B3FFF9F, 0x916C90C8F323F517}, // e =  -767, k = -277
        {0x08DD01FAD907FFC3, 0xBAE3DA7D97F6792F}, // e =  -764, k = -276
        {0x0B1442798F49FFB4, 0xA99CD11CFDF4177A}, // e =  -762, k = -275
        {0x0DD95317F31C7FA1, 0xD40405643D711D59}, // e =  -760, k = -274
        {0x08A7D3EEF7F1CFC5, 0x2482835EA666B258}, // e =  -757, k = -273
        {0x0AD1C8EAB5EE43B6, 0x6DA3243650005EED}, // e =  -755, k = -272
        {0x0D863B256369D4A4, 0x090BED43E40076A9}, // e =  -753, k = -271
        {0x0873E4F75E2224E6, 0x85A7744A6E804A2A}, // e =  -750, k = -270
        {0x0A90DE3535AAAE20, 0x2711515D0A205CB4}, // e =  -748, k = -269
        {0x0D3515C2831559A8, 0x30D5A5B44CA873E1}, // e =  -746, k = -268
        {0x08412D9991ED5809, 0x1E858790AFE9486D}, // e =  -743, k = -267
        {0x0A5178FFF668AE0B, 0x6626E974DBE39A88}, // e =  -741, k = -266
        {0x0CE5D73FF402D98E, 0x3FB0A3D212DC8129}, // e =  -739, k = -265
        {0x080FA687F881C7F8, 0xE7CE66634BC9D0BA}, // e =  -736, k = -264
        {0x0A139029F6A239F7, 0x21C1FFFC1EBC44E9}, // e =  -734, k = -263
        {0x0C987434744AC874, 0xEA327FFB266B5623}, // e =  -732, k = -262
        {0x0FBE9141915D7A92, 0x24BF1FF9F0062BAB}, // e =  -730, k = -261
        {0x09D71AC8FADA6C9B, 0x56F773FC3603DB4B}, // e =  -727, k = -260
        {0x0C4CE17B399107C2, 0x2CB550FB4384D21E}, // e =  -725, k = -259
        {0x0F6019DA07F549B2, 0xB7E2A53A146606A5}, // e =  -723, k = -258
        {0x099C102844F94E0F, 0xB2EDA7444CBFC427}, // e =  -720, k = -257
        {0x0C0314325637A193, 0x9FA911155FEFB531}, // e =  -718, k = -256
        {0x0F03D93EEBC589F8, 0x8793555AB7EBA27D}, // e =  -716, k = -255
        {0x096267C7535B763B, 0x54BC1558B2F3458E}, // e =  -713, k = -254
        {0x0BBB01B9283253CA, 0x29EB1AAEDFB016F2}, // e =  -711, k = -253
        {0x0EA9C227723EE8BC, 0xB465E15A979C1CAE}, // e =  -709, k = -252
        {0x092A1958A7675175, 0xF0BFACD89EC191ED}, // e =  -706, k = -251
        {0x0B749FAED14125D3, 0x6CEF980EC671F668}, // e =  -704, k = -250
        {0x0E51C79A85916F48, 0x482B7E12780E7402}, // e =  -702, k = -249
        {0x08F31CC0937AE58D, 0x2D1B2ECB8B090882}, // e =  -699, k = -248
        {0x0B2FE3F0B8599EF0, 0x7861FA7E6DCB4AA2}, // e =  -697, k = -247
        {0x0DFBDCECE67006AC, 0x967A791E093E1D4A}, // e =  -695, k = -246
        {0x08BD6A141006042B, 0xDE0C8BB2C5C6D24F}, // e =  -692, k = -245
        {0x0AECC49914078536, 0xD58FAE9F773886E2}, // e =  -690, k = -244
        {0x0DA7F5BF59096684, 0x8AF39A475506A89A}, // e =  -688, k = -243
        {0x0888F99797A5E012, 0xD6D8406C95242961}, // e =  -685, k = -242
        {0x0AAB37FD7D8F5817, 0x8C8E5087BA6D33B9}, // e =  -683, k = -241
        {0x0D5605FCDCF32E1D, 0x6FB1E4A9A90880A7}, // e =  -681, k = -240
        {0x0855C3BE0A17FCD2, 0x65CF2EEA09A55068}, // e =  -678, k = -239
        {0x0A6B34AD8C9DFC06, 0xFF42FAA48C0EA482}, // e =  -676, k = -238
        {0x0D0601D8EFC57B08, 0xBF13B94DAF124DA3}, // e =  -674, k = -237
        {0x0823C12795DB6CE5, 0x776C53D08D6B7086}, // e =  -671, k = -236
        {0x0A2CB1717B52481E, 0xD54768C4B0C64CA7}, // e =  -669, k = -235
        {0x0CB7DDCDDA26DA26, 0x8A9942F5DCF7DFD1}, // e =  -667, k = -234
        {0x0FE5D54150B090B0, 0x2D3F93B35435D7C5}, // e =  -665, k = -233
        {0x09EFA548D26E5A6E, 0x1C47BC5014A1A6DB}, // e =  -662, k = -232
        {0x0C6B8E9B0709F109, 0xA359AB6419CA1092}, // e =  -660, k = -231
        {0x0F867241C8CC6D4C, 0x0C30163D203C94B7}, // e =  -658, k = -230
        {0x09B407691D7FC44F, 0x879E0DE63425DCF2}, // e =  -655, k = -229
        {0x0C21094364DFB563, 0x6985915FC12F542F}, // e =  -653, k = -228
        {0x0F294B943E17A2BC, 0x43E6F5B7B17B293A}, // e =  -651, k = -227
        {0x0979CF3CA6CEC5B5, 0xAA705992CEECF9C5}, // e =  -648, k = -226
        {0x0BD8430BD0827723, 0x150C6FF782A83836}, // e =  -646, k = -225
        {0x0ECE53CEC4A314EB, 0xDA4F8BF563524643}, // e =  -644, k = -224
        {0x0940F4613AE5ED13, 0x6871B7795E136BEA}, // e =  -641, k = -223
        {0x0B913179899F6858, 0x428E2557B59846E4}, // e =  -639, k = -222
        {0x0E757DD7EC07426E, 0x5331AEADA2FE589D}, // e =  -637, k = -221
        {0x09096EA6F3848984, 0xF3FF0D2C85DEF763}, // e =  -634, k = -220
        {0x0B4BCA50B065ABE6, 0x30FED077A756B53B}, // e =  -632, k = -219
        {0x0E1EBCE4DC7F16DF, 0xBD3E8495912C628A}, // e =  -630, k = -218
        {0x08D3360F09CF6E4B, 0xD64712DD7ABBBD96}, // e =  -627, k = -217
        {0x0B080392CC4349DE, 0xCBD8D794D96AACFC}, // e =  -625, k = -216
        {0x0DCA04777F541C56, 0x7ECF0D7A0FC5583B}, // e =  -623, k = -215
        {0x089E42CAAF9491B6, 0x0F41686C49DB5725}, // e =  -620, k = -214
        {0x0AC5D37D5B79B623, 0x9311C2875C522CEE}, // e =  -618, k = -213
        {0x0D77485CB25823AC, 0x77D633293366B829}, // e =  -616, k = -212
        {0x086A8D39EF77164B, 0xCAE5DFF9C020331A}, // e =  -613, k = -211
        {0x0A8530886B54DBDE, 0xBD9F57F830283FE0}, // e =  -611, k = -210
        {0x0D267CAA862A12D6, 0x6D072DF63C324FD8}, // e =  -609, k = -209
        {0x08380DEA93DA4BC6, 0x04247CB9E59F71E7}, // e =  -606, k = -208
        {0x0A46116538D0DEB7, 0x852D9BE85F074E61}, // e =  -604, k = -207
        {0x0CD795BE87051665, 0x667902E276C921F9}, // e =  -602, k = -206
        {0x0806BD9714632DFF, 0x600BA1CD8A3DB53C}, // e =  -599, k = -205
        {0x0A086CFCD97BF97F, 0x380E8A40ECCD228B}, // e =  -597, k = -204
        {0x0C8A883C0FDAF7DF, 0x06122CD128006B2D}, // e =  -595, k = -203
        {0x0FAD2A4B13D1B5D6, 0xC796B805720085F9}, // e =  -593, k = -202
        {0x09CC3A6EEC6311A6, 0x3CBE3303674053BC}, // e =  -590, k = -201
        {0x0C3F490AA77BD60F, 0xCBEDBFC4411068AA}, // e =  -588, k = -200
        {0x0F4F1B4D515ACB93, 0xBEE92FB5515482D5}, // e =  -586, k = -199
        {0x0991711052D8BF3C, 0x5751BDD152D4D1C5}, // e =  -583, k = -198
        {0x0BF5CD54678EEF0B, 0x6D262D45A78A0636}, // e =  -581, k = -197
        {0x0EF340A98172AACE, 0x486FB897116C87C4}, // e =  -579, k = -196
        {0x09580869F0E7AAC0, 0xED45D35E6AE3D4DB}, // e =  -576, k = -195
        {0x0BAE0A846D219571, 0x28974836059CCA11}, // e =  -574, k = -194
        {0x0E998D258869FACD, 0x72BD1A438703FC95}, // e =  -572, k = -193
        {0x091FF83775423CC0, 0x67B6306A34627DDD}, // e =  -569, k = -192
        {0x0B67F6455292CBF0, 0x81A3BC84C17B1D55}, // e =  -567, k = -191
        {0x0E41F3D6A7377EEC, 0xA20CABA5F1D9E4AA}, // e =  -565, k = -190
        {0x08E938662882AF53, 0xE547EB47B7282EEA}, // e =  -562, k = -189
        {0x0B23867FB2A35B28, 0xDE99E619A4F23AA5}, // e =  -560, k = -188
        {0x0DEC681F9F4C31F3, 0x16405FA00E2EC94E}, // e =  -558, k = -187
        {0x08B3C113C38F9F37, 0xEDE83BC408DD3DD1}, // e =  -555, k = -186
        {0x0AE0B158B4738705, 0xE9624AB50B148D45}, // e =  -553, k = -185
        {0x0D98DDAEE19068C7, 0x63BADD624DD9B096}, // e =  -551, k = -184
        {0x087F8A8D4CFA417C, 0x9E54CA5D70A80E5E}, // e =  -548, k = -183
        {0x0A9F6D30A038D1DB, 0xC5E9FCF4CCD211F5}, // e =  -546, k = -182
        {0x0D47487CC8470652, 0xB7647C3200069672}, // e =  -544, k = -181
        {0x084C8D4DFD2C63F3, 0xB29ECD9F40041E08}, // e =  -541, k = -180
        {0x0A5FB0A17C777CF0, 0x9F4681071005258A}, // e =  -539, k = -179
        {0x0CF79CC9DB955C2C, 0xC7182148D4066EEC}, // e =  -537, k = -178
        {0x081AC1FE293D599B, 0xFC6F14CD84840554}, // e =  -534, k = -177
        {0x0A21727DB38CB002, 0xFB8ADA00E5A506A8}, // e =  -532, k = -176
        {0x0CA9CF1D206FDC03, 0xBA6D90811F0E4852}, // e =  -530, k = -175
        {0x0FD442E4688BD304, 0xA908F4A166D1DA67}, // e =  -528, k = -174
        {0x09E4A9CEC15763E2, 0xE9A598E4E0432880}, // e =  -525, k = -173
        {0x0C5DD44271AD3CDB, 0xA40EFF1E1853F2A0}, // e =  -523, k = -172
        {0x0F7549530E188C12, 0x8D12BEE59E68EF48}, // e =  -521, k = -171
        {0x09A94DD3E8CF578B, 0x982BB74F8301958D}, // e =  -518, k = -170
        {0x0C13A148E3032D6E, 0x7E36A52363C1FAF1}, // e =  -516, k = -169
        {0x0F18899B1BC3F8CA, 0x1DC44E6C3CB279AD}, // e =  -514, k = -168
        {0x096F5600F15A7B7E, 0x529AB103A5EF8C0C}, // e =  -511, k = -167
        {0x0BCB2B812DB11A5D, 0xE7415D448F6B6F0F}, // e =  -509, k = -166
        {0x0EBDF661791D60F5, 0x6111B495B3464AD3}, // e =  -507, k = -165
        {0x0936B9FCEBB25C99, 0x5CAB10DD900BEEC4}, // e =  -504, k = -164
        {0x0B84687C269EF3BF, 0xB3D5D514F40EEA75}, // e =  -502, k = -163
        {0x0E65829B3046B0AF, 0xA0CB4A5A3112A512}, // e =  -500, k = -162
        {0x08FF71A0FE2C2E6D, 0xC47F0E785EABA72B}, // e =  -497, k = -161
        {0x0B3F4E093DB73A09, 0x359ED216765690F6}, // e =  -495, k = -160
        {0x0E0F218B8D25088B, 0x8306869C13EC3533}, // e =  -493, k = -159
        {0x08C974F738372557, 0x31E414218C73A140}, // e =  -490, k = -158
        {0x0AFBD2350644EEAC, 0xFE5D1929EF908990}, // e =  -488, k = -157
        {0x0DBAC6C247D62A58, 0x3DF45F746B74ABF4}, // e =  -486, k = -156
        {0x0894BC396CE5DA77, 0x26B8BBA8C328EB79}, // e =  -483, k = -155
        {0x0AB9EB47C81F5114, 0xF066EA92F3F32657}, // e =  -481, k = -154
        {0x0D686619BA27255A, 0x2C80A537B0EFEFEC}, // e =  -479, k = -153
        {0x08613FD014587758, 0x5BD06742CE95F5F4}, // e =  -476, k = -152
        {0x0A798FC4196E952E, 0x72C48113823B7371}, // e =  -474, k = -151
        {0x0D17F3B51FCA3A7A, 0x0F75A15862CA504D}, // e =  -472, k = -150
        {0x082EF85133DE648C, 0x49A984D73DBE7230}, // e =  -469, k = -149
        {0x0A3AB66580D5FDAF, 0x5C13E60D0D2E0EBC}, // e =  -467, k = -148
        {0x0CC963FEE10B7D1B, 0x3318DF905079926B}, // e =  -465, k = -147
        {0x0FFBBCFE994E5C61, 0xFFDF17746497F706}, // e =  -463, k = -146
        {0x09FD561F1FD0F9BD, 0x3FEB6EA8BEDEFA64}, // e =  -460, k = -145
        {0x0C7CABA6E7C5382C, 0x8FE64A52EE96B8FD}, // e =  -458, k = -144
        {0x0F9BD690A1B68637, 0xB3DFDCE7AA3C673C}, // e =  -456, k = -143
        {0x09C1661A651213E2, 0xD06BEA10CA65C085}, // e =  -453, k = -142
        {0x0C31BFA0FE5698DB, 0x8486E494FCFF30A7}, // e =  -451, k = -141
        {0x0F3E2F893DEC3F12, 0x65A89DBA3C3EFCD0}, // e =  -449, k = -140
        {0x0986DDB5C6B3A76B, 0x7F89629465A75E02}, // e =  -446, k = -139
        {0x0BE8952338609146, 0x5F6BBB397F113583}, // e =  -444, k = -138
        {0x0EE2BA6C0678B597, 0xF746AA07DED582E3}, // e =  -442, k = -137
        {0x094DB483840B717E, 0xFA8C2A44EB4571CE}, // e =  -439, k = -136
        {0x0BA121A4650E4DDE, 0xB92F34D62616CE42}, // e =  -437, k = -135
        {0x0E896A0D7E51E156, 0x677B020BAF9C81D2}, // e =  -435, k = -134
        {0x0915E2486EF32CD6, 0x00ACE1474DC1D123}, // e =  -432, k = -133
        {0x0B5B5ADA8AAFF80B, 0x80D819992132456C}, // e =  -430, k = -132
        {0x0E3231912D5BF60E, 0x610E1FFF697ED6C7}, // e =  -428, k = -131
        {0x08DF5EFABC5979C8, 0xFCA8D3FFA1EF463D}, // e =  -425, k = -130
        {0x0B1736B96B6FD83B, 0x3BD308FF8A6B17CC}, // e =  -423, k = -129
        {0x0DDD0467C64BCE4A, 0x0AC7CB3F6D05DDBE}, // e =  -421, k = -128
        {0x08AA22C0DBEF60EE, 0x46BCDF07A423AA97}, // e =  -418, k = -127
        {0x0AD4AB7112EB3929, 0xD86C16C98D2C953D}, // e =  -416, k = -126
        {0x0D89D64D57A60774, 0x4E871C7BF077BA8C}, // e =  -414, k = -125
        {0x087625F056C7C4A8, 0xB11471CD764AD498}, // e =  -411, k = -124
        {0x0A93AF6C6C79B5D2, 0xDD598E40D3DD89BD}, // e =  -409, k = -123
        {0x0D389B4787982347, 0x94AFF1D108D4EC2D}, // e =  -407, k = -122
        {0x0843610CB4BF160C, 0xBCEDF722A585139C}, // e =  -404, k = -121
        {0x0A54394FE1EEDB8F, 0xEC2974EB4EE65883}, // e =  -402, k = -120
        {0x0CE947A3DA6A9273, 0xE733D226229FEEA4}, // e =  -400, k = -119
        {0x0811CCC668829B88, 0x70806357D5A3F526}, // e =  -397, k = -118
        {0x0A163FF802A3426A, 0x8CA07C2DCB0CF270}, // e =  -395, k = -117
        {0x0C9BCFF6034C1305, 0x2FC89B393DD02F0C}, // e =  -393, k = -116
        {0x0FC2C3F3841F17C6, 0x7BBAC2078D443ACF}, // e =  -391, k = -115
        {0x09D9BA7832936EDC, 0x0D54B944B84AA4C1}, // e =  -388, k = -114
        {0x0C5029163F384A93, 0x10A9E795E65D4DF2}, // e =  -386, k = -113
        {0x0F64335BCF065D37, 0xD4D4617B5FF4A16E}, // e =  -384, k = -112
        {0x099EA0196163FA42, 0xE504BCED1BF8E4E5}, // e =  -381, k = -111
        {0x0C06481FB9BCF8D3, 0x9E45EC2862F71E1E}, // e =  -379, k = -110
        {0x0F07DA27A82C3708, 0x85D767327BB4E5A5}, // e =  -377, k = -109
        {0x0964E858C91BA265, 0x53A6A07F8D510F87}, // e =  -374, k = -108
        {0x0BBE226EFB628AFE, 0xA890489F70A55369}, // e =  -372, k = -107
        {0x0EADAB0ABA3B2DBE, 0x52B45AC74CCEA843}, // e =  -370, k = -106
        {0x092C8AE6B464FC96, 0xF3B0B8BC9001292A}, // e =  -367, k = -105
        {0x0B77ADA0617E3BBC, 0xB09CE6EBB4017375}, // e =  -365, k = -104
        {0x0E55990879DDCAAB, 0xDCC420A6A101D052}, // e =  -363, k = -103
        {0x08F57FA54C2A9EAB, 0x69FA946824A12233}, // e =  -360, k = -102
        {0x0B32DF8E9F354656, 0x447939822DC96AC0}, // e =  -358, k = -101
        {0x0DFF9772470297EB, 0xD59787E2B93BC570}, // e =  -356, k = -100
        {0x08BFBEA76C619EF3, 0x657EB4EDB3C55B66}, // e =  -353, k =  -99
        {0x0AEFAE51477A06B0, 0x3EDE622920B6B240}, // e =  -351, k =  -98
        {0x0DAB99E59958885C, 0x4E95FAB368E45ECF}, // e =  -349, k =  -97
        {0x088B402F7FD75539, 0xB11DBCB0218EBB42}, // e =  -346, k =  -96
        {0x0AAE103B5FCD2A88, 0x1D652BDC29F26A12}, // e =  -344, k =  -95
        {0x0D59944A37C0752A, 0x24BE76D3346F0496}, // e =  -342, k =  -94
        {0x0857FCAE62D8493A, 0x56F70A4400C562DE}, // e =  -339, k =  -93
        {0x0A6DFBD9FB8E5B88, 0xECB4CCD500F6BB96}, // e =  -337, k =  -92
        {0x0D097AD07A71F26B, 0x27E2000A41346A7B}, // e =  -335, k =  -91
        {0x0825ECC24C873782, 0xF8ED400668C0C28D}, // e =  -332, k =  -90
        {0x0A2F67F2DFA90563, 0xB728900802F0F330}, // e =  -330, k =  -89
        {0x0CBB41EF979346BC, 0xA4F2B40A03AD2FFC}, // e =  -328, k =  -88
        {0x0FEA126B7D78186B, 0xCE2F610C84987BFB}, // e =  -326, k =  -87
        {0x09F24B832E6B0F43, 0x60DD9CA7D2DF4D7D}, // e =  -323, k =  -86
        {0x0C6EDE63FA05D314, 0x391503D1C79720DC}, // e =  -321, k =  -85
        {0x0F8A95FCF88747D9, 0x475A44C6397CE913}, // e =  -319, k =  -84
        {0x09B69DBE1B548CE7, 0xCC986AFBE3EE11AC}, // e =  -316, k =  -83
        {0x0C24452DA229B021, 0xBFBE85BADCE99617}, // e =  -314, k =  -82
        {0x0F2D56790AB41C2A, 0x2FAE27299423FB9D}, // e =  -312, k =  -81
        {0x097C560BA6B0919A, 0x5DCCD879FC967D42}, // e =  -309, k =  -80
        {0x0BDB6B8E905CB600, 0xF5400E987BBC1C93}, // e =  -307, k =  -79
        {0x0ED246723473E381, 0x3290123E9AAB23B7}, // e =  -305, k =  -78
        {0x09436C0760C86E30, 0xBF9A0B6720AAF653}, // e =  -302, k =  -77
        {0x0B94470938FA89BC, 0xEF808E40E8D5B3E7}, // e =  -300, k =  -76
        {0x0E7958CB87392C2C, 0x2B60B1D1230B20E1}, // e =  -298, k =  -75
        {0x090BD77F3483BB9B, 0x9B1C6F22B5E6F48D}, // e =  -295, k =  -74
        {0x0B4ECD5F01A4AA82, 0x81E38AEB6360B1B0}, // e =  -293, k =  -73
        {0x0E2280B6C20DD523, 0x225C6DA63C38DE1C}, // e =  -291, k =  -72
        {0x08D590723948A535, 0xF579C487E5A38AD1}, // e =  -288, k =  -71
        {0x0B0AF48EC79ACE83, 0x72D835A9DF0C6D86}, // e =  -286, k =  -70
        {0x0DCDB1B279818224, 0x4F8E431456CF88E7}, // e =  -284, k =  -69
        {0x08A08F0F8BF0F156, 0xB1B8E9ECB641B590}, // e =  -281, k =  -68
        {0x0AC8B2D36EED2DAC, 0x5E272467E3D222F4}, // e =  -279, k =  -67
        {0x0D7ADF884AA87917, 0x75B0ED81DCC6ABB1}, // e =  -277, k =  -66
        {0x086CCBB52EA94BAE, 0xA98E947129FC2B4F}, // e =  -274, k =  -65
        {0x0A87FEA27A539E9A, 0x53F2398D747B3623}, // e =  -272, k =  -64
        {0x0D29FE4B18E88640, 0xE8EEC7F0D19A03AB}, // e =  -270, k =  -63
        {0x083A3EEEEF9153E8, 0x91953CF68300424B}, // e =  -267, k =  -62
        {0x0A48CEAAAB75A8E2, 0xB5FA8C3423C052DE}, // e =  -265, k =  -61
        {0x0CDB02555653131B, 0x63792F412CB06795}, // e =  -263, k =  -60
        {0x0808E17555F3EBF1, 0x1E2BBD88BBEE40BE}, // e =  -260, k =  -59
        {0x0A0B19D2AB70E6ED, 0x65B6ACEAEAE9D0ED}, // e =  -258, k =  -58
        {0x0C8DE047564D20A8, 0xBF245825A5A44528}, // e =  -256, k =  -57
        {0x0FB158592BE068D2, 0xEEED6E2F0F0D5672}, // e =  -254, k =  -56
        {0x09CED737BB6C4183, 0xD55464DD69685607}, // e =  -251, k =  -55
        {0x0C428D05AA4751E4, 0xCAA97E14C3C26B89}, // e =  -249, k =  -54
        {0x0F53304714D9265D, 0xFD53DD99F4B3066B}, // e =  -247, k =  -53
        {0x0993FE2C6D07B7FA, 0xBE546A8038EFE403}, // e =  -244, k =  -52
        {0x0BF8FDB78849A5F9, 0x6DE98520472BDD04}, // e =  -242, k =  -51
        {0x0EF73D256A5C0F77, 0xC963E66858F6D445}, // e =  -240, k =  -50
        {0x095A8637627989AA, 0xDDDE7001379A44AB}, // e =  -237, k =  -49
        {0x0BB127C53B17EC15, 0x95560C018580D5D6}, // e =  -235, k =  -48
        {0x0E9D71B689DDE71A, 0xFAAB8F01E6E10B4B}, // e =  -233, k =  -47
        {0x09226712162AB070, 0xDCAB3961304CA70F}, // e =  -230, k =  -46
        {0x0B6B00D69BB55C8D, 0x13D607B97C5FD0D3}, // e =  -228, k =  -45
        {0x0E45C10C42A2B3B0, 0x58CB89A7DB77C507}, // e =  -226, k =  -44
        {0x08EB98A7A9A5B04E, 0x377F3608E92ADB25}, // e =  -223, k =  -43
        {0x0B267ED1940F1C61, 0xC55F038B237591EE}, // e =  -221, k =  -42
        {0x0DF01E85F912E37A, 0x36B6C46DEC52F669}, // e =  -219, k =  -41
        {0x08B61313BBABCE2C, 0x62323AC4B3B3DA02}, // e =  -216, k =  -40
        {0x0AE397D8AA96C1B7, 0x7ABEC975E0A0D082}, // e =  -214, k =  -39
        {0x0D9C7DCED53C7225, 0x596E7BD358C904A3}, // e =  -212, k =  -38
        {0x0881CEA14545C757, 0x57E50D64177DA2E6}, // e =  -209, k =  -37
        {0x0AA242499697392D, 0x2DDE50BD1D5D0B9F}, // e =  -207, k =  -36
        {0x0D4AD2DBFC3D0778, 0x7955E4EC64B44E87}, // e =  -205, k =  -35
        {0x084EC3C97DA624AB, 0x4BD5AF13BEF0B114}, // e =  -202, k =  -34
        {0x0A6274BBDD0FADD6, 0x1ECB1AD8AEACDD59}, // e =  -200, k =  -33
        {0x0CFB11EAD453994B, 0xA67DE18EDA5814B0}, // e =  -198, k =  -32
        {0x081CEB32C4B43FCF, 0x480EACF948770CEE}, // e =  -195, k =  -31
        {0x0A2425FF75E14FC3, 0x1A1258379A94D029}, // e =  -193, k =  -30
        {0x0CAD2F7F5359A3B3, 0xE096EE45813A0434}, // e =  -191, k =  -29
        {0x0FD87B5F28300CA0, 0xD8BCA9D6E1888540}, // e =  -189, k =  -28
        {0x09E74D1B791E07E4, 0x8775EA264CF55348}, // e =  -186, k =  -27
        {0x0C612062576589DD, 0xA95364AFE032A81A}, // e =  -184, k =  -26
        {0x0F79687AED3EEC55, 0x13A83DDBD83F5221}, // e =  -182, k =  -25
        {0x09ABE14CD44753B5, 0x2C4926A967279355}, // e =  -179, k =  -24
        {0x0C16D9A0095928A2, 0x775B7053C0F1782A}, // e =  -177, k =  -23
        {0x0F1C90080BAF72CB, 0x15324C68B12DD634}, // e =  -175, k =  -22
        {0x0971DA05074DA7BE, 0xED3F6FC16EBCA5E1}, // e =  -172, k =  -21
        {0x0BCE5086492111AE, 0xA88F4BB1CA6BCF59}, // e =  -170, k =  -20
        {0x0EC1E4A7DB69561A, 0x52B31E9E3D06C32F}, // e =  -168, k =  -19
        {0x09392EE8E921D5D0, 0x73AFF322E62439FD}, // e =  -165, k =  -18
        {0x0B877AA3236A4B44, 0x909BEFEB9FAD487D}, // e =  -163, k =  -17
        {0x0E69594BEC44DE15, 0xB4C2EBE687989A9C}, // e =  -161, k =  -16
        {0x0901D7CF73AB0ACD, 0x90F9D37014BF60A2}, // e =  -158, k =  -15
        {0x0B424DC35095CD80, 0xF538484C19EF38CA}, // e =  -156, k =  -14
        {0x0E12E13424BB40E1, 0x32865A5F206B06FC}, // e =  -154, k =  -13
        {0x08CBCCC096F5088C, 0xBF93F87B7442E45E}, // e =  -151, k =  -12
        {0x0AFEBFF0BCB24AAF, 0xEF78F69A51539D75}, // e =  -149, k =  -11
        {0x0DBE6FECEBDEDD5B, 0xEB573440E5A884D2}, // e =  -147, k =  -10
        {0x089705F4136B4A59, 0x731680A88F895304}, // e =  -144, k =   -9
        {0x0ABCC77118461CEF, 0xCFDC20D2B36BA7C4}, // e =  -142, k =   -8
        {0x0D6BF94D5E57A42B, 0xC3D32907604691B5}, // e =  -140, k =   -7
        {0x08637BD05AF6C69B, 0x5A63F9A49C2C1B11}, // e =  -137, k =   -6
        {0x0A7C5AC471B47842, 0x30FCF80DC33721D6}, // e =  -135, k =   -5
        {0x0D1B71758E219652, 0xBD3C36113404EA4B}, // e =  -133, k =   -4
        {0x083126E978D4FDF3, 0xB645A1CAC083126F}, // e =  -130, k =   -3
        {0x0A3D70A3D70A3D70, 0xA3D70A3D70A3D70B}, // e =  -128, k =   -2
        {0x0CCCCCCCCCCCCCCC, 0xCCCCCCCCCCCCCCCD}, // e =  -126, k =   -1
        {0x0800000000000000, 0x0000000000000000}, // e =  -123, k =    0
        {0x0A00000000000000, 0x0000000000000000}, // e =  -121, k =    1
        {0x0C80000000000000, 0x0000000000000000}, // e =  -119, k =    2
        {0x0FA0000000000000, 0x0000000000000000}, // e =  -117, k =    3
        {0x09C4000000000000, 0x0000000000000000}, // e =  -114, k =    4
        {0x0C35000000000000, 0x0000000000000000}, // e =  -112, k =    5
        {0x0F42400000000000, 0x0000000000000000}, // e =  -110, k =    6
        {0x0989680000000000, 0x0000000000000000}, // e =  -107, k =    7
        {0x0BEBC20000000000, 0x0000000000000000}, // e =  -105, k =    8
        {0x0EE6B28000000000, 0x0000000000000000}, // e =  -103, k =    9
        {0x09502F9000000000, 0x0000000000000000}, // e =  -100, k =   10
        {0x0BA43B7400000000, 0x0000000000000000}, // e =   -98, k =   11
        {0x0E8D4A5100000000, 0x0000000000000000}, // e =   -96, k =   12
        {0x09184E72A0000000, 0x0000000000000000}, // e =   -93, k =   13
        {0x0B5E620F48000000, 0x0000000000000000}, // e =   -91, k =   14
        {0x0E35FA931A000000, 0x0000000000000000}, // e =   -89, k =   15
        {0x08E1BC9BF0400000, 0x0000000000000000}, // e =   -86, k =   16
        {0x0B1A2BC2EC500000, 0x0000000000000000}, // e =   -84, k =   17
        {0x0DE0B6B3A7640000, 0x0000000000000000}, // e =   -82, k =   18
        {0x08AC7230489E8000, 0x0000000000000000}, // e =   -79, k =   19
        {0x0AD78EBC5AC62000, 0x0000000000000000}, // e =   -77, k =   20
        {0x0D8D726B7177A800, 0x0000000000000000}, // e =   -75, k =   21
        {0x0878678326EAC900, 0x0000000000000000}, // e =   -72, k =   22
        {0x0A968163F0A57B40, 0x0000000000000000}, // e =   -70, k =   23
        {0x0D3C21BCECCEDA10, 0x0000000000000000}, // e =   -68, k =   24
        {0x084595161401484A, 0x0000000000000000}, // e =   -65, k =   25
        {0x0A56FA5B99019A5C, 0x8000000000000000}, // e =   -63, k =   26
        {0x0CECB8F27F4200F3, 0xA000000000000000}, // e =   -61, k =   27
        {0x0813F3978F894098, 0x4400000000000000}, // e =   -58, k =   28
        {0x0A18F07D736B90BE, 0x5500000000000000}, // e =   -56, k =   29
        {0x0C9F2C9CD04674ED, 0xEA40000000000000}, // e =   -54, k =   30
        {0x0FC6F7C404581229, 0x64D0000000000000}, // e =   -52, k =   31
        {0x09DC5ADA82B70B59, 0xDF02000000000000}, // e =   -49, k =   32
        {0x0C5371912364CE30, 0x56C2800000000000}, // e =   -47, k =   33
        {0x0F684DF56C3E01BC, 0x6C73200000000000}, // e =   -45, k =   34
        {0x09A130B963A6C115, 0xC3C7F40000000000}, // e =   -42, k =   35
        {0x0C097CE7BC90715B, 0x34B9F10000000000}, // e =   -40, k =   36
        {0x0F0BDC21ABB48DB2, 0x01E86D4000000000}, // e =   -38, k =   37
        {0x096769950B50D88F, 0x4131444800000000}, // e =   -35, k =   38
        {0x0BC143FA4E250EB3, 0x117D955A00000000}, // e =   -33, k =   39
        {0x0EB194F8E1AE525F, 0xD5DCFAB080000000}, // e =   -31, k =   40
        {0x092EFD1B8D0CF37B, 0xE5AA1CAE50000000}, // e =   -28, k =   41
        {0x0B7ABC627050305A, 0xDF14A3D9E4000000}, // e =   -26, k =   42
        {0x0E596B7B0C643C71, 0x96D9CCD05D000000}, // e =   -24, k =   43
        {0x08F7E32CE7BEA5C6, 0xFE4820023A200000}, // e =   -21, k =   44
        {0x0B35DBF821AE4F38, 0xBDDA2802C8A80000}, // e =   -19, k =   45
        {0x0E0352F62A19E306, 0xED50B2037AD20000}, // e =   -17, k =   46
        {0x08C213D9DA502DE4, 0x54526F422CC34000}, // e =   -14, k =   47
        {0x0AF298D050E4395D, 0x69670B12B7F41000}, // e =   -12, k =   48
        {0x0DAF3F04651D47B4, 0xC3C0CDD765F11400}, // e =   -10, k =   49
        {0x088D8762BF324CD0, 0xFA5880A69FB6AC80}, // e =    -7, k =   50
        {0x0AB0E93B6EFEE005, 0x38EEA0D047A457A0}, // e =    -5, k =   51
        {0x0D5D238A4ABE9806, 0x872A4904598D6D88}, // e =    -3, k =   52
        {0x085A36366EB71F04, 0x147A6DA2B7F86475}, // e =     0, k =   53
        {0x0A70C3C40A64E6C5, 0x1999090B65F67D93}, // e =     2, k =   54
        {0x0D0CF4B50CFE2076, 0x5FFF4B4E3F741CF7}, // e =     4, k =   55
        {0x082818F1281ED449, 0xFBFF8F10E7A8921B}, // e =     7, k =   56
        {0x0A321F2D7226895C, 0x7AFF72D52192B6A1}, // e =     9, k =   57
        {0x0CBEA6F8CEB02BB3, 0x99BF4F8A69F7644A}, // e =    11, k =   58
        {0x0FEE50B7025C36A0, 0x802F236D04753D5C}, // e =    13, k =   59
        {0x09F4F2726179A224, 0x501D762422C9465A}, // e =    16, k =   60
        {0x0C722F0EF9D80AAD, 0x6424D3AD2B7B97F0}, // e =    18, k =   61
        {0x0F8EBAD2B84E0D58, 0xBD2E0898765A7DEC}, // e =    20, k =   62
        {0x09B934C3B330C857, 0x763CC55F49F88EB3}, // e =    23, k =   63
        {0x0C2781F49FFCFA6D, 0x53CBF6B71C76B260}, // e =    25, k =   64
        {0x0F316271C7FC3908, 0xA8BEF464E3945EF8}, // e =    27, k =   65
        {0x097EDD871CFDA3A5, 0x697758BF0E3CBB5B}, // e =    30, k =   66
        {0x0BDE94E8E43D0C8E, 0xC3D52EEED1CBEA32}, // e =    32, k =   67
        {0x0ED63A231D4C4FB2, 0x74CA7AAA863EE4BE}, // e =    34, k =   68
        {0x0945E455F24FB1CF, 0x88FE8CAA93E74EF7}, // e =    37, k =   69
        {0x0B975D6B6EE39E43, 0x6B3E2FD538E122B5}, // e =    39, k =   70
        {0x0E7D34C64A9C85D4, 0x460DBBCA87196B62}, // e =    41, k =   71
        {0x090E40FBEEA1D3A4, 0xABC8955E946FE31D}, // e =    44, k =   72
        {0x0B51D13AEA4A488D, 0xD6BABAB6398BDBE5}, // e =    46, k =   73
        {0x0E264589A4DCDAB1, 0x4C696963C7EED2DE}, // e =    48, k =   74
        {0x08D7EB76070A08AE, 0xCFC1E1DE5CF543CB}, // e =    51, k =   75
        {0x0B0DE65388CC8ADA, 0x83B25A55F43294BD}, // e =    53, k =   76
        {0x0DD15FE86AFFAD91, 0x249EF0EB713F39EC}, // e =    55, k =   77
        {0x08A2DBF142DFCC7A, 0xB6E3569326C78434}, // e =    58, k =   78
        {0x0ACB92ED9397BF99, 0x649C2C37F0796541}, // e =    60, k =   79
        {0x0D7E77A8F87DAF7F, 0xBDC33745EC97BE91}, // e =    62, k =   80
        {0x086F0AC99B4E8DAF, 0xD69A028BB3DED71B}, // e =    65, k =   81
        {0x0A8ACD7C0222311B, 0xCC40832EA0D68CE1}, // e =    67, k =   82
        {0x0D2D80DB02AABD62, 0xBF50A3FA490C301A}, // e =    69, k =   83
        {0x083C7088E1AAB65D, 0xB792667C6DA79E10}, // e =    72, k =   84
        {0x0A4B8CAB1A1563F5, 0x2577001B89118594}, // e =    74, k =   85
        {0x0CDE6FD5E09ABCF2, 0x6ED4C0226B55E6F9}, // e =    76, k =   86
        {0x080B05E5AC60B617, 0x8544F8158315B05C}, // e =    79, k =   87
        {0x0A0DC75F1778E39D, 0x6696361AE3DB1C73}, // e =    81, k =   88
        {0x0C913936DD571C84, 0xC03BC3A19CD1E38F}, // e =    83, k =   89
        {0x0FB5878494ACE3A5, 0xF04AB48A04065C73}, // e =    85, k =   90
        {0x09D174B2DCEC0E47, 0xB62EB0D64283F9C8}, // e =    88, k =   91
        {0x0C45D1DF942711D9, 0xA3BA5D0BD324F83A}, // e =    90, k =   92
        {0x0F5746577930D650, 0x0CA8F44EC7EE3648}, // e =    92, k =   93
        {0x09968BF6ABBE85F2, 0x07E998B13CF4E1ED}, // e =    95, k =   94
        {0x0BFC2EF456AE276E, 0x89E3FEDD8C321A68}, // e =    97, k =   95
        {0x0EFB3AB16C59B14A, 0x2C5CFE94EF3EA102}, // e =    99, k =   96
        {0x095D04AEE3B80ECE, 0x5BBA1F1D158724A2}, // e =   102, k =   97
        {0x0BB445DA9CA61281, 0xF2A8A6E45AE8EDCA}, // e =   104, k =   98
        {0x0EA1575143CF9722, 0x6F52D09D71A3293C}, // e =   106, k =   99
        {0x0924D692CA61BE75, 0x8593C2626705F9C6}, // e =   109, k =  100
        {0x0B6E0C377CFA2E12, 0xE6F8B2FB00C77837}, // e =   111, k =  101
        {0x0E498F455C38B997, 0xA0B6DFB9C0F95645}, // e =   113, k =  102
        {0x08EDF98B59A373FE, 0xC4724BD4189BD5EB}, // e =   116, k =  103
        {0x0B2977EE300C50FE, 0x758EDEC91EC2CB66}, // e =   118, k =  104
        {0x0DF3D5E9BC0F653E, 0x12F2967B66737E3F}, // e =   120, k =  105
        {0x08B865B215899F46, 0xCBD79E0D20082EE8}, // e =   123, k =  106
        {0x0AE67F1E9AEC0718, 0x7ECD8590680A3AA2}, // e =   125, k =  107
        {0x0DA01EE641A708DE, 0x9E80E6F4820CC94A}, // e =   127, k =  108
        {0x0884134FE908658B, 0x23109058D147FDCE}, // e =   130, k =  109
        {0x0AA51823E34A7EED, 0xEBD4B46F0599FD42}, // e =   132, k =  110
        {0x0D4E5E2CDC1D1EA9, 0x66C9E18AC7007C92}, // e =   134, k =  111
        {0x0850FADC09923329, 0xE03E2CF6BC604DDC}, // e =   137, k =  112
        {0x0A6539930BF6BFF4, 0x584DB8346B786152}, // e =   139, k =  113
        {0x0CFE87F7CEF46FF1, 0x6E612641865679A7}, // e =   141, k =  114
        {0x081F14FAE158C5F6, 0xE4FCB7E8F3F60C08}, // e =   144, k =  115
        {0x0A26DA3999AEF774, 0x9E3BE5E330F38F0A}, // e =   146, k =  116
        {0x0CB090C8001AB551, 0xC5CADF5BFD3072CD}, // e =   148, k =  117
        {0x0FDCB4FA002162A6, 0x373D9732FC7C8F80}, // e =   150, k =  118
        {0x09E9F11C4014DDA7, 0xE2867E7FDDCDD9B0}, // e =   153, k =  119
        {0x0C646D63501A1511, 0xDB281E1FD541501C}, // e =   155, k =  120
        {0x0F7D88BC24209A56, 0x51F225A7CA91A423}, // e =   157, k =  121
        {0x09AE757596946075, 0xF3375788DE9B0696}, // e =   160, k =  122
        {0x0C1A12D2FC397893, 0x70052D6B1641C83B}, // e =   162, k =  123
        {0x0F209787BB47D6B8, 0x4C0678C5DBD23A4A}, // e =   164, k =  124
        {0x09745EB4D50CE633, 0x2F840B7BA963646F}, // e =   167, k =  125
        {0x0BD176620A501FBF, 0xFB650E5A93BC3D8A}, // e =   169, k =  126
        {0x0EC5D3FA8CE427AF, 0xFA3E51F138AB4CEC}, // e =   171, k =  127
        {0x093BA47C980E98CD, 0xFC66F336C36B1014}, // e =   174, k =  128
        {0x0B8A8D9BBE123F01, 0x7B80B0047445D419}, // e =   176, k =  129
        {0x0E6D3102AD96CEC1, 0xDA60DC059157491F}, // e =   178, k =  130
        {0x09043EA1AC7E4139, 0x287C89837AD68DB3}, // e =   181, k =  131
        {0x0B454E4A179DD187, 0x729BABE4598C3120}, // e =   183, k =  132
        {0x0E16A1DC9D8545E9, 0x4F4296DD6FEF3D68}, // e =   185, k =  133
        {0x08CE2529E2734BB1, 0xD1899E4A65F58661}, // e =   188, k =  134
        {0x0B01AE745B101E9E, 0x45EC05DCFF72E7F9}, // e =   190, k =  135
        {0x0DC21A1171D42645, 0xD76707543F4FA1F8}, // e =   192, k =  136
        {0x0899504AE72497EB, 0xA6A06494A791C53B}, // e =   195, k =  137
        {0x0ABFA45DA0EDBDE6, 0x90487DB9D176368A}, // e =   197, k =  138
        {0x0D6F8D7509292D60, 0x345A9D2845D3C42C}, // e =   199, k =  139
        {0x0865B86925B9BC5C, 0x20B8A2392BA45A9C}, // e =   202, k =  140
        {0x0A7F26836F282B73, 0x28E6CAC7768D7142}, // e =   204, k =  141
        {0x0D1EF0244AF2364F, 0xF3207D795430CD93}, // e =   206, k =  142
        {0x08335616AED761F1, 0xF7F44E6BD49E807C}, // e =   209, k =  143
        {0x0A402B9C5A8D3A6E, 0x75F16206C9C6209B}, // e =   211, k =  144
        {0x0CD036837130890A, 0x136DBA887C37A8C1}, // e =   213, k =  145
        {0x0802221226BE55A6, 0x4C2494954DA2C979}, // e =   216, k =  146
        {0x0A02AA96B06DEB0F, 0xDF2DB9BAA10B7BD7}, // e =   218, k =  147
        {0x0C83553C5C8965D3, 0xD6F92829494E5ACD}, // e =   220, k =  148
        {0x0FA42A8B73ABBF48, 0xCCB772339BA1F180}, // e =   222, k =  149
        {0x09C69A97284B578D, 0x7FF2A760414536F0}, // e =   225, k =  150
        {0x0C38413CF25E2D70, 0xDFEF5138519684AC}, // e =   227, k =  151
        {0x0F46518C2EF5B8CD, 0x17EB258665FC25D7}, // e =   229, k =  152
        {0x098BF2F79D599380, 0x2EF2F773FFBD97A7}, // e =   232, k =  153
        {0x0BEEEFB584AFF860, 0x3AAFB550FFACFD90}, // e =   234, k =  154
        {0x0EEAABA2E5DBF678, 0x495BA2A53F983CF4}, // e =   236, k =  155
        {0x0952AB45CFA97A0B, 0x2DD945A747BF2619}, // e =   239, k =  156
        {0x0BA756174393D88D, 0xF94F971119AEEF9F}, // e =   241, k =  157
        {0x0E912B9D1478CEB1, 0x77A37CD5601AAB86}, // e =   243, k =  158
        {0x091ABB422CCB812E, 0xEAC62E055C10AB34}, // e =   246, k =  159
        {0x0B616A12B7FE617A, 0xA577B986B314D601}, // e =   248, k =  160
        {0x0E39C49765FDF9D9, 0x4ED5A7E85FDA0B81}, // e =   250, k =  161
        {0x08E41ADE9FBEBC27, 0xD14588F13BE84731}, // e =   253, k =  162
        {0x0B1D219647AE6B31, 0xC596EB2D8AE258FD}, // e =   255, k =  163
        {0x0DE469FBD99A05FE, 0x36FCA5F8ED9AEF3C}, // e =   257, k =  164
        {0x08AEC23D680043BE, 0xE25DE7BB9480D586}, // e =   260, k =  165
        {0x0ADA72CCC20054AE, 0x9AF561AA79A10AE7}, // e =   262, k =  166
        {0x0D910F7FF28069DA, 0x41B2BA1518094DA1}, // e =   264, k =  167
        {0x087AA9AFF7904228, 0x690FB44D2F05D085}, // e =   267, k =  168
        {0x0A99541BF57452B2, 0x8353A1607AC744A6}, // e =   269, k =  169
        {0x0D3FA922F2D1675F, 0x242889B8997915CF}, // e =   271, k =  170
        {0x0847C9B5D7C2E09B, 0x769956135FEBADA2}, // e =   274, k =  171
        {0x0A59BC234DB398C2, 0x543FAB9837E6990A}, // e =   276, k =  172
        {0x0CF02B2C21207EF2, 0xE94F967E45E03F4C}, // e =   278, k =  173
        {0x08161AFB94B44F57, 0xD1D1BE0EEBAC2790}, // e =   281, k =  174
        {0x0A1BA1BA79E1632D, 0xC6462D92A6973174}, // e =   283, k =  175
        {0x0CA28A291859BBF9, 0x37D7B8F7503CFDD0}, // e =   285, k =  176
        {0x0FCB2CB35E702AF7, 0x85CDA735244C3D44}, // e =   287, k =  177
        {0x09DEFBF01B061ADA, 0xB3A0888136AFA64B}, // e =   290, k =  178
        {0x0C56BAEC21C7A191, 0x6088AAA1845B8FDE}, // e =   292, k =  179
        {0x0F6C69A72A3989F5, 0xB8AAD549E57273D5}, // e =   294, k =  180
        {0x09A3C2087A63F639, 0x936AC54E2F678865}, // e =   297, k =  181
        {0x0C0CB28A98FCF3C7, 0xF84576A1BB416A7E}, // e =   299, k =  182
        {0x0F0FDF2D3F3C30B9, 0xF656D44A2A11C51E}, // e =   301, k =  183
        {0x0969EB7C47859E74, 0x39F644AE5A4B1B33}, // e =   304, k =  184
        {0x0BC4665B59670611, 0x4873D5D9F0DDE1FF}, // e =   306, k =  185
        {0x0EB57FF22FC0C795, 0x9A90CB506D155A7F}, // e =   308, k =  186
        {0x09316FF75DD87CBD, 0x809A7F12442D5890}, // e =   311, k =  187
        {0x0B7DCBF5354E9BEC, 0xE0C11ED6D538AEB3}, // e =   313, k =  188
        {0x0E5D3EF282A242E8, 0x18F1668C8A86DA60}, // e =   315, k =  189
        {0x08FA475791A569D1, 0x0F96E017D694487C}, // e =   318, k =  190
        {0x0B38D92D760EC445, 0x537C981DCC395A9B}, // e =   320, k =  191
        {0x0E070F78D3927556, 0xA85BBE253F47B142}, // e =   322, k =  192
        {0x08C469AB843B8956, 0x293956D7478CCEC9}, // e =   325, k =  193
        {0x0AF58416654A6BAB, 0xB387AC8D1970027C}, // e =   327, k =  194
        {0x0DB2E51BFE9D0696, 0xA06997B05FCC031A}, // e =   329, k =  195
        {0x088FCF317F22241E, 0x2441FECE3BDF81F1}, // e =   332, k =  196
        {0x0AB3C2FDDEEAAD25, 0xAD527E81CAD7626D}, // e =   334, k =  197
        {0x0D60B3BD56A5586F, 0x18A71E223D8D3B08}, // e =   336, k =  198
        {0x085C705656275745, 0x6F6872D5667844E5}, // e =   339, k =  199
        {0x0A738C6BEBB12D16, 0xCB428F8AC016561E}, // e =   341, k =  200
        {0x0D106F86E69D785C, 0x7E13336D701BEBA6}, // e =   343, k =  201
        {0x082A45B450226B39, 0xCECC002466117348}, // e =   346, k =  202
        {0x0A34D721642B0608, 0x427F002D7F95D01A}, // e =   348, k =  203
        {0x0CC20CE9BD35C78A, 0x531EC038DF7B4420}, // e =   350, k =  204
        {0x0FF290242C83396C, 0xE7E67047175A1528}, // e =   352, k =  205
        {0x09F79A169BD203E4, 0x10F0062C6E984D39}, // e =   355, k =  206
        {0x0C75809C42C684DD, 0x152C07B78A3E6087}, // e =   357, k =  207
        {0x0F92E0C353782614, 0x5A7709A56CCDF8A9}, // e =   359, k =  208
        {0x09BBCC7A142B17CC, 0xB88A66076400BB6A}, // e =   362, k =  209
        {0x0C2ABF989935DDBF, 0xE6ACFF893D00EA44}, // e =   364, k =  210
        {0x0F356F7EBF83552F, 0xE0583F6B8C4124D5}, // e =   366, k =  211
        {0x098165AF37B2153D, 0xEC3727A337A8B705}, // e =   369, k =  212
        {0x0BE1BF1B059E9A8D, 0x6744F18C0592E4C6}, // e =   371, k =  213
        {0x0EDA2EE1C7064130, 0xC1162DEF06F79DF8}, // e =   373, k =  214
        {0x09485D4D1C63E8BE, 0x78ADDCB5645AC2BB}, // e =   376, k =  215
        {0x0B9A74A0637CE2EE, 0x16D953E2BD71736A}, // e =   378, k =  216
        {0x0E8111C87C5C1BA9, 0x9C8FA8DB6CCDD044}, // e =   380, k =  217
        {0x0910AB1D4DB9914A, 0x01D9C9892400A22B}, // e =   383, k =  218
        {0x0B54D5E4A127F59C, 0x82503BEB6D00CAB5}, // e =   385, k =  219
        {0x0E2A0B5DC971F303, 0xA2E44AE64840FD62}, // e =   387, k =  220
        {0x08DA471A9DE737E2, 0x45CEAECFED289E5E}, // e =   390, k =  221
        {0x0B10D8E1456105DA, 0xD7425A83E872C5F5}, // e =   392, k =  222
        {0x0DD50F1996B94751, 0x8D12F124E28F7772}, // e =   394, k =  223
        {0x08A5296FFE33CC92, 0xF82BD6B70D99AAA7}, // e =   397, k =  224
        {0x0ACE73CBFDC0BFB7, 0xB636CC64D1001551}, // e =   399, k =  225
        {0x0D8210BEFD30EFA5, 0xA3C47F7E05401AA5}, // e =   401, k =  226
        {0x08714A775E3E95C7, 0x865ACFAEC34810A8}, // e =   404, k =  227
        {0x0A8D9D1535CE3B39, 0x67F1839A741A14D1}, // e =   406, k =  228
        {0x0D31045A8341CA07, 0xC1EDE48111209A06}, // e =   408, k =  229
        {0x083EA2B892091E44, 0xD934AED0AAB46044}, // e =   411, k =  230
        {0x0A4E4B66B68B65D6, 0x0F81DA84D5617854}, // e =   413, k =  231
        {0x0CE1DE40642E3F4B, 0x936251260AB9D669}, // e =   415, k =  232
        {0x080D2AE83E9CE78F, 0x3C1D72B7C6B42602}, // e =   418, k =  233
        {0x0A1075A24E442173, 0x0B24CF65B8612F82}, // e =   420, k =  234
        {0x0C94930AE1D529CF, 0xCDEE033F26797B63}, // e =   422, k =  235
        {0x0FB9B7CD9A4A7443, 0xC169840EF017DA3C}, // e =   424, k =  236
        {0x09D412E0806E88AA, 0x58E1F289560EE865}, // e =   427, k =  237
        {0x0C491798A08A2AD4, 0xEF1A6F2BAB92A27F}, // e =   429, k =  238
        {0x0F5B5D7EC8ACB58A, 0x2AE10AF696774B1E}, // e =   431, k =  239
        {0x09991A6F3D6BF176, 0x5ACCA6DA1E0A8EF3}, // e =   434, k =  240
        {0x0BFF610B0CC6EDD3, 0xF17FD090A58D32B0}, // e =   436, k =  241
        {0x0EFF394DCFF8A948, 0xEDDFC4B4CEF07F5C}, // e =   438, k =  242
        {0x095F83D0A1FB69CD, 0x94ABDAF101564F99}, // e =   441, k =  243
        {0x0BB764C4CA7A4440, 0xF9D6D1AD41ABE380}, // e =   443, k =  244
        {0x0EA53DF5FD18D551, 0x384C86189216DC5F}, // e =   445, k =  245
        {0x092746B9BE2F8552, 0xC32FD3CF5B4E49BC}, // e =   448, k =  246
        {0x0B7118682DBB66A7, 0x73FBC8C33221DC2B}, // e =   450, k =  247
        {0x0E4D5E82392A4051, 0x50FABAF3FEAA5335}, // e =   452, k =  248
        {0x08F05B1163BA6832, 0xD29CB4D87F2A7401}, // e =   455, k =  249
        {0x0B2C71D5BCA9023F, 0x8743E20E9EF51102}, // e =   457, k =  250
        {0x0DF78E4B2BD342CF, 0x6914DA9246B25542}, // e =   459, k =  251
        {0x08BAB8EEFB6409C1, 0xA1AD089B6C2F7549}, // e =   462, k =  252
        {0x0AE9672ABA3D0C32, 0x0A184AC2473B529C}, // e =   464, k =  253
        {0x0DA3C0F568CC4F3E, 0x8C9E5D72D90A2742}, // e =   466, k =  254
        {0x08865899617FB187, 0x17E2FA67C7A6588A}, // e =   469, k =  255
        {0x0AA7EEBFB9DF9DE8, 0xDDDBB901B98FEEAC}, // e =   471, k =  256
        {0x0D51EA6FA8578563, 0x1552A74227F3EA57}, // e =   473, k =  257
        {0x08533285C936B35D, 0xED53A88958F87276}, // e =   476, k =  258
        {0x0A67FF273B846035, 0x68A892ABAF368F14}, // e =   478, k =  259
        {0x0D01FEF10A657842, 0xC2D2B7569B0432D9}, // e =   480, k =  260
        {0x08213F56A67F6B29, 0xB9C3B29620E29FC8}, // e =   483, k =  261
        {0x0A298F2C501F45F4, 0x28349F3BA91B47B9}, // e =   485, k =  262
        {0x0CB3F2F764271771, 0x3241C70A936219A8}, // e =   487, k =  263
        {0x0FE0EFB53D30DD4D, 0x7ED238CD383AA012}, // e =   489, k =  264
        {0x09EC95D1463E8A50, 0x6F4363804324A40B}, // e =   492, k =  265
        {0x0C67BB4597CE2CE4, 0x8B143C6053EDCD0E}, // e =   494, k =  266
        {0x0F81AA16FDC1B81D, 0xADD94B7868E94051}, // e =   496, k =  267
        {0x09B10A4E5E991312, 0x8CA7CF2B4191C833}, // e =   499, k =  268
        {0x0C1D4CE1F63F57D7, 0x2FD1C2F611F63A40}, // e =   501, k =  269
        {0x0F24A01A73CF2DCC, 0xFBC633B39673C8CF}, // e =   503, k =  270
        {0x0976E41088617CA0, 0x1D5BE0503E085D82}, // e =   506, k =  271
        {0x0BD49D14AA79DBC8, 0x24B2D8644D8A74E2}, // e =   508, k =  272
        {0x0EC9C459D51852BA, 0x2DDF8E7D60ED121A}, // e =   510, k =  273
        {0x093E1AB8252F33B4, 0x5CABB90E5C942B51}, // e =   513, k =  274
        {0x0B8DA1662E7B00A1, 0x73D6A751F3B93625}, // e =   515, k =  275
        {0x0E7109BFBA19C0C9, 0xD0CC512670A783AE}, // e =   517, k =  276
        {0x0906A617D450187E, 0x227FB2B80668B24D}, // e =   520, k =  277
        {0x0B484F9DC9641E9D, 0xAB1F9F660802DEE0}, // e =   522, k =  278
        {0x0E1A63853BBD2645, 0x15E7873F8A039698}, // e =   524, k =  279
        {0x08D07E33455637EB, 0x2DB0B487B6423E1F}, // e =   527, k =  280
        {0x0B049DC016ABC5E5, 0xF91CE1A9A3D2CDA7}, // e =   529, k =  281
        {0x0DC5C5301C56B75F, 0x77641A140CC78110}, // e =   531, k =  282
        {0x089B9B3E11B6329B, 0xAA9E904C87FCB0AA}, // e =   534, k =  283
        {0x0AC2820D9623BF42, 0x9546345FA9FBDCD5}, // e =   536, k =  284
        {0x0D732290FBACAF13, 0x3A97C177947AD40A}, // e =   538, k =  285
        {0x0867F59A9D4BED6C, 0x049ED8EABCCCC486}, // e =   541, k =  286
        {0x0A81F301449EE8C7, 0x05C68F256BFFF5A8}, // e =   543, k =  287
        {0x0D226FC195C6A2F8, 0xC73832EEC6FFF312}, // e =   545, k =  288
        {0x083585D8FD9C25DB, 0x7C831FD53C5FF7EB}, // e =   548, k =  289
        {0x0A42E74F3D032F52, 0x5BA3E7CA8B77F5E6}, // e =   550, k =  290
        {0x0CD3A1230C43FB26, 0xF28CE1BD2E55F35F}, // e =   552, k =  291
        {0x080444B5E7AA7CF8, 0x57980D163CF5B81C}, // e =   555, k =  292
        {0x0A0555E361951C36, 0x6D7E105BCC332622}, // e =   557, k =  293
        {0x0C86AB5C39FA6344, 0x08DD9472BF3FEFAB}, // e =   559, k =  294
        {0x0FA856334878FC15, 0x0B14F98F6F0FEB96}, // e =   561, k =  295
        {0x09C935E00D4B9D8D, 0x26ED1BF9A569F33E}, // e =   564, k =  296
        {0x0C3B8358109E84F0, 0x70A862F80EC4700D}, // e =   566, k =  297
        {0x0F4A642E14C6262C, 0x8CD27BB612758C10}, // e =   568, k =  298
        {0x098E7E9CCCFBD7DB, 0xD8038D51CB89778A}, // e =   571, k =  299
        {0x0BF21E44003ACDD2, 0xCE0470A63E6BD56D}, // e =   573, k =  300
        {0x0EEEA5D500498147, 0x81858CCFCE06CAC8}, // e =   575, k =  301
        {0x095527A5202DF0CC, 0xB0F37801E0C43EBD}, // e =   578, k =  302
        {0x0BAA718E68396CFF, 0xDD30560258F54E6C}, // e =   580, k =  303
        {0x0E950DF20247C83F, 0xD47C6B82EF32A207}, // e =   582, k =  304
        {0x091D28B7416CDD27, 0xE4CDC331D57FA545}, // e =   585, k =  305
        {0x0B6472E511C81471, 0xDE0133FE4ADF8E96}, // e =   587, k =  306
        {0x0E3D8F9E563A198E, 0x558180FDDD97723B}, // e =   589, k =  307
        {0x08E679C2F5E44FF8, 0xF570F09EAA7EA765}, // e =   592, k =  308
        {0x0B201833B35D63F7, 0x32CD2CC6551E513E}, // e =   594, k =  309
        {0x0DE81E40A034BCF4, 0xFF8077F7EA65E58E}, // e =   596, k =  310
        {0x08B112E86420F619, 0x1FB04AFAF27FAF79}, // e =   599, k =  311
        {0x0ADD57A27D29339F, 0x679C5DB9AF1F9B57}, // e =   601, k =  312
        {0x0D94AD8B1C738087, 0x418375281AE7822C}, // e =   603, k =  313
        {0x087CEC76F1C83054, 0x88F2293910D0B15C}, // e =   606, k =  314
        {0x0A9C2794AE3A3C69, 0xAB2EB3875504DDB3}, // e =   608, k =  315
        {0x0D433179D9C8CB84, 0x15FA60692A46151F}, // e =   610, k =  316
        {0x0849FEEC281D7F32, 0x8DBC7C41BA6BCD34}, // e =   613, k =  317
        {0x0A5C7EA73224DEFF, 0x312B9B522906C081}, // e =   615, k =  318
        {0x0CF39E50FEAE16BE, 0xFD768226B34870A1}, // e =   617, k =  319
        {0x081842F29F2CCE37, 0x5E6A1158300D4665}, // e =   620, k =  320
        {0x0A1E53AF46F801C5, 0x360495AE3C1097FE}, // e =   622, k =  321
        {0x0CA5E89B18B60236, 0x8385BB19CB14BDFD}, // e =   624, k =  322
        {0x0FCF62C1DEE382C4, 0x246729E03DD9ED7C}, // e =   626, k =  323
        {0x09E19DB92B4E31BA, 0x96C07A2C26A8346E}, // e =   629, k =  324
        {0x0C5A05277621BE29, 0x3C7098B730524189}, // e =   631, k =  325
    };

    RYU_ASSERT(k >= MinDecExp);
    RYU_ASSERT(k <= MaxDecExp);
    return Pow5[static_cast<unsigned>(k - MinDecExp)];
}

#if RYU_USE_INTRINSICS() && defined(__SIZEOF_INT128__)

static inline uint64_t MulShift(uint64_t m, const Uint64x2* mul, int j)
{
    __extension__ using uint128_t = unsigned __int128;

    RYU_ASSERT(j >= 65);
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

#elif RYU_USE_INTRINSICS() && ( defined(_MSC_VER) && defined(_M_X64) )

static inline uint64_t MulShift(uint64_t m, const Uint64x2* mul, int j)
{
    RYU_ASSERT(j >= 65);
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
    // j >= 117 and m has at most 53 + 2 = 55 bits.
    // The product along with the subsequent shift therefore requires
    // 55 + 124 - 117 = 62 bits.

    const auto k = FloorLog2Pow5(e5) + 1 - BitsPerPow5_Double;
    const auto j = e2 - k;
    RYU_ASSERT(j >= BitsPerPow5_Double - 7); // 117 - 64 = 53
    RYU_ASSERT(j <= BitsPerPow5_Double - 1); // 123 - 64 = 59

    const auto pow5 = ComputePow5_Double(e5);

    a = MulShift(u, &pow5, j);
    b = MulShift(v, &pow5, j);
    c = MulShift(w, &pow5, j);
}

// Returns whether value is divisible by 5^e5
static inline bool MultipleOfPow5(uint64_t value, int e5)
{
    RYU_ASSERT(e5 >= 0);
    RYU_ASSERT(e5 <= 22);

    struct MulCmp {
        uint64_t mul;
        uint64_t cmp;
    };

    static constexpr MulCmp Mod5[] = { // 368 bytes
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
    };

    return value * Mod5[e5].mul <= Mod5[e5].cmp;
}

// Returns whether value is divisible by 2^e2
static inline bool MultipleOfPow2(uint64_t value, int e2)
{
    RYU_ASSERT(e2 >= 0);
    RYU_ASSERT(e2 <= 63);

    return (value & ((uint64_t{1} << e2) - 1)) == 0;
}

struct ToDecimalResultDouble {
    uint64_t digits; // num_digits <= 17
    int exponent;
};

#if !RYU_KEEP_TRAILING_ZEROS_IN_SMALL_INT()
static inline ToDecimalResultDouble RemoveTrailingZeros64(uint64_t m2)
{
    // m2 < 2^53, which has 16 decimal digits.
    // We therefore remove at most 15 digits.
    RYU_ASSERT(m2 < 9007199254740992);

    int k = 0;
    for (;;)
    {
        const uint64_t q = m2 / 10;
        const uint32_t r = Lo32(m2) - 10 * Lo32(q);
        if (r != 0)
            break;
        m2 = q;
        ++k;
    }

    return {m2, k};
}
#endif

static inline ToDecimalResultDouble ToDecimal(double value)
{
    using Double = IEEE<double>;

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
#if RYU_KEEP_TRAILING_ZEROS_IN_SMALL_INT()
            return {m2 >> -e2, 0};
#else
            return RemoveTrailingZeros64(m2 >> -e2);
#endif
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
            // NB: q <= 22 implies e2 <= 79.
            // NB: Since w - u <= 4, only one of u, v, and w can be a multiple of 5, if any.

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
        //          = (u,v,w) * 5^(-e10) / 2^(e10 - e2),

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
            // NB: q <= 1 implies -4 <= e2 <= -1
            // NB: 2 <= q <= 55 implies -81 <= e2 <= -5

            za = MultipleOfPow2(u, q);
            zb = MultipleOfPow2(v, q);
            zc = MultipleOfPow2(w, q);
        }
    }

    uint64_t a;
    uint64_t b;
    uint64_t c;
    MulPow5DivPow2_Double(u, v, w, -e10, e10 - e2, a, b, c);

    //
    // Step 4:
    // Find the shortest decimal representation in the interval of valid
    // representations.
    //

    //
    // TODO:
    //  Use Theorem 6.2 from Floitsch's Grisu paper here?
    //  Might be at least beneficial for 32-bit platforms...
    //

    c -= !accept_upper && zc;

    const uint64_t aq = a;
    const uint64_t bq = b;

    uint64_t mask = 1;
    // mask = 10^(number of digits removed),
    // i.e., (bq % mask) contains the actual digits removed from bq.

    while (a / 10000 < c / 10000)
    {
        mask *= 10000;
        a /= 10000;
        b /= 10000;
        c /= 10000;
        e10 += 4;
    }

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
// ToDecimal
//
// Single-precision implementation
//==================================================================================================
// Constant data: 624 (+ 88) bytes

static constexpr int BitsPerPow5_Single = 64;

static inline uint64_t ComputePow5_Single(int k)
{
    // Let e = FloorLog2Pow5(k) + 1 - 64
    // For k >= 0, stores 5^k in the form: ceil( 5^k / 2^e )
    // For k <= 0, stores 5^k in the form: ceil(2^-e / 5^-k)
    static constexpr int MinDecExp = -29;
    static constexpr int MaxDecExp =  47;
    static constexpr uint64_t Pow5[MaxDecExp - MinDecExp + 1] = {
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
        0x813F3978F8940985, // e =     2, k =   28
        0xA18F07D736B90BE6, // e =     4, k =   29
        0xC9F2C9CD04674EDF, // e =     6, k =   30
        0xFC6F7C4045812297, // e =     8, k =   31
        0x9DC5ADA82B70B59E, // e =    11, k =   32
        0xC5371912364CE306, // e =    13, k =   33
        0xF684DF56C3E01BC7, // e =    15, k =   34
        0x9A130B963A6C115D, // e =    18, k =   35
        0xC097CE7BC90715B4, // e =    20, k =   36
        0xF0BDC21ABB48DB21, // e =    22, k =   37
        0x96769950B50D88F5, // e =    25, k =   38
        0xBC143FA4E250EB32, // e =    27, k =   39
        0xEB194F8E1AE525FE, // e =    29, k =   40
        0x92EFD1B8D0CF37BF, // e =    32, k =   41
        0xB7ABC627050305AE, // e =    34, k =   42
        0xE596B7B0C643C71A, // e =    36, k =   43
        0x8F7E32CE7BEA5C70, // e =    39, k =   44
        0xB35DBF821AE4F38C, // e =    41, k =   45
        0xE0352F62A19E306F, // e =    43, k =   46
        0x8C213D9DA502DE46, // e =    46, k =   47
    };

    RYU_ASSERT(k >= MinDecExp);
    RYU_ASSERT(k <= MaxDecExp);
    return Pow5[static_cast<unsigned>(k - MinDecExp)];
}

static inline uint64_t MulShift(uint32_t m, uint64_t mul, int j)
{
    RYU_ASSERT(j >= 0);
    RYU_ASSERT(j <= 63);

#if RYU_USE_INTRINSICS() && defined(__SIZEOF_INT128__)
    __extension__ using uint128_t = unsigned __int128;
    const uint64_t shifted_sum = static_cast<uint64_t>((uint128_t{mul} * m) >> (j & 63));
#elif RYU_USE_INTRINSICS() && ( defined(_MSC_VER) && defined(_M_X64) )
    uint64_t hi;
    uint64_t lo = _umul128(m, mul, &hi);
    const uint64_t shifted_sum = __shiftright128(lo, hi, static_cast<unsigned char>(j));
#else
    const uint64_t bits0 = uint64_t{m} * Lo32(mul);
    const uint64_t bits1 = uint64_t{m} * Hi32(mul);
    const uint64_t sum = bits1 + Hi32(bits0);
#if RYU_USE_INTRINSICS() && ( defined(_MSC_VER) && defined(_M_IX86) && !defined(__clang__) )
    const uint64_t shifted_sum = __ull_rshift(sum, j);
#else
    const int shift = j & 31;
    const uint64_t shifted_sum = sum >> shift;
#endif
#endif

    return shifted_sum;
}

static inline void MulPow5DivPow2_Single(uint32_t u, uint32_t v, uint32_t w, int e5, int e2, uint64_t& a, uint64_t& b, uint64_t& c)
{
    static constexpr int BitsPerPow5 = 64;

    // j >= 57 and m has at most 24 + 2 = 26 bits.
    // The product along with the subsequent shift therefore requires
    // 26 + 64 - 57 = 33 bits.

    const auto k = FloorLog2Pow5(e5) + 1 - BitsPerPow5;
    const auto j = e2 - k;
    RYU_ASSERT(j >= BitsPerPow5_Single - 7); // 57
    RYU_ASSERT(j <= BitsPerPow5_Single - 1); // 63

    const auto pow5 = ComputePow5_Single(e5);

    a = MulShift(u, pow5, j);
    b = MulShift(v, pow5, j);
    c = MulShift(w, pow5, j);
}

static inline bool MultipleOfPow5(uint32_t value, int e5)
{
    RYU_ASSERT(e5 >= 0);
    RYU_ASSERT(e5 <= 10);

    struct MulCmp {
        uint32_t mul;
        uint32_t cmp;
    };

    static constexpr MulCmp Mod5[] = { // 88 bytes
        {0x00000001u, 0xFFFFFFFFu}, // 5^0
        {0xCCCCCCCDu, 0x33333333u}, // 5^1
        {0xC28F5C29u, 0x0A3D70A3u}, // 5^2
        {0x26E978D5u, 0x020C49BAu}, // 5^3
        {0x3AFB7E91u, 0x0068DB8Bu}, // 5^4
        {0x0BCBE61Du, 0x0014F8B5u}, // 5^5
        {0x68C26139u, 0x000431BDu}, // 5^6
        {0xAE8D46A5u, 0x0000D6BFu}, // 5^7
        {0x22E90E21u, 0x00002AF3u}, // 5^8
        {0x3A2E9C6Du, 0x00000897u}, // 5^9
        {0x3ED61F49u, 0x000001B7u}, // 5^10
    };

    return value * Mod5[e5].mul <= Mod5[e5].cmp;
}

static inline bool MultipleOfPow2(uint32_t value, int e2)
{
    RYU_ASSERT(e2 >= 0);
    RYU_ASSERT(e2 <= 31);

    return (value & ((uint32_t{1} << e2) - 1)) == 0;
}

struct ToDecimalResultSingle {
    uint32_t digits; // num_digits <= 9
    int exponent;
};

#if !RYU_KEEP_TRAILING_ZEROS_IN_SMALL_INT()
static inline ToDecimalResultSingle RemoveTrailingZeros32(uint32_t m2)
{
    // m2 < 2^24, which has 8 decimal digits.
    // We therefore remove at most 7 digits.
    RYU_ASSERT(m2 < 16777216);

    int k = 0;
    for (;;)
    {
        const uint32_t q = m2 / 10;
        const uint32_t r = m2 - 10 * q;
        if (r != 0)
            break;
        m2 = q;
        ++k;
    }

    return {m2, k};
}
#endif

static inline ToDecimalResultSingle ToDecimal(float value)
{
    using Single = IEEE<float>;

    RYU_ASSERT(Single(value).IsFinite());
    RYU_ASSERT(value > 0);

    //
    // Step 1:
    // Decode the floating point number, and unify normalized and subnormal cases.
    //

    const Single ieee_value(value);

    const uint32_t ieee_mantissa = ieee_value.PhysicalSignificand();
    const uint32_t ieee_exponent = ieee_value.PhysicalExponent();

    uint32_t m2;
    int e2;
    if (ieee_exponent == 0)
    {
        m2 = ieee_mantissa;
        e2 = 1 - Single::ExponentBias;
    }
    else
    {
        m2 = Single::HiddenBit | ieee_mantissa;
        e2 = static_cast<int>(ieee_exponent) - Single::ExponentBias;

        if /*unlikely*/ ((0 <= -e2 && -e2 < Single::SignificandSize) && MultipleOfPow2(m2, -e2))
        {
#if RYU_KEEP_TRAILING_ZEROS_IN_SMALL_INT()
            return {m2 >> -e2, 0};
#else
            return RemoveTrailingZeros32(m2 >> -e2);
#endif
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
    const uint32_t u = 4 * m2 - 2 + lower_boundary_is_closer;
    const uint32_t v = 4 * m2;
    const uint32_t w = 4 * m2 + 2;

    //
    // Step 3:
    // Convert to a decimal power base.
    //

    int e10;

    bool za = false;
    bool zb = false;
    bool zc = false;

    if (e2 >= 0)
    {
        const int q = FloorLog10Pow2(e2) - (e2 > 3); // == max(0, q' - 1)
        RYU_ASSERT(q >= 0);

        e10 = q;
        RYU_ASSERT(e10 >= 0);
        RYU_ASSERT(e10 - e2 <= 0);

        if (q <= 10) // 10 = floor(log_5(2^24))
        {
            za = MultipleOfPow5(u, q);
            zb = MultipleOfPow5(v, q);
            zc = MultipleOfPow5(w, q);
        }
    }
    else
    {
        const int q = FloorLog10Pow5(-e2) - (-e2 > 1); // == max(0, q' - 1)
        RYU_ASSERT(q >= 0);

        e10 = q + e2;
        RYU_ASSERT(e10 < 0);
        RYU_ASSERT(e10 - e2 >= 0);

        if (q <= Single::SignificandSize + 2)
        {
            za = MultipleOfPow2(u, q);
            zb = MultipleOfPow2(v, q);
            zc = MultipleOfPow2(w, q);
        }
    }

    uint64_t aq;
    uint64_t bq;
    uint64_t cq;
    MulPow5DivPow2_Single(u, v, w, -e10, e10 - e2, aq, bq, cq);

    //
    // Step 4:
    // Find the shortest decimal representation in the interval of legal representations.
    //

    cq -= !accept_upper && zc;

    // c < 2^33 = 8'589'934'592,
    // and we will therefore remove at most 9 decimal digits, i.e. mask fits into an uint32_t.
    uint32_t mask = 1;

    // aq,bq,cq sometimes have 33 bits and we want to use 32-bit operations as much as
    // possible. In this case, we remove the first decimal digit and then use 32-bit
    // integers.
    //
    // TODO:
    //  Do this only for 32-bit platforms?!
    //

    uint32_t a = Lo32(aq);
    uint32_t b = Lo32(bq);
    uint32_t c = Lo32(cq);

    if (Hi32(cq) != 0)
    {
        RYU_ASSERT(aq / 10 < cq / 10);
        RYU_ASSERT(Hi32(aq / 2) == 0);
        RYU_ASSERT(Hi32(bq / 2) == 0);
        RYU_ASSERT(Hi32(cq / 2) == 0);

        mask = 10;
        a = Lo32(aq / 2) / 5; // = aq / 10
        b = Lo32(bq / 2) / 5; // = bq / 10
        c = Lo32(cq / 2) / 5; // = cq / 10
        ++e10;
    }

    while (a / 100 < c / 100)
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
        const uint32_t br = Lo32(bq) - b * mask; // Digits removed from bq
        const uint32_t half = mask / 2;

        b += (a == b || br >= half);
    }
    else
    {
        const bool can_use_lower = accept_lower && za && (Lo32(aq) - a * mask == 0);
        if (can_use_lower)
        {
            RYU_ASSERT(a != 0);
            for (;;)
            {
                const uint32_t q = a / 10;
                const uint32_t r = a - 10 * q;
                if (r != 0)
                    break;
                mask *= 10;
                a = q;
                b = q;
//              c = q;
                ++e10;
            }
        }

        const uint32_t br = Lo32(bq) - b * mask; // Digits removed from bq
        const uint32_t half = mask / 2;

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

static inline int DecimalLength(uint32_t v)
{
    RYU_ASSERT(v >= 1);
    RYU_ASSERT(v <= 999999999);

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

static inline void PrintDecimalDigits(char* buf, uint32_t output, int output_length)
{
    while (output >= 10000)
    {
        RYU_ASSERT(output_length > 4);
        const uint32_t q = output / 10000;
        const uint32_t r = output % 10000;
        output = q;
        output_length -= 4;
        Utoa_4Digits(buf + output_length, r);
    }

    if (output >= 100)
    {
        RYU_ASSERT(output_length > 2);
        const uint32_t q = output / 100;
        const uint32_t r = output % 100;
        output = q;
        output_length -= 2;
        Utoa_2Digits(buf + output_length, r);
    }

    if (output >= 10)
    {
        RYU_ASSERT(output_length == 2);
        Utoa_2Digits(buf, output);
    }
    else
    {
        RYU_ASSERT(output_length == 1);
        buf[0] = static_cast<char>('0' + output);
    }
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

template <typename UnsignedInt>
static inline char* FormatDigits(char* buffer, UnsignedInt digits, int decimal_exponent, bool force_trailing_dot_zero = false)
{
    //
    // TODO:
    // Specialized version for single-precision???
    //
    // constexpr bool is_single = sizeof(UnsignedInt) == sizeof(uint32_t);

    RYU_ASSERT(digits >= 1);
    RYU_ASSERT(digits <= 99999999999999999ull);
    RYU_ASSERT(decimal_exponent >= -999);
    RYU_ASSERT(decimal_exponent <=  999);

    const int num_digits = DecimalLength(digits);
    const int decimal_point = num_digits + decimal_exponent;

#if !RYU_SCIENTIFIC_NOTATION_ONLY()
    // single-precision: MaxDigits10 = 9, MaxIntLength = 8. And MaxAdditionalZeros = 4?
    constexpr int MaxIntLength = 16; // 2^53 = 9'007'199'254'740'992
#if 1
    constexpr int MaxAdditionalZeros = 1;
#else
    constexpr int MaxAdditionalZeros = 5;
#endif
    constexpr int MaxFixedDecimalPoint = MaxIntLength + MaxAdditionalZeros; //   digits[000]
    constexpr int MinFixedDecimalPoint = -MaxAdditionalZeros;               // 0.[000]digits

    const bool use_fixed = MinFixedDecimalPoint <= decimal_point && decimal_point <= MaxFixedDecimalPoint;
#endif

    // Prepare the buffer.
    // Avoid calling memset/memcpy with variable arguments below...

    int decimal_digits_position;
#if !RYU_SCIENTIFIC_NOTATION_ONLY()
    if (use_fixed)
    {
        if (decimal_point <= 0)
        {
            // 0.[000]digits
            // -5 <= decimal_point <= 0
            //  ==> 2 <= 2 + -decimal_point <= 7
            // Pre-filling the buffer with 7 '0's is therefore sufficient.
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
            // 1 <= num_digits <= 17 <= decimal_point <= 21.
            // Pre-filling buffer with 21 '0's is therefore sufficient.
            std::memset(buffer, '0', 24); // sp: 12
            decimal_digits_position = 0;
        }
    }
    else
#endif
    {
        // dE+123 or d.igitsE+123
        // We only need to copy the first digit one position to the left.
        decimal_digits_position = 1;
    }

    PrintDecimalDigits(buffer + decimal_digits_position, digits, num_digits);

#if !RYU_SCIENTIFIC_NOTATION_ONLY()
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
            std::memmove(buffer + (decimal_point + 1), buffer + decimal_point, 16); // sp: 8
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
#endif
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
        else if (/* sp || */ k < 100)
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

template <typename Float>
static inline char* ToChars(char* buffer, Float value, bool force_trailing_dot_zero = false)
{
    using Fp = IEEE<Float>;
    const Fp v(value);

    if (!v.IsFinite())
    {
        if (v.IsNaN())
        {
            std::memcpy(buffer, "NaN ", 4);
            return buffer + 3;
        }
        if (v.SignBit())
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

char* RyuDtoa(char* buffer, double value)
{
    return ToChars(buffer, value);
}

char* RyuFtoa(char* buffer, float value)
{
    return ToChars(buffer, value);
}
