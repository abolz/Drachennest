#include "benchmark/benchmark.h"

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <random>
#include <vector>
#include <chrono>

#include <math.h>

//#define BENCH_CHARCONV 1
//#define BENCH_GRISU2 1
//#define BENCH_GRISU3 1
#define BENCH_RYU 1
//#define BENCH_RYU_UPSTREAM 1

#define BENCH_SINGLE 1
#define BENCH_DOUBLE 1
#define BENCH_TO_DECIMAL 0

//==================================================================================================
//
//==================================================================================================

#if BENCH_GRISU2
#include "grisu2.h"
struct D2S
{
    static char const* Name() { return "Grisu2"; }
    char* operator()(char* buf, int /*buflen*/, float f) const { return grisu2::ToChars(buf, f); }
    char* operator()(char* buf, int /*buflen*/, double f) const { return grisu2::ToChars(buf, f); }
};
#endif
#if BENCH_GRISU3
#include "grisu3.h"
struct D2S
{
    static char const* Name() { return "Grisu3 (Dragon4)"; }
    char* operator()(char* buf, int /*buflen*/, float f) const { return grisu3::ToChars(buf, f); }
    char* operator()(char* buf, int /*buflen*/, double f) const { return grisu3::ToChars(buf, f); }
};
#endif
#if BENCH_RYU
#if 1
#include "ryu_charconv.h"
struct D2S
{
    static char const* Name() { return "Ryu"; }
    char* operator()(char* buf, int /*buflen*/, float f) const { return RyuFtoa(buf, f); }
    char* operator()(char* buf, int /*buflen*/, double f) const { return RyuDtoa(buf, f); }
};
#else
#include "ryu.h"
struct D2S
{
    static char const* Name() { return "Ryu"; }
    char* operator()(char* buf, int /*buflen*/, float f) const { return ryu::ToChars(buf, f); }
    char* operator()(char* buf, int /*buflen*/, double f) const { return ryu::ToChars(buf, f); }
};
#endif
#endif
#if BENCH_RYU_UPSTREAM
#include "ryu/ryu/ryu.h"
struct D2S
{
    static char const* Name() { return "Ryu (upstream)"; }
    char* operator()(char* buf, int /*buflen*/, float f) const { return buf + f2s_buffered_n(f, buf); }
    char* operator()(char* buf, int /*buflen*/, double f) const { return buf + d2s_buffered_n(f, buf); }
};
#endif
#if BENCH_CHARCONV
#include <charconv>
struct D2S
{
    static char const* Name() { return "charconv"; }
#if 0
    char* operator()(char* buf, int buflen, float f) const { return std::to_chars(buf, buf + buflen, f, std::chars_format::general).ptr; }
    char* operator()(char* buf, int buflen, double f) const { return std::to_chars(buf, buf + buflen, f, std::chars_format::general).ptr; }
#else
    char* operator()(char* buf, int buflen, float f) const { return std::to_chars(buf, buf + buflen, f).ptr; }
    char* operator()(char* buf, int buflen, double f) const { return std::to_chars(buf, buf + buflen, f).ptr; }
#endif
};
#endif

//==================================================================================================
//
//==================================================================================================

template <typename Target, typename Source>
static Target ReinterpretBits(Source const& source)
{
    static_assert(sizeof(Target) == sizeof(Source), "size mismatch");

    Target target;
    std::memcpy(&target, &source, sizeof(Source));
    return target;
}

class JenkinsRandom
{
    // A small noncryptographic PRNG
    // http://burtleburtle.net/bob/rand/smallprng.html

    uint32_t a;
    uint32_t b;
    uint32_t c;
    uint32_t d;

    static uint32_t Rotate(uint32_t value, int n) {
        return (value << n) | (value >> (32 - n));
    }

    uint32_t Gen() {
        const uint32_t e = a - Rotate(b, 27);
        a = b ^ Rotate(c, 17);
        b = c + d;
        c = d + e;
        d = e + a;
        return d;
    }

public:
    using result_type = uint32_t;

    static constexpr uint32_t min() { return 0; }
    static constexpr uint32_t max() { return UINT32_MAX; }

    explicit JenkinsRandom(uint32_t seed = 0) {
        a = 0xF1EA5EED;
        b = seed;
        c = seed;
        d = seed;
        for (int i = 0; i < 20; ++i) {
            static_cast<void>(Gen());
        }
    }

    uint32_t operator()() { return Gen(); }
};

static JenkinsRandom random;

//==================================================================================================
//
//==================================================================================================

static constexpr int BufSize = 64;
static constexpr int NumFloats = 1 << 13;

//template <typename Float>
//static inline void PrintFloat(Float v)
//{
//    char buf[32];
//    char* end = D2S_DoubleConversion{}(buf, 32, v);
//    end[0] = '\0';
//    printf("%s\n", buf);
//}

template <typename ...Args>
static inline char const* StrPrintf(char const* format, Args&&... args)
{
    char buf[1024];
    snprintf(buf, 1024, format, std::forward<Args>(args)...);
#ifdef _MSC_VER
    return _strdup(buf); // leak...
#else
    return strdup(buf); // leak...
#endif
}

#if BENCH_TO_DECIMAL
#if BENCH_GRISU2
inline uint32_t ToDecimal(int& exponent, float  value) { const auto dec = grisu2::ToDecimal(value); /*exponent = dec.exponent;*/ return dec.digits; }
inline uint64_t ToDecimal(int& exponent, double value) { const auto dec = grisu2::ToDecimal(value); /*exponent = dec.exponent;*/ return dec.digits; }
#endif
#if BENCH_GRISU3
inline uint32_t ToDecimal(int& exponent, float  value) { const auto dec = grisu3::ToDecimal(value); /*exponent = dec.exponent;*/ return dec.digits; }
inline uint64_t ToDecimal(int& exponent, double value) { const auto dec = grisu3::ToDecimal(value); /*exponent = dec.exponent;*/ return dec.digits; }
#endif
#if BENCH_RYU
inline uint32_t ToDecimal(int& exponent, float  value) { const auto dec = ryu::ToDecimal(value); /*exponent = dec.exponent;*/ return dec.digits; }
inline uint64_t ToDecimal(int& exponent, double value) { const auto dec = ryu::ToDecimal(value); /*exponent = dec.exponent;*/ return dec.digits; }
#endif

template <typename, typename Float>
static inline void BenchIt(benchmark::State& state, std::vector<Float> const& numbers)
{
    int index = 0;

    uint64_t sum = 0;
    for (auto _ : state)
    {
        int exponent;
        sum += ToDecimal(exponent, numbers[index]) & 0xFF;
        index = (index + 1) & (NumFloats - 1);
    }

    if (sum == UINT64_MAX)
        abort();
}
#else
template <typename D2S, typename Float>
static inline void BenchIt(benchmark::State& state, std::vector<Float> const& numbers)
{
    D2S d2s;

    int index = 0;

    uint64_t sum = 0;
    for (auto _ : state)
    {
        char buffer[BufSize];
        d2s(buffer, BufSize, numbers[index]);
        sum += static_cast<unsigned char>(buffer[0]);
        index = (index + 1) & (NumFloats - 1);
    }

    if (sum == UINT64_MAX)
        abort();
}
#endif

template <typename Float>
static inline void RegisterBenchmarks(char const* name, std::vector<Float> const& numbers)
{
    const char* float_name = sizeof(Float) == 4 ? "single" : "double";
    auto* bench = benchmark::RegisterBenchmark(StrPrintf("%s - %s   ", float_name, name), BenchIt<D2S, Float>, numbers);

    bench->ComputeStatistics("min", [](const std::vector<double>& v) -> double {
        return *(std::min_element(std::begin(v), std::end(v)));
    });
    bench->Repetitions(3);
    bench->ReportAggregatesOnly();
}

//----------------------------------------------------------------------------------------------------------
//
//----------------------------------------------------------------------------------------------------------
static inline void Register_RandomBits_double()
{
    std::vector<double> numbers(NumFloats);

    std::uniform_int_distribution<uint64_t> gen(1, 0x7FF0000000000000ull - 1);
    std::generate(numbers.begin(), numbers.end(), [&] { return ReinterpretBits<double>(gen(random)); });

    RegisterBenchmarks("Random-bits", numbers);
}

static inline void Register_RandomBits_single()
{
    std::vector<float> numbers(NumFloats);

    std::uniform_int_distribution<uint32_t> gen(1, 0x7F800000u - 1);
    std::generate(numbers.begin(), numbers.end(), [&] { return ReinterpretBits<float>(gen(random)); });

    RegisterBenchmarks("Random-bits", numbers);
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------
template <typename Float>
static inline void Register_Uniform(Float low, Float high)
{
    std::vector<Float> numbers(NumFloats);

    std::uniform_real_distribution<Float> gen(low, high);
    std::generate(numbers.begin(), numbers.end(), [&] { return gen(random); });

    RegisterBenchmarks(StrPrintf("Uniform %.1g/%.1g", low, high), numbers);
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------
static constexpr int64_t kPow10_i64[] = {
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
    1000000000000000000,
};
static constexpr double kPow10_f64[] = {
    1.0e+00,
    1.0e+01,
    1.0e+02,
    1.0e+03,
    1.0e+04,
    1.0e+05,
    1.0e+06,
    1.0e+07,
    1.0e+08,
    1.0e+09,
    1.0e+10,
    1.0e+11,
    1.0e+12,
    1.0e+13,
    1.0e+14,
    1.0e+15,
    1.0e+16,
    1.0e+17,
    1.0e+18,
    1.0e+19,
    1.0e+20,
    1.0e+21,
    1.0e+22,
};

static inline void Register_Digits_double(const char* name, int digits, int e10)
{
    assert(digits >= 1);
    assert(digits <= 18);
    assert(e10 >= -22);
    assert(e10 <= 22);

    std::vector<double> numbers(NumFloats);

    std::uniform_int_distribution<int64_t> gen(kPow10_i64[digits - 1], kPow10_i64[digits] - 1);

    std::generate(numbers.begin(), numbers.end(), [&] {
        int64_t n = gen(random);
        if (n % 10 == 0)
            n |= 1;
        double v = static_cast<double>(n);
        if (e10 < 0)
            v /= kPow10_f64[-e10];
        else
            v *= kPow10_f64[e10];
        //PrintFloat(v);
        return v;
    });

    RegisterBenchmarks(name, numbers);
}

static constexpr int32_t kPow10_i32[] = {
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
};
static constexpr float kPow10_f32[] = {
    1.0e+00f,
    1.0e+01f,
    1.0e+02f,
    1.0e+03f,
    1.0e+04f,
    1.0e+05f,
    1.0e+06f,
    1.0e+07f,
    1.0e+08f,
    1.0e+09f,
    1.0e+10f,
};

static inline void Register_Digits_single(const char* name, int digits, int e10)
{
    assert(digits >= 1);
    assert(digits <= 9);
    assert(e10 >= -10);
    assert(e10 <= 10);

    std::vector<float> numbers(NumFloats);

    std::uniform_int_distribution<int32_t> gen(kPow10_i32[digits - 1], kPow10_i32[digits] - 1);

    std::generate(numbers.begin(), numbers.end(), [&] {
        int32_t n = gen(random);
        if (n % 10 == 0)
            n |= 1;
        float v = static_cast<float>(n);
        if (e10 < 0)
            v /= kPow10_f32[-e10];
        else
            v *= kPow10_f32[e10];
        //PrintFloat(v);
        return v;
    });

    RegisterBenchmarks(name, numbers);
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

int main(int argc, char** argv)
{
#if defined(__clang__)
    printf("Clang %d.%d\n", __clang_major__, __clang_minor__);
#elif defined(__GNUC__)
    printf("Gcc %s\n", __VERSION__);
#elif defined(_MSC_VER)
    printf("Msc %d\n", _MSC_FULL_VER);
#endif

    printf("Benchmarking %s\n", D2S::Name());

#if 1
#if BENCH_DOUBLE
    Register_RandomBits_double();
    Register_Uniform(0.0, 1.0);
    Register_Uniform(0.0, 1.0e+308);

    //for (int e = 10; e < 20; ++e) {
    //    Register_Uniform(std::pow(10.0, e), std::pow(10.0, e+1));
    //}

#if 1
#if 0
    for (int d = 1; d <= 15; ++d) {
//  for (int d = 15; d >= 1; --d) {
        Register_Digits_double(StrPrintf("1.%d-digits", d - 1), d, -(d - 1));
    }
    for (int d = 1; d <= 15; ++d) {
        Register_Digits_double(StrPrintf("%d.1-digits", d - 1), d, -1);
    }
    for (int d = 1; d <= 15; ++d) {
        Register_Digits_double(StrPrintf("%d-digits / 10^22", d), d, -22);
    }
    for (int d = 1; d <= 15; ++d) {
        Register_Digits_double(StrPrintf("%d-digits * 10^22", d), d, 22);
    }
#else
    for (int d = 1; d <= 18; ++d) {
		if (d != 1 && d != 18) continue;
    //for (int d = 18; d >= 1; --d) {
        //for (int e = 22; e >= -22; e -= 1) {
        for (int e = -22; e <= 22; e += 1) {
            Register_Digits_double(StrPrintf("%2d,%3d", d, e), d, e);
        }
    }
#endif
#endif
#endif

#if BENCH_SINGLE
    Register_RandomBits_single();
    Register_Uniform(0.0f, 1.0f);
    Register_Uniform(0.0f, 1.0e+38f);

#if 0
#if 0
    for (int d = 1; d <= 7; ++d) {
        Register_Digits_single(StrPrintf("1.%d-digits", d - 1), d, -(d - 1));
    }
    for (int d = 1; d <= 7; ++d) {
        Register_Digits_single(StrPrintf("%d.1-digits", d - 1), d, -1);
    }
    for (int d = 1; d <= 7; ++d) {
        Register_Digits_single(StrPrintf("%d-digits / 10^10", d), d, -10);
    }
    for (int d = 1; d <= 7; ++d) {
        Register_Digits_single(StrPrintf("%d-digits * 10^10", d), d, 10);
    }
#else
    for (int d = 1; d <= 9; ++d) {
        for (int e = -10; e <= 10; e += 1) {
            Register_Digits_single(StrPrintf("%2d,%3d", d, e), d, e);
        }
    }
#endif
#endif
#endif
#endif

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}
