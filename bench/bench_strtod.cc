#include "benchmark/benchmark.h"

#include <cstring>

#include <algorithm>
#include <random>
#include <string>

#include "ryu_64.h"

#define BENCH_RYU()                 1
#define BENCH_STD_STRTOD()          0
#define BENCH_STD_CHARCONV()        0
#define BENCH_DOUBLE_CONVERSION()   0

static constexpr int NumFloats = 1 << 14;

#if BENCH_RYU()
struct S2DRyu
{
    using value_type = double;

    value_type operator()(std::string const& str) const
    {
        value_type flt = 0;
        const auto res = ryu::Strtod(str.data(), str.data() + str.size(), flt);
        assert(res.status != ryu::StrtodStatus::invalid);
        return flt;
    }
};
#endif

#if BENCH_STD_STRTOD()
struct S2DStdStrtod
{
    using value_type = double;

    value_type operator()(std::string const& str) const
    {
        value_type flt = std::strtod(str.c_str(), nullptr);
        return flt;
    }
};
#endif

#if BENCH_STD_CHARCONV()
#include <charconv>
struct S2DStdCharconv
{
    using value_type = double;

    value_type operator()(std::string const& str) const
    {
        value_type flt = 0;
        const bool ok = std::from_chars(str.data(), str.data() + str.size(), flt).ec == std::errc{};
        assert(ok);
        return flt;
    }
};
#endif

#if BENCH_DOUBLE_CONVERSION()
#include "double-conversion/double-conversion.h"
struct S2DDoubleConversion
{
    using value_type = double;

    value_type operator()(std::string const& str) const
    {
        double_conversion::StringToDoubleConverter s2d(0, 0.0, 1.0, "inf", "nan");
        int processed_characters_count = 0;
        return s2d.StringToDouble(str.data(), static_cast<int>(str.size()), &processed_characters_count);
    }
};
#endif

template <typename Converter>
static void BenchIt(benchmark::State& state, std::vector<std::string> const& numbers)
{
    Converter convert;

    size_t index = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize( convert(numbers[index]) );
        index = (index + 1) & (NumFloats - 1);
    }
}

template <typename Converter>
static void RegisterBenchmarks(char const* name, std::vector<std::string> const& numbers)
{
    auto* bench = benchmark::RegisterBenchmark(name, BenchIt<Converter>, numbers);

    bench->ComputeStatistics("min", [](std::vector<double> const& v) -> double {
        return *(std::min_element(std::begin(v), std::end(v)));
    });
    //bench->Repetitions(3);
    bench->ReportAggregatesOnly();
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

static inline void RegisterUniform_double(char const* name, double min, double max)
{
    std::vector<std::string> numbers(NumFloats);

    std::uniform_real_distribution<double> gen(min, max);

    std::generate(numbers.begin(), numbers.end(), [&] {
        char buf[128];

        char* const end = ryu::Dtoa(buf, gen(random));
        //char* const end = buf + std::snprintf(buf, 128, "%.17g", gen(random));
        //char* const end = buf + std::snprintf(buf, 128, "%.19g", gen(random));
        //char* const end = buf + std::snprintf(buf, 128, "%.20g", gen(random));

        return std::string(buf, end);
    });

#if BENCH_RYU()
    RegisterBenchmarks<S2DRyu             >(StrPrintf("%s Ryu               ", name), numbers);
#endif
#if BENCH_STD_STRTOD()
    RegisterBenchmarks<S2DStdStrtod       >(StrPrintf("%s std::strtod       ", name), numbers);
#endif
#if BENCH_STD_CHARCONV()
    RegisterBenchmarks<S2DStdCharconv     >(StrPrintf("%s std::charconv     ", name), numbers);
#endif
#if BENCH_DOUBLE_CONVERSION()
    RegisterBenchmarks<S2DDoubleConversion>(StrPrintf("%s double_conversion ", name), numbers);
#endif
}

int main(int argc, char** argv)
{
#if defined(__clang__)
    printf("clang %d.%d\n", __clang_major__, __clang_minor__);
#elif defined(__GNUC__)
    printf("gcc %s\n", __VERSION__);
#elif defined(_MSC_VER)
    printf("msc %d\n", _MSC_FULL_VER);
#endif

    RegisterUniform_double("warm up", 0, 1);
    RegisterUniform_double("warm up", 0, 1);
    RegisterUniform_double("warm up", 0, 1);

    RegisterUniform_double("uniform [0,1/2]", 0.0, 0.5);
    RegisterUniform_double("uniform [1/4,1/2]", 0.25, 0.5);
    RegisterUniform_double("uniform [1/2,1]", 0.5, 1.0);
    RegisterUniform_double("uniform [0,1]", 0.0, 1.0);
    RegisterUniform_double("uniform [1,2]", 1.0, 2.0);
    RegisterUniform_double("uniform [2,4]", 2.0, 4.0);
    RegisterUniform_double("uniform [4,8]", 4.0, 8.0);
    RegisterUniform_double("uniform [8,2^10]", 8.0, 1ll << 10);
    RegisterUniform_double("uniform [2^10,2^20]", 1ll << 10, 1ll << 20);
    RegisterUniform_double("uniform [2^20,2^50]", 1ll << 20, 1ll << 50);
    RegisterUniform_double("uniform [0,max]", 0.0, std::numeric_limits<double>::max());

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    return 0;
}
