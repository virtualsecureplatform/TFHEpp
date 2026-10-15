#include <modular_rns32.hpp>
#include <ckks/narrow_matrix.hpp>
#include <chrono>
#include <iostream>
#include <random>

using namespace TFHEpp;
using Clock = std::chrono::steady_clock;

int main()
{
    std::mt19937_64 random(20262119);
    std::size_t checked = 0;
    for (const auto prime : modular_rns32::primes) {
        const auto p = prime.value;
        for (std::uint64_t divisor = 2; divisor*divisor <= p; ++divisor)
            if (p % divisor == 0) throw std::runtime_error("small NTT modulus is composite");
        modular_ntt32::Montgomery arithmetic(p);
        for (std::size_t i = 0; i < 100000; ++i) {
            const auto wide = i == 0 ? (std::uint64_t{1} << 52)-1 :
                i == 1 ? (std::uint64_t{1} << 60)-1 :
                i == 2 ? (p << 32)-1 : random() % (std::uint64_t{1} << 60);
            const auto encoded = arithmetic.encode(wide);
            if (arithmetic.decode(encoded) != wide % p) throw std::runtime_error("wide Montgomery encoding mismatch");
            const auto input = i == 0 ? p-1 : random() % p;
            if (arithmetic.multiply(input, encoded) != static_cast<unsigned __int128>(input)*wide % p)
                throw std::runtime_error("Montgomery multiplication mismatch");
            ++checked;
        }
        for (std::size_t columns : {1, 2, 3, 4, 5, 8, 9, 12, 17}) {
            std::vector<std::uint32_t> left(columns), right(columns), first(columns), second(columns);
            std::vector<std::uint64_t> a(columns), b(columns);
            for (std::size_t repetition = 0; repetition < 1024; ++repetition) {
                auto selector = repetition ? random() % p : p-1;
                for (std::size_t i = 0; i < columns; ++i) {
                    left[i] = repetition ? random() % p : p-1;
                    right[i] = repetition ? random() % p : 0;
                    a[i] = (a[i]+static_cast<unsigned __int128>(left[i])*selector) % p;
                    b[i] = (b[i]+static_cast<unsigned __int128>(right[i])*selector) % p;
                }
                arithmetic.addPair(first.data(), second.data(), left.data(), right.data(), columns, arithmetic.encode(selector));
                for (std::size_t i = 0; i < columns; ++i) {
                    if (first[i] >= 2*p || second[i] >= 2*p || first[i] % p != a[i] || second[i] % p != b[i])
                        throw std::runtime_error("lazy paired matrix differs from integer oracle");
                    ++checked;
                }
            }
        }
        for(std::size_t columns:{1,2,3,4,5,8,9,12,17}){
            std::array<std::vector<std::uint32_t>,4> left,right;
            std::array<const std::uint32_t*,4> lp,rp;
            for(std::size_t t=0;t<4;++t){left[t].resize(columns);right[t].resize(columns);lp[t]=left[t].data();rp[t]=right[t].data();}
            std::vector<std::uint32_t> first(columns),second(columns);
            std::vector<std::uint64_t> a(columns),b(columns);
            for(std::size_t repeat=0;repeat<1024;++repeat){
                std::array<std::uint32_t,4> selectors;
                for(std::size_t t=0;t<4;++t){
                    // First trial forces the largest raw products, including
                    // a Montgomery residue p-1 rather than its encoding.
                    selectors[t]=repeat?random()%p:p-1;
                    for(std::size_t c=0;c<columns;++c){
                        left[t][c]=repeat?random()%p:p-1;right[t][c]=repeat?random()%p:0;
                        const auto ordinary=arithmetic.decode(selectors[t]);
                        a[c]=(a[c]+static_cast<unsigned __int128>(left[t][c])*ordinary)%p;
                        b[c]=(b[c]+static_cast<unsigned __int128>(right[t][c])*ordinary)%p;
                    }
                }
                arithmetic.addDot4Pair(first.data(),second.data(),lp,rp,selectors,columns);
                for(std::size_t c=0;c<columns;++c){
                    if(first[c]>=2*p || second[c]>=2*p || first[c]%p!=a[c] || second[c]%p!=b[c])throw std::runtime_error("four-term Montgomery dot differs from integer oracle");
                    ++checked;
                }
            }
        }
        for (std::size_t n : {2, 8, 512, 16384, 65536}) {
            modular_ntt32::NegacyclicNTTPlan small(n, prime);
            modular_ntt::NegacyclicNTTPlan reference(n, prime);
            std::vector<std::uint32_t> data(n), encoded(n);
            std::vector<std::uint64_t> expected(n);
            for (std::size_t i = 0; i < n; ++i) expected[i] = data[i] = random() % p;
            const auto original = data;
            for (std::size_t i = 0; i < n; ++i) encoded[i] = arithmetic.encode(data[i]);
            small.forward(data); small.forward(encoded); reference.forward(expected);
            for (std::size_t i = 0; i < n; ++i)
                if (data[i] != expected[i] || arithmetic.decode(encoded[i]) != expected[i])
                    throw std::runtime_error("small-prime NTT natural order differs from reference");
            auto reversed=original, reversed_encoded=original;
            for(auto&v:reversed_encoded)v=arithmetic.encode(v);
            small.forwardBitReversed(reversed);
            small.forwardBitReversed(reversed_encoded);
            for(std::size_t i=0;i<n;++i) {
                std::size_t index=0,bits=i;
                for(std::size_t width=n;width>1;width>>=1){index=(index<<1)|(bits&1);bits>>=1;}
                if(reversed[index]!=expected[i] || arithmetic.decode(reversed_encoded[index])!=expected[i])
                    throw std::runtime_error("bit-reversed NTT differs from independent natural-order reference");
            }
            small.inverseBitReversed(reversed);
            small.inverseBitReversed(reversed_encoded);
            if(reversed!=original)throw std::runtime_error("bit-reversed NTT roundtrip mismatch");
            for(std::size_t i=0;i<n;++i)if(arithmetic.decode(reversed_encoded[i])!=original[i])
                throw std::runtime_error("Montgomery bit-reversed roundtrip mismatch");
            // Exercise zero and extremal canonical residues as well as random inputs.
            for(auto edge:{std::uint32_t{0},std::uint32_t(p-1)}) {
                std::vector<std::uint32_t> values(n,edge),natural=values;
                small.forward(natural);small.forwardBitReversed(values);
                for(std::size_t i=0;i<n;++i) {
                    std::size_t index=0,bits=i;
                    for(std::size_t width=n;width>1;width>>=1){index=(index<<1)|(bits&1);bits>>=1;}
                    if(values[index]!=natural[i])throw std::runtime_error("bit-reversed endpoint spectrum differs from reference");
                }
                small.inverseBitReversed(values);
                if(values!=std::vector<std::uint32_t>(n,edge))throw std::runtime_error("bit-reversed endpoint roundtrip mismatch");
            }
            std::vector<std::uint32_t> short_input(n-1);
            for(bool inverse:{false,true}) {
                bool rejected=false;
                try {
                    if(inverse)small.inverseBitReversed(short_input);
                    else small.forwardBitReversed(short_input);
                } catch(const std::invalid_argument&) {rejected=true;}
                if(!rejected)throw std::runtime_error("bit-reversed NTT accepted invalid length");
            }
            small.inverse(data);
            if (data != original) throw std::runtime_error("small-prime NTT roundtrip mismatch");
            checked += n;
        }
    }
    const modular_ntt32::FourPrimeCRT crt;
    for (std::size_t i = 0; i < 100000; ++i) {
        const auto magnitude = i == 0 ? modular_ntt32::crt_product/2 :
            ((static_cast<unsigned __int128>(random()) << 64) | random()) % (modular_ntt32::crt_product/2);
        const auto value = i & 1 ? -static_cast<__int128>(magnitude) : static_cast<__int128>(magnitude);
        std::array<std::uint32_t, 4> residues;
        for (std::size_t j = 0; j < 4; ++j) {
            const auto reduced = magnitude % modular_ntt32::primes[j].value;
            residues[j] = value < 0 && reduced ? modular_ntt32::primes[j].value-reduced : reduced;
        }
        if (crt.reconstructSigned(residues) != value) throw std::runtime_error("four-prime signed CRT mismatch");
        ++checked;
    }
    // Exact signed negacyclic convolution, independently accumulated in int128.
    constexpr std::size_t degree = 64;
    std::array<std::int64_t, degree> left, right;
    std::array<__int128, degree> expected{};
    for (std::size_t i = 0; i < degree; ++i) {
        left[i] = std::int64_t(random() % (std::uint64_t{1} << 52)) - (std::int64_t{1} << 51);
        right[i] = std::int64_t(random() % (std::uint64_t{1} << 42)) - (std::int64_t{1} << 41);
    }
    for (std::size_t i = 0; i < degree; ++i)
        for (std::size_t j = 0; j < degree; ++j)
            expected[(i+j)%degree] += (i+j < degree ? 1 : -1) * static_cast<__int128>(left[i])*right[j];
    std::array<std::vector<std::uint32_t>, 4> products;
    for (std::size_t prime = 0; prime < 4; ++prime) {
        const auto p = modular_ntt32::primes[prime].value;
        modular_ntt32::NegacyclicNTTPlan plan(degree, modular_ntt32::primes[prime]);
        auto residue = [&](std::int64_t x) { return std::uint32_t(x < 0 ? (p-std::uint64_t(-x)%p)%p : std::uint64_t(x)%p); };
        std::vector<std::uint32_t> a(degree), b(degree);
        for (std::size_t i = 0; i < degree; ++i) { a[i] = plan.arithmetic().encode(residue(left[i])); b[i] = residue(right[i]); }
        auto reversed_a=a,reversed_b=b;
        plan.forwardBitReversed(reversed_a);plan.forwardBitReversed(reversed_b);
        for(std::size_t i=0;i<degree;++i)reversed_b[i]=plan.arithmetic().multiply(reversed_a[i],reversed_b[i]);
        plan.inverseBitReversed(reversed_b);
        plan.forward(a); plan.forward(b);
        for (std::size_t i = 0; i < degree; ++i) b[i] = plan.arithmetic().multiply(a[i], b[i]);
        plan.inverse(b);
        if(b!=reversed_b)throw std::runtime_error("permutation-free convolution differs from natural-order convolution");
        products[prime] = std::move(b);
    }
    for (std::size_t i = 0; i < degree; ++i)
        if (crt.reconstructSigned({products[0][i],products[1][i],products[2][i],products[3][i]}) != expected[i])
            throw std::runtime_error("signed NTT convolution differs from integer oracle");
    std::cout << "{\"oracle_checks\":" << checked+degree << ",\"passed\":true}" << std::endl;

    std::vector<modular_ntt::NegacyclicNTTPlan> wide_plans;
    std::vector<modular_ntt32::NegacyclicNTTPlan> small_plans;
    std::array<std::vector<std::uint64_t>, 2> wide_data;
    std::array<std::vector<std::uint32_t>, 4> small_data;
    for (std::size_t i = 0; i < 4; ++i) {
        small_plans.emplace_back(16384,modular_ntt32::primes[i]);
        small_data[i].resize(16384);
        for (auto &v : small_data[i]) v = random() % modular_ntt32::primes[i].value;
        if (i < 2) {
            wide_plans.emplace_back(16384,modular_ntt::degree65536_primes[i]);
            wide_data[i].resize(16384);
            for (auto &v : wide_data[i]) v = random() % modular_ntt::degree65536_primes[i].value;
        }
    }
    const auto original_small = small_data;
    const auto original_wide = wide_data;
    for (bool small : {false, true, true, false}) {
        const auto start = Clock::now();
        for (std::size_t repetition = 0; repetition < 64; ++repetition)
            if (small)
                for (std::size_t i = 0; i < 4; ++i) { small_plans[i].forward(small_data[i]); small_plans[i].inverse(small_data[i]); }
            else
                for (std::size_t i = 0; i < 2; ++i) { wide_plans[i].forward(wide_data[i]); wide_plans[i].inverse(wide_data[i]); }
        if (small_data != original_small || wide_data != original_wide) throw std::runtime_error("timed NTT roundtrip failed");
        std::cout << "{\"small_primes\":" << (small ? "true" : "false")
                  << ",\"transform_seconds\":" << std::chrono::duration<double>(Clock::now()-start).count() << "}" << std::endl;
    }

    constexpr std::size_t n = 16384, columns = 4, selectors = 32, rows = 16;
    std::vector<std::uint64_t> source(rows*selectors*n), tables(2*n*columns);
    for (auto &v : source) v = random() % (std::uint64_t{1} << 52);
    for (auto &v : tables) v = random() % (std::uint64_t{1} << 42);
    std::vector<std::uint64_t> wide_source(2*source.size()), wide_tables(2*tables.size());
    for (std::size_t prime = 0; prime < 2; ++prime) {
        std::copy(source.begin(),source.end(),wide_source.begin()+prime*source.size());
        std::copy(tables.begin(),tables.end(),wide_tables.begin()+prime*tables.size());
    }
    std::vector<std::uint32_t> narrow_source(4*source.size()), narrow_tables(4*tables.size());
    for (std::size_t prime = 0; prime < 4; ++prime) {
        modular_ntt32::Montgomery arithmetic(modular_ntt32::primes[prime].value);
        for (std::size_t i = 0; i < source.size(); ++i) narrow_source[prime*source.size()+i] = arithmetic.encode(source[i]);
        for (std::size_t i = 0; i < tables.size(); ++i) narrow_tables[prime*tables.size()+i] = tables[i] % arithmetic.modulus();
    }
    std::vector<std::uint64_t> wide_sums(2*2*rows*n*columns);
    std::vector<std::uint32_t> narrow_sums(4*2*rows*n*columns);
    for (bool narrow : {false, true, true, false}) {
        if (narrow) std::fill(narrow_sums.begin(), narrow_sums.end(), 0);
        else std::fill(wide_sums.begin(), wide_sums.end(), 0);
        const auto start = Clock::now();
        for (std::size_t repetition = 0; repetition < 4; ++repetition) {
            if (narrow) {
#pragma omp parallel for collapse(2) schedule(static)
                for (std::size_t prime = 0; prime < 4; ++prime)
                    for (std::size_t row = 0; row < rows; ++row) {
                        const auto *table = narrow_tables.data()+prime*tables.size();
                        modular_ntt32::Montgomery arithmetic(modular_ntt32::primes[prime].value);
                        auto *first = narrow_sums.data()+(prime*2*rows+row)*n*columns;
                        auto *second = first+rows*n*columns;
                        for (std::size_t first_point = 0; first_point < n; first_point += 256)
                            for (std::size_t y = 0; y < selectors; ++y) {
                                const auto exponent = modular_ntt::power(5,(n/2-y)%(n/2),2*n);
                                for (std::size_t i = first_point; i < first_point+256; ++i) {
                                    const auto index = narrow_detail::narrowNTTPermutationIndex<n>(exponent,i)*columns;
                                    arithmetic.addPair(first+i*columns,second+i*columns,table+index,table+n*columns+index,columns,
                                        narrow_source[prime*source.size()+(row*selectors+y)*n+i]);
                                }
                            }
                    }
            }
            else {
#pragma omp parallel for collapse(2) schedule(static)
                for (std::size_t prime = 0; prime < 2; ++prime)
                    for (std::size_t row = 0; row < rows; ++row) {
                        const auto modulus = modular_ntt::degree65536_primes[prime].value;
                        narrow_detail::SelectorMultiplierFactory factory(modulus);
                        narrow_detail::addSelectorMatrixRowsPair<n>(
                            wide_sums.data()+(prime*2*rows+row)*n*columns,
                            wide_sums.data()+(prime*2*rows+row+rows)*n*columns,
                            wide_tables.data()+prime*tables.size(),wide_tables.data()+prime*tables.size()+n*columns,columns,selectors,
                            [&](std::size_t y) { return narrow_detail::SelectorMatrixRow{wide_source.data()+prime*source.size()+(row*selectors+y)*n,
                                modular_ntt::power(5,(n/2-y)%(n/2),2*n)}; }, factory,modulus);
                    }
            }
        }
        std::cout << "{\"small_primes\":" << (narrow ? "true" : "false")
                  << ",\"matrix_seconds\":" << std::chrono::duration<double>(Clock::now()-start).count() << "}" << std::endl;
    }
    const modular_ntt::TwoPrimeCRT old_crt(modular_ntt::degree65536_primes[0],modular_ntt::degree65536_primes[1]);
    for (std::size_t i = 0; i < 2*rows*n*columns; ++i) {
        const auto old_value = old_crt.reconstructSigned(wide_sums[i],wide_sums[2*rows*n*columns+i]);
        const auto new_value = crt.reconstructSigned({narrow_sums[i]%std::uint32_t(modular_ntt32::primes[0].value),
            narrow_sums[2*rows*n*columns+i]%std::uint32_t(modular_ntt32::primes[1].value),
            narrow_sums[4*rows*n*columns+i]%std::uint32_t(modular_ntt32::primes[2].value),
            narrow_sums[6*rows*n*columns+i]%std::uint32_t(modular_ntt32::primes[3].value)});
        if (new_value != old_value) throw std::runtime_error("small-prime full matrix differs from wide CRT");
    }
    std::cout << "{\"matrix_crt_coefficients\":" << 2*rows*n*columns << ",\"passed\":true}" << std::endl;
}
