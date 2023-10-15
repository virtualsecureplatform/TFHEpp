#include <chrono>
#include <tfhe++.hpp>

using namespace TFHEpp;

const size_t NUM_TEST = 1000;

void perform(const SecretKey &sk, std::vector<TRGSWFFT<lvl1param>> &cs)
{
    Polynomial<lvl1param> poly = {};
    for (auto &&t : cs) {
        t = trgswfftSymEncrypt<lvl1param>(poly, sk.key.lvl1);
    }
}

int main()
{
    std::unique_ptr<SecretKey> sk = std::make_unique<SecretKey>();
    std::vector<TRGSWFFT<lvl1param>> cs(NUM_TEST);

    using namespace std::chrono;
    system_clock::time_point start = chrono::system_clock::now();
    perform(*sk, cs);
    system_clock::time_point end = chrono::system_clock::now();

    double elapsed = duration_cast<microseconds>(end - start).count();
    std::cout << "Total: " << elapsed << "(us)\n"
              << "Mean: " << elapsed / NUM_TEST << "(us)\n";
}