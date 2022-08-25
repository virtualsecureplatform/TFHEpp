#include <array>
#include <tfhe++.hpp>

using namespace TFHEpp;

template <class P>
void HomExpandStep2(TRGSW<P>& out, const TRGSWFFT<P>& A,
                    const std::array<TRLWE<P>, P::l>& c)
{
    static_assert(P::k == 1);
    for (size_t k = 0; k < P::l; k++) {
        TFHEpp::trgswfftExternalProduct<P>(out[k], c[k], A);
        out[k + P::l] = c[k];
    }
}

int main()
{
    const size_t num_tests = 100;
    using P = CCR2019;

    for (size_t i = 0; i < num_tests; i++) {
        Key<P> key;
        bool ans = true;

        const TRGSWFFT<P> A = [&] {
            Polynomial<P> mu = {};
            for (size_t i = 0; i < P::n; i++) mu[i] = -static_cast<int>(key[i]);
            return trgswfftSymEncrypt<P>(mu, P::α, key);
        }();

        std::array<TRLWE<P>, P::l> c;
        for (size_t i = 0; i < P::l; i++) {
            Polynomial<P> mu = {};
            c[i] = trlweSymEncrypt<P>(mu, P::α, key);
            c[i][1][0] += (ans ? 1 : 0) * (1u << (32 - (i + 1) * P::Bgbit));
        }

        TRGSW<P> out;
        HomExpandStep2<P>(out, A, c);
        TRGSWFFT<P> out_fft = ApplyFFT2trgsw<P>(out);

        TRLWE<P> c0 = trlweSymEncrypt<P>(Polynomial<P>{-P::μ}, P::α, key);
        TRLWE<P> c1 = trlweSymEncrypt<P>(Polynomial<P>{P::μ}, P::α, key);

        TRLWE<P> res;
        CMUXFFT<P>(res, out_fft, c1, c0);
        bool res_bool = trlweSymDecrypt<P>(res, key)[0];
        if (res_bool == ans) {
            std::cerr << ".";
        }
        else {
            std::cerr << "ERROR\t" << ans << "\t" << res_bool << "\n";
            return 1;
        }
    }

    return 0;
}
