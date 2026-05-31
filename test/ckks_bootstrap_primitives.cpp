#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <tfhe++.hpp>

namespace {

using P = TFHEpp::ckkslvl3param;

void fill_key(TFHEpp::Key<P> &key)
{
    for (std::size_t i = 0; i < P::n; i++) {
        const int v = static_cast<int>(i % 3) - 1;
        key[i] = static_cast<typename P::T>(v);
    }
}

double max_error(const TFHEpp::CKKSSlotVector<P> &got,
                 const TFHEpp::CKKSSlotVector<P> &want)
{
    double err = 0.0;
    for (std::size_t i = 0; i < P::n / 2; i++)
        err = std::max(err, std::abs(got[i] - want[i]));
    return err;
}

void require_close(const TFHEpp::CKKSSlotVector<P> &got,
                   const TFHEpp::CKKSSlotVector<P> &want, double tol,
                   const char *label)
{
    const double err = max_error(got, want);
    std::cout << label << " max_error=" << err << std::endl;
    if (err > tol) {
        std::cerr << label << " exceeded tolerance " << tol << std::endl;
        std::exit(1);
    }
}

}  // namespace

int main()
{
    constexpr std::uint32_t low_log_q = 82;
    constexpr std::uint32_t boot_log_q = 110;
    constexpr std::uint32_t log_delta = 40;
    constexpr std::uint32_t plain_log_delta = 20;
    constexpr double tol = 0.05;

    using LowCt = TFHEpp::CKKSCiphertext<P, low_log_q, log_delta>;
    using BootCt = TFHEpp::CKKSCiphertext<P, boot_log_q, log_delta>;
    using LinearOut =
        TFHEpp::CKKSPlainMulResult<P, boot_log_q, log_delta, plain_log_delta>;

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_key(*key);

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    slots->fill({0.0, 0.0});
    for (std::size_t i = 0; i < 32; i++) {
        const double re = static_cast<double>(static_cast<int>(i % 9) - 4) /
                          128.0;
        const double im = static_cast<double>(static_cast<int>(i % 5) - 2) /
                          128.0;
        (*slots)[i] = {re, im};
    }
    (*slots)[P::n / 2 - 1] = {-0.0625, 0.03125};

    auto low_ct = std::make_unique<LowCt>();
    TFHEpp::ckksSlotEncrypt<P, low_log_q, log_delta>(*low_ct, *slots, *key,
                                                     {0.0, 0});

    auto raised = std::make_unique<BootCt>();
    TFHEpp::CKKSModRaise<P, low_log_q, boot_log_q, log_delta>(*raised,
                                                              *low_ct);

    auto low_phase = std::make_unique<TFHEpp::Polynomial<P>>(
        TFHEpp::trlwePhase<P>(low_ct->ct, *key));
    auto raised_phase = std::make_unique<TFHEpp::Polynomial<P>>(
        TFHEpp::trlwePhase<P>(raised->ct, *key));
    for (std::uint32_t i = 0; i < P::n; i++) {
        const auto low =
            TFHEpp::ckks_detail::reduceToLevel<P, low_log_q>((*low_phase)[i]);
        const auto raised_low = TFHEpp::ckks_detail::reduceToLevel<P, low_log_q>(
            (*raised_phase)[i]);
        if (low != raised_low) {
            std::cerr << "mod raise phase mismatch at coeff " << i << std::endl;
            return 1;
        }
    }
    std::cout << "CKKS mod raise preserves low-level phase" << std::endl;

    auto boot_ct = std::make_unique<BootCt>();
    TFHEpp::ckksSlotEncrypt<P, boot_log_q, log_delta>(*boot_ct, *slots, *key,
                                                      {0.0, 0});

    auto boot_gk = std::make_unique<TFHEpp::CKKSGaloisKey<P, boot_log_q>>();
    TFHEpp::CKKSGaloisKeyGen<P, boot_log_q>(*boot_gk, *key, {0.0, 0});
    auto out_gk =
        std::make_unique<TFHEpp::CKKSGaloisKey<P, LinearOut::log_q>>();
    TFHEpp::CKKSGaloisKeyGen<P, LinearOut::log_q>(*out_gk, *key, {0.0, 0});

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    {
        auto conjugated = std::make_unique<BootCt>();
        TFHEpp::CKKSConjugateSlots<P, boot_log_q>(conjugated->ct, boot_ct->ct,
                                                  *boot_gk);
        TFHEpp::ckksSlotDecrypt<P, boot_log_q, log_delta>(*decoded,
                                                          *conjugated, *key);
        for (std::size_t i = 0; i < P::n / 2; i++)
            (*expected)[i] = std::conj((*slots)[i]);
        require_close(*decoded, *expected, tol, "CKKS conjugation");
    }

    std::vector<TFHEpp::CKKSSlotVector<P>> diagonals(4);
    for (auto &d : diagonals) d.fill({0.0, 0.0});
    for (std::size_t i = 0; i < P::n / 2; i++) {
        diagonals[0][i] = {0.25, 0.0};
        diagonals[1][i] = {0.125, 0.0625};
        diagonals[2][i] = {-0.0625, 0.03125};
        diagonals[3][i] = {i % 2 == 0 ? 0.03125 : -0.03125, -0.015625};
    }
    const std::vector<int> offsets{0, 1, 5, -3};
    constexpr int k_step = 4;
    TFHEpp::CKKSLinearTransformPlan<P, boot_log_q, log_delta, plain_log_delta>
        plan;
    TFHEpp::CKKSBuildLinearTransformBSGSPlan<P, boot_log_q, log_delta,
                                             plain_log_delta>(
        plan, diagonals, offsets, k_step);

    auto transformed = std::make_unique<LinearOut>();
    TFHEpp::CKKSLinearTransformBSGS<P, boot_log_q, log_delta,
                                    plain_log_delta>(
        *transformed, *boot_ct, plan, *boot_gk, *out_gk);
    TFHEpp::ckksSlotDecrypt<P, LinearOut::log_q, LinearOut::log_delta>(
        *decoded, *transformed, *key);

    constexpr int half = static_cast<int>(P::n) / 2;
    for (int i = 0; i < half; i++) {
        (*expected)[i] =
            diagonals[0][i] * (*slots)[i] +
            diagonals[1][i] * (*slots)[(i + 1) % half] +
            diagonals[2][i] * (*slots)[(i + 5) % half] +
            diagonals[3][i] * (*slots)[(i - 3 + half) % half];
    }
    require_close(*decoded, *expected, tol, "CKKS precomputed BSGS");

    std::cout << "Passed" << std::endl;
}
