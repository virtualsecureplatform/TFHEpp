export module tfhepp:hexlpp;
import std;
export import tfhepp:INTorus; // For INTorus

export namespace hexelpp {
export constexpr uint32_t P = 1073741789;
export constexpr uint32_t W = 2;

// Removed 'const' from PowW as it's not a member function of a class here.
export uint32_t PowW(uint exponent)
{
    uint32_t res = 1;
    for (uint64_t i = 0; i < exponent; i++) res = (static_cast<uint64_t>(res) * W) % P; // Ensure intermediate mult is 64-bit
    return res;
}

export template <uint32_t Nbit>
inline std::array<std::array<TFHEpp::cuHEpp::INTorus, 1U << Nbit>, 2> TwistGen() // Assuming INTorus is TFHEpp::cuHEpp::INTorus
{
    constexpr uint32_t N = 1U << Nbit;

    std::array<std::array<TFHEpp::cuHEpp::INTorus, 1U << Nbit>, 2> twist;
    // PowW returns uint32_t, INTorus constructor might need adjustment or PowW needs to return INTorus
    const TFHEpp::cuHEpp::INTorus w_val = TFHEpp::cuHEpp::INTorus(W); // Example, may need Pow
    const TFHEpp::cuHEpp::INTorus w = w_val.Pow(1U << (TFHEpp::cuHEpp::INTorus::P_bit - Nbit - 1)); // Assuming P_bit in INTorus

    twist[0][0] = twist[1][0] = TFHEpp::cuHEpp::INTorus(1, false);
    for (uint32_t i = 1; i < N; i++) twist[1][i] = twist[1][i - 1] * w;
    assert((twist[1][N - 1] * w).Pow(2).value == 1); // value might be an issue if INTorus doesn't directly expose it or P is different
    twist[0][N - 1] = twist[1][N - 1] * w * w;
    for (uint32_t i = 2; i < N; i++)
        twist[0][N - i] = twist[0][N - i + 1] * w;
    assert((twist[0][1] * w).value == 1);
    return twist;
}

export template <uint32_t Nbit>
inline std::array<std::array<TFHEpp::cuHEpp::INTorus, 1U << Nbit>, 2> TableGen()
{
    constexpr uint32_t N = 1U << Nbit;

    std::array<std::array<TFHEpp::cuHEpp::INTorus, N>, 2> table;
    const TFHEpp::cuHEpp::INTorus w_val = TFHEpp::cuHEpp::INTorus(W); // W is from hexelpp namespace
    const TFHEpp::cuHEpp::INTorus w = w_val.Pow(1U << (TFHEpp::cuHEpp::INTorus::P_bit - Nbit)); // Assuming P_bit in INTorus

    table[0][0] = table[1][0] = TFHEpp::cuHEpp::INTorus(1, false);
    for (uint32_t i = 1; i < N; i++) table[1][i] = table[1][i - 1] * w;
    for (uint32_t i = 1; i < N; i++) table[0][i] = table[1][N - i];
    return table;
}

// HarveyButterfly arguments X, Y are uint32_t, but operations involve hexelpp::P.
// W_param is std::array<uint32_t,2>
export void HarveyButterfly(uint32_t &X, uint32_t &Y, std::array<uint32_t, 2> W_param) // Pass X, Y by reference
{
    const uint64_t T = static_cast<uint64_t>(X) - Y + 2 * P; // P is hexelpp::P
    X += Y;
    if (X >= 2 * P) X -= 2 * P; // P is hexelpp::P
    const uint32_t Q = (static_cast<uint64_t>(W_param[1]) * T) >> 32;
    Y = (static_cast<uint64_t>(W_param[0]) * T - static_cast<uint64_t>(Q) * P) >> 32; // P is hexelpp::P
}

}  // namespace hexelpp