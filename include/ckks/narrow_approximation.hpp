#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>
#include <cereal/types/complex.hpp>
#include <cereal/types/optional.hpp>
#include <cereal/types/vector.hpp>

namespace TFHEpp {

// Section 7's offline approximation. Positive coefficients use
// p(x)=c0+2*Re(sum_{k=1}^d c_k exp(2*pi*i*k*x/gap)).
struct CKKSNarrowTrigonometricApproximation {
    std::vector<std::complex<double>> coefficients;
    double message_bound = 0;
    double gap = 0;
    double max_sampled_error = 0;
    std::optional<double> uniform_error_bound;
    std::size_t iterations = 0;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(coefficients, message_bound, gap, max_sampled_error,
                uniform_error_bound, iterations);
    }
};

inline double CKKSEvaluateNarrowApproximation(
    const CKKSNarrowTrigonometricApproximation &approximation, double message)
{
    if (approximation.coefficients.empty() || !std::isfinite(message) ||
        !std::isfinite(approximation.gap) || approximation.gap <= 0)
        throw std::invalid_argument("invalid narrow trigonometric approximation");
    const long double angle = 2 * std::acos(-1.0L) * message / approximation.gap;
    long double result = approximation.coefficients[0].real();
    for (std::size_t k = 1; k < approximation.coefficients.size(); k++) {
        const auto c = approximation.coefficients[k];
        result += 2 * (static_cast<long double>(c.real()) * std::cos(k * angle) -
                       static_cast<long double>(c.imag()) * std::sin(k * angle));
    }
    return static_cast<double>(result);
}

namespace narrow_approximation_detail {

inline std::vector<long double> basis(long double x, std::size_t degree, double gap)
{
    std::vector<long double> result(2 * degree + 1);
    result[0] = 1;
    const long double angle = 2 * std::acos(-1.0L) * x / gap;
    for (std::size_t k = 1; k <= degree; k++) {
        result[2 * k - 1] = std::cos(k * angle);
        result[2 * k] = std::sin(k * angle);
    }
    return result;
}

inline std::vector<long double> solve(std::vector<std::vector<long double>> matrix)
{
    const auto n = matrix.size();
    for (std::size_t col = 0; col < n; col++) {
        auto pivot = col;
        for (std::size_t row = col + 1; row < n; row++)
            if (std::abs(matrix[row][col]) > std::abs(matrix[pivot][col])) pivot = row;
        if (std::abs(matrix[pivot][col]) < 1e-28L)
            throw std::invalid_argument("ill-conditioned trigonometric Remez domain or degree");
        std::swap(matrix[col], matrix[pivot]);
        for (std::size_t row = col + 1; row < n; row++) {
            const long double multiplier = matrix[row][col] / matrix[col][col];
            for (std::size_t k = col + 1; k <= n; k++)
                matrix[row][k] -= multiplier * matrix[col][k];
        }
    }
    std::vector<long double> result(n);
    for (std::size_t row = n; row-- > 0;) {
        long double value = matrix[row][n];
        for (std::size_t col = row + 1; col < n; col++) value -= matrix[row][col] * result[col];
        result[row] = value / matrix[row][row];
    }
    return result;
}

struct Extremum {
    long double x;
    long double error;
};

}  // namespace narrow_approximation_detail

// Discrete trigonometric Remez exchange on a strict subinterval of a period.
// Validate the resulting double coefficients on an independent shifted grid.
// If a true Lipschitz bound for the target is supplied, the reported uniform
// bound additionally covers between-grid points using |p'| <= sum 4*pi*k|c_k|/gap.
// An optional target curvature bound provides a quadratic-spacing bound using
// interpolation of the residual. Bounds are conditional on the caller's target
// derivative bounds and include a floating-point allowance. A second-derivative
// callback plus a target fourth-derivative bound additionally bounds residual
// curvature on the declared interval, preserving cancellation of large terms.
inline CKKSNarrowTrigonometricApproximation CKKSBuildNarrowTrigonometricApproximation(
    const std::function<double(double)> &target, double message_bound, double gap,
    std::size_t degree, std::optional<double> target_lipschitz_bound = std::nullopt,
    std::size_t grid_intervals = 8192, std::size_t max_iterations = 40,
    std::optional<double> target_second_derivative_bound = std::nullopt,
    const std::function<double(double)> &target_second_derivative = {},
    std::optional<double> target_fourth_derivative_bound = std::nullopt)
{
    using namespace narrow_approximation_detail;
    if (!target || !std::isfinite(message_bound) || message_bound <= 0 ||
        !std::isfinite(gap) || gap <= 2 * message_bound || degree > 128 ||
        grid_intervals < 16 * (2 * degree + 2) || max_iterations == 0 ||
        (target_lipschitz_bound && (!std::isfinite(*target_lipschitz_bound) ||
                                   *target_lipschitz_bound < 0)))
        throw std::invalid_argument("invalid narrow Remez approximation request");
    if (target_second_derivative_bound && (!std::isfinite(*target_second_derivative_bound) ||
                                          *target_second_derivative_bound < 0))
        throw std::invalid_argument("invalid target curvature bound");
    if (static_cast<bool>(target_second_derivative) != target_fourth_derivative_bound.has_value() ||
        (target_fourth_derivative_bound && (!std::isfinite(*target_fourth_derivative_bound) ||
                                           *target_fourth_derivative_bound < 0)))
        throw std::invalid_argument("residual curvature requires a derivative callback and fourth-derivative bound");
    const std::size_t dimension = 2 * degree + 1;
    const std::size_t count = dimension + 1;
    auto evaluate_target = [&](long double x) {
        const double value = target(static_cast<double>(x));
        if (!std::isfinite(value)) throw std::invalid_argument("nonfinite Remez target value");
        return static_cast<long double>(value);
    };
    std::vector<long double> reference(count);
    for (std::size_t i = 0; i < count; i++)
        reference[i] = -message_bound * std::cos(std::acos(-1.0L) * i / (count - 1));
    std::vector<long double> solution;
    CKKSNarrowTrigonometricApproximation result;
    result.gap = gap;
    result.message_bound = message_bound;
    bool converged = false;
    for (std::size_t iteration = 0; iteration < max_iterations; iteration++) {
        std::vector<std::vector<long double>> matrix(count);
        for (std::size_t i = 0; i < count; i++) {
            matrix[i] = basis(reference[i], degree, gap);
            matrix[i].push_back(i % 2 == 0 ? 1 : -1);
            matrix[i].push_back(evaluate_target(reference[i]));
        }
        solution = solve(std::move(matrix));
        auto error_at = [&](long double x) {
            const auto values = basis(x, degree, gap);
            long double polynomial = 0;
            for (std::size_t j = 0; j < dimension; j++) polynomial += solution[j] * values[j];
            return evaluate_target(x) - polynomial;
        };
        // Include the old alternation points to retain alternation even when
        // extrema are narrower than the uniform search-grid spacing.
        std::vector<long double> samples = reference;
        for (std::size_t i = 0; i <= grid_intervals; i++)
            samples.push_back(-message_bound + 2.0L * message_bound * i / grid_intervals);
        std::sort(samples.begin(), samples.end());
        samples.erase(std::unique(samples.begin(), samples.end()), samples.end());
        std::vector<long double> errors(samples.size());
        long double maximum = 0;
        for (std::size_t i = 0; i < samples.size(); i++) {
            errors[i] = error_at(samples[i]);
            maximum = std::max(maximum, std::abs(errors[i]));
        }
        result.iterations = iteration + 1;
        if (maximum - std::abs(solution.back()) <= 1e-10L * std::max(1.0L, maximum)) {
            converged = true;
            break;
        }
        std::vector<Extremum> extrema;
        for (std::size_t i = 0; i < samples.size(); i++) {
            const bool endpoint = i == 0 || i + 1 == samples.size();
            const bool local = !endpoint &&
                ((errors[i] >= errors[i - 1] && errors[i] >= errors[i + 1]) ||
                 (errors[i] <= errors[i - 1] && errors[i] <= errors[i + 1]));
            if (!endpoint && !local) continue;
            if (!extrema.empty() && std::signbit(errors[i]) == std::signbit(extrema.back().error)) {
                if (std::abs(errors[i]) > std::abs(extrema.back().error)) extrema.back() = {samples[i], errors[i]};
            }
            else extrema.push_back({samples[i], errors[i]});
        }
        if (extrema.size() < count)
            throw std::runtime_error("trigonometric Remez lost alternation; increase grid resolution");
        while (extrema.size() > count) {
            if (extrema.size() == count + 1) {
                if (std::abs(extrema.front().error) < std::abs(extrema.back().error)) extrema.erase(extrema.begin());
                else extrema.pop_back();
                continue;
            }
            const auto weakest = std::min_element(extrema.begin(), extrema.end(),
                [](const auto &left, const auto &right) { return std::abs(left.error) < std::abs(right.error); });
            const auto index = static_cast<std::size_t>(weakest - extrema.begin());
            if (index == 0 || index + 1 == extrema.size()) extrema.erase(weakest);
            else {
                const auto keep = std::abs(extrema[index - 1].error) > std::abs(extrema[index + 1].error)
                    ? extrema[index - 1] : extrema[index + 1];
                extrema[index - 1] = keep;
                extrema.erase(extrema.begin() + index, extrema.begin() + index + 2);
            }
        }
        for (std::size_t i = 0; i < count; i++) reference[i] = extrema[i].x;
    }
    if (!converged) throw std::runtime_error("trigonometric Remez did not converge");
    result.coefficients.resize(degree + 1);
    result.coefficients[0] = {static_cast<double>(solution[0]), 0};
    for (std::size_t k = 1; k <= degree; k++)
        result.coefficients[k] = {static_cast<double>(solution[2 * k - 1] / 2),
                                  static_cast<double>(-solution[2 * k] / 2)};
    double sampled_residual_curvature = 0;
    auto sample_curvature = [&](double x) {
        if (!target_second_derivative) return;
        const double target_curvature = target_second_derivative(x);
        if (!std::isfinite(target_curvature))
            throw std::invalid_argument("nonfinite target second derivative");
        const long double omega = 2 * std::acos(-1.0L) / gap;
        long double polynomial_curvature = 0;
        for (std::size_t k = 1; k <= degree; k++) {
            const auto c = result.coefficients[k];
            const long double angle = omega * k * x;
            polynomial_curvature -= 2 * omega * omega * k * k *
                (static_cast<long double>(c.real()) * std::cos(angle) -
                 static_cast<long double>(c.imag()) * std::sin(angle));
        }
        sampled_residual_curvature = std::max(sampled_residual_curvature,
            static_cast<double>(std::abs(target_curvature - polynomial_curvature)));
    };
    // This grid is twice as fine and offset from the exchange grid.
    const std::size_t validation_intervals = 2 * grid_intervals;
    for (std::size_t i = 0; i <= validation_intervals; i++) {
        const double x = -message_bound + 2 * message_bound *
            std::min(static_cast<double>(i) + 0.5, static_cast<double>(validation_intervals)) /
            validation_intervals;
        sample_curvature(x);
        result.max_sampled_error = std::max(result.max_sampled_error,
            std::abs(static_cast<double>(evaluate_target(x)) - CKKSEvaluateNarrowApproximation(result, x)));
    }
    sample_curvature(-message_bound);
    result.max_sampled_error = std::max(result.max_sampled_error,
        std::abs(static_cast<double>(evaluate_target(-message_bound)) -
                 CKKSEvaluateNarrowApproximation(result, -message_bound)));
    if (target_lipschitz_bound) {
        double derivative_bound = *target_lipschitz_bound;
        for (std::size_t k = 1; k <= degree; k++)
            derivative_bound += 4 * std::acos(-1.0) * k * std::abs(result.coefficients[k]) / gap;
        // The shifted grid plus both endpoints has covering radius at most
        // one grid spacing. Include a conservative floating-point allowance.
        result.uniform_error_bound = result.max_sampled_error + derivative_bound *
            2 * message_bound / validation_intervals +
            64 * std::numeric_limits<double>::epsilon() * (1 + derivative_bound * gap);
    }
    if (target_second_derivative_bound) {
        double curvature_bound = *target_second_derivative_bound;
        double coefficient_norm = std::abs(result.coefficients[0]);
        const double omega = 2 * std::acos(-1.0) / gap;
        for (std::size_t k = 1; k <= degree; k++) {
            curvature_bound += 2 * omega * omega * k * k * std::abs(result.coefficients[k]);
            coefficient_norm += 2 * std::abs(result.coefficients[k]);
        }
        // Linear interpolation of the residual between adjacent validation
        // points has error <= sup|residual''| * spacing^2 / 8. The shifted
        // grid includes both endpoints and no gap exceeds one spacing.
        const double spacing = 2 * message_bound / validation_intervals;
        const double bound = result.max_sampled_error + curvature_bound * spacing * spacing / 8 +
            64 * std::numeric_limits<double>::epsilon() * (1 + coefficient_norm);
        if (!result.uniform_error_bound || bound < *result.uniform_error_bound)
            result.uniform_error_bound = bound;
    }
    if (target_second_derivative) {
        const double omega = 2 * std::acos(-1.0) / gap;
        double fourth_bound = *target_fourth_derivative_bound;
        double curvature_norm = 0, coefficient_norm = std::abs(result.coefficients[0]);
        for (std::size_t k = 1; k <= degree; k++) {
            const double norm = 2 * std::abs(result.coefficients[k]);
            const double frequency_squared = omega * omega * k * k;
            coefficient_norm += norm;
            curvature_norm += frequency_squared * norm;
            fourth_bound += frequency_squared * frequency_squared * norm;
        }
        const double spacing = 2 * message_bound / validation_intervals;
        const double interpolation_factor = spacing * spacing / 8;
        const double roundoff = 64 * std::numeric_limits<double>::epsilon();
        // Apply interpolation first to residual curvature, using the fourth
        // derivative bound, then to the residual itself. This retains
        // cancellation inside the interval even if its extension is large.
        const double residual_curvature_bound = sampled_residual_curvature +
            fourth_bound * interpolation_factor + roundoff * (1 + curvature_norm);
        const double bound = result.max_sampled_error + residual_curvature_bound * interpolation_factor +
            roundoff * (1 + coefficient_norm);
        if (!result.uniform_error_bound || bound < *result.uniform_error_bound)
            result.uniform_error_bound = bound;
    }
    return result;
}

// C8 taper used only outside the protected message interval. Bernstein form
// avoids cancellation in the degree-17 smoothstep near its endpoints.
namespace narrow_approximation_detail {
inline std::array<double, 3> window(double x, double inner, double outer)
{
    const double width = outer - inner;
    const double u = (std::abs(x) - inner) / width;
    if (u <= 0) return {1, 0, 0};
    if (u >= 1) return {0, 0, 0};
    double weight = 0, binomial = 1;
    for (int j = 0; j <= 8; j++) {
        if (j) binomial *= static_cast<double>(18 - j) / j;
        weight += binomial * std::pow(u, j) * std::pow(1 - u, 17 - j);
    }
    const double g = u * (1 - u);
    const double first = -218790 * std::pow(g, 8) / width;
    const double second = -218790 * 8 * std::pow(g, 7) * (1 - 2 * u) / (width * width);
    return {weight, std::copysign(1.0, x) * first, second};
}
}  // namespace narrow_approximation_detail

// Fit a smooth periodic extension to control coefficient growth while keeping
// the target exactly unchanged throughout [-message_bound,message_bound].
// The caller supplies first/second derivatives and bounds for derivative
// orders 0..4 throughout [-gap/2,gap/2]. The C8 taper is zero at both period
// endpoints. Remez fits a wider interval; its uniform bound also covers the
// protected interval. No encrypted approximation or additional level is used.
inline CKKSNarrowTrigonometricApproximation CKKSBuildNarrowWindowedTrigonometricApproximation(
    const std::function<double(double)> &target,
    const std::function<double(double)> &target_first_derivative,
    const std::function<double(double)> &target_second_derivative,
    const std::array<double, 5> &target_derivative_bounds,
    double message_bound, double gap, std::size_t degree,
    std::size_t grid_intervals = 8192, std::size_t max_iterations = 40)
{
    if (!target || !target_first_derivative || !target_second_derivative ||
        !std::isfinite(message_bound) || message_bound <= 0 ||
        !std::isfinite(gap) || gap <= 2 * message_bound)
        throw std::invalid_argument("invalid windowed trigonometric approximation request");
    for (const double bound : target_derivative_bounds)
        if (!std::isfinite(bound) || bound < 0)
            throw std::invalid_argument("invalid windowed target derivative bound");
    const double outer = gap / 2, width = outer - message_bound;
    const double fitting_bound = outer - width / 8;
    if (!(message_bound < fitting_bound && fitting_bound < outer))
        throw std::invalid_argument("window transition cannot be represented");
    const double w1 = 218790 * std::pow(0.25, 8) / width;
    const double w2 = 218790 * 8 * std::pow(0.25, 7) / std::pow(width, 2);
    const double w3 = 218790 * 8 * (7 * std::pow(0.25, 6) + 2 * std::pow(0.25, 7)) /
        std::pow(width, 3);
    const double w4 = 218790 * 8 * (7 * 6 * std::pow(0.25, 5) + 6 * 7 * std::pow(0.25, 6)) /
        std::pow(width, 4);
    const auto &b = target_derivative_bounds;
    const double first_bound = b[1] + b[0] * w1;
    const double second_bound = b[2] + 2 * b[1] * w1 + b[0] * w2;
    const double fourth_bound = b[4] + 4 * b[3] * w1 + 6 * b[2] * w2 + 4 * b[1] * w3 + b[0] * w4;
    const auto extended = [&](double x) {
        return target(x) * narrow_approximation_detail::window(x, message_bound, outer)[0];
    };
    const auto curvature = [&](double x) {
        const auto w = narrow_approximation_detail::window(x, message_bound, outer);
        return target_second_derivative(x) * w[0] + 2 * target_first_derivative(x) * w[1] + target(x) * w[2];
    };
    auto result = CKKSBuildNarrowTrigonometricApproximation(extended, fitting_bound, gap, degree,
        first_bound, grid_intervals, max_iterations, second_bound, curvature, fourth_bound);
    // The target is unchanged on this interval; the wider fit's bound remains
    // valid here. Keep the public metadata about the actual message domain.
    result.message_bound = message_bound;
    return result;
}

}  // namespace TFHEpp
