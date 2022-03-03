#pragma once

#include <functional>
#include <algorithm>

#include "Eigen/Dense"



using uint = unsigned int;
using sint = int;

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

using ScalarFunction = std::function<double(double)>;
using VectorFunction = std::function<Vector(const Vector&)>;

using ODERightSideFunction = std::function<Vector(double, const Vector&)>;

constexpr double INF = std::numeric_limits<double>::infinity();
constexpr double PI = 3.14159265358979323846;

// Formatting
constexpr auto SEPARATOR = "\n\n----------------------------------------";
constexpr auto STARTER = "\n>>> ";

Eigen::IOFormat DEFAULT(6, 0, ", ", "\n", "[", "]");
Eigen::IOFormat INLINE(6, 0, ", ", ", ", "", "", "{", "}");
Eigen::IOFormat NONE(6, 0, " ", " ", "", "", "", "");

template<typename T>
constexpr bool is_zero(T value) {
	constexpr T EPSILON = static_cast<T>(1e-16);

	return std::abs(value) < EPSILON;
}

// sign(x) = { -1, 0, 1 } signum
template<typename T>
constexpr sint sign(T value) {
	return (T(0) < value) - (value < T(0));
}

// sqr(x) = x^2
template<typename T>
constexpr T sqr(T value) {
	return value * value;
}

// cube(x) = x^3
template<typename T>
constexpr T cube(T value) {
	return value * value * value;
}

constexpr double middle(double a, double b) {
	return 0.5 * (a + b);
}

// Differentiation
constexpr double DIFF_EPS = 2e-8;
constexpr double DIFF_INVERSE_EPS = 1. / DIFF_EPS;

inline double derivative(ScalarFunction f, double x) {
	return (f(x + DIFF_EPS) - f(x)) * DIFF_INVERSE_EPS;
}

inline Matrix jacobian(VectorFunction F, Vector X) {
	const auto N = X.size();

	Matrix J = Matrix::Zero(N, N);

	Vector dXj = Vector::Zero(N);

	for (int j = 0; j < N; ++j) {
		// Differentiate F with respect to x_j
		dXj(j) += DIFF_EPS;
		const auto dFj = (F(X + 0.5 * dXj) - F(X - 0.5 * dXj)) * DIFF_INVERSE_EPS;
		dXj(j) = 0.;

		// As a result we get an entire column of jacobian
		J.col(j) = dFj;
	}

	return J;
}