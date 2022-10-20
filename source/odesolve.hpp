#pragma once

#include <tuple>
#include <vector>
#include <fstream>
#include <string>

#include "core_math.hpp"
#include "nonlinear_solver.hpp"



// ODE Problem
// - Represents u' = f(t, u) ODE problem
struct ODEProblem {
	ODERightSideFunction f;
	Vector y0;
	double timeInterval;
	double tau; // time step
};

// ODE Solution
using dvector = std::vector<double>;

inline void save_as(const dvector &vector, const std::string &filename) {
	std::ofstream file(filename);
	for (const auto &elem : vector) file << elem << '\n';
}

using vvector = std::vector<Vector>;

inline void save_as(const vvector &vector, const std::string &filename) {
	std::ofstream file(filename);
	for (const auto &elem : vector) file << elem.format(NONE) << '\n';
}

using ODESolution = std::tuple<dvector, vvector>;



// ODE Solvers
// Explicit Euler's method
ODESolution odesolve_explicit_euler(const ODEProblem &problem) {
	const auto f = problem.f;
	auto y0 = problem.y0;
	const auto timeInterval = problem.timeInterval;
	const auto tau = problem.tau;

	dvector out_t;
	vvector out_y;
	out_t.reserve(static_cast<size_t>(timeInterval / tau));
	out_y.reserve(static_cast<size_t>(timeInterval / tau));

	Vector y = y0;

	for (double t = 0; t < problem.timeInterval; t += tau) {
		out_t.push_back(t);
		out_y.push_back(y);

		y0 = y;

		y = y0 + tau * f(t, y0);
	}

	return { out_t, out_y };
}


// Implicit Euler's method (backwards Euler)
ODESolution odesolve_implicit_euler(const ODEProblem &problem) {
	const auto f = problem.f;
	auto y0 = problem.y0;
	const auto timeInterval = problem.timeInterval;
	const auto tau = problem.tau;

	dvector out_t;
	vvector out_y;
	out_t.reserve(static_cast<size_t>(timeInterval / tau));
	out_y.reserve(static_cast<size_t>(timeInterval / tau));

	Vector y = y0;

	for (double t = 0; t < problem.timeInterval; t += tau) {
		out_t.push_back(t);
		out_y.push_back(y);

		y0 = y;

		VectorFunction implicitEquation = [&](const Vector &yn) -> Vector {
			return yn - y0 - tau * f(t + tau, yn);
		};

		y = nonlinear_solve(implicitEquation, y0, 1e-6);
	}

	return { out_t, out_y };
}


// Symmetric Euler's method
ODESolution odesolve_symmetric_euler(const ODEProblem &problem) {
	const auto f = problem.f;
	auto y0 = problem.y0;
	const auto timeInterval = problem.timeInterval;
	const auto tau = problem.tau;

	dvector out_t;
	vvector out_y;
	out_t.reserve(static_cast<size_t>(timeInterval / tau));
	out_y.reserve(static_cast<size_t>(timeInterval / tau));

	Vector y = y0;

	for (double t = 0; t < problem.timeInterval; t += tau) {
		out_t.push_back(t);
		out_y.push_back(y);

		y0 = y;

		VectorFunction implicitEquation = [&](const Vector &yn) -> Vector {
			return yn - y0 - tau * 0.5 * (f(t + tau, yn) + f(t, y0));
		};

		y = nonlinear_solve(implicitEquation, y0, 1e-8);
	}

	return { out_t, out_y };
}


// 4th order Runge-Kutta method with adaptive step
Vector _RK4_iteration(ODERightSideFunction f, const Vector &y0, double t, double tau) {
	const auto k1 = f(t, y0);
	const auto k2 = f(t + 0.5 * tau, y0 + 0.5 * tau * k1);
	const auto k3 = f(t + 0.5 * tau, y0 + 0.5 * tau * k2);
	const auto k4 = f(t + tau, y0 + tau * k3);

	return y0 + tau / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
}

double error_estimate(ODERightSideFunction f, const Vector &y0, double t, double tau) {
	const auto fullStepY = _RK4_iteration(f, y0, t, tau);
	const auto halfStepY = _RK4_iteration(f, _RK4_iteration(f, y0, t, 0.5 * tau), t + 0.5 * tau, 0.5 * tau);

	return (fullStepY - halfStepY).norm() / 15.;
}

std::tuple<dvector, vvector, dvector, dvector> odesolve_adaptiveRK4(const ODEProblem &problem, double epsilon) {
	const auto f = problem.f;
	auto y0 = problem.y0;
	const auto timeInterval = problem.timeInterval;
	auto tau = problem.tau; // tau not const since it's changed adaptively

	const auto h = tau;

	dvector out_t;
	vvector out_y;
	dvector out_err;
	dvector out_tau;
	out_t.reserve(static_cast<size_t>(timeInterval / tau));
	out_y.reserve(static_cast<size_t>(timeInterval / tau));

	Vector y = y0;

	for (double t = 0; t < timeInterval; t += tau) {
		out_t.push_back(t);
		out_y.push_back(y);
		
		y0 = y;

		// Adjust step using Runge's law
		while (error_estimate(f, y0, t, tau) < 1e-1 * epsilon) tau *= 2.;
		while (error_estimate(f, y0, t, tau) > epsilon) tau *= 0.5;
		
		out_tau.push_back(tau);
		out_err.push_back(error_estimate(f, y0, t, tau));

		y = _RK4_iteration(f, y0, t, tau);
	}

	return { out_t, out_y, out_err, out_tau };
}

// 4th order Runge-Kutta method
ODESolution odesolve_RK4(const ODEProblem &problem) {
	const auto f = problem.f;
	auto y0 = problem.y0;
	const auto timeInterval = problem.timeInterval;
	const auto tau = problem.tau;

	dvector out_t;
	vvector out_y;
	out_t.reserve(static_cast<size_t>(timeInterval / tau));
	out_y.reserve(static_cast<size_t>(timeInterval / tau));

	Vector y = y0;

	for (double t = 0; t < timeInterval; t += tau) {
		out_t.push_back(t);
		out_y.push_back(y);

		y0 = y;

		const auto k1 = f(t, y0);
		const auto k2 = f(t + 0.5 * tau, y0 + 0.5 * tau * k1);
		const auto k3 = f(t + 0.5 * tau, y0 + 0.5 * tau * k2);
		const auto k4 = f(t + tau, y0 + tau * k3);

		y = y0 + tau / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
	}

	return { out_t, out_y };
}


// 2nd order Runge-Kutta method
ODESolution odesolve_RK2(const ODEProblem &problem) {
	const auto f = problem.f;
	auto y0 = problem.y0;
	const auto timeInterval = problem.timeInterval;
	const auto tau = problem.tau;

	dvector out_t;
	vvector out_y;
	out_t.reserve(static_cast<size_t>(timeInterval / tau));
	out_y.reserve(static_cast<size_t>(timeInterval / tau));

	Vector y = y0;

	for (double t = 0; t < timeInterval; t += tau) {
		out_t.push_back(t);
		out_y.push_back(y);

		y0 = y;

		const auto k1 = f(t, y0);
		const auto k2 = f(t + tau, y0 + tau * k1);

		y = y0 + tau * (0.5 * k1 + 0.5 * k2);
	}

	return { out_t, out_y };
}


// 4th order Predictor–corrector method
ODESolution odesolve_PredictorCorrector4(const ODEProblem &problem) {
	const auto f = problem.f;
	auto y0 = problem.y0;
	const auto timeInterval = problem.timeInterval;
	auto tau = problem.tau; // tau not const since it's changed adaptively

	dvector out_t;
	vvector out_y;
	out_t.reserve(static_cast<size_t>(timeInterval / tau));
	out_y.reserve(static_cast<size_t>(timeInterval / tau));

	vvector fvals;
	fvals.reserve(static_cast<size_t>(timeInterval / tau));

	Vector y = y0;

	// y0, y1, y2, y3 are obtained through RK4
	double t = 0;

	for (uint i = 0; i < 3; ++i, t += tau) {
		out_t.push_back(t);
		out_y.push_back(y);
		fvals.push_back(f(t, y));

		y0 = y;

		const auto k1 = f(t, y0);
		const auto k2 = f(t + 0.5 * tau, y0 + 0.5 * tau * k1);
		const auto k3 = f(t + 0.5 * tau, y0 + 0.5 * tau * k2);
		const auto k4 = f(t + tau, y0 + tau * k3);

		y = y0 + tau / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
	}

	// Starting from 5th approximation we can use intended multistep method
	for (; t < timeInterval; t += tau) {
		out_t.push_back(t);
		out_y.push_back(y);
		fvals.push_back(f(t, y));

		y0 = y;

		// Prediction
		const auto py = y0 + tau / 24. * (55. * fvals.end()[-1] - 59. * fvals.end()[-2]
			+ 37. * fvals.end()[-3] - 9. * fvals.end()[-4]);
		
		const auto pf = f(t + tau, py);

		// Correction
		y = y0 + tau / 24. * (9. * pf + 19. * fvals.end()[-1] - 5. * fvals.end()[-2] + fvals.end()[-3]);
	}

	return { out_t, out_y };
}


// 4th order Adams method
ODESolution odesolve_Adams4(const ODEProblem &problem) {
	const auto f = problem.f;
	auto y0 = problem.y0;
	const auto timeInterval = problem.timeInterval;
	auto tau = problem.tau; // tau not const since it's changed adaptively

	dvector out_t;
	vvector out_y;
	out_t.reserve(static_cast<size_t>(timeInterval / tau));
	out_y.reserve(static_cast<size_t>(timeInterval / tau));

	vvector fvals;
	fvals.reserve(static_cast<size_t>(timeInterval / tau));

	Vector y = y0;

	// y0, y1, y2, y3 are obtained through RK4
	double t = 0;

	for (uint i = 0; i < 3; ++i, t += tau) {
		out_t.push_back(t);
		out_y.push_back(y);
		fvals.push_back(f(t, y));

		y0 = y;

		const auto k1 = f(t, y0);
		const auto k2 = f(t + 0.5 * tau, y0 + 0.5 * tau * k1);
		const auto k3 = f(t + 0.5 * tau, y0 + 0.5 * tau * k2);
		const auto k4 = f(t + tau, y0 + tau * k3);

		y = y0 + tau / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
	}

	// Starting from 5th approximation we can use intended multistep method
	for (; t < timeInterval; t += tau) {
		out_t.push_back(t);
		out_y.push_back(y);
		fvals.push_back(f(t, y));

		y0 = y;

		y = y0 + tau / 24. * (55. * fvals.end()[-1] - 59. * fvals.end()[-2]
			+ 37. * fvals.end()[-3] - 9. * fvals.end()[-4]);
	}

	return { out_t, out_y };
}