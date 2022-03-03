#pragma once

#include "odesolve.hpp"



// ¹1 - variant (2) system
namespace variant_equation {
	constexpr double a = 1.5;
	constexpr double b = -0.3;
	constexpr double delta = 0.01;
	constexpr double F = 0.05;
	inline double v = std::sqrt(a) + 0.15; // BECAUSE std::sqrt() IS STILL NOT CONSTEXPR

	inline Vector f(double t, const Vector &u) {
		return Vector{ { u(1), F * cos(v * t) - 2. * delta * u(1) - a * u(0) - b * cube(u(0)) } };
	}

	inline auto y0 = Vector{ { 0.1, 0.1 } };

	constexpr double timeInterval = 10;
	constexpr double tau1 = 1e-2;
	constexpr double tau2 = 2e-2;
	constexpr double tau3 = 4e-2;

	inline ODEProblem odeproblem1{ &f, y0, timeInterval, tau1 };
	inline ODEProblem odeproblem2{ &f, y0, timeInterval, tau2 };
	inline ODEProblem odeproblem3{ &f, y0, timeInterval, tau3 };
}

namespace variant13_equation {
	constexpr double r1 = 0.4;
	constexpr double r2 = 0.2;
	constexpr double b11 = 0.005;
	constexpr double b12 = 0.08;
	constexpr double b21 = 0.04;
	constexpr double b22 = 0.003;

	inline Vector f(double t, const Vector &u) {
		return Vector{ {
				r1 * u(0) - b11 * sqr(u(0)) - b12 * u(0) * u(1),
				-r2 * u(1) - b22 * sqr(u(1)) + b21 * u(0) * u(1)
			} };
	}

	inline auto y0 = Vector{ { 0.1, 1.0 } };

	constexpr double timeInterval = 150;
	constexpr double tau1 = 1e-2;
	constexpr double tau2 = 2e-2;
	constexpr double tau3 = 4e-2;

	inline ODEProblem odeproblem1{ &f, y0, timeInterval, tau1 };
	inline ODEProblem odeproblem2{ &f, y0, timeInterval, tau2 };
	inline ODEProblem odeproblem3{ &f, y0, timeInterval, tau3 };

	constexpr double epsilon = 1e-9;
}

namespace oscillating_equation {
	inline Vector f(double t, const Vector &u) {
		// sin * sin
		/*return Vector{ {
				1.,
				40. * cos(40. * t) * std::pow(sin(0.4 * t), 20) + 1.2 * cos(0.4 * t) * std::pow(sin(0.4 * t), 19) * sin(40. * t)
			} };*/
		// (5 - t)^10
		/*return Vector{ {
				1.,
				10. * std::pow(t - 5, 9)
			} };*/
		// sin
		return Vector{ {
				1.,
				0.1 * cos(0.1 * t)
		} };
	}

	inline auto y0 = Vector{ { 0., 0. } };

	constexpr double timeInterval = 20 * PI;
	constexpr double tau = 1e-2;

	inline ODEProblem odeproblem{ &f, y0, timeInterval, tau };

	constexpr double epsilon = 1e-9;
}

// ¹2 - Spring equation
namespace spring_equation {
	constexpr double w = 1.; // w = sqrt(k/m)

	inline Vector f(double t, const Vector &u) {
		return Vector{ { u(1), -sqr(w) * u(0) } };
	}

	inline auto y0 = Vector{ { 0., 1. } };

	constexpr double timeInterval = 100;
	constexpr double tau = 1e-2;
	
	inline ODEProblem odeproblem{ &f, y0, timeInterval, tau };
}

// ¹3 - Tests
namespace test1_equation {
	inline Vector f(double t, const Vector &u) {
		return Vector{ { 2. * u(0) + sqr(u(1)) - 1., 6. * u(0) - sqr(u(1)) + 1. } };
	}

	inline auto y0 = Vector{ { 0., 0. } };

	constexpr double timeInterval = 0.3;
	constexpr double tau = 1e-2;

	inline ODEProblem odeproblem{ &f, y0, timeInterval, tau };

	// Phase portrait
	constexpr double xmin = -2;
	constexpr double xmax = 2;
	constexpr double ymin = -2;
	constexpr double ymax = 2;

	constexpr uint Nx = 21;
	constexpr uint Ny = 21;

	constexpr double hx = (xmax - xmin) / (Nx - 1);
	constexpr double hy = (ymax - ymin) / (Ny - 1);
}

namespace test2_equation {
	inline Vector f(double t, const Vector &u) {
		return Vector{ { 1. - sqr(u(0)) - sqr(u(1)), 2. * u(0) } };
	}

	inline auto y0 = Vector{ { 0., 0. } };

	constexpr double timeInterval = 0.3;
	constexpr double tau = 1e-2;

	inline ODEProblem odeproblem{ &f, y0, timeInterval, tau };

	// Phase portrait
	constexpr double xmin = -2;
	constexpr double xmax = 2;
	constexpr double ymin = -2;
	constexpr double ymax = 2;

	constexpr uint Nx = 21;
	constexpr uint Ny = 21;

	constexpr double hx = (xmax - xmin) / (Nx - 1);
	constexpr double hy = (ymax - ymin) / (Ny - 1);
}

namespace test3_equation {
	constexpr double sigma = 10.;
	constexpr double r = 28.;
	constexpr double b = 8. / 3.;

	inline Vector f(double t, const Vector &u) {
		return Vector{ { sigma * (u(1) - u(0)), u(0) * (r - u(2)) - u(1), u(0) * u(1) - b * u(2) } };
	}

	inline auto y0 = Vector{ { 0., 0., 0. } };

	constexpr double timeInterval = 0.1;
	constexpr double tau = 1e-2;

	inline ODEProblem odeproblem{ &f, y0, timeInterval, tau };

	// Phase portrait
	constexpr double xmin = -2;
	constexpr double xmax = 2;
	constexpr double ymin = -2;
	constexpr double ymax = 2;
	constexpr double zmin = -2;
	constexpr double zmax = 2;

	constexpr uint Nx = 7;
	constexpr uint Ny = 7;
	constexpr uint Nz = 7;

	constexpr double hx = (xmax - xmin) / (Nx - 1);
	constexpr double hy = (ymax - ymin) / (Ny - 1);
	constexpr double hz = (zmax - zmin) / (Ny - 1);
}