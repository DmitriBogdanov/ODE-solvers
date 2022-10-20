#include <iostream>

#include <limits>

#include "static_timer.hpp"

#include "odesolve.hpp"
#include "odeproblems.hpp"


const std::string PATH = "[output]/";

// ¹1
void No1_adaptiveRK4() {
	std::cout
		<< SEPARATOR
		<< STARTER << "Running adaptive 4th order Runge-Kutta method...\n";

	// Solve with 3 different time steps
	StaticTimer::start();
	const auto [t, y, err, tau] = odesolve_adaptiveRK4(variant_equation::odeproblem1, variant_equation::epsilon);
	StaticTimer::end();

	// Save
	save_as(t, PATH + "[adaptiveRK4][variant_equation](t).dat");
	save_as(y, PATH + "[adaptiveRK4][variant_equation](y).dat");
	save_as(err, PATH + "[adaptiveRK4][variant_equation](err).dat");
	save_as(tau, PATH + "[adaptiveRK4][variant_equation](tau).dat");

	std::ofstream file(PATH + "[adaptiveRK4][variant_equation](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Adaptive test equation"
		<< '\n' << "Method" << '\n' << "Adaptive 4th order Runge-Kutta method"
		<< '\n' << "Time interval" << '\n' << variant_equation::odeproblem1.timeInterval
		<< '\n' << "Initial time step" << '\n' << variant_equation::odeproblem1.tau
		<< '\n' << "y0" << '\n' << variant_equation::odeproblem1.y0.format(INLINE)
		<< '\n' << "Epsilon" << '\n' << variant_equation::epsilon;
}

void No1_RK4() {
	std::cout
		<< SEPARATOR
		<< STARTER << "Running 4th order Runge-Kutta method...\n";

	// Solve with 3 different time steps
	StaticTimer::start();
	const auto [t1, y1] = odesolve_RK4(variant_equation::odeproblem1);
	const auto [t2, y2] = odesolve_RK4(variant_equation::odeproblem2);
	const auto [t3, y3] = odesolve_RK4(variant_equation::odeproblem3);
	StaticTimer::end();

	// Find convergence order
	const double p = std::log((y3.back() - y2.back()).norm() / (y2.back() - y1.back()).norm()) / std::log(2);

	// Save
	save_as(t1, PATH + "[RK4][variant_equation](t).dat");
	save_as(y1, PATH + "[RK4][variant_equation](y).dat");

	std::ofstream file(PATH + "[RK4][variant_equation](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Variant 2 equation"
		<< '\n' << "Method" << '\n' << "4th order Runge-Kutta method"
		<< '\n' << "Time interval" << '\n' << variant_equation::odeproblem1.timeInterval
		<< '\n' << "Time step 1" << '\n' << variant_equation::odeproblem1.tau
		<< '\n' << "Time step 2" << '\n' << variant_equation::odeproblem2.tau
		<< '\n' << "Time step 3" << '\n' << variant_equation::odeproblem3.tau
		<< '\n' << "y0" << '\n' << variant_equation::odeproblem1.y0.format(INLINE)
		<< '\n' << "Convergence order" << '\n' << p;
}

void No1_RK2() {
	std::cout
		<< SEPARATOR
		<< STARTER << "Running 2nd order Runge-Kutta method...\n";

	// Solve with 3 different time steps
	using namespace variant_equation;

	StaticTimer::start();
	const auto [t1, y1] = odesolve_RK2({ &f, variant_equation::y0, timeInterval, tau1 });
	const auto [t2, y2] = odesolve_RK2(variant_equation::odeproblem2);
	const auto [t3, y3] = odesolve_RK2(variant_equation::odeproblem3);
	StaticTimer::end();

	// Find convergence order
	const double p = std::log((y3.back() - y2.back()).norm() / (y2.back() - y1.back()).norm()) / std::log(2);

	// Save
	save_as(t1, PATH + "[RK2][variant_equation](t).dat");
	save_as(y1, PATH + "[RK2][variant_equation](y).dat");

	std::ofstream file(PATH + "[RK2][variant_equation](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Variant 2 equation"
		<< '\n' << "Method" << '\n' << "2nd order Runge-Kutta method"
		<< '\n' << "Time interval" << '\n' << variant_equation::odeproblem1.timeInterval
		<< '\n' << "Time step 1" << '\n' << variant_equation::odeproblem1.tau
		<< '\n' << "Time step 2" << '\n' << variant_equation::odeproblem2.tau
		<< '\n' << "Time step 3" << '\n' << variant_equation::odeproblem3.tau
		<< '\n' << "y0" << '\n' << variant_equation::odeproblem1.y0.format(INLINE)
		<< '\n' << "Convergence order" << '\n' << p;
}

void No1_PredictorCorrector4() {
	std::cout
		<< SEPARATOR
		<< STARTER << "Running 4th order Predictor-corrector method...\n";

	// Solve with 3 different time steps
	StaticTimer::start();
	const auto[t1, y1] = odesolve_PredictorCorrector4(variant_equation::odeproblem1);
	const auto[t2, y2] = odesolve_PredictorCorrector4(variant_equation::odeproblem2);
	const auto[t3, y3] = odesolve_PredictorCorrector4(variant_equation::odeproblem3);
	StaticTimer::end();

	// Find convergence order
	const double p = std::log((y3.back() - y2.back()).norm() / (y2.back() - y1.back()).norm()) / std::log(2);

	// Save
	save_as(t1, PATH + "[PredictorCorrector4][variant_equation](t).dat");
	save_as(y1, PATH + "[PredictorCorrector4][variant_equation](y).dat");

	std::ofstream file(PATH + "[PredictorCorrector4][variant_equation](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Variant 2 equation"
		<< '\n' << "Method" << '\n' << "4th order Predictor-corrector method"
		<< '\n' << "Time interval" << '\n' << variant_equation::odeproblem1.timeInterval
		<< '\n' << "Time step 1" << '\n' << variant_equation::odeproblem1.tau
		<< '\n' << "Time step 2" << '\n' << variant_equation::odeproblem2.tau
		<< '\n' << "Time step 3" << '\n' << variant_equation::odeproblem3.tau
		<< '\n' << "y0" << '\n' << variant_equation::odeproblem1.y0.format(INLINE)
		<< '\n' << "Convergence order" << '\n' << p;
}

void No1_Adams4() {
	std::cout
		<< SEPARATOR
		<< STARTER << "Running 4th order Adams-Bashforth method...\n";

	// Solve with 3 different time steps
	StaticTimer::start();
	const auto[t1, y1] = odesolve_Adams4(variant_equation::odeproblem1);
	const auto[t2, y2] = odesolve_Adams4(variant_equation::odeproblem2);
	const auto[t3, y3] = odesolve_Adams4(variant_equation::odeproblem3);
	StaticTimer::end();

	// Find convergence order
	const double p = std::log((y3.back() - y2.back()).norm() / (y2.back() - y1.back()).norm()) / std::log(2);

	// Save
	save_as(t1, PATH + "[Adams4][variant_equation](t).dat");
	save_as(y1, PATH + "[Adams4][variant_equation](y).dat");

	std::ofstream file(PATH + "[Adams4][variant_equation](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Variant 2 equation"
		<< '\n' << "Method" << '\n' << "4th order Adams-Bashforth method"
		<< '\n' << "Time interval" << '\n' << variant_equation::odeproblem1.timeInterval
		<< '\n' << "Time step 1" << '\n' << variant_equation::odeproblem1.tau
		<< '\n' << "Time step 2" << '\n' << variant_equation::odeproblem2.tau
		<< '\n' << "Time step 3" << '\n' << variant_equation::odeproblem3.tau
		<< '\n' << "y0" << '\n' << variant_equation::odeproblem1.y0.format(INLINE)
		<< '\n' << "Convergence order" << '\n' << p;
}

void No1_ExplicitEuler() {
	std::cout
		<< SEPARATOR
		<< STARTER << "Running explicit Euler's method...\n";

	// Solve with 3 different time steps
	StaticTimer::start();
	const auto[t1, y1] = odesolve_explicit_euler(variant_equation::odeproblem1);
	const auto[t2, y2] = odesolve_explicit_euler(variant_equation::odeproblem2);
	const auto[t3, y3] = odesolve_explicit_euler(variant_equation::odeproblem3);
	StaticTimer::end();

	// Find convergence order
	const double p = std::log((y3.back() - y2.back()).norm() / (y2.back() - y1.back()).norm()) / std::log(2);

	// Save
	save_as(t1, PATH + "[explicit_euler][variant_equation](t).dat");
	save_as(y1, PATH + "[explicit_euler][variant_equation](y).dat");

	std::ofstream file(PATH + "[explicit_euler][variant_equation](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Variant 2 equation"
		<< '\n' << "Method" << '\n' << "Explicit Euler's method"
		<< '\n' << "Time interval" << '\n' << variant_equation::odeproblem1.timeInterval
		<< '\n' << "Time step 1" << '\n' << variant_equation::odeproblem1.tau
		<< '\n' << "Time step 2" << '\n' << variant_equation::odeproblem2.tau
		<< '\n' << "Time step 3" << '\n' << variant_equation::odeproblem3.tau
		<< '\n' << "y0" << '\n' << variant_equation::odeproblem1.y0.format(INLINE)
		<< '\n' << "Convergence order" << '\n' << p;
}

void No1_ImplicitEuler() {
	std::cout
		<< SEPARATOR
		<< STARTER << "Running implicit Euler's method...\n";

	// Solve with 3 different time steps
	StaticTimer::start();
	const auto[t1, y1] = odesolve_implicit_euler(variant_equation::odeproblem1);
	const auto[t2, y2] = odesolve_implicit_euler(variant_equation::odeproblem2);
	const auto[t3, y3] = odesolve_implicit_euler(variant_equation::odeproblem3);
	StaticTimer::end();

	// Find convergence order
	const double p = std::log((y3.back() - y2.back()).norm() / (y2.back() - y1.back()).norm()) / std::log(2);

	// Save
	save_as(t1, PATH + "[implicit_euler][variant_equation](t).dat");
	save_as(y1, PATH + "[implicit_euler][variant_equation](y).dat");

	std::ofstream file(PATH + "[implicit_euler][variant_equation](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Variant 2 equation"
		<< '\n' << "Method" << '\n' << "Implicit Euler's method"
		<< '\n' << "Time interval" << '\n' << variant_equation::odeproblem1.timeInterval
		<< '\n' << "Time step 1" << '\n' << variant_equation::odeproblem1.tau
		<< '\n' << "Time step 2" << '\n' << variant_equation::odeproblem2.tau
		<< '\n' << "Time step 3" << '\n' << variant_equation::odeproblem3.tau
		<< '\n' << "y0" << '\n' << variant_equation::odeproblem1.y0.format(INLINE)
		<< '\n' << "Convergence order" << '\n' << p;
}

void No1_SymmetricEuler() {
	std::cout
		<< SEPARATOR
		<< STARTER << "Running symmetric Euler's method...\n";

	// Solve with 3 different time steps
	StaticTimer::start();
	const auto [t1, y1] = odesolve_symmetric_euler(variant_equation::odeproblem1);
	const auto [t2, y2] = odesolve_symmetric_euler(variant_equation::odeproblem2);
	const auto [t3, y3] = odesolve_symmetric_euler(variant_equation::odeproblem3);
	StaticTimer::end();

	// Find convergence order
	const double p = std::log((y3.back() - y2.back()).norm() / (y2.back() - y1.back()).norm()) / std::log(2);

	// Save
	save_as(t1, PATH + "[symmetric_euler][variant_equation](t).dat");
	save_as(y1, PATH + "[symmetric_euler][variant_equation](y).dat");

	std::ofstream file(PATH + "[symmetric_euler][variant_equation](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Variant 2 equation"
		<< '\n' << "Method" << '\n' << "Symmetric Euler's method"
		<< '\n' << "Time interval" << '\n' << variant_equation::odeproblem1.timeInterval
		<< '\n' << "Time step 1" << '\n' << variant_equation::odeproblem1.tau
		<< '\n' << "Time step 2" << '\n' << variant_equation::odeproblem2.tau
		<< '\n' << "Time step 3" << '\n' << variant_equation::odeproblem3.tau
		<< '\n' << "y0" << '\n' << variant_equation::odeproblem1.y0.format(INLINE)
		<< '\n' << "Convergence order" << '\n' << p;
}

// ¹2
void No2_explicit_euler() {
	std::cout
		<< SEPARATOR
		<< STARTER << "Running explicit euler method...\n";

	// Solve
	StaticTimer::start();
	const auto [t, y] = odesolve_explicit_euler(spring_equation::odeproblem);
	StaticTimer::end();

	// Save
	save_as(t, PATH + "[explicit_euler][spring_equation](t).dat");
	save_as(y, PATH + "[explicit_euler][spring_equation](y).dat");
	
	std::ofstream file(PATH + "[explicit_euler][spring_equation](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Harmonic oscillator equation"
		<< '\n' << "Method" << '\n' << "Explicit Euler's method"
		<< '\n' << "Time interval" << '\n' << spring_equation::odeproblem.timeInterval
		<< '\n' << "Time step" << '\n' << spring_equation::odeproblem.tau
		<< '\n' << "y0" << '\n' << spring_equation::odeproblem.y0.format(INLINE);
}

void No2_implicit_euler() {
	std::cout
		<< SEPARATOR
		<< STARTER << "Running implicit euler method...\n";

	// Solve
	StaticTimer::start();
	const auto [t, y] = odesolve_implicit_euler(spring_equation::odeproblem);
	StaticTimer::end();

	// Save
	save_as(t, PATH + "[implicit_euler][spring_equation](t).dat");
	save_as(y, PATH + "[implicit_euler][spring_equation](y).dat");

	std::ofstream file(PATH + "[implicit_euler][spring_equation](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Harmonic oscillator equation"
		<< '\n' << "Method" << '\n' << "Implicit Euler's method"
		<< '\n' << "Time interval" << '\n' << spring_equation::odeproblem.timeInterval
		<< '\n' << "Time step" << '\n' << spring_equation::odeproblem.tau
		<< '\n' << "y0" << '\n' << spring_equation::odeproblem.y0.format(INLINE);
}

void No2_symmetric_euler() {
	std::cout
		<< SEPARATOR
		<< STARTER << "2) Harmonic oscillator"
		<< STARTER << "Running symmetric euler method...\n";

	// Solve
	StaticTimer::start();
	const auto[t, y] = odesolve_symmetric_euler(spring_equation::odeproblem);
	StaticTimer::end();

	// Save
	save_as(t, PATH + "[symmetric_euler][spring_equation](t).dat");
	save_as(y, PATH + "[symmetric_euler][spring_equation](y).dat");

	std::ofstream file(PATH + "[symmetric_euler][spring_equation](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Harmonic oscillator equation"
		<< '\n' << "Method" << '\n' << "Symmetric Euler's method"
		<< '\n' << "Time interval" << '\n' << spring_equation::odeproblem.timeInterval
		<< '\n' << "Time step" << '\n' << spring_equation::odeproblem.tau
		<< '\n' << "y0" << '\n' << spring_equation::odeproblem.y0.format(INLINE);
}

void No2_RK4() {
	std::cout
		<< SEPARATOR
		<< STARTER << "2) Harmonic oscillator"
		<< STARTER << "Running RK4 method...\n";

	// Solve
	StaticTimer::start();
	const auto [t, y] = odesolve_RK4(spring_equation::odeproblem);
	StaticTimer::end();

	// Save
	save_as(t, PATH + "[RK4][spring_equation](t).dat");
	save_as(y, PATH + "[RK4][spring_equation](y).dat");

	std::ofstream file(PATH + "[RK4][spring_equation](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Harmonic oscillator equation"
		<< '\n' << "Method" << '\n' << "4th order Runge-Kutta"
		<< '\n' << "Time interval" << '\n' << spring_equation::odeproblem.timeInterval
		<< '\n' << "Time step" << '\n' << spring_equation::odeproblem.tau
		<< '\n' << "y0" << '\n' << spring_equation::odeproblem.y0.format(INLINE);
}

// ¹3
void No3_phase_portrait_1() {
	std::cout
		<< SEPARATOR
		<< STARTER << "3) Phase portrairs"
		<< STARTER << "Running RK4 for system 1...\n";

	// Solve
	StaticTimer::start();

	namespace eq = test1_equation;
	auto problem = eq::odeproblem;

	dvector all_t;
	vvector all_y;
	const auto size = static_cast<size_t>(eq::timeInterval / eq::tau * eq::Nx * eq::Ny);
	all_t.reserve(size);
	all_y.reserve(size);

	// Go over all points in a grid
	for (uint i = 0; i < eq::Ny; ++i)
		for (uint j = 0; j < eq::Nx; ++j) {
			// Simulate short time interval starting from point (i, j)
			problem.y0 = Vector{ { eq::xmin + j * eq::hx, eq::ymin + i * eq::hy } };

			const auto [t, y] = odesolve_RK4(problem);

			all_t.insert(all_t.end(), t.begin(), t.end());
			all_y.insert(all_y.end(), y.begin(), y.end());
		}
	
	StaticTimer::end();

	// Save
	save_as(all_t, PATH + "[phase_portrait_1](t).dat");
	save_as(all_y, PATH + "[phase_portrait_1](y).dat");

	std::ofstream file(PATH + "[phase_portrait_1](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Test equation 1"
		<< '\n' << "Method" << '\n' << "4th order Runge-Kutta"
		<< '\n' << "Grid: x min" << '\n' << eq::xmin
		<< '\n' << "Grid: x max" << '\n' << eq::xmax
		<< '\n' << "Grid: y min" << '\n' << eq::ymin
		<< '\n' << "Grid: y max" << '\n' << eq::ymax
		<< '\n' << "Grid: Nx" << '\n' << eq::Nx
		<< '\n' << "Grid: Ny" << '\n' << eq::Ny
		<< '\n' << "Time interval" << '\n' << problem.timeInterval
		<< '\n' << "Time step" << '\n' << problem.tau;
}

void No3_phase_portrait_2() {
	std::cout
		<< SEPARATOR
		<< STARTER << "3) Phase portrairs"
		<< STARTER << "Running RK4 for system 2...\n";

	// Solve
	StaticTimer::start();

	namespace eq = test2_equation;
	auto problem = eq::odeproblem;

	dvector all_t;
	vvector all_y;
	const auto size = static_cast<size_t>(eq::timeInterval / eq::tau * eq::Nx * eq::Ny);
	all_t.reserve(size);
	all_y.reserve(size);

	// Go over all points in a grid
	for (uint i = 0; i < eq::Ny; ++i)
		for (uint j = 0; j < eq::Nx; ++j) {
			// Simulate short time interval starting from point (i, j)
			problem.y0 = Vector{ { eq::xmin + j * eq::hx, eq::ymin + i * eq::hy } };

			const auto [t, y] = odesolve_RK4(problem);

			all_t.insert(all_t.end(), t.begin(), t.end());
			all_y.insert(all_y.end(), y.begin(), y.end());
		}

	StaticTimer::end();

	// Save
	save_as(all_t, PATH + "[phase_portrait_2](t).dat");
	save_as(all_y, PATH + "[phase_portrait_2](y).dat");

	std::ofstream file(PATH + "[phase_portrait_2](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Test equation 2"
		<< '\n' << "Method" << '\n' << "4th order Runge-Kutta"
		<< '\n' << "Grid: x min" << '\n' << eq::xmin
		<< '\n' << "Grid: x max" << '\n' << eq::xmax
		<< '\n' << "Grid: y min" << '\n' << eq::ymin
		<< '\n' << "Grid: y max" << '\n' << eq::ymax
		<< '\n' << "Grid: Nx" << '\n' << eq::Nx
		<< '\n' << "Grid: Ny" << '\n' << eq::Ny
		<< '\n' << "Time interval" << '\n' << problem.timeInterval
		<< '\n' << "Time step" << '\n' << problem.tau;
}

void No3_phase_portrait_3() {
	std::cout
		<< SEPARATOR
		<< STARTER << "3) Phase portrairs"
		<< STARTER << "Running RK4 for system 3...\n";

	// Solve
	StaticTimer::start();

	namespace eq = test3_equation;
	auto problem = eq::odeproblem;

	dvector all_t;
	vvector all_y;
	const auto size = static_cast<size_t>(eq::timeInterval / eq::tau * eq::Nx * eq::Ny);
	all_t.reserve(size);
	all_y.reserve(size);

	// Go over all points in a grid
	for (uint i = 0; i < eq::Ny; ++i)
		for (uint j = 0; j < eq::Nx; ++j)
			for (uint k = 0; k < eq::Nz; ++k) {
				// Simulate short time interval starting from point (i, j)
				problem.y0 = Vector{ { eq::xmin + j * eq::hx, eq::ymin + i * eq::hy,  eq::zmin + k * eq::hz } };

				const auto [t, y] = odesolve_RK4(problem);

				all_t.insert(all_t.end(), t.begin(), t.end());
				all_y.insert(all_y.end(), y.begin(), y.end());
			}

	StaticTimer::end();

	// Save
	save_as(all_t, PATH + "[phase_portrait_3](t).dat");
	save_as(all_y, PATH + "[phase_portrait_3](y).dat");

	std::ofstream file(PATH + "[phase_portrait_3](problem).dat");
	file
		<< '\n' << "Problem" << '\n' << "Test equation 3"
		<< '\n' << "Method" << '\n' << "4th order Runge-Kutta"
		<< '\n' << "Grid: x min" << '\n' << eq::xmin
		<< '\n' << "Grid: x max" << '\n' << eq::xmax
		<< '\n' << "Grid: y min" << '\n' << eq::ymin
		<< '\n' << "Grid: y max" << '\n' << eq::ymax
		<< '\n' << "Grid: z min" << '\n' << eq::zmin
		<< '\n' << "Grid: z max" << '\n' << eq::zmax
		<< '\n' << "Grid: Nx" << '\n' << eq::Nx
		<< '\n' << "Grid: Ny" << '\n' << eq::Ny
		<< '\n' << "Grid: Nz" << '\n' << eq::Ny
		<< '\n' << "Time interval" << '\n' << problem.timeInterval
		<< '\n' << "Time step" << '\n' << problem.tau;
}

int main(int argc, char *argv[]) {
	// ¹1
	No1_adaptiveRK4();
	No1_RK4();
	No1_RK2();
	No1_PredictorCorrector4();
	No1_Adams4();
	No1_ExplicitEuler();
	No1_ImplicitEuler();
	No1_SymmetricEuler();

	// ¹2
	No2_explicit_euler();
	No2_implicit_euler();
	No2_symmetric_euler();
	No2_RK4();

	// ¹3
	No3_phase_portrait_1();
	No3_phase_portrait_2();
	No3_phase_portrait_3();

	return 0;
}