#pragma once

#include "core_math.hpp"


inline Vector nonlinear_solve(VectorFunction F, Vector X0, double precision, uint maxIterations = 100) {
	Vector X = X0;

	// Fill approximations
	for (uint iterations = 0; iterations < maxIterations; ++iterations) {
		X = X0 - jacobian(F, X0).inverse() * F(X0);

		if ((X - X0).norm() < precision) break;

		X0 = X;
	}

	return X;
}