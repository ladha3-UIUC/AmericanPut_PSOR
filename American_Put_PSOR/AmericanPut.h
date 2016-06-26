#pragma once

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

class American {
public:
	double ** putPrice();

	// parameters
	int N, M;
	double r, sigma, K, S_max, T;
	double omega;
	const double epsilon_stop = 0.0000001;

	// constructors
	American() {
		N = 10;
		M = 100;
		r = 0.01;
		sigma = 0.2;
		K = 100;
		S_max = 100;
		T = 1;
		omega = 1.0;
	}

	American(int iN, int iM, double ir, double isigma, double iK, double iS, double iT, double iOmega) {
		N = iN;
		M = iM;

		r = ir;
		sigma = isigma;
		K = iK;
		S_max = iS;
		T = iT;

		omega = iOmega;
	}

private:
	void projectedSOR(gsl_vector * f, gsl_vector * g, gsl_vector * v, const gsl_matrix * A);
};
