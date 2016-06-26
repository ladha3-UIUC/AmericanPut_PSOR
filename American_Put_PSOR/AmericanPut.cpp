#define SOR_HARD_STOP 500

#include <iostream>
#include <algorithm>
#include "AmericanPut.h"

double ** American::putPrice() {
	double theta = 0.5; // weight parameter, 0.5 for Crank-Nicolson
	double dS = S_max / (N - 1);
	double dt = T / M;

	double ** result;
	result = new double *[M + 1];
	for (int i = 0; i < M + 1; i++) {
		result[i] = new double[N];
		for (int j = 0; j < N; j++)
			result[i][j] = 0;
	}

	gsl_matrix * M1 = gsl_matrix_calloc(N, N);
	gsl_matrix * M2 = gsl_matrix_calloc(N, N);

	// set up tridiagonal matrices
	double d, l, u;
	double alpha = pow(sigma, 2) * dt;
	double beta = r * dt;

	for (int i = 0; i < N; i++) {
		d = 1 + theta * (alpha * pow(i + 1, 2) + beta);
		gsl_matrix_set(M1, i, i, d);
		double d = 1 - (1 - theta) * (alpha * pow(i + 1, 2) + beta);
		gsl_matrix_set(M2, i, i, d);
	}

	for (int i = 0; i < N - 1; i++) {
		l = -0.5 * theta * (alpha * pow(i + 2, 2) - beta * (i + 2));
		gsl_matrix_set(M1, i + 1, i, l);
		l = 0.5 * (1 - theta) * (alpha * pow(i + 2, 2) - beta * (i + 2));
		gsl_matrix_set(M2, i + 1, i, l);
	}

	for (int i = 0; i < N - 1; i++) {
		u = -0.5 * theta * (alpha * pow(i + 1, 2) + beta * (i + 1));
		gsl_matrix_set(M1, i, i + 1, u);
		u = 0.5 * (1 - theta) * (alpha * pow(i + 1, 2) + beta * (i + 1));
		gsl_matrix_set(M2, i, i + 1, u);
	}

	// boundary conditions
	gsl_matrix_set(M1, N - 1, N - 1, 1);
	gsl_matrix_set(M1, N - 1, N - 2, 0);
	gsl_matrix_set(M2, N - 1, N - 1, 1);
	gsl_matrix_set(M2, N - 1, N - 2, 0);

	gsl_vector * w = gsl_vector_alloc(N);
	gsl_vector * v = gsl_vector_alloc(N);
	gsl_vector * g = gsl_vector_alloc(N);
	gsl_vector * b = gsl_vector_alloc(N);

	// initial value for American put
	for (int i = 0; i < N; i++) {
		if (i * dS >= K) {
			gsl_vector_set(g, i, 0);
		}
		else {
			gsl_vector_set(g, i, K - i * dS);
		}
	}
	gsl_vector_memcpy(w, g);

	// step through time
	for (int t = 0; t < M + 1; t++) {
		memcpy(result[t], w->data, w->size * sizeof(double));

		// solve Cryer's problem
		gsl_blas_dgemv(CblasNoTrans, 1.0, M2, w, 0.0, b);
		projectedSOR(b, g, v, M1);
		gsl_vector_swap(w, v);
	}

	gsl_matrix_free(M1);
	gsl_matrix_free(M2);
	gsl_vector_free(w);
	gsl_vector_free(v);
	gsl_vector_free(g);
	gsl_vector_free(b);

	return result;
};

void American::projectedSOR(gsl_vector * f, gsl_vector * g, gsl_vector * v, const gsl_matrix * A) {
	int N = (int)v->size;

	gsl_matrix * cpyA = gsl_matrix_alloc(N, N);
	gsl_matrix * invA = gsl_matrix_alloc(N, N);
	gsl_matrix_memcpy(cpyA, A);

	// find A-inverse
	int s;
	gsl_permutation * p = gsl_permutation_alloc(N);
	gsl_linalg_LU_decomp(cpyA, p, &s);
	gsl_linalg_LU_invert(cpyA, p, invA);

	// initialize v
	gsl_blas_dgemv(CblasNoTrans, 1.0, invA, f, 0.0, v);
	for (int i = 0; i < N; i++) {
		double v_i = gsl_vector_get(v, i);
		double g_i = gsl_vector_get(g, i);
		gsl_vector_set(v, i, std::max(v_i, g_i));
	}

	gsl_vector * v_new = gsl_vector_calloc(N);
	gsl_vector * v_epsilon = gsl_vector_alloc(N);

	bool converged = false;
	for (int k = 0; k <= SOR_HARD_STOP; k++) {
		// forward substitution
		for (int i = 1; i <= N; i++) {
			double sum1 = 0;
			double sum2 = 0;

			for (int j = 1; j <= i - 1; j++) {
				sum1 += gsl_matrix_get(A, i - 1, j - 1) * gsl_vector_get(v_new, j - 1);
			}
			for (int j = i + 1; j <= N; j++) {
				sum2 += gsl_matrix_get(A, i - 1, j - 1) * gsl_vector_get(v, j - 1);
			}

			double a_ii = gsl_matrix_get(A, i - 1, i - 1);
			double f_i = gsl_vector_get(f, i - 1);
			double t_i = (f_i - sum1 - sum2) / a_ii;
			double v_i = gsl_vector_get(v, i - 1);
			double g_i = gsl_vector_get(g, i - 1);
			gsl_vector_set(v_new, i - 1, std::max(v_i + omega * (t_i - v_i), g_i));
		}

		gsl_vector_memcpy(v, v_new);

		// measure how much v has changed
		gsl_vector_memcpy(v_epsilon, v_new);
		gsl_vector_sub(v_epsilon, v); // subtracts. v_epsilon - v stores in x_epsilon
		gsl_vector_mul(v_epsilon, v_epsilon); // make everything positive

		//std::cout << "max epsilon " << gsl_vector_max(v_new) << std::endl;

		// decide whether to stop the iteration.
		if (gsl_vector_max(v_epsilon) < epsilon_stop * epsilon_stop) {
			converged = true;
			break;
		}
	}

	if (!converged)
		std::cout << "Warning: projected SOR DID NOT CONVERGE!!!" << std::endl;

	gsl_matrix_free(cpyA);
	gsl_matrix_free(invA);
	gsl_vector_free(v_new);
	gsl_vector_free(v_epsilon);
};
