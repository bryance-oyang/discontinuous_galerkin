/*
 * To solve D_t u + D_x (u^2 / 2) = 0
 *
 * D_t c_m = (2m + 1) / dx * (-[wm f] + 1/2 Sum_n Sum_l c_n c_p DPPP_{mnl})
 * where
 * DPPP_{mnl} = integral of dP_m/dx P_n P_l (Legendre poly)
 *
 * Shape functions are w_m = P_m(2(x - a) / dx) where a is the cell center coord
 */

#ifdef DEBUG
#define _GNU_SOURCE 1
#include <fenv.h>
#endif /* DEBUG */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <math.h>

#include "ws_ctube.h"

#define NGHOST 1
#define ORDER 2
#define CFL (0.5 * (1.0 / (2*(3) + 1)))
#define WAVE_AMP 1e-6
#define MAX_TIME 1

/* Legendre P at x = 0 */
double Pcc[3];
/* integral of dP_m/dx P_n P_l */
double DPPP[3][3][3];
double ssprk3_w[3][3];

struct ws_ctube *ctube;

double hlle(double ul, double ur, double fl, double fr, double ca, double cb)
{
	double cl = fmin(ca, cb);
	double cr = fmax(ca, cb);

	if (cl >= 0) {
		return fl;
	} else if (cr <= 0) {
		return fr;
	} else {
		return (cr*cl*(ur - ul) + cr*fl - cl*fr) / (cr - cl);
	}
}

struct sim {
	int ncell;
	double dx;

	double *c0;
	double *c;

	double *ul;
	double *ur;

	double *dcdt;
	double *J;

	double t;
	double dt;
};

void init_weights()
{
	ssprk3_w[0][0] = 1;
	ssprk3_w[0][1] = 0;
	ssprk3_w[0][2] = 1;

	ssprk3_w[1][0] = 0.75;
	ssprk3_w[1][1] = 0.25;
	ssprk3_w[1][2] = 0.25;

	ssprk3_w[2][0] = 1.0/3.0;
	ssprk3_w[2][1] = 2.0/3.0;
	ssprk3_w[2][2] = 2.0/3.0;

	Pcc[0] = 1;
	Pcc[2] = -0.5;

	for (int n = 0; n < ORDER; n++) {
		DPPP[1][n][n] = 2.0 / (2*n + 1);
	}

	DPPP[2][0][1] = 2;
	DPPP[2][1][0] = 2;

	DPPP[2][1][2] = 0.8;
	DPPP[2][2][1] = 0.8;
}

/* (2m+1) / 2 * P_m(2(x-a)/dx) sin(2pi(x - t)) */
double integral_Psin(double m, double a, double dx, double t)
{
	double answer = 0;
	if (m == 0) {
		answer = sin(2*M_PI*(a - t)) * sin(dx*M_PI) / M_PI;
	} else if (m == 1) {
		answer = cos(2*M_PI*(a-t)) * (sin(dx*M_PI) - dx*M_PI*cos(dx*M_PI)) / (dx * M_PI * M_PI);
	} else if (m == 2) {
		answer = sin(2*M_PI*(a-t)) * ((dx*dx*M_PI*M_PI - 3)*sin(dx*M_PI) + 3*dx*M_PI*cos(dx*M_PI)) / (dx*dx * M_PI*M_PI*M_PI);
	}

	return answer * (2*m + 1) / dx;
}

void init_cond(struct sim *restrict sim)
{
	sim->t = 0;

	for (int i = 0; i < sim->ncell + 2*NGHOST; i++) {
		double x = (i - NGHOST) * sim->dx + 0.5*sim->dx;

		sim->c[ORDER*i + 0] = 1.0 + WAVE_AMP * integral_Psin(0, x, sim->dx, 0);
		for (int m = 1; m < ORDER; m++) {
			sim->c[ORDER*i + m] = WAVE_AMP * integral_Psin(m, x, sim->dx, 0);
		}
	}
}

void boundary(struct sim *restrict sim)
{
	int ncell = sim->ncell;

	/* periodic */
	for (int k = 0; k < NGHOST; k++) {
		for (int m = 0; m < ORDER; m++) {
			sim->c[ORDER*(k) + m] = sim->c[ORDER*(ncell + k) + m];
			sim->c[ORDER*(NGHOST + ncell + k) + m] = sim->c[ORDER*(NGHOST + k) + m];
		}
	}
}

void compute_ulur(struct sim *restrict sim)
{
	for (int i = 0; i < sim->ncell + 2*NGHOST; i++) {
		sim->ul[i] = 0;
		sim->ur[i] = 0;

		for (int m = 0, parity = 1; m < ORDER; m++, parity *= -1) {
			sim->ul[i] += sim->c[ORDER*i + m] * parity;
			sim->ur[i] += sim->c[ORDER*i + m];
		}
	}
}

void determine_dt(struct sim *restrict sim)
{
	double max_abs_wavespeed = 0;
	compute_ulur(sim);
	for (int i = 0; i < sim->ncell + 2*NGHOST; i++) {
		max_abs_wavespeed = fmax(fabs(sim->ul[i]), max_abs_wavespeed);
		max_abs_wavespeed = fmax(fabs(sim->ur[i]), max_abs_wavespeed);
	}

	sim->dt = sim->dx / max_abs_wavespeed;
	sim->dt *= CFL;
}

void compute_J(struct sim *restrict sim)
{
	compute_ulur(sim);

	for (int i = 0; i < sim->ncell + 1; i++) {
		double ul = sim->ur[NGHOST + i - 1];
		double ur = sim->ul[NGHOST + i];
		double fl = 0.5 * ul * ul;
		double fr = 0.5 * ur * ur;
		sim->J[NGHOST + i] = hlle(ul, ur, fl, fr, ul, ur);
	}
}

void compute_dcdt(struct sim *restrict sim)
{
	compute_J(sim);

	for (int i = 0; i < sim->ncell; i++) {
		for (int m = 0, parity = 1; m < ORDER; m++, parity *= -1) {
			double flux = parity * sim->J[NGHOST + i] - sim->J[NGHOST + i + 1];

			double tensor = 0;
			for (int n = 0; n < ORDER; n++) {
				double cn = sim->c[ORDER*(NGHOST + i) + n];
				for (int l = 0; l < ORDER; l++) {
					double cl = sim->c[ORDER*(NGHOST + i) + l];
					tensor += 0.5 * cn * cl
						* DPPP[m][n][l];
				}
			}

			sim->dcdt[ORDER*(NGHOST + i) + m] = flux + tensor;
			sim->dcdt[ORDER*(NGHOST + i) + m] *= (2*m + 1) / sim->dx;
		}
	}
}

void take_big_timestep(struct sim *restrict sim)
{
	memcpy(sim->c0, sim->c, ORDER * (sim->ncell + 2*NGHOST) * sizeof(sim->c));

	for (int s = 0; s < 3; s++) {
		compute_dcdt(sim);
		for (int i = 0; i < sim->ncell; i++) {
			for (int m = 0; m < ORDER; m++) {
				sim->c[ORDER*(NGHOST + i) + m] =
					ssprk3_w[s][0] * sim->c0[ORDER*(NGHOST + i) + m] +
					ssprk3_w[s][1] * sim->c[ORDER*(NGHOST + i) + m] +
					ssprk3_w[s][2] * sim->dcdt[ORDER*(NGHOST + i) + m] * sim->dt;
			}
		}
		boundary(sim);
	}

	sim->t += sim->dt;
}

void l1_error(struct sim *restrict sim)
{
	double error = 0;
	for (int i = 0; i < sim->ncell; i++) {
		double x = i*sim->dx + 0.5*sim->dx;

		double numerical = sim->c[ORDER*(NGHOST + i) + 0];
		double actual = 1 + WAVE_AMP * sin(sim->dx * M_PI) * sin(2*M_PI*(x - sim->t)) / (sim->dx * M_PI);
		error += fabs(numerical - actual) * sim->dx;
	}

	printf("ncell: %d\tlog error: %lf\tlog dx: %lf\n", sim->ncell, log2(error), log2(sim->dx));
}

void print_c(struct sim *restrict sim)
{
	for (int i = 0; i < sim->ncell + 2*NGHOST; i++) {
		for (int m = 0; m < ORDER; m++) {
			printf("c[%02d][%d]=%lf\t", i, m, sim->c[ORDER*i + m]);
		}
		printf("\n");
	}
}

void run_sim(int ncell)
{
	struct sim sim;

	sim.ncell = ncell;
	sim.dx = 1.0 / ncell;

	sim.c0 = malloc(ORDER * (ncell + 2*NGHOST) * sizeof(*sim.c0));
	sim.c = malloc(ORDER * (ncell + 2*NGHOST) * sizeof(*sim.c));

	sim.ul = malloc((ncell + 2*NGHOST) * sizeof(*sim.c));
	sim.ur = malloc((ncell + 2*NGHOST) * sizeof(*sim.c));

	sim.dcdt = malloc(ORDER * (ncell + 2*NGHOST) * sizeof(*sim.dcdt));
	sim.J = malloc((ncell + 2*NGHOST + 1) * sizeof(*sim.J));

	if (sim.c0 == NULL || sim.c == NULL || sim.dcdt == NULL || sim.J == NULL) {
		raise(SIGSEGV);
	}

	init_cond(&sim);
	boundary(&sim);

	for (;;) {
		//print_c(&sim);

		determine_dt(&sim);
		//printf("t = %lf\tdt = %lf\n", sim.t, sim.dt);

		take_big_timestep(&sim);

		if (sim.t >= MAX_TIME) {
			break;
		}
	}

	l1_error(&sim);

	free(sim.c0);
	free(sim.c);
	free(sim.ul);
	free(sim.ur);
	free(sim.dcdt);
	free(sim.J);
}

int main()
{
#ifdef DEBUG
	feenableexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
#endif /* DEBUG */

	init_weights();

	ctube = ws_ctube_open(9743, 2, 0, 24);

	for (int ncell = 2; ncell <= 2048; ncell *= 2) {
		run_sim(ncell);
	}

	ws_ctube_close(ctube);

	return 0;
}
