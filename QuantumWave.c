
/*
***************************
 * MADE BY X-3306 
 * NOTES
 *
 *  This program uses the 4th-order Runge-Kutta method to
 *  solve the Schrodinger equation for a given potential.
 *  This is done by computing the derivatives of the wavefunction
 *  at each step and using these to calculate the next step.
 *
 *  The wavefunction is a complex array, i.e. it contains an
 *  array of real parts and an array of imaginary parts.
 *
 *  The program can take either a fixed potential, or a file
 *  of potential values as input.
 *
 *  The output is the wavefunction at each step, which can
 *  then be used to calculate expectation values such as
 *  energy, momentum, etc.
 ***************************
 */

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <complex.h>

// define constants
#define h 6.626e-34
#define pi 3.14159

// define global variables
double *pot; // potential array
int N; // number of grid points
int nsteps; // number of timesteps
double dx; // grid spacing
double dt; // timestep

// define complex arrays
double complex *psi_real;
double complex *psi_imag;
double complex *psi_real_temp;
double complex *psi_imag_temp;

// define functions
void init(char* filename);
void rk4();
void print();

int main (int argc, char *argv[]) {
  if (argc != 3) {
    printf("Usage: %s <potential_file> <timesteps>\n", argv[0]);
    return 1;
  }

  // read in arguments
  char *filename = argv[1];
  nsteps = atoi(argv[2]);

  // initialize arrays
  init(filename);

  // timestep
  for (int i = 0; i < nsteps; i++) {
    rk4();
  }

  // output results
  print();

  // free memory
  free(pot);
  free(psi_real);
  free(psi_imag);
  free(psi_real_temp);
  free(psi_imag_temp);

  return 0;
}

void init (char* filename) {
  // read in potential
  FILE *file;
  file = fopen(filename, "r");
  if (file == NULL) {
    fprintf(stderr, "Error opening file: %s\n", strerror(errno));
    exit(EXIT_FAILURE);
  }

  // get number of grid points
  fscanf(file, "%d", &N);

  // allocate memory for potential
  pot = (double*) malloc(N * sizeof(double));

  // read in potential
  for (int i = 0; i < N; i++) {
    fscanf(file, "%lf", &pot[i]);
  }

  // close file
  fclose(file);

  // set grid spacing
  dx = 1.0 / (N - 1);

  // set timestep
  dt = h / (4 * pi * pot[N/2]);

  // allocate memory for wavefunction
  psi_real = (double complex*) malloc(N * sizeof(double complex));
  psi_imag = (double complex*) malloc(N * sizeof(double complex));
  psi_real_temp = (double complex*) malloc(N * sizeof(double complex));
  psi_imag_temp = (double complex*) malloc(N * sizeof(double complex));

  // set initial conditions
  for (int i = 0; i < N; i++) {
    psi_real[i] = 0.0 + 0.0*I;
    psi_imag[i] = 0.0 + 0.0*I;
  }
  psi_real[N/2] = 1.0 + 0.0*I;
  psi_imag[N/2] = 0.0 + 0.0*I;
}

void rk4 () {
  // calculate dpsi/dx
  double complex *dpsidx_real = (double complex*) malloc(N * sizeof(double complex));
  double complex *dpsidx_imag = (double complex*) malloc(N * sizeof(double complex));

  // central difference for first and last points
  dpsidx_real[0] = (psi_real[1] - psi_real[0]) / (2 * dx);
  dpsidx_imag[0] = (psi_imag[1] - psi_imag[0]) / (2 * dx);
  dpsidx_real[N-1] = (psi_real[N-1] - psi_real[N-2]) / (2 * dx);
  dpsidx_imag[N-1] = (psi_imag[N-1] - psi_imag[N-2]) / (2 * dx);

  // central difference for interior points
  for (int i = 1; i < N-1; i++) {
    dpsidx_real[i] = (psi_real[i+1] - psi_real[i-1]) / (2 * dx);
    dpsidx_imag[i] = (psi_imag[i+1] - psi_imag[i-1]) / (2 * dx);
  }

  // calculate k1
  double complex *k1_real = (double complex*) malloc(N * sizeof(double complex));
  double complex *k1_imag = (double complex*) malloc(N * sizeof(double complex));
  for (int i = 0; i < N; i++) {
    k1_real[i] = -I * (pot[i] * psi_real[i] + dpsidx_imag[i]) * dt;
    k1_imag[i] = -I * (pot[i] * psi_imag[i] - dpsidx_real[i]) * dt;
  }

  // calculate k2
  double complex *k2_real = (double complex*) malloc(N * sizeof(double complex));
  double complex *k2_imag = (double complex*) malloc(N * sizeof(double complex));
  for (int i = 0; i < N; i++) {
    psi_real_temp[i] = psi_real[i] + 0.5 * k1_real[i];
    psi_imag_temp[i] = psi_imag[i] + 0.5 * k1_imag[i];
  }

  // calculate dpsi_temp/dx
  double complex *dpsi_temp_real = (double complex*) malloc(N * sizeof(double complex));
  double complex *dpsi_temp_imag = (double complex*) malloc(N * sizeof(double complex));

  // central difference for first and last points
  dpsi_temp_real[0] = (psi_real_temp[1] - psi_real_temp[0]) / (2 * dx);
  dpsi_temp_imag[0] = (psi_imag_temp[1] - psi_imag_temp[0]) / (2 * dx);
  dpsi_temp_real[N-1] = (psi_real_temp[N-1] - psi_real_temp[N-2]) / (2 * dx);
  dpsi_temp_imag[N-1] = (psi_imag_temp[N-1] - psi_imag_temp[N-2]) / (2 * dx);

  // central difference for interior points
  for (int i = 1; i < N-1; i++) {
    dpsi_temp_real[i] = (psi_real_temp[i+1] - psi_real_temp[i-1]) / (2 * dx);
    dpsi_temp_imag[i] = (psi_imag_temp[i+1] - psi_imag_temp[i-1]) / (2 * dx);
  }

  for (int i = 0; i < N; i++) {
    k2_real[i] = -I * (pot[i] * psi_real_temp[i] + dpsi_temp_imag[i]) * dt;
    k2_imag[i] = -I * (pot[i] * psi_imag_temp[i] - dpsi_temp_real[i]) * dt;
  }

  // calculate k3
  double complex *k3_real = (double complex*) malloc(N * sizeof(double complex));
  double complex *k3_imag = (double complex*) malloc(N * sizeof(double complex));
  for (int i = 0; i < N; i++) {
    psi_real_temp[i] = psi_real[i] + 0.5 * k2_real[i];
    psi_imag_temp[i] = psi_imag[i] + 0.5 * k2_imag[i];
  }

  // calculate dpsi_temp/dx
  for (int i = 0; i < N; i++) {
    dpsi_temp_real[i] = (psi_real_temp[i+1] - psi_real_temp[i-1]) / (2 * dx);
    dpsi_temp_imag[i] = (psi_imag_temp[i+1] - psi_imag_temp[i-1]) / (2 * dx);
  }

  for (int i = 0; i < N; i++) {
    k3_real[i] = -I * (pot[i] * psi_real_temp[i] + dpsi_temp_imag[i]) * dt;
    k3_imag[i] = -I * (pot[i] * psi_imag_temp[i] - dpsi_temp_real[i]) * dt;
  }

  // calculate k4
  double complex *k4_real = (double complex*) malloc(N * sizeof(double complex));
  double complex *k4_imag = (double complex*) malloc(N * sizeof(double complex));
  for (int i = 0; i < N; i++) {
    psi_real_temp[i] = psi_real[i] + k3_real[i];
    psi_imag_temp[i] = psi_imag[i] + k3_imag[i];
  }

  // calculate dpsi_temp/dx
  for (int i = 0; i < N; i++) {
    dpsi_temp_real[i] = (psi_real_temp[i+1] - psi_real_temp[i-1]) / (2 * dx);
    dpsi_temp_imag[i] = (psi_imag_temp[i+1] - psi_imag_temp[i-1]) / (2 * dx);
  }

  for (int i = 0; i < N; i++) {
    k4_real[i] = -I * (pot[i] * psi_real_temp[i] + dpsi_temp_imag[i]) * dt;
    k4_imag[i] = -I * (pot[i] * psi_imag_temp[i] - dpsi_temp_real[i]) * dt;
  }

  // calculate psi at next time step
  for (int i = 0; i < N; i++) {
    psi_real[i] = psi_real[i] + (1.0/6.0) * (k1_real[i] + 2*k2_real[i] + 2*k3_real[i] + k4_real[i]);
    psi_imag[i] = psi_imag[i] + (1.0/6.0) * (k1_imag[i] + 2*k2_imag[i] + 2*k3_imag[i] + k4_imag[i]);
  }

  // free memory
  free(dpsidx_real);
  free(dpsidx_imag);
  free(k1_real);
  free(k1_imag);
  free(k2_real);
  free(k2_imag);
  free(k3_real);
  free(k3_imag);
  free(k4_real);
  free(k4_imag);
  free(dpsi_temp_real);
  free(dpsi_temp_imag);
}

void print () {
  // print wavefunction
  printf("# x\treal\timag\n");
  for (int i = 0; i < N; i++) {
    double x = i * dx;
    printf("%lf\t%lf\t%lf\n", x, creal(psi_real[i]), cimag(psi_imag[i]));
  }
}
