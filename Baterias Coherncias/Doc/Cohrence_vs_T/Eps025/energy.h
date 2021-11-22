#ifndef ENERGY_H
#define ENERGY_H
double calc_Energy(double complex *state, double complex *hamiltonian, int dim);
double calc_Heat(double heat_old, double complex *state, double complex *bath_state, double complex *bath_hamiltonian, double complex *interaction, int dim, double dt);
double calc_Work(double energy, double heat);
#endif