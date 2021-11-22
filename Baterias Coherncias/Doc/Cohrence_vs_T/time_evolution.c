#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "matrixoperations.h"
#include "RK4.h"
#include "propagators.h"
#include "energy.h"
#include "dephasing.h"
#include "pauli.h"
#include "hamiltonian.h"
void Open_evolution(double complex *state, double complex *hamiltonian, double complex *bath_state, double complex *bath_hamiltonian, double complex *interaction, int dim, double tf, double dt, double *E, double *Q, double *W, int prev_steps, double complex *hot_bath_state)
{
    int N = (int)(tf/dt), step;
    double complex *propagator, *dissipator;
    double E_0;
    propagator = (double complex*) calloc(dim*dim, sizeof(double complex));
    dissipator = (double complex*) calloc(dim*dim, sizeof(double complex));
    if(prev_steps == 0)
    {
        *E  = calc_Energy(state, hamiltonian, dim);
        *Q = 0;
        *W = 0;
    }
    E_0 = calc_Energy(hot_bath_state, hamiltonian, dim);
    for(step = prev_steps + 1; step < N; step++)
    {
        RK4_open(propagator, dissipator, state, hamiltonian, bath_state, bath_hamiltonian, interaction, dt, dim);
        *(E + step) = calc_Energy(state, hamiltonian, dim) - E_0;
        *(Q + step) = calc_Heat(*(Q + step - 1), state, bath_state, bath_hamiltonian, interaction, dim, dt);
        *(W + step) = calc_Work(*(E + step), *(Q + step));
    }
    if(prev_steps == 0)
    {
        *E = 0;
    }
    free(propagator);
    free(dissipator);
}
void Driven_evolution(double complex *state, double complex *hamiltonian, double complex *bath_state, double complex *bath_hamiltonian, double complex *interaction, int dim, double tf, double dt, double *E, double *Q, double *W)
{
    int N = (int)(tf/dt), step;
    double complex *propagator, *dissipator, *hamiltonian_z;
    propagator = (double complex*) calloc(dim*dim, sizeof(double complex));
    dissipator = (double complex*) calloc(dim*dim, sizeof(double complex));
    hamiltonian_z = (double complex*) calloc(dim*dim, sizeof(double complex));
    Pauli_z(hamiltonian_z, *(hamiltonian + 1));
    *E  = calc_Energy(state, hamiltonian_z, dim);
    *Q = 0;
    *W = 0;
    for(step = 1; step < N; step+=2)
    {
        RK4_open(propagator, dissipator, state, hamiltonian, bath_state, bath_hamiltonian, interaction, dt, dim);
        *(E + step) = calc_Energy(state, hamiltonian, dim) - *E;
        *(Q + step) = calc_Heat(*(Q + step - 1), state, bath_state, bath_hamiltonian, interaction, dim, dt);
        *(W + step) = calc_Work(*(E + step), *(Q + step));
        RK4_closed(propagator, state, hamiltonian, dt, dim);
        *(E + step + 1) = calc_Energy(state, hamiltonian, dim) - *E;
        *(Q + step + 1) = *(Q + step);
        *(W + step + 1) = calc_Work(*(E + step + 1), *(Q + step + 1));
    }
    *E = 0;
    free(propagator);
    free(dissipator);
    free(hamiltonian_z);
}
void Driven_evolution_dephasing(double complex *state, double complex *hamiltonian, double complex *bath_state, double complex *bath_hamiltonian, double complex *interaction, int dim, double tf, double dt, double *E, double *Q, double *W, double p)
{
   int N = (int)(tf/dt), step;
    double complex *propagator, *dissipator, *hamiltonian_z, *change_basis_x;
    propagator = (double complex*) calloc(dim*dim, sizeof(double complex));
    dissipator = (double complex*) calloc(dim*dim, sizeof(double complex));
    hamiltonian_z = (double complex*) calloc(dim*dim, sizeof(double complex));
    change_basis_x = (double complex*) calloc(dim*dim, sizeof(double complex));
    Change_basis_x(change_basis_x);
    Pauli_z(hamiltonian_z, *(hamiltonian + 1));
    *E  = calc_Energy(state, hamiltonian_z, dim);
    *Q = 0;
    *W = 0;
    for(step = 1; step < N; step+=2)
    {
        RK4_open(propagator, dissipator, state, hamiltonian, bath_state, bath_hamiltonian, interaction, dt, dim);
        *(E + step) = calc_Energy(state, hamiltonian, dim) - *E;
        *(Q + step) = calc_Heat(*(Q + step - 1), state, bath_state, bath_hamiltonian, interaction, dim, dt);
        *(W + step) = calc_Work(*(E + step), *(Q + step));
        RK4_closed(propagator, state, hamiltonian, dt, dim);
        *(E + step + 1) = calc_Energy(state, hamiltonian, dim) - *E;
        *(Q + step + 1) = *(Q + step);
        *(W + step + 1) = calc_Work(*(E + step + 1), *(Q + step + 1));
        Dephasing_noise(state, change_basis_x, p, dim);
    }
    *E = 0;
    free(propagator);
    free(dissipator);
    free(hamiltonian_z);
    free(change_basis_x);
}
void Closed_evolution(double complex *state, double complex *hamiltonian, int dim, double tf, double dt)
{
    int N = (int)(tf/dt), step;
    double complex *propagator;
    propagator = (double complex*) calloc(dim*dim, sizeof(double complex));
    for(step = 0; step < N; step++)
    {
        RK4_closed(propagator, state, hamiltonian, dt, dim);
    }
    free(propagator);
}
