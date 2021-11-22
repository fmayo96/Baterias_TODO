#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "matrixoperations.h"

double calc_Energy(double complex *state, double complex *hamiltonian, int dim)
{
    double E = 0;
    double complex *state_hamiltonian;
    state_hamiltonian = (double complex *) calloc(dim*dim, sizeof(double complex));
    dot(state_hamiltonian, state, hamiltonian, dim);
    E = trace(state_hamiltonian, dim);
    free(state_hamiltonian);
    return E;
}
double calc_Heat(double heat_old, double complex *state, double complex *bath_state, double complex *bath_hamiltonian, double complex *interaction, int dim, double dt)
{
    int i;
    double complex *product_state, *product_hamiltonian, *eye, *comm, *double_comm, *product_double_comm_state;
    double heat;
    product_state = (double complex*) calloc(dim*dim*dim*dim, sizeof(double complex));
    product_hamiltonian = (double complex*) calloc(dim*dim*dim*dim, sizeof(double complex));
    comm = (double complex*) calloc(dim*dim*dim*dim, sizeof(double complex));
    double_comm = (double complex*) calloc(dim*dim*dim*dim, sizeof(double complex));
    product_double_comm_state = (double complex*) calloc(dim*dim*dim*dim, sizeof(double complex));
    eye = (double complex*) calloc(dim*dim, sizeof(double complex));
    for(i = 0; i < dim; i++)
    {
        *(eye + dim*i + i) = 1;
    }
    kron(product_hamiltonian, eye, bath_hamiltonian, dim);   
    kron(product_state, state, bath_state, dim);
    commutator(comm, interaction, product_hamiltonian, dim*dim);
    commutator(double_comm, interaction, comm, dim*dim);
    dot(product_double_comm_state, double_comm, product_state, dim*dim);
    heat = heat_old + 0.5 * dt * creal(trace(product_double_comm_state, dim*dim));
    free(product_state);
    free(product_hamiltonian);
    free(eye);
    free(comm);
    free(double_comm);
    free(product_double_comm_state);
    return heat;
} 
double calc_Work(double energy, double heat)
{
    double work = heat - energy;
    return work;
}
