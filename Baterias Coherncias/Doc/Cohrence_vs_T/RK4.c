#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "matrixoperations.h"
#include "propagators.h"
void RK4_open(double complex *propagator, double complex *dissipator, double complex *state, double complex *hamiltonian, double complex *bath_state, double complex *bath_hamiltinian, double complex *interaction, double dt, int dim)
{
    int i;
    double complex *K1, *K2, *K3, *K4, *aux_state;
    K1 = (double complex*) calloc(dim*dim, sizeof(double complex));
    K2 = (double complex*) calloc(dim*dim, sizeof(double complex));
    K3 = (double complex*) calloc(dim*dim, sizeof(double complex));
    K4 = (double complex*) calloc(dim*dim, sizeof(double complex));
    aux_state = (double complex*) calloc(dim*dim, sizeof(double complex));
    Dissipator(dissipator, state, bath_state, interaction, dim);
    Propagator(propagator, state, hamiltonian, dissipator, dim);
    for(i = 0; i < dim*dim; i++)
    {
        *(K1 + i) = *(propagator + i) * dt;
        *(aux_state + i) = *(state + i) + 0.5 * *(K1 + i);
    }
    Dissipator(dissipator, aux_state, bath_state, interaction, dim);
    Propagator(propagator, aux_state, hamiltonian, dissipator, dim);
    for(i = 0; i < dim*dim; i++)
    {
        *(K2 + i) = *(propagator + i) * dt;
        *(aux_state + i) = *(state + i) + 0.5 * *(K2 + i);
    }
    Dissipator(dissipator, aux_state, bath_state, interaction, dim);
    Propagator(propagator, aux_state, hamiltonian, dissipator, dim);
    for(i = 0; i < dim*dim; i++)
    {
        *(K3 + i) = *(propagator + i) * dt;
        *(aux_state + i) = *(state + i) + *(K3 + i);
    }
    Dissipator(dissipator, aux_state, bath_state, interaction, dim);
    Propagator(propagator, aux_state, hamiltonian, dissipator, dim);
    for(i = 0; i < dim*dim; i++)
    {
        *(K4 + i) = *(propagator + i) * dt;
    }
    for(i = 0; i < dim * dim; i++)
    {
        *(state + i) = *(state + i) + (double)1/6 * (*(K1 + i) + 2* *(K2 + i) + 2* *(K3 + i) + *(K4 + i));
    }
    free(K1);
    free(K2);
    free(K3);
    free(K4);
    free(aux_state);
}
void RK4_closed(double complex *propagator, double complex *state, double complex *hamiltonian, double dt, int dim)
{
    int i;
    double complex *K1, *K2, *K3, *K4,  *aux_state;
    K1 = (double complex*) calloc(dim*dim, sizeof(double complex));
    K2 = (double complex*) calloc(dim*dim, sizeof(double complex));
    K3 = (double complex*) calloc(dim*dim, sizeof(double complex));
    K4 = (double complex*) calloc(dim*dim, sizeof(double complex));
    aux_state = (double complex*) calloc(dim*dim, sizeof(double complex));
    Unitary(propagator, state, hamiltonian, dim);
    for(i = 0; i < dim*dim; i++)
    {
        *(K1 + i) = *(propagator + i) * dt;
        *(aux_state + i) = *(state + i) + 0.5 * *(K1 + i);
    }
    Unitary(propagator, aux_state, hamiltonian, dim);
    for(i = 0; i < dim*dim; i++)
    {
        *(K2 + i) = *(propagator + i) * dt;
        *(aux_state + i) = *(state + i) + 0.5 * *(K2 + i);
    }
    Unitary(propagator, aux_state, hamiltonian, dim);
    for(i = 0; i < dim*dim; i++)
    {
        *(K3 + i) = *(propagator + i) * dt;
        *(aux_state + i) = *(state + i) + *(K3 + i);
    }
    Unitary(propagator, aux_state, hamiltonian, dim);
    for(i = 0; i < dim*dim; i++)
    {
        *(K4 + i) = *(propagator + i) * dt;
    }
    for(i = 0; i < dim * dim; i++)
    {
        *(state + i) = *(state + i) + (double)1/6 * (*(K1 + i) + 2* *(K2 + i) + 2* *(K3 + i) + *(K4 + i));
    }
    free(K1);
    free(K2);
    free(K3);
    free(K4);
    free(aux_state);
}


