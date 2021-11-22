#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "matrixoperations.h"

void Change_basis_x(double complex *C)
{
    *C = (double)1/sqrt(2);
    *(C + 1) = (double)1/sqrt(2);
    *(C + 2) = (double)1/sqrt(2);
    *(C + 3) = -(double)1/sqrt(2);
}
void Dephasing_noise(double complex *state, double complex *C, double p, int dim)
{
    double complex *state_C, *C_state_C;
    state_C = (double complex *) calloc(dim*dim, sizeof(double complex));
    C_state_C = (double complex *) calloc(dim*dim, sizeof(double complex));
    dot(state_C, state, C, dim);
    dot(C_state_C, C, state_C, dim);
    *(C_state_C + 1) = (1-(double)p/2)* *(C_state_C + 1);
    *(C_state_C + 2) = (1-(double)p/2)* *(C_state_C + 2);
    dot(state_C, C, C_state_C, dim);
    dot(state, state_C, C, dim);
    free(state_C);
    free(C_state_C);
}