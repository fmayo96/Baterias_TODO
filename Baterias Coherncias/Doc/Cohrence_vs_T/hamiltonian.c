#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "matrixoperations.h"
#include "pauli.h"
void Decimal_to_binary(int * binario, int x, int N)
{
    int *binario_inv, i = 0, z;
    binario_inv = (int *) calloc(2*N, sizeof(int));
    z = x;
    while(z > 0)
    {   
        *(binario_inv + i) = z%2;
        z = z/2;
        i++;
    }
    for(i = 0; i < 2*N; i++)
    {
        *(binario + i) = *(binario_inv + 2*N - 1 - i);
    }
    free(binario_inv);
}
int Correct_term(int x, int N)
{
    int *binario, i, sum = 0;
    binario = (int *) calloc(2*N, sizeof(int));
    Decimal_to_binary(binario, x, N);
    for(i = 0; i < 2*N; i++)
    {
        if(*(binario + i) == 0)
        {
            *(binario + i) = -1;
        }
    }
    for(i = 0; i < 2*N; i++)
    {
        sum = sum + *(binario + i);
    }
    free(binario);
    return sum;
}
int Factorial(int N)
{
    int i, factorial = 1;
    for(i = 0; i < N; i++)
    {
        factorial = factorial * (N - i);
    }
    return factorial;
}
int Combinatorio(int N, int p)
{
    int N_factorial, p_factorial, resta_factorial, combinatorio;
    N_factorial = Factorial(N);
    p_factorial = Factorial(p);
    resta_factorial = Factorial(N - p);
    combinatorio = N_factorial /(resta_factorial*p_factorial);
    return combinatorio;
}
void Hamiltonian(double complex *hamiltonian, double gap, int N)
{
    int i, dim = pow(2,N), cant_spin_down, cant_term, component = 0;
    for(i = 0; i < dim*dim; i++)
    {
        *(hamiltonian + i) = 0;
    }
    for(cant_spin_down = 0; cant_spin_down < N + 1; cant_spin_down++)
    {
        cant_term = Combinatorio(N, cant_spin_down);
        for(i = 0; i < cant_term; i++)
        {
            *(hamiltonian + (component + i)*dim + (component + i)) = (N - 2*cant_spin_down)*(double)gap/2;
        }
        component = component + cant_term;
    }
    
}
void Interaction(double complex *interaction, double eps, int N)
{
    int i, dim = pow(2, N), correct_term;
    for(i = 0; i < pow(dim, 4); i++)
    {
        *(interaction + i) = 0;
    }
    for(i = 0; i < dim; i++)
    {
        correct_term = Correct_term((dim + 1)*i, N);
        if(correct_term != 0)
        {
            *(interaction + dim*dim*(dim + 1)*i + (dim*dim - 1 - (dim+1)*i)) = eps;
        }
    }
}
void Thermal_state(double complex *state, double complex *hamiltonian, double beta, int dim)
{
    int i;
    double Z = 0;
    for(i = 0; i < dim*dim; i++)
    {
        *(state + i) = 0;
    }
    for(i = 0; i < dim; i++)
    {
        *(state + dim*i + i) = exp(-(double)beta* *(hamiltonian + dim*i + i));
    }
    for(i = 0; i < dim; i++)
    {
        Z += *(state + dim*i + i);
    }
    for(i = 0; i < dim; i++)
    {
        *(state + dim*i + i) = (double)*(state + dim*i + i)/Z; 
    }
}
void Active_state(double complex *state, double complex *hamiltonian, double beta, int dim)
{
    int i;
    double Z = 0;
    for(i = 0; i < dim*dim; i++)
    {
        *(state + i) = 0;
    }
    for(i = 0; i < dim; i++)
    {
        *(state + dim*i + i) = exp((double)beta* *(hamiltonian + dim*i + i));
    }
    for(i = 0; i < dim; i++)
    {
        Z += *(state + dim*i + i);
    }
    for(i = 0; i < dim; i++)
    {
        *(state + dim*i + i) = (double)*(state + dim*i + i)/Z; 
    }
}