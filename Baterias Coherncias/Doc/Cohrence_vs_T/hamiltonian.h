#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
void Hamiltonian(double complex *hamiltonian, double gap, int N);
void Interaction(double complex *interaction, double eps, int N);
void Decimal_to_binary(int * binario, int x, int N);
int Correct_term(int x, int N);
int Factorial(int N);
void Thermal_state(double complex *state, double complex *hamiltonian, double beta, int dim);
void Active_state(double complex *state, double complex *hamiltonian, double beta, int dim);
#endif 