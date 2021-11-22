#ifndef DEPHASING_H
#define DEPHASING_H
void Change_basis_x(double complex *C);
void Dephasing_noise(double complex *state, double complex *C, double p, int dim);
#endif