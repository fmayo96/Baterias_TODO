#ifndef PAULI_H
#define PAULI_H
void Pauli_z(double complex *matrix, double complex gap);
void Pauli_x(double complex *matrix, double complex gap);
void Pauli_y(double complex *matrix, double complex gap);
void Pauli_pl(double complex *matrix, double complex gap);
void Pauli_mn(double complex *matrix, double complex gap);
void Eye(double complex *matrix, double complex gap);
#endif