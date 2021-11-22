#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
void Pauli_z(double complex *matrix, double complex gap)
{
    int i;
    for(i = 0; i < 4; i++)
    {
        *(matrix + i) = 0;
    }
    *(matrix) = gap;
    *(matrix + 3) = -gap;
}
void Pauli_x(double complex *matrix, double complex gap)
{
    int i;
    for(i = 0; i < 4; i++)
    {
        *(matrix + i) = 0;
    }
    *(matrix + 1) = gap;
    *(matrix + 2) = gap;
}
void Pauli_y(double complex *matrix, double complex gap)
{
    int i;
    for(i = 0; i < 4; i++)
    {
        *(matrix + i) = 0;
    }
    *(matrix + 1) = -I*gap;
    *(matrix + 2) = I*gap;
}
void Pauli_pl(double complex *matrix, double complex gap)
{
    int i;
    for(i = 0; i < 4; i++)
    {
        *(matrix + i) = 0;
    }
    *(matrix + 1) = gap;
}
void Pauli_mn(double complex *matrix, double complex gap)
{
    int i;
    for(i = 0; i < 4; i++)
    {
        *(matrix + i) = 0;
    }
    *(matrix + 2) = gap;
}
void Eye(double complex *matrix, double complex gap)
{
    int i;
    for(i = 0; i < 2; i++)
    {
        *(matrix + i*2 + i) = 1;
    }
}
