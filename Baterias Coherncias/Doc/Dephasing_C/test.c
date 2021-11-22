#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "matrixoperations.h"
#include "propagators.h"
#include "RK4.h"
#include "time_evolution.h"
#include "energy.h"
#include "pauli.h"
#include "hamiltonian.h"
#include "dephasing.h"
int main()
{
//-----------Declaración de variblaes--------------------------------
double t_f, td = 2.0, dt = 0.001, h_hot = 6, h_cold = 1, beta_hot = 0.2, beta_cold = 10, eps = sqrt(0.5);
int dim, N_qubits = 1, N_tf = 200, i, *N, N_step, N_d = (int)(td/dt);
dim = pow(2, N_qubits);
double *tf, tf_min = 2.1, tf_max = 18;
double complex *hamiltonian_x, *cold_bath_hamiltonian, *hot_bath_hamiltonian, *interaction;
double complex *state, *state_deph, *cold_bath_state, *hot_bath_state, *active_hot, *active_cold;
double *eff, *eff_deph, W1, W2, W4, W4_deph, W5, W5_deph, Q_hot, Q_hot_deph;
tf = (double *) calloc(N_tf, sizeof(double));
N = (int *) calloc(N_tf, sizeof(int));
eff = (double *) calloc(N_tf, sizeof(double));
eff_deph = (double *) calloc(N_tf, sizeof(double));
state = (double complex *) calloc(dim*dim, sizeof(double complex));
state_deph = (double complex *) calloc(dim*dim, sizeof(double complex));
cold_bath_state = (double complex *) calloc(dim*dim, sizeof(double complex));
hot_bath_state = (double complex *) calloc(dim*dim, sizeof(double complex));
active_cold = (double complex *) calloc(dim*dim, sizeof(double complex));
active_hot = (double complex *) calloc(dim*dim, sizeof(double complex));
hamiltonian_x = (double complex *) calloc(dim*dim, sizeof(double complex));
cold_bath_hamiltonian = (double complex *) calloc(dim*dim, sizeof(double complex));
hot_bath_hamiltonian = (double complex *) calloc(dim*dim, sizeof(double complex));
interaction = (double complex *) calloc(dim*dim*dim*dim, sizeof(double complex));
//-----------Tiempos----------------------------------------------------
for(i = 0; i < N_tf; i++)
{
    *(tf + i) = tf_min + i*(double)(tf_max - tf_min) / (N_tf - 1);
    *(N + i) = (int)(*(tf + i)/ dt);
}
double *E, *Q, *W, *t, *Pow, *E_deph, *Q_deph, *W_deph, *t_deph, *Pow_deph; 
E = (double*) calloc(*(N + N_tf - 1), sizeof(double));
Q = (double*) calloc(*(N + N_tf - 1), sizeof(double));
W = (double*) calloc(*(N + N_tf - 1), sizeof(double));
t = (double*) calloc(*(N + N_tf - 1), sizeof(double));
Pow = (double*) calloc(*(N + N_tf - 1), sizeof(double));
E_deph = (double*) calloc(*(N + N_tf - 1), sizeof(double));
Q_deph = (double*) calloc(*(N + N_tf - 1), sizeof(double));
W_deph = (double*) calloc(*(N + N_tf - 1), sizeof(double));
t_deph = (double*) calloc(*(N + N_tf - 1), sizeof(double));
Pow_deph = (double*) calloc(*(N + N_tf - 1), sizeof(double));
//----------Hamiltonianos y Condiciones iniciales-----------------------
Pauli_z(hot_bath_hamiltonian, (double)h_hot/2);
Pauli_z(cold_bath_hamiltonian, (double)h_cold/2);
Pauli_x(hamiltonian_x, (double)h_cold/2);
Interaction(interaction, eps, N_qubits);
Thermal_state(hot_bath_state, hot_bath_hamiltonian, beta_hot, dim);
Thermal_state(cold_bath_state, cold_bath_hamiltonian, beta_cold, dim);
Active_state(active_hot, hot_bath_hamiltonian, beta_hot, dim);
Active_state(active_cold, cold_bath_hamiltonian, beta_cold, dim);

W1 = calc_Energy(active_hot, hot_bath_hamiltonian, dim) - calc_Energy(hot_bath_state, hot_bath_hamiltonian, dim);
W2 = calc_Energy(hot_bath_state, hot_bath_hamiltonian, dim) - calc_Energy(hot_bath_state, cold_bath_hamiltonian, dim);
for(i = 0; i < N_tf; i++)
{
    t_f = *(tf + i);
    N_step = *(N + i) - 1;
    Thermal_state(state, hot_bath_hamiltonian, beta_hot, dim);
    Thermal_state(state_deph, hot_bath_hamiltonian, beta_hot, dim);
    //------------Evolucion--------------------------------------------------------------------
    Driven_evolution_dephasing(state_deph, hamiltonian_x, cold_bath_state, cold_bath_hamiltonian, interaction, dim, td, dt, E_deph, Q_deph, W_deph, 1);    
    Open_evolution(state_deph, cold_bath_hamiltonian, cold_bath_state, cold_bath_hamiltonian, interaction, dim, t_f, dt, E_deph, Q_deph, W_deph, N_d, hot_bath_state);
    Driven_evolution(state, hamiltonian_x, cold_bath_state, cold_bath_hamiltonian, interaction, dim, td, dt, E, Q, W);
    Open_evolution(state, cold_bath_hamiltonian, cold_bath_state, cold_bath_hamiltonian, interaction, dim, t_f, dt, E, Q, W, N_d, hot_bath_state);
    W4 = calc_Energy(state, cold_bath_hamiltonian, dim) - calc_Energy(state, hot_bath_hamiltonian, dim);
    W4_deph = calc_Energy(state_deph, cold_bath_hamiltonian, dim) - calc_Energy(state_deph, hot_bath_hamiltonian, dim);
    W5 = 2.0*(calc_Energy(state, hot_bath_hamiltonian, dim) - calc_Energy(active_hot, hot_bath_hamiltonian, dim));
    W5_deph = 2.0*(calc_Energy(state_deph, hot_bath_hamiltonian, dim) - calc_Energy(active_hot, hot_bath_hamiltonian, dim));
    Q_hot = fabs((double)W5/2);
    Q_hot_deph = fabs((double)W5_deph/2);
    *(eff + i) = (double)((W1 + W2 + W4 + W5 + *(W + N_step)) / Q_hot);
    *(eff_deph + i) = (double)(W1 + W2 + W4_deph + W5_deph + *(W_deph + N_step)) / Q_hot_deph;
    *(Pow + i) = (double)((W1 + W2 + W4 + W5 + *(W + N_step)) / *(tf + i));
    *(Pow_deph + i) = (double)(W1 + W2 + W4_deph + W5_deph + *(W_deph + N_step)) / *(tf + i);
}
//-------------Guardar Datos---------------------------------------------------------------
char filename[255];
sprintf(filename,"eff_eff_deph_vs_tf.txt");

char filename_2[255];
sprintf(filename_2,"Test_constant.txt");


FILE *fp=fopen(filename,"w");

for(i = 0; i < N_tf; i++)
{
    fprintf(fp, "%lf %lf %lf %lf %lf \n", *(tf + i),*(eff + i), *(eff_deph + i), *(Pow + i), *(Pow_deph + i));
}
FILE *fp2=fopen(filename_2,"w");

for(i = 0; i < *(N + N_tf - 1); i++)
{
    fprintf(fp2, "%lf %lf %lf %lf\n", *(t + i),*(E + i), *(Q + i), *(W + i));
}
    fclose(fp);
    fclose(fp2);
    free(state);
    free(state_deph);
    free(cold_bath_state);
    free(hot_bath_state);
    free(active_hot);
    free(active_cold);
    free(cold_bath_hamiltonian);
    free(hot_bath_hamiltonian);
    free(hamiltonian_x);
    free(interaction);
    free(tf);
    free(N);
    free(eff);
    free(eff_deph);
    free(E);
    free(Q);
    free(W);
    free(t);
    free(E_deph);
    free(Q_deph);
    free(W_deph);
    free(t_deph);
    return 0;
}


