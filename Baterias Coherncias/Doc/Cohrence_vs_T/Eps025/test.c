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
double tf = 10, td = 2, dt = 0.0001, *h_hot, h_cold = 1, *beta_hot, beta_cold = 10.0, eps = sqrt(0.5);
double h_start = 2, h_finish = 20;
int dim, N_qubits = 1, i, j, N = (int)(tf/dt), N_d = (int)(td/dt), N_betas = 4, N_h_hot = 40;
dim = pow(2, N_qubits);
double complex *hamiltonian_x, *cold_bath_hamiltonian, *hot_bath_hamiltonian, *interaction;
double complex *state, *state_deph, *cold_bath_state, *hot_bath_state, *active_hot, *active_cold;
double *eff, *eff_deph, W1, W2, W4, W4_deph, W5, W5_deph, Q_hot, Q_hot_deph;
h_hot = (double *) calloc(N_h_hot, sizeof(double));
beta_hot = (double *) calloc(N_betas, sizeof(double));
eff = (double *) calloc(N_h_hot*N_betas, sizeof(double));
eff_deph = (double *) calloc(N_h_hot*N_betas, sizeof(double));
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
*(beta_hot) = 0.01;
*(beta_hot + 1) = 0.05;
*(beta_hot + 2) = 0.1;
*(beta_hot + 3) = 0.2;
for(i = 0; i < N_h_hot; i++)
{
    *(h_hot + i) = h_start + (h_finish - h_start)*(double)i/(N_h_hot - 1);
}
double *E, *Q, *W, *t, *Pow, *E_deph, *Q_deph, *W_deph, *t_deph, *Pow_deph; 
E = (double*) calloc(N, sizeof(double));
Q = (double*) calloc(N, sizeof(double));
W = (double*) calloc(N, sizeof(double));
t = (double*) calloc(N, sizeof(double));
Pow = (double*) calloc(N_h_hot*N_betas, sizeof(double));
E_deph = (double*) calloc(N, sizeof(double));
Q_deph = (double*) calloc(N, sizeof(double));
W_deph = (double*) calloc(N, sizeof(double));
t_deph = (double*) calloc(N, sizeof(double));
Pow_deph = (double*) calloc(N_h_hot*N_betas, sizeof(double));
//----------Hamiltonianos y Condiciones iniciales-----------------------
Pauli_z(cold_bath_hamiltonian, (double)h_cold/2.0);
Pauli_x(hamiltonian_x, (double)h_cold/2.0);
Interaction(interaction, eps, N_qubits);
for(i = 0; i < N_betas; i++)
{
    for(j = 0; j < N_h_hot; j++)
    {
        Pauli_z(hot_bath_hamiltonian, (double)*(h_hot + j)/2.0);
        Thermal_state(hot_bath_state, hot_bath_hamiltonian, *(beta_hot + i), dim);
        Thermal_state(cold_bath_state, cold_bath_hamiltonian, beta_cold, dim);
        Active_state(active_hot, hot_bath_hamiltonian, *(beta_hot + i), dim);
        Active_state(active_cold, cold_bath_hamiltonian, beta_cold, dim);

        W1 = calc_Energy(active_hot, hot_bath_hamiltonian, dim) - calc_Energy(hot_bath_state, hot_bath_hamiltonian, dim);
        W2 = calc_Energy(hot_bath_state, hot_bath_hamiltonian, dim) - calc_Energy(hot_bath_state, cold_bath_hamiltonian, dim);
        Thermal_state(state, hot_bath_hamiltonian, *(beta_hot + i), dim);
        Thermal_state(state_deph, hot_bath_hamiltonian, *(beta_hot + i), dim);
        //------------Evolucion--------------------------------------------------------------------
        Driven_evolution_dephasing(state_deph, hamiltonian_x, cold_bath_state, cold_bath_hamiltonian, interaction, dim, td, dt, E_deph, Q_deph, W_deph, 1);    
        Open_evolution(state_deph, cold_bath_hamiltonian, cold_bath_state, cold_bath_hamiltonian, interaction, dim, tf, dt, E_deph, Q_deph, W_deph, N_d, hot_bath_state);
        Driven_evolution(state, hamiltonian_x, cold_bath_state, cold_bath_hamiltonian, interaction, dim, td, dt, E, Q, W);
        Open_evolution(state, cold_bath_hamiltonian, cold_bath_state, cold_bath_hamiltonian, interaction, dim, tf, dt, E, Q, W, N_d, hot_bath_state);
        W4 = calc_Energy(state, cold_bath_hamiltonian, dim) - calc_Energy(state, hot_bath_hamiltonian, dim);
        W4_deph = calc_Energy(state_deph, cold_bath_hamiltonian, dim) - calc_Energy(state_deph, hot_bath_hamiltonian, dim);
        W5 = 2.0*(calc_Energy(state, hot_bath_hamiltonian, dim) - calc_Energy(active_hot, hot_bath_hamiltonian, dim));
        W5_deph = 2.0*(calc_Energy(state_deph, hot_bath_hamiltonian, dim) - calc_Energy(active_hot, hot_bath_hamiltonian, dim));
        Q_hot = ((double)W5/2.0);
        Q_hot_deph = ((double)W5_deph/2.0);
        *(eff + N_h_hot*i + j) = (double)(W1 + W2 + W4 + W5 + *(W + N - 1)) / Q_hot;
        *(eff_deph + N_h_hot*i + j) = (double)(W1 + W2 + W4_deph + W5_deph + *(W_deph + N - 1)) / Q_hot_deph;
        *(Pow + N_h_hot*i + j) = (double)(W1 + W2 + W4 + W5 + *(W + N - 1)) / (2.0 * tf);
        *(Pow_deph + N_h_hot*i + j) = (double)(W1 + W2 + W4_deph + W5_deph + *(W_deph + N - 1)) / (2.0 * tf);
    }
}
//-------------Guardar Datos---------------------------------------------------------------
char filename[255];
sprintf(filename,"beta=0.01.txt");
char filename_2[255];
sprintf(filename_2,"beta=0.05.txt");
char filename_3[255];
sprintf(filename_3,"beta=0.1.txt");
char filename_4[255];
sprintf(filename_4,"beta=0.2.txt");

FILE *fp=fopen(filename,"w");

    for(j = 0; j < N_h_hot; j++)
    {
        fprintf(fp, "%lf %lf %lf %lf %lf \n", *(h_hot + j),*(eff + j), *(eff_deph + j), *(Pow + j), *(Pow_deph + j));
    }
FILE *fp2=fopen(filename_2,"w");

    for(j = 0; j <N_h_hot; j++)
    {
        fprintf(fp2, "%lf %lf %lf %lf %lf \n", *(h_hot + j),*(eff + j + N_h_hot), *(eff_deph + j + N_h_hot), *(Pow + j + N_h_hot), *(Pow_deph + j + N_h_hot));
    }
FILE *fp3=fopen(filename_3,"w");

    for(j = 0; j < N_h_hot; j++)
    {
        fprintf(fp3, "%lf %lf %lf %lf %lf \n", *(h_hot + j),*(eff + j + 2*N_h_hot), *(eff_deph + j + 2*N_h_hot), *(Pow + j + 2*N_h_hot), *(Pow_deph + j + 2*N_h_hot));
    }
FILE *fp4=fopen(filename_4,"w");

    for(j = 0; j < N_h_hot; j++)
    {
        fprintf(fp4, "%lf %lf %lf %lf %lf \n", *(h_hot + j),*(eff + j + 3*N_h_hot), *(eff_deph + j + 3*N_h_hot), *(Pow + j + 3*N_h_hot), *(Pow_deph + j + 3*N_h_hot));
    }
fclose(fp);
fclose(fp2);
fclose(fp3);
fclose(fp4);
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
free(h_hot);
free(beta_hot);
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


