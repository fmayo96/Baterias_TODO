import numpy as np
import matplotlib.pyplot as plt 
from scipy.linalg import expm 
from repeated_interactions import *
from pauli import *

#----------Create Pauli matrix instance--------
pauli = Pauli()
#----------------------------------------------
#----------Define Hamiltonians-----------------
h = 1.5 

beta = 1
tf = 5
td = 2
dt = 0.001
N = int(tf/dt)
Nd = int(td/dt)
H = np.zeros([N, 2, 2])
for i in range(Nd):
    H[i] = h/2.0 * pauli.z
for i in range(Nd, N):
    H[i] = h/2.0 * pauli.z
H_bath = h/2.0 * pauli.z
eps = np.sqrt(10)
V = eps*(np.kron(pauli.pl, pauli.pl) + np.kron(pauli.mn, pauli.mn))
#-------------------------------------------------
#----------Initialize systems---------------------
bath = Reservoir(H_bath, beta)
qubit = Qubit(bath.Thermal_state(), H)
qubit.Open_evolution(bath.Thermal_state(), bath.hamiltonian, V, dt, tf)

total_ent_prod = np.zeros(N)
for i in range(N-1):
    total_ent_prod[i+1] = total_ent_prod[i] + dt * qubit.ent_prod[i]
#-----------Plot----------------------------------
t = np.linspace(0, tf, N)
plt.figure()
plt.plot(t, total_ent_prod, linewidth = 2)
plt.show()