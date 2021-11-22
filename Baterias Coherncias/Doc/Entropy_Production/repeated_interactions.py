import numpy as np 
from scipy.linalg import expm, logm 
from pauli import *

pauli = Pauli()

def Relative_Entropy(rho, sigma):
        log_rho = logm(rho)
        log_sigma = logm(sigma)
        return np.trace(np.dot(rho, log_rho)) - np.trace(np.dot(rho, log_sigma))

class Reservoir():
    def __init__(self, hamiltonian, beta):
        self.hamiltonian = hamiltonian
        self.beta = beta
    def Thermal_state(self):
        return np.diag(np.diag(np.exp(-self.beta * self.hamiltonian)))/np.sum(np.diag(np.diag(np.exp(-self.beta * self.hamiltonian)))) 
    def Active_state(self): 
        return np.diag(np.diag(np.exp(self.beta * self.hamiltonian)))/np.sum(np.diag(np.diag(np.exp(self.beta * self.hamiltonian)))) 


class Qubit():
    def __init__(self, state, hamiltonian):
        self.state = state
        self.hamiltonian = hamiltonian 
        self.energy = []
        self.heat = []
        self.work = []
        self.ent_prod = []
    def Open_evolution(self, bath_state, bath_hamiltonian, V, dt, tf):
        N = int(tf/dt)
        rho_tot = np.kron(self.state, bath_state) 
        def L(x):
            Comm_v_rho = np.dot(V, x) - np.dot(x, V)
            Comm_v_v_rho = np.dot(V, Comm_v_rho) - np.dot(Comm_v_rho, V)
            
            return -0.5*Comm_v_v_rho
        def K_1(x):
            return dt * L(x)
        def K_2(x):
            return dt*L(x + K_1(x)*0.5)
        def K_3(x): 
            return dt*L(x + K_2(x)*0.5)
        def K_4(x):
            return dt*L(x + K_3(x))    
        for i in range(N):
            rho_tot = rho_tot + 1/6.0 * (K_1(rho_tot) + 2.0 * K_2(rho_tot) + 2.0 * K_3(rho_tot) + K_4(rho_tot)) 
            self.state = np.trace(rho_tot.reshape([2,2,2,2]), axis1= 1, axis2= 3)
            self.energy.append(np.trace(np.dot(self.state, self.hamiltonian[i])))
            self.ent_prod.append(Relative_Entropy(rho_tot, np.kron(self.state, bath_state)))
            rho_tot = np.kron(self.state, bath_state)
            