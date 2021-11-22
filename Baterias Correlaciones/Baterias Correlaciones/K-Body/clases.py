import numpy as np 

class Reservoir():
    def __init__(self, hamiltonian, beta):
        self.hamiltonian = hamiltonian
        self.beta = beta
    def thermal_state(self):        
        return np.array(np.diag(np.diag(np.exp(-self.beta * self.hamiltonian))) / np.sum(np.diag(np.diag(np.exp(-self.beta * self.hamiltonian)))), dtype = np.complex)
    def active_state(self):        
        return np.array(np.diag(np.diag(np.exp(self.beta * self.hamiltonian))) / np.sum(np.diag(np.diag(np.exp(self.beta * self.hamiltonian)))), dtype = np.complex)

class System():
    def __init__(self, hamiltonian, state):
        self.hamiltonian = hamiltonian
        self.state = state
        self.energy = [] 
    def dissipator(self, state, bath, V):
            dim = np.shape(self.hamiltonian)[0]
            rho = np.kron(state, bath)
            Conmutator = np.dot(V,np.dot(V,rho))-np.dot(V,np.dot(rho,V))-np.dot(V,np.dot(rho,V))+np.dot(rho,np.dot(V,V))
            return -0.5*np.trace(Conmutator.reshape([dim, dim, dim, dim]), axis1 = 1, axis2 = 3)
    def time_evol(self, bath, V, tf, dt):
        self.energy.append(np.trace(np.dot(self.state, self.hamiltonian)))
        N = int(tf/dt)
        for i in range(N-1):
            K1 = dt*self.dissipator(self.state, bath, V)
            K2 = dt*self.dissipator(self.state + 0.5*K1, bath, V)
            K3 = dt*self.dissipator(self.state + 0.5*K2, bath, V)
            K4 = dt*self.dissipator(self.state + K3, bath, V)
            self.state = self.state + (1/6) * (K1 + 2*K2 + 2*K3 + K4)
            self.energy.append(np.trace(np.dot(self.state, self.hamiltonian)))
        for i in range(1,N):
            self.energy[i] -= self.energy[0]
        self.energy[0] = 0

def K1(f, H, x, y, dt, V):
    return dt * f(H, x, y, V)
def K2(f, H, x, y, dt, V):
    return dt * (f(H, x + 0.5 * K1(f, H, x, y, dt, V), y, V))
def K3(f, H, x, y, dt, V):
    return dt * (f(H, x + 0.5 * K2(f, H, x, y, dt, V), y, V))
def K4(f, H, x, y, dt, V):
    return dt * (f(H, x + K3(f, H, x, y, dt, V), y, V))

def RK4(f, H, x, y, dt, V):
    return x + (1.0/6) * (K1(f, H, x, y, dt, V)+2*K2(f, H, x, y, dt, V)+2*K3(f, H, x, y, dt, V)+K4(f, H, x, y, dt, V))

