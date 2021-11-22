import numpy as np 
import matplotlib.pyplot as plt 
from clases import *

#=======Pauli Matrices==============#
sz = np.array([[1,0],[0, -1]])
sx = np.array([[0,1],[1, 0]])
sy = np.array([[0,-1j],[1j, 0]])
sd = np.array([[0,0],[1, 0]])
su = np.array([[0,1],[0, 0]])
eye = np.array([[1,0],[0, 1]])
#===================================#
tf = 20
dt = 1e-3
N = int(tf/dt)
beta = 20
h = 1.5
eps = np.sqrt(0.5)
H = h/2*sz
H2 = np.kron(H, eye) + np.kron(eye, H)
H3 = np.kron(H, np.kron(eye, eye)) + np.kron(eye, np.kron(H, eye)) + np.kron(eye, np.kron(eye, H))
H4 = np.kron(H, np.kron(eye, np.kron(eye,eye))) + np.kron(eye, np.kron(H, np.kron(eye,eye))) + np.kron(eye, np.kron(eye, np.kron(H,eye))) + np.kron(eye, np.kron(eye, np.kron(eye,H)))
V_21 = eps*(np.kron(su, np.kron(eye, np.kron(su, eye))) + np.kron(sd, np.kron(eye, np.kron(sd, eye))) + np.kron(eye, np.kron(su, np.kron(eye, su))) + np.kron(eye, np.kron(sd, np.kron(eye, sd))))
V_22 = 2*eps*(np.kron(su, np.kron(su, np.kron(su, su))) + np.kron(sd, np.kron(sd, np.kron(sd, sd))))
#V_31 = eps*(np.kron(su, np.kron(eye, np.kron(eye, np.kron(su, np.kron(eye, eye))))) + np.kron(sd, np.kron(eye, np.kron(eye, np.kron(sd, np.kron(eye, eye))))) + np.kron(eye, np.kron(su, np.kron(eye, np.kron(eye, np.kron(su, eye))))) + np.kron(eye, np.kron(sd, np.kron(eye, np.kron(eye, np.kron(sd, eye))))) + np.kron(eye, np.kron(eye, np.kron(su, np.kron(eye, np.kron(eye, su))))) + np.kron(eye, np.kron(eye, np.kron(sd, np.kron(eye, np.kron(eye, sd))))))
#V_32 = 0.65138*eps*(np.kron(su, np.kron(su, np.kron(eye, np.kron(su, np.kron(su, eye))))) + np.kron(sd, np.kron(sd, np.kron(eye, np.kron(sd, np.kron(sd, eye))))) + np.kron(su, np.kron(eye, np.kron(su, np.kron(su, np.kron(eye, su))))) + np.kron(sd, np.kron(eye, np.kron(sd, np.kron(sd, np.kron(eye, sd))))) + np.kron(eye, np.kron(su, np.kron(su, np.kron(eye, np.kron(su, su))))) + np.kron(eye, np.kron(sd, np.kron(sd, np.kron(eye, np.kron(sd, sd))))) + np.kron(su, np.kron(eye, np.kron(eye, np.kron(su, np.kron(eye,eye))))) + np.kron(sd, np.kron(eye, np.kron(eye, np.kron(sd, np.kron(eye,eye))))) + np.kron(eye, np.kron(su, np.kron(eye, np.kron(eye, np.kron(su,eye))))) + np.kron(eye, np.kron(sd, np.kron(eye, np.kron(eye, np.kron(sd,eye))))) + np.kron(eye, np.kron(eye, np.kron(su, np.kron(eye, np.kron(eye,su))))) + np.kron(eye, np.kron(eye, np.kron(sd, np.kron(eye, np.kron(eye,sd))))))
V_33 = 3*eps*(np.kron(su,np.kron(su,np.kron(su,np.kron(su,np.kron(su,su))))) + np.kron(sd,np.kron(sd,np.kron(sd,np.kron(sd,np.kron(sd,sd))))) + np.kron(su,np.kron(sd,np.kron(su,np.kron(su,np.kron(sd,su))))) + np.kron(sd,np.kron(su,np.kron(sd,np.kron(sd,np.kron(su,sd)))))+ np.kron(su,np.kron(su,np.kron(sd,np.kron(su,np.kron(su,sd))))) + np.kron(sd,np.kron(sd,np.kron(su,np.kron(sd,np.kron(sd,su)))))+ np.kron(sd,np.kron(su,np.kron(su,np.kron(sd,np.kron(su,su))))) + np.kron(su,np.kron(sd,np.kron(sd,np.kron(su,np.kron(sd,sd))))))
V_44 = 4*eps*(np.kron(su,np.kron(su,np.kron(su,np.kron(su,np.kron(su,np.kron(su,np.kron(su,su))))))) + np.kron(sd,np.kron(sd,np.kron(sd,np.kron(sd,np.kron(sd,np.kron(sd,np.kron(sd,sd))))))) + np.kron(su,np.kron(su,np.kron(su,np.kron(sd,np.kron(su,np.kron(su,np.kron(su,sd))))))) + np.kron(sd,np.kron(sd,np.kron(sd,np.kron(su,np.kron(sd,np.kron(sd,np.kron(sd,su)))))))+ np.kron(su,np.kron(su,np.kron(sd,np.kron(su,np.kron(su,np.kron(su,np.kron(sd,su))))))) + np.kron(sd,np.kron(sd,np.kron(su,np.kron(sd,np.kron(sd,np.kron(sd,np.kron(su,sd))))))) + np.kron(su,np.kron(sd,np.kron(su,np.kron(su,np.kron(su,np.kron(sd,np.kron(su,su))))))) + np.kron(sd,np.kron(su,np.kron(sd,np.kron(sd,np.kron(sd,np.kron(su,np.kron(sd,sd)))))))+ np.kron(su,np.kron(sd,np.kron(sd,np.kron(sd,np.kron(su,np.kron(sd,np.kron(sd,sd))))))) + np.kron(sd,np.kron(su,np.kron(su,np.kron(su,np.kron(sd,np.kron(su,np.kron(su,su))))))))
V_42 = eps*(np.kron(su, np.kron(su, np.kron(eye, np.kron(eye, np.kron(su, np.kron(su, np.kron(eye, eye) )))))) + np.kron(su, np.kron(su, np.kron(eye, np.kr8 on(eye, np.kron(sd, np.kron(sd, np.kron(eye, eye))))))))
eival_21, eivec_21 = np.linalg.eig(V_21)
print(max(abs(eival_21)))
eival_22, eivec_22 = np.linalg.eig(V_22)
print(max(abs(eival_22)))
eival_33, eivec_33 = np.linalg.eig(V_33)
print(max(abs(eival_33)))
eival_44, eivec_44 = np.linalg.eig(V_44)

bath_2 = Reservoir(H2, beta)
bath_3 = Reservoir(H3, beta)
bath_4 = Reservoir(H4, beta)
system_par = System(H2, bath_2.thermal_state())
system_par.time_evol(bath_2.thermal_state(), V_21, tf, dt)
E_21 = np.array(system_par.energy)
system_col = System(H2, bath_2.thermal_state())
system_col.time_evol(bath_2.thermal_state(), V_22, tf, dt)
E_22 = np.array(system_col.energy)
#system_31 = System(H3, bath_3.thermal_state())
#system_31.time_evol(bath_3.thermal_state(), V_31, tf, dt)
#E_31 = np.array(system_31.energy)
#system_32 = System(H3, bath_3.thermal_state())
#system_32.time_evol(bath_3.thermal_state(), V_32, tf, dt)
#E_32 = np.array(system_32.energy)
system_33 = System(H3, bath_3.thermal_state())
system_33.time_evol(bath_3.thermal_state(), V_33, tf, dt)
E_33 = np.array(system_33.energy)
system_44 = System(H4, bath_4.thermal_state())
system_44.time_evol(bath_4.thermal_state(), V_44, tf, dt)
E_44 = np.array(system_44.energy)


t = np.linspace(0,tf,N)
plt.figure()
plt.plot(t, E_21/E_21[-1], linewidth = 2, label = r'$N = 2, k = 1$')
plt.plot(t, E_22/E_22[-1], linewidth = 2, label = r'$N = 2, k = 2$')
#plt.plot(t, E_31/E_33[-1], '--',linewidth = 2, label = r'$N = 3, k = 1$')
#plt.plot(t, E_32/E_33[-1],linewidth = 2, label = r'$N = 3, k = 2$')
plt.plot(t, E_33/E_33[-1], '--',linewidth = 2, label = r'$N = 3, k = 3$')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.legend()
plt.grid()
plt.show()
