import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import sqrtm

tf=1
dt = 0.0001
N = int(tf/dt)
h = 1.5
beta = 5
bh = h*beta/2.0
eps = np.sqrt(10)
Z = np.exp(bh) + np.exp(-bh)
rho1 = np.zeros([N,2,2],dtype=np.complex)
rho2 = np.zeros([N,4,4],dtype=np.complex)
rho3 = np.zeros((N,8,8),dtype=np.complex)
rho4 = np.zeros((N,16,16),dtype=np.complex)
rho5 = np.zeros((N,32,32),dtype=np.complex)
thermal = [[np.exp(bh)/Z,0],[0,np.exp(-bh)/Z]]
thermal_posta = [[np.exp(-bh)/Z,0],[0,np.exp(bh)/Z]]
rho1[0] = thermal_posta
rho2[0] = np.kron(thermal_posta,thermal_posta)
rho3[0] = np.kron(np.kron(thermal_posta,thermal_posta),thermal_posta)
rho4[0] = np.kron(thermal_posta,np.kron(thermal_posta,np.kron(thermal_posta,thermal_posta)))
rho5[0] = np.kron(thermal_posta,np.kron(thermal_posta,np.kron(thermal_posta,np.kron(thermal_posta,thermal_posta))))
rho1E = thermal
rho2E = np.kron(thermal,thermal)
rho3E = np.kron(np.kron(thermal,thermal),thermal)
rho4E = np.kron(thermal,np.kron(thermal,np.kron(thermal,thermal)))
rho5E = np.kron(thermal,np.kron(thermal,np.kron(thermal,np.kron(thermal,thermal))))
Id = [[1,0],[0,1]]
sz = [[1,0],[0,-1]]
su = [[0,0],[1,0]]
sd = [[0,1],[0,0]]

H1 = (np.eye(2)*(h/2))* sz
H2 = np.kron(H1,Id)+np.kron(Id,H1)
H3 = np.kron(H1,np.kron(Id,Id))+np.kron(Id,np.kron(H1,Id))+np.kron(Id,np.kron(Id,H1))
H4 = np.kron(H1,np.kron(Id,np.kron(Id,Id)))+np.kron(Id,np.kron(H1,np.kron(Id,Id)))+np.kron(Id,np.kron(Id,np.kron(H1,Id)))+np.kron(Id,np.kron(Id,np.kron(Id,H1)))
H5 = np.kron(H1,np.kron(Id,np.kron(Id,np.kron(Id,Id))))+np.kron(Id,np.kron(H1,np.kron(Id,np.kron(Id,Id))))+np.kron(Id,np.kron(Id,np.kron(H1,np.kron(Id,Id))))+np.kron(Id,np.kron(Id,np.kron(Id,np.kron(H1,Id)))) + np.kron(Id,np.kron(Id,np.kron(Id,np.kron(Id,H1))))
V1 = eps*(np.kron(su,sd) + np.kron(sd,su))
V2 = 2*eps*(np.kron(su,np.kron(su,np.kron(sd,sd))) + np.kron(sd,np.kron(sd,np.kron(su,su)))+np.kron(su,np.kron(sd,np.kron(sd,su))) + np.kron(sd,np.kron(su,np.kron(su,sd))))
V3 = 3*eps*(np.kron(su,np.kron(su,np.kron(su,np.kron(sd,np.kron(sd,sd))))) + np.kron(sd,np.kron(sd,np.kron(sd,np.kron(su,np.kron(su,su))))) + np.kron(su,np.kron(sd,np.kron(su,np.kron(sd,np.kron(su,sd))))) + np.kron(sd,np.kron(su,np.kron(sd,np.kron(su,np.kron(sd,su)))))+ np.kron(su,np.kron(su,np.kron(sd,np.kron(sd,np.kron(sd,su))))) + np.kron(sd,np.kron(sd,np.kron(su,np.kron(su,np.kron(su,sd)))))+ np.kron(sd,np.kron(su,np.kron(su,np.kron(su,np.kron(sd,sd))))) + np.kron(su,np.kron(sd,np.kron(sd,np.kron(sd,np.kron(su,su))))))
V4 = 4*eps*(np.kron(su,np.kron(su,np.kron(su,np.kron(su,np.kron(sd,np.kron(sd,np.kron(sd,sd))))))) + np.kron(sd,np.kron(sd,np.kron(sd,np.kron(sd,np.kron(su,np.kron(su,np.kron(su,su))))))) + np.kron(su,np.kron(su,np.kron(su,np.kron(sd,np.kron(sd,np.kron(sd,np.kron(sd,su))))))) + np.kron(sd,np.kron(sd,np.kron(sd,np.kron(su,np.kron(su,np.kron(su,np.kron(su,sd)))))))+ np.kron(su,np.kron(su,np.kron(sd,np.kron(su,np.kron(sd,np.kron(sd,np.kron(su,sd))))))) + np.kron(sd,np.kron(sd,np.kron(su,np.kron(sd,np.kron(su,np.kron(su,np.kron(sd,su))))))) + np.kron(su,np.kron(sd,np.kron(su,np.kron(su,np.kron(sd,np.kron(su,np.kron(sd,sd))))))) + np.kron(sd,np.kron(su,np.kron(sd,np.kron(sd,np.kron(su,np.kron(sd,np.kron(su,su)))))))+ np.kron(su,np.kron(sd,np.kron(sd,np.kron(sd,np.kron(sd,np.kron(su,np.kron(su,su))))))) + np.kron(sd,np.kron(su,np.kron(su,np.kron(su,np.kron(su,np.kron(sd,np.kron(sd,sd))))))))
V5 = np.zeros([1024,1024])
for i in range(1024):
    V5[i,1024-i-1] = 5*eps
print(V2)
eigval1,eigvec1 = LA.eig(V1)
eigval2,eigvec2 = LA.eig(V2)
eigval3,eigvec3 = LA.eig(V3)
eigval4,eigvec4 = LA.eig(V4)
eigval5, eigvec5 = LA.eig(V5)
print(np.max(eigval1))
print(np.max(eigval2)/2)
print(np.max(eigval3)/3)
print(np.max(eigval4)/4)
print(np.max(eigval5)/5)
def D_1(x):
        rho = np.kron(x,rho1E)
        Conmutator = np.dot(V1,np.dot(V1,rho))-np.dot(V1,np.dot(rho,V1))-np.dot(V1,np.dot(rho,V1))+np.dot(rho,np.dot(V1,V1))
        return -0.5*np.trace(np.reshape(Conmutator,[2,2,2,2]), axis1 = 1, axis2 = 3)
def D_2(x):
        rho = np.kron(x,rho2E)
        Conmutator = np.dot(V2,np.dot(V2,rho))-np.dot(V2,np.dot(rho,V2))-np.dot(V2,np.dot(rho,V2))+np.dot(rho,np.dot(V2,V2))
        return -0.5*np.trace(np.reshape(Conmutator,[4,4,4,4]), axis1 = 1, axis2 = 3)

def D_3(x):
        rho = np.kron(x,rho3E)
        Conmutator = np.dot(V3,np.dot(V3,rho))-np.dot(V3,np.dot(rho,V3))-np.dot(V3,np.dot(rho,V3))+np.dot(rho,np.dot(V3,V3))
        return -0.5*np.trace(np.reshape(Conmutator,[8,8,8,8]), axis1 = 1, axis2 = 3)

def D_4(x):
        rho = np.kron(x,rho4E)
        Conmutator = np.dot(V4,np.dot(V4,rho))-np.dot(V4,np.dot(rho,V4))-np.dot(V4,np.dot(rho,V4))+np.dot(rho,np.dot(V4,V4))
        return -0.5*np.trace(np.reshape(Conmutator,[16,16,16,16]), axis1 = 1, axis2 = 3)

def D_5(x):
        rho = np.kron(x,rho5E)
        Conmutator = np.dot(V5,np.dot(V5,rho))-np.dot(V5,np.dot(rho,V5))-np.dot(V5,np.dot(rho,V5))+np.dot(rho,np.dot(V5,V5))
        return -0.5*np.trace(np.reshape(Conmutator,[32,32,32,32]), axis1 = 1, axis2 = 3)



def L_1(x):
    return (-1j*(np.dot(H1,rho1[i]) - np.dot(rho1[i],H1)) + D_1(rho1[i]))
def K_1_1(x):
    return dt * L_1(x)
def K_1_2(x):
    return dt*L_1(x + K_1_1(x)*0.5)
def K_1_3(x): 
    return dt*L_1(x + K_1_2(x)*0.5)
def K_1_4(x):
    return dt*L_1(x + K_1_3(x))    

def L_2(x):
    return (-1j*(np.dot(H2,rho2[i]) - np.dot(rho2[i],H2)) + D_2(rho2[i]))
def K_2_1(x):
    return dt * L_2(x)
def K_2_2(x):
    return dt*L_2(x + K_2_1(x)*0.5)
def K_2_3(x): 
    return dt*L_2(x + K_2_2(x)*0.5)
def K_2_4(x):
    return dt*L_2(x + K_2_3(x))    


def L_3(x):
    return (-1j*(np.dot(H3,rho3[i]) - np.dot(rho3[i],H3)) + D_3(rho3[i]))
def K_3_1(x):
    return dt * L_3(x)
def K_3_2(x):
    return dt*L_3(x + K_3_1(x)*0.5)
def K_3_3(x): 
    return dt*L_3(x + K_3_2(x)*0.5)
def K_3_4(x):
    return dt*L_3(x + K_3_3(x))    


def L_4(x):
    return (-1j*(np.dot(H4,rho4[i]) - np.dot(rho4[i],H4)) + D_4(rho4[i]))
def K_4_1(x):
    return dt * L_4(x)
def K_4_2(x):
    return dt*L_4(x + K_4_1(x)*0.5)
def K_4_3(x): 
    return dt*L_4(x + K_4_2(x)*0.5)
def K_4_4(x):
    return dt*L_4(x + K_4_3(x))    


def L_5(x):
    return (-1j*(np.dot(H5,rho5[i]) - np.dot(rho5[i],H5)) + D_5(rho5[i]))

for i in range(0,N-1):
    rho1[i+1] = rho1[i] + K_1_1(rho1[i])/6.0 + K_1_2(rho1[i])/3.0 + K_1_3(rho1[i])/3.0 + K_1_4(rho1[i])/6.0

for i in range(0,N-1):
    rho2[i+1] = rho2[i] + K_2_1(rho2[i])/6.0 + K_2_2(rho2[i])/3.0 + K_2_3(rho2[i])/3.0 + K_2_4(rho2[i])/6.0

for i in range(0,N-1):
    rho3[i+1] = rho3[i] + K_3_1(rho3[i])/6.0 + K_3_2(rho3[i])/3.0 + K_3_3(rho3[i])/3.0 + K_3_4(rho3[i])/6.0

for i in range(0,N-1):
    rho4[i+1] = rho4[i] + K_4_1(rho4[i])/6.0 + K_4_2(rho4[i])/3.0 + K_4_3(rho4[i])/3.0 + K_4_4(rho4[i])/6.0


"""
for i in range(0,N-1):
    rho5[i+1] = rho5[i] + dt * L_5(rho5[i] + (dt/2) * L_5(rho5[i]))

"""

E1 = np.zeros(N)
E2 = np.zeros(N)
E3 = np.zeros(N)
E4 = np.zeros(N)
E5 = np.zeros(N)
for i in range(0,N):
    E1[i] = np.trace(np.dot(H1,rho1[i]))
    E2[i] = np.trace(np.dot(H2,rho2[i]))
    E3[i] = np.trace(np.dot(H3,rho3[i]))
    E4[i] = np.trace(np.dot(H4,rho4[i]))
    #E5[i] = np.trace(np.dot(H5,rho5[i]))
for i in range(1,N):
    E1[i] -= E1[0]
    E2[i] -= E2[0]
    E3[i] -= E3[0]
    E4[i] -= E4[0]
    #E5[i] -= E5[0]
E1[0] = 0
E2[0] = 0
E3[0] = 0
E4[0] = 0
#E5[0] = 0

t = np.linspace(0,tf,N)

plt.figure()
plt.plot(t,E1,linewidth = 2)
plt.plot(t,E2/2.0,linewidth = 2)
plt.plot(t,E3/3.0,linewidth = 2)
plt.plot(t,E4/4.0,linewidth = 2)
plt.legend(["1 qubit","2 qubits", "3 qubits","4 qubits"])
plt.ylabel("Energy per qubit")
plt.xlabel("Time")
plt.grid()
plt.show()


tcarga1 = 0
tcarga2 = 0
tcarga3 = 0
tcarga4 = 0


for i in range(0,N):
    if (E1[i]<=E1[N-1]*0.999):
        tcarga1 = i
    else:
         break
print(tcarga1)

for i in range(0,N):
    if (E2[i]<=E2[N-1]*0.999):
        tcarga2 = i
    else:
         break
print(tcarga2)
for i in range(0,N):
    if (E3[i]<=E3[N-1]*0.999):
        tcarga3 = i
    else:
         break
print(tcarga3)
for i in range(0,N):
    if (E4[i]<=E4[N-1]*0.999):
        tcarga4 = i
    else:
         break


area_1 = 0
area_2 = 0
area_3 = 0
area_4 = 0

for i in range(tcarga1):
    area_1 += E1[i]/tcarga1

for i in range(tcarga2):
    area_2 += E2[i]/tcarga2

for i in range(tcarga3):
    area_3 += E3[i]/tcarga3

for i in range(tcarga4):
    area_4 += E4[i]/tcarga4

print('Area_1 = {}'.format(area_1))
print('Area_2 = {}'.format(area_2))
print('Area_3 = {}'.format(area_3))
print('Area_4 = {}'.format(area_4))

print(tcarga4)
print("Energias:")
print(E1[-1])
print(E2[-1]/2.0)
print(E3[-1]/3.0)
print(E4[-1]/4.0)

P1 = E1[int(tcarga1)]/(tcarga1*dt)
P2 = E2[int(tcarga2)]/(tcarga2*dt)
P3 = E3[int(tcarga3)]/(tcarga3*dt)
P4 = E4[int(tcarga4)]/(tcarga4*dt)
Pot = [P1,P2,P3,P4]
print(Pot)
N = [1,2,3,4]
x = np.linspace(0,np.log(4),1000)
x2 = np.linspace(0,np.log(4),1000)
plt.figure()
plt.plot(x,x*(2),'r--', linewidth = 2.5)
plt.plot(x,x,'r',linewidth = 2.5)
plt.plot(np.log(N),np.log(Pot/P1),"o")
plt.plot(x, 3*x, 'k--', linewidth = 2)
plt.fill_between(x,x,x*2,facecolor ='salmon' )
plt.legend(["N^2","N","P(N)/P(1)"], fontsize = 14)
plt.xlabel("N", fontsize = 16)
plt.ylabel("Average Power", fontsize = 16)
plt.show()


