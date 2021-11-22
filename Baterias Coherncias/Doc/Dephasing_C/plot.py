import numpy as np 
import matplotlib.pyplot as plt 

#data_dr = np.loadtxt('Test_driven.txt')
data = np.loadtxt('Test_constant.txt')

t = data[:,0]
E = data[:,1]
Q = data[:,2]
W = data[:,3]
"""
t_dr = data_dr[:,0]
E_dr = data_dr[:,1]
Q_dr = data_dr[:,2]
W_dr = data_dr[:,3]
"""

print(np.argmax(E))

plt.figure()
#plt.plot(t, E_dr, linewidth = 2)
plt.plot(E,linewidth = 2)
plt.plot(Q,linewidth = 2)
plt.plot(W,linewidth = 2)
plt.show()
"""
plt.figure()
plt.plot(t, E, linewidth = 2)
plt.plot(t, Q, linewidth = 2)
plt.plot(t, W, linewidth = 2)
plt.plot(t, E_dr, linewidth = 2)
plt.plot(t, Q_dr, linewidth = 2)
plt.plot(t, W_dr, linewidth = 2)
plt.xlabel('Time', fontsize = 12)
plt.ylabel('Energy', fontsize = 12)
plt.legend(['E', 'Q', 'W','E_dr', 'Q_dr', 'W_dr'], fontsize = 11)
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.show()"""