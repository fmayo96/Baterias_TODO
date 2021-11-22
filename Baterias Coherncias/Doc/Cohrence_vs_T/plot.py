import numpy as np 
import matplotlib.pyplot as plt 

data1 = np.loadtxt("beta=0.01.txt")
data2 = np.loadtxt("beta=0.05.txt")
data3 = np.loadtxt("beta=0.1.txt")
data4 = np.loadtxt("beta=0.2.txt")

h_hot = data1[:,0]
eff_1 = data1[:,1]
eff_deph_1 = data1[:,2]
Pow_1 = data1[:,3]
Pow_deph_1 = data1[:,4]

eff_2 = data2[:,1]
eff_deph_2 = data2[:,2]
Pow_2 = data2[:,3]
Pow_deph_2 = data2[:,4]

eff_3 = data3[:,1]
eff_deph_3 = data3[:,2]
Pow_3 = data3[:,3]
Pow_deph_3 = data3[:,4]

eff_4 = data4[:,1]
eff_deph_4 = data4[:,2]
Pow_4 = data4[:,3]
Pow_deph_4 = data4[:,4]

plt.figure()
plt.plot(h_hot, eff_1, 'C0', linewidth = 2)
plt.plot(h_hot, eff_2, 'C1', linewidth = 2)
plt.plot(h_hot, eff_3, 'C2', linewidth = 2)
plt.plot(h_hot, eff_4, 'C3', linewidth = 2)
plt.plot(h_hot, 1-(1/h_hot), '--k', linewidth = 2)
plt.plot(h_hot, eff_deph_1, '--C0', linewidth = 2)
plt.plot(h_hot, eff_deph_2, '--C1', linewidth = 2)
plt.plot(h_hot, eff_deph_3, '--C2', linewidth = 2)
plt.plot(h_hot, eff_deph_4, '--C3', linewidth = 2)
plt.xlabel(r"$h_{hot}$", fontsize = 14)
plt.ylabel(r"$\eta$", fontsize = 14)
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.legend([r'$\beta_h = 0.01$', r'$\beta_h = 0.05$', r'$\beta_h = 0.1$', r'$\beta_h = 0.2$', r'$\eta_{Otto}$'], fontsize = 11)
plt.show()

plt.figure()
plt.plot(h_hot, Pow_1, 'C0', linewidth = 2)
plt.plot(h_hot, Pow_2, 'C1', linewidth = 2)
plt.plot(h_hot, Pow_3, 'C2', linewidth = 2)
plt.plot(h_hot, Pow_4, 'C3', linewidth = 2)
plt.plot(h_hot, Pow_deph_1, '--C0', linewidth = 2)
plt.plot(h_hot, Pow_deph_2, '--C1', linewidth = 2)
plt.plot(h_hot, Pow_deph_3, '--C2', linewidth = 2)
plt.plot(h_hot, Pow_deph_4, '--C3', linewidth = 2)
plt.xlabel(r"$h_{hot}$", fontsize = 14)
plt.ylabel(r"$P$", fontsize = 14)
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.legend([r'$\beta_h = 0.01$', r'$\beta_h = 0.05$', r'$\beta_h = 0.1$', r'$\beta_h = 0.2$'], fontsize = 11)
plt.show()