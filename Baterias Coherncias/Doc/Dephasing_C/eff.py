import numpy as np 
import matplotlib.pyplot as plt 

data = np.loadtxt('eff_eff_deph_vs_tf.txt')

tf = data[:,0]
eff = data[:,1]
eff_deph = data[:,2]
Pow = data[:,3]
Pow_deph = data[:,4]
eff_otto = np.ones(len(eff))*(1 - 1/6.0)

x = np.linspace(2,18,1000)
zero = np.zeros(len(x))

plt.figure()
plt.plot(tf, eff, linewidth = 2)
plt.plot(tf, eff_deph, 'C2',linewidth = 2)
plt.plot(tf, eff_otto, '--C7',linewidth = 2)
plt.plot(x, zero, '--r',linewidth = 2)
plt.fill_between(x,zero,-0.5*np.ones(len(x)),facecolor ='salmon' )
plt.xlabel(r'$t_f$', fontsize = 14)
plt.ylabel(r'$\eta$', fontsize = 14)
plt.legend([r'$\eta$', r'$\eta_{deph}$', r'$\eta_{Otto}$'], fontsize = 11)
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.ylim([-0.1,1])
plt.xlim([2,18])
plt.show()

plt.figure()
plt.plot(tf, Pow, linewidth = 2)
plt.plot(tf, Pow_deph, 'C2',linewidth = 2)
plt.plot(x, zero, '--r',linewidth = 2)
plt.fill_between(x,zero,-0.1*np.ones(len(x)),facecolor ='salmon' )
plt.xlabel(r'$t_f$', fontsize = 14)
plt.ylabel(r'$P$', fontsize = 14)
plt.legend([r'$P$', r'$P_{deph}$'], fontsize = 11)
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.ylim([-0.05,0.05])
plt.xlim([2,18])
plt.show()

plt.figure()
plt.plot(Pow, eff, linewidth = 2)
plt.plot(Pow_deph, eff_deph, 'C2',linewidth = 2)
#plt.plot(x, zero, '--r',linewidth = 2)
#plt.fill_between(x,zero,-0.1*np.ones(len(x)),facecolor ='salmon' )
plt.xlabel(r'$P$', fontsize = 14)
plt.ylabel(r'$\eta$', fontsize = 14)
plt.legend([r'$\eta$', r'$\eta_{deph}$'], fontsize = 11)
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.ylim([-0.05,1])
plt.xlim([0,0.6])
plt.show()