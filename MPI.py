import matplotlib.pylab as plt;
import numpy as np
n = np.array([1,4,9,16,25,36])
t1 = np.array([2.5474, 5.5087, 9.0467, 13.9604, 22.7355, 33.7416])
ti1 = np.array([t1[0], 2*t1[0], 3*t1[0], 4*t1[0], 5*t1[0], 6*t1[0]])
t2 = np.array([2.5474, 2.7543, 3.0155, 3.4901, 4.5471, 5.6236])
ti2 = np.array([2.5474,2.5474,2.5474,2.5474,2.5474,2.5474])
t3 = np.array([571.185, 155.9973, 79.2225, 53.1466, 38.7138, 33.7416])
ti3 = np.array([t3[0]/n[0], t3[0]/n[1], t3[0]/n[2], t3[0]/n[3], t3[0]/n[4], t3[0]/n[5]])

plt.figure(1)
plt.clf()
plt.grid()
plt.plot(n,t1,'-o', n, ti1,'r-o');
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.title(r'Weak Scaling', fontsize = 12);
plt.xlabel(r'Number of Processors', fontsize = 12);
plt.ylabel(r'Run Time (mins)', fontsize = 12);
plt.legend(['Simulated','Ideal'])
plt.show(); 

plt.figure(2)
plt.clf()
plt.grid()
plt.plot(n,t2,'-o',n,ti2,'r-o');
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.title(r'Normalized Weak Scaling', fontsize = 12);
plt.xlabel(r'Number of Processors', fontsize = 12);
plt.ylabel(r'Run Time (mins) / Sqrt of number of Processors', fontsize = 12);
plt.legend(['Simulated','Ideal'])
plt.show(); 

plt.figure(3)
plt.clf()
plt.grid()
plt.loglog(n,t3,'-o',n,ti3,'r-o')
plt.xticks(n,[1,4,9,16,25,36],fontsize = 12)
plt.yticks(fontsize = 12)
plt.title(r'Strong Scaling', fontsize = 12);
plt.xlabel(r'Number of Processors', fontsize = 12);
plt.ylabel(r'Run Time (mins)', fontsize = 12);
plt.legend(['Simulated','Ideal'])
plt.show(); 


plt.figure(4)
plt.clf()
plt.grid()
plt.plot(n,t3,'-o',n,ti3,'r-o');
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.title(r'Strong Scaling', fontsize = 12);
plt.xlabel(r'Number of Processors', fontsize = 12);
plt.ylabel(r'Run Time (mins)', fontsize = 12);
plt.legend(['Simulated','Ideal'])
plt.show(); 