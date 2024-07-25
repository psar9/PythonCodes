import matplotlib.pylab as plt;
import numpy as np
n = np.array([1,4,9,16,25,36])
t1 = np.array([0.712520, 1.39322, 2.20624,  2.89251, 4.47319, 6.08303])
ti1 = np.array([t1[0], 2*t1[0], 3*t1[0], 4*t1[0], 5*t1[0], 6*t1[0]])
t2 = np.array([t1[0]/np.sqrt(n[0]),t1[1]/np.sqrt(n[1]),t1[2]/np.sqrt(n[2]),t1[3]/np.sqrt(n[3]),\
               t1[4]/np.sqrt(n[4]), t1[5]/np.sqrt(n[5])])
ti2 = np.array([t1[0], t1[0], t1[0], t1[0], t1[0], t1[0]])
t3 = np.array([571.185, 72.2931, 26.6366, 14.5929, 8.23379, 6.08303])
ti3 = np.array([t3[0]/n[0], t3[0]/n[1], t3[0]/n[2], t3[0]/n[3], t3[0]/n[4], t3[0]/n[5]])
p = np.array([16,25,36,64,81,100])
t4 = np.array([68.4821, 38.1650, 26.2019, 13.3805, 28.6558, 13.8148])
ti4 = np.array([t4[0]/p[0], t4[0]/p[1], t4[0]/p[2], t4[0]/p[3], t4[0]/p[4], t4[0]/p[5]])

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

plt.figure(5)
plt.clf()
plt.grid()
plt.loglog(p,t4,'-o',p,ti4,'r-o')
plt.xticks(n,[1,4,9,16,25,36],fontsize = 12)
plt.yticks(fontsize = 12)
plt.title(r'Strong Scaling', fontsize = 12);
plt.xlabel(r'Number of Processors', fontsize = 12);
plt.ylabel(r'Run Time (mins)', fontsize = 12);
plt.legend(['Simulated','Ideal'])
plt.show(); 


plt.figure(6)
plt.clf()
plt.grid()
plt.plot(p,t4,'-o',p,ti4,'r-o');
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.title(r'Strong Scaling', fontsize = 12);
plt.xlabel(r'Number of Processors', fontsize = 12);
plt.ylabel(r'Run Time (mins)', fontsize = 12);
plt.legend(['Simulated','Ideal'])
plt.show(); 