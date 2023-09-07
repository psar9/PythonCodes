#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 12:04:24 2022

@author: PreetiSar2
"""

# import critical libraries
import numpy as np;
import matplotlib;
import matplotlib.pylab as plt;

# set some plotting options
matplotlib.rcParams.update({'font.size':14, 'font.family':'serif'});
matplotlib.rcParams.update({'text.usetex':'true'});
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
  
'''X = -1
Nx = 80
Nv = 80
t = 0.04
dx = abs(2*X/Nx)'''

X = 0
Nx = 200
Nv = 200#81
t = 0.16
dx = abs(1/Nx)
#tau = np.inf 
V = -10
dv = abs(2*V/(Nv-1))     
dt = abs(0.9*dx/V)
nt = int(np.ceil(t/dt))
dt = t/nt

C = 1
rho0 = 1
T0 = 1
e = 5
Kn = 2**-e
tau = Kn/C

a = ([0.5,0,0],[-0.5,0.5,0],[0,0.5,0.5])
w = [0,0.5,0.5]
at = ([0,0,0],[0,0,0],[0,1,0])
wt = [0,0.5,0.5]
nu = len(a)

R = 1
N = 1
sigma0 = 10
sigma = np.zeros((Nx+1,Nv))
x = np.zeros(Nx+1)
v = np.zeros(Nv)
f = np.zeros((Nx+1,Nv))
ft = np.zeros((Nx+1,Nv))
rho = rho0*np.ones(Nx+1)
rhot = np.zeros(Nx+1)
T = T0*np.ones(Nx+1)
u = np.zeros(Nx+1)
ut = np.zeros(Nx+1)
m = np.zeros(Nx+1)
E = np.zeros(Nx+1)
F = np.zeros((Nx+1,Nv))
f1 = np.zeros((Nx+1,Nv))
f2 = np.zeros((Nx+1,Nv))
f3 = np.zeros((Nx+1,Nv))
F1 = np.zeros((Nx+1,Nv))
F2 = np.zeros((Nx+1,Nv))
F3 = np.zeros((Nx+1,Nv))
fM1 = np.zeros((Nx+1,Nv))
fM2 = np.zeros((Nx+1,Nv))
fM3 = np.zeros((Nx+1,Nv))
B = np.zeros((Nx+1,Nv))

'''for j in range(0,Nx+1):
    x[j] = X + j*dx
    #u0[j] = 0
    u[j] = (1/sigma0)*(np.exp(-(sigma0*x[j]-1)**2)-2*np.exp(-(sigma0*x[j]+3)**2))    
    m[j] = rho[j]*u[j]
    E[j] = 0.5*(N*R*rho[j]*T[j] + (rho[j]*u[j]**2))
    
    for k in range(0,Nv):
        v[k] = V + k*dv        
        f[j,k] = (rho[j]/((2*np.pi*R*T[j])**0.5))*np.exp((-(v[k]-u[j])**2)/(2*R*T[j]))
        #ut[j] = (1/sigma0)*(np.exp(-(sigma0*(x[j] - v[k]*t)-1)**2)-2*np.exp(-(sigma0*(x[j] - v[k]*t)+3)**2))
        #ft[j,k] = (rho[j]/((2*np.pi*R*T[j])**0.5))*np.exp((-(v[k]-ut[j])**2)/(2*R*T[j]))'''

for j in range(0,Nx+1):
    x[j] = X + j*dx
    if(x[j] <= 0.5):
        rho[j] = 1#2.25
        u[j] = 0
        T[j] = 1#1.125
        
    else:
        rho[j] = 0.125#(3/7)
        u[j] = 0
        T[j] = 0.8#(1/6)
    
    m[j] = rho[j]*u[j]
    E[j] = 0.5*(N*R*rho[j]*T[j] + (rho[j]*u[j]**2))
    for k in range(0,Nv):
        v[k] = V + k*dv        
        f[j,k] = (rho[j]/((2*np.pi*R*T[j])**0.5))*np.exp(-((v[k]-u[j])**2)/(2*R*T[j])) 

def trap(v,F,phi,dv):
    sum = np.zeros(Nx+1)
    sum = sum + 0.5*F[:,0]*phi(v[0]) + 0.5*F[:,-1]*phi(v[-1])
    for k in range(1,len(v)-1):
        sum = sum + F[:,k]*phi(v[k])          
    return dv*sum  

'''def trap(v,F,phi,dv):
    sum = np.zeros(Nx+1)
    for j in range(0,Nx+1):
        sum[j] = sum[j] + 0.5*F[j,0]*phi(v[0]) + 0.5*F[j,-1]*phi(v[-1])
    for j in range(0,Nx+1):
        for k in range(1,len(v)-1):
            sum[j] = sum[j] + F[j,k]*phi(v[k])    
    return dv*sum'''

def phi1(y):
    return 1

def phi2(y):
    return y

def phi3(y):
    return 0.5*y**2

#rhot = trap(v,ft,phi1,dv)

#time = 0.0
for n in range(0,nt): 
   
    st = np.zeros((Nx+1,Nv))
    sc = np.zeros((Nx+1,Nv))
    
    for i in range(0,nu):
            
        smrho = np.zeros(Nx+1)
        smm = np.zeros(Nx+1)
        smE = np.zeros(Nx+1)
        
        if(i == 0):  
            
            for j in range(0,Nx+1):                
                  
                for k in range(0,Nv):
                                                                                  
                    fM1[j,k] = (rho[j]/((2*np.pi*R*T[j])**(N/2)))*np.exp((-(v[k]-u[j])**2)/(2*R*T[j])) 
            
            #f1 = f
            f1 = (tau*f + dt*a[0][0]*fM1)/(tau + dt*a[0][0])
            fs = f1
                    
        else:
            
            for l in range(0,i):
                if(l == 0):
                    F = F1
                                
                if(l == 1):
                    F = F2
                            
                smrho = np.add(smrho,at[i][l]*trap(v,F,phi1,dv))
                smm = np.add(smm,at[i][l]*trap(v,F,phi2,dv))
                smE = np.add(smE,at[i][l]*trap(v,F,phi3,dv))  
                    
            rhos = rho - (dt/dx)*smrho
            ms = m - (dt/dx)*smm
            Es = E - (dt/dx)*smE
            us = np.divide(ms,rhos)
            Ts = np.divide((2*Es-np.multiply(rhos,np.square(us))),(rhos*N*R))
                    
            for j in range(0,Nx+1):
                      
                for k in range(0,Nv):
                    
                    if(i == 1):
                        fM2[j,k] = (rhos[j]/((2*np.pi*R*Ts[j])**(N/2)))*np.exp((-(v[k]-us[j])**2)/(2*R*Ts[j]))
                    if(i == 2):
                        fM3[j,k] = (rhos[j]/((2*np.pi*R*Ts[j])**(N/2)))*np.exp((-(v[k]-us[j])**2)/(2*R*Ts[j]))
            
            sm = np.zeros((Nx+1,Nv))
            sf = np.zeros((Nx+1,Nv))        
            
            for l in range(0,i):
                
                if(l == 0):
                    sm = sm + (a[i][l]/tau)*(fM1 - f1) 
                    sf = sf + at[i][l]*F1
                                
                if(l == 1):
                    sm = sm + (a[i][l]/tau)*(fM2 - f2) 
                    sf = sf + at[i][l]*F2
                        
            B = f - (dt/dx)*sf + dt*sm
           
            if(i == 1):
                #f2 = B
                f2 = (tau*B + dt*a[i][i]*fM2)/(tau + dt*a[i][i])
                fs = f2
            if(i == 2):
                #f3 = B
                f3 = (tau*B + dt*a[i][i]*fM3)/(tau + dt*a[i][i])
                fs = f3
                
        for j in range(1,Nx):
            for k in range(0,Nv):
    
               #sigma[j,k] = (fs[j+1,k] - fs[j-1,k])/2
               s1 = fs[j+1,k] - fs[j,k]
               s2 = fs[j,k] - fs[j-1,k]
               if (s1>0 and s2>0):
                   sigma[j,k] = min(s1,s2)
               elif (s1<0 and s2<0):
                   sigma[j,k] = max(s1,s2)
               else:
                   sigma[j,k] = 0
               '''
               r = 1
               if ((fs[j,k] - fs[j-1,k]) != 0):
                   r = (fs[j+1,k] - fs[j,k])/(fs[j,k] - fs[j-1,k])
               sigma[j,k] = max(0,min(1,r))'''
        
        for j in range(0,Nx+1):
                  
            for k in range(0,Nv):

                if(j == 0):
                    Fp = max(v[k],0)*fs[j,k] + min(v[k],0)*fs[j+1,k]
                    Fn = v[k]*fs[j,k]           
        
                elif(j == Nx):  
                    Fp = v[k]*fs[j,k]
                    Fn = max(v[k],0)*fs[j-1,k] + min(v[k],0)*fs[j,k]
                    
                else:
                    Fp = max(v[k],0)*(fs[j,k]+0.5*sigma[j,k]) + min(v[k],0)*(fs[j+1,k]-0.5*sigma[j+1,k])
                    Fn = max(v[k],0)*(fs[j-1,k]+0.5*sigma[j-1,k]) + min(v[k],0)*(fs[j,k]-0.5*sigma[j,k])
                    #Fn = max(v[k],0)*(fs[j-1,k]-0.5*sigma[j-1]) + min(v[k],0)*(fs[j,k]+0.5*sigma[j])
                
                if (i == 0):
                    F1[j,k] = Fp - Fn
                    
                elif(i == 1):
                    F2[j,k] = Fp - Fn
                    
                elif(i == 2):
                    F3[j,k] = Fp - Fn
        
        if(i == 0):
            st = st + wt[i]*F1
            sc = sc + (w[i]/tau)*(fM1 - f1)
        elif(i == 1):
            st = st + wt[i]*F2
            sc = sc + (w[i]/tau)*(fM2 - f2)
        elif(i == 2):
            st = st + wt[i]*F3
            sc = sc + (w[i]/tau)*(fM3 - f3)
        
    f = f - (dt/dx)*st + dt*sc
    
    rho = trap(v,f,phi1,dv)
    m = trap(v,f,phi2,dv)
    E = trap(v,f,phi3,dv)
        
    u = np.divide(m,rho)
    T = np.divide((2*E-np.multiply(rho,np.square(u))),(rho*N*R))
    #time = time + dt
    #print(time)

plt.figure(1)
#plt.clf()
plt.grid()
plt.plot(x,rho,label = 'Numerical');
#plt.plot(x,rhot,label = 'Exact');
#plt.plot(x,rhoc,label = 'Exact_C');
plt.title(r'Plot of rho vs x');
plt.xlabel(r'x');
plt.ylabel(r'rho');
plt.legend();
#plt.show();  

plt.figure(2)
#plt.clf()
plt.grid()
plt.plot(x,T);
plt.title(r'Plot of T vs x');
plt.xlabel(r'x');
plt.ylabel(r'T');
#plt.show();   

plt.figure(3)
#plt.clf()
plt.grid()
plt.plot(x,u);
plt.title(r'Plot of u vs x');
plt.xlabel(r'x');
plt.ylabel(r'u');
#plt.show();

plt.figure(4)
plt.clf()
plt.pcolor(np.transpose(f))
plt.colorbar()
plt.xlabel(r'x');
plt.ylabel(r'v');
plt.show(); 



#q = list(map(float,open("rho.txt").read().split()))
#q = rhot
#diff = q - rho
#z = int((len(q)-1)/Nx)
#for i in range(0,Nx+1):
#    diff[i] = rho[i] - q[z*i]    

#print(np.linalg.norm(diff, 1)/np.linalg.norm(q, 1))
#print(np.linalg.norm(diff, 1))
#print(np.linalg.norm(diff, np.inf)/np.linalg.norm(q, np.inf))
#print(np.linalg.norm(diff, np.inf))

#np.savetxt('rho.txt',rho)
        


