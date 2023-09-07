#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 10:09:57 2022

@author: sar
"""

# import critical libraries
import numpy as np;
import matplotlib;
import matplotlib.pylab as plt;

# set some plotting options
matplotlib.rcParams.update({'font.size':14, 'font.family':'serif'});
matplotlib.rcParams.update({'text.usetex':'true'});
matplotlib.rcParams.update(matplotlib.rcParamsDefault)

X = 0
Y = 0
Nx = 10
Ny = 10
Vx = -10
Vy = -10
Nvx = 10
Nvy = 10
t = 0.14
dx = abs(1/Nx)
dy = abs(1/Ny)
dvx = abs((2*Vx)/(Nvx-1))
dvy = abs((2*Vy)/(Nvy-1))
e = 5
Kn = 10**-e
dtx = abs(0.9*dx/Vx)
dty = abs(0.9*dy/Vy)
nt = int(np.ceil(t/dtx))
dt = t/nt

C = 1
#tau = np.inf
tau = Kn/C

a = ([0.5,0,0],[-0.5,0.5,0],[0,0.5,0.5])
w = [0,0.5,0.5]
at = ([0,0,0],[0,0,0],[0,1,0])
wt = [0,0.5,0.5]
nu = len(a)

rho0 = 1
T0 = 1
R = 1
N = 2
sigma0 = 10
sigmax = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
sigmay = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
x = np.zeros(Nx+1)
y = np.zeros(Ny+1)
vx = np.zeros(Nvx)
vy = np.zeros(Nvy)
f = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
rho = rho0*np.ones((Nx+1,Ny+1))
T = T0*np.ones((Nx+1,Ny+1))
ux = np.zeros((Nx+1,Ny+1))
uy = np.zeros((Nx+1,Ny+1))
mx = np.zeros((Nx+1,Ny+1))
my = np.zeros((Nx+1,Ny+1))
#E11 = np.zeros((Nx+1,Ny+1))
#E12 = np.zeros((Nx+1,Ny+1))
#E22 = np.zeros((Nx+1,Ny+1))
E = np.zeros((Nx+1,Ny+1))
fs = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
f1 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
f2 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
f3 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
F1x = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
F2x = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
F3x = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
F1y = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
F2y = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
F3y = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
fM1 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
fM2 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
fM3 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
B = np.zeros((Nx+1,Ny+1,Nvx,Nvy))

for j in range(0,Nx+1):   
    x[j] = X + j*dx
    '''
    if(x[j] <= 0.5):
        rho[j,:] = 2.25#1
        ux[j,:] = 0
        uy[j,:] = 0
        T[j,:] = 1.125#1
            
    else:
        rho[j,:] = 3.0/7.0# 0.125
        ux[j,:] = 0
        uy[j,:] = 0
        T[j,:] = (1.0/6.0)#0.8
    mx[j,:] = rho[j,:]*ux[j,:]
    my[j,:] = rho[j,:]*uy[j,:]
    E[j,:] = 0.5*(N*R*rho[j,:]*T[j,:] + (rho[j,:]*(ux[j,:]**2+uy[j,:]**2)))'''
    
    for k in range(0,Ny+1):
        y[k] = Y + k*dy
        if(y[k] <= 0.5):
            rho[:,k] = 1
            ux[:,k] = 0
            uy[:,k] = 0
            T[:,k] = 1
                
        else:
            rho[:,k] = 0.125
            ux[:,k] = 0
            uy[:,k] = 0
            T[:,k] = 0.8
        mx[:,k] = rho[:,k]*ux[:,k]
        my[:,k] = rho[:,k]*uy[:,k]
        E[:,k] = 0.5*(N*R*rho[:,k]*T[:,k] + (rho[:,k]*(ux[:,k]**2+uy[:,k]**2)))
             
        for l in range(0,Nvx):
            vx[l] = Vx + l*dvx 
            
            for m in range(0,Nvy):
                vy[m] = Vy + m*dvy     
                f[j,k,l,m] = (rho[j,k]/((2*np.pi*R*T[j,k])**(N/2)))*np.exp(-((vx[l]-ux[j,k])**2+(vy[m]-uy[j,k])**2)/(2*R*T[j,k]))


def trap2(vx,vy,F,phix,phiy,dvx,dvy):
    sum = np.zeros((Nx+1,Ny+1))
    
    for j in range(0,Nx+1):                 
        for k in range(0,Ny+1):
            for l in range(1,len(vx)-1):
                for m in range(1,len(vy)-1):
                    sum[j,k] = sum[j,k] + 4*F[j,k,l,m]*phix(vx[l])*phiy(vy[m])
    for j in range(0,Nx+1):                 
        for k in range(0,Ny+1):
            for l in range(1,len(vx)-1):
                sum[j,k] = sum[j,k] + 2*(F[j,k,l,0]*phix(vx[l])*phiy(vy[0])+F[j,k,l,-1]*phix(vx[l])*phiy(vy[-1]))
    for j in range(0,Nx+1):                 
        for k in range(0,Ny+1):
            for m in range(1,len(vy)-1):
                sum[j,k] = sum[j,k] + 2*(F[j,k,0,m]*phix(vx[0])*phiy(vy[m])+F[j,k,-1,m]*phix(vx[-1])*phiy(vy[m]))
    for j in range(0,Nx+1):                 
        for k in range(0,Ny+1):
            sum[j,k] = sum[j,k] + F[j,k,0,0]*phix(vx[0])*phiy(vy[0]) + F[j,k,0,-1]*phix(vx[0])*phiy(vy[-1])\
            +F[j,k,-1,0]*phix(vx[-1])*phiy(vy[0]) + F[j,k,-1,-1]*phix(vx[-1])*phiy(vy[-1])
    
    return dvx*dvy*sum/4


def phi1(p):
    return 1

def phi2(p):
    return p

def phi3(p):
    return 0.5*p**2

#rhot = trap(v,ft,phi1,dv)
#time = 0.0

for n in range(0,nt): 
   
    stx = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
    sty = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
    sc = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
    
    for i in range(0,nu):
            
        smrhox = np.zeros((Nx+1,Ny+1))
        smmxx = np.zeros((Nx+1,Ny+1))
        smmxy = np.zeros((Nx+1,Ny+1))
        smEx11 = np.zeros((Nx+1,Ny+1))
        smEx22 = np.zeros((Nx+1,Ny+1))
        smrhoy = np.zeros((Nx+1,Ny+1))
        smmyx = np.zeros((Nx+1,Ny+1))
        smmyy = np.zeros((Nx+1,Ny+1))
        smEy11 = np.zeros((Nx+1,Ny+1))
        smEy22 = np.zeros((Nx+1,Ny+1))
        
        if(i == 0):  
            
            for j in range(0,Nx+1):                 
                for k in range(0,Ny+1):
                    for l in range(0,Nvx):                     
                        for m in range(0,Nvy):
                                                                                              
                            fM1[j,k,l,m] = (rho[j,k]/((2*np.pi*R*T[j,k])**(N/2)))*np.exp(-((vx[l]-ux[j,k])**2+(vy[m]-uy[j,k])**2)/(2*R*T[j,k]))            
            #f1 = f
            f1 = (tau*f + dt*a[0][0]*fM1)/(tau + dt*a[0][0])
            fs = f1       
        else:
            
            for q in range(0,i):
                
                if(q == 0):
                    Fx = F1x
                    Fy = F1y
                                
                if(q == 1):
                    Fx = F2x
                    Fy = F2y
                        
                smrhox = np.add(smrhox,at[i][q]*trap2(vx,vy,Fx,phi1,phi1,dvx,dvy))
                smrhoy = np.add(smrhoy,at[i][q]*trap2(vx,vy,Fy,phi1,phi1,dvx,dvy))
                
                smmxx = np.add(smmxx,at[i][q]*trap2(vx,vx,Fx,phi2,phi1,dvx,dvy))
                smmxy = np.add(smmxy,at[i][q]*trap2(vx,vx,Fy,phi2,phi1,dvx,dvy))
                smmyx = np.add(smmyx,at[i][q]*trap2(vy,vy,Fx,phi1,phi2,dvx,dvy))
                smmyy = np.add(smmyy,at[i][q]*trap2(vy,vy,Fy,phi1,phi2,dvx,dvy))
                
                smEx11 = np.add(smEx11,at[i][q]*trap2(vx,vx,Fx,phi3,phi1,dvx,dvy))  
                smEy11 = np.add(smEy11,at[i][q]*trap2(vx,vx,Fy,phi3,phi1,dvx,dvy))
                smEx22 = np.add(smEx22,at[i][q]*trap2(vy,vy,Fx,phi1,phi3,dvx,dvy))  
                smEy22 = np.add(smEy22,at[i][q]*trap2(vy,vy,Fy,phi1,phi3,dvx,dvy))
                    
            rhos = rho - (dt/dx)*smrhox - (dt/dy)*smrhoy
            msx = mx - (dt/dx)*smmxx - (dt/dy)*smmxy
            msy = my - (dt/dx)*smmyx - (dt/dy)*smmyy
            #Es11 = E11 - (dt/dx)*smEx11 - (dt/dy)*smEy11
            #Es22 = E22 - (dt/dx)*smEx22 - (dt/dy)*smEy22 
            Es = E - (dt/dx)*smEx11 - (dt/dy)*smEy11 - (dt/dx)*smEx22 - (dt/dy)*smEy22 
            usx = np.divide(msx,rhos)
            usy = np.divide(msy,rhos)
            Ts = np.divide((2*(Es)-np.multiply(rhos,np.square(usx)+np.square(usy))),(rhos*N*R))
                    
            for j in range(0,Nx+1):                    
                for k in range(0,Ny+1):
                    for l in range(0,Nvx):                              
                        for m in range(0,Nvy):
                    
                            if(i == 1):
                                fM2[j,k,l,m] = (rhos[j,k]/((2*np.pi*R*Ts[j,k])**(N/2)))*np.exp(-((vx[l]-usx[j,k])**2+(vy[m]-usy[j,k])**2)/(2*R*Ts[j,k]))
                            if(i == 2):
                                fM3[j,k,l,m] = (rhos[j,k]/((2*np.pi*R*Ts[j,k])**(N/2)))*np.exp(-((vx[l]-usx[j,k])**2+(vy[m]-usy[j,k])**2)/(2*R*Ts[j,k]))
            
            sm = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
            sfx = np.zeros((Nx+1,Ny+1,Nvx,Nvy))   
            sfy = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
            
            for q in range(0,i):
                
                if(q == 0): 
                    sm = sm + (a[i][q]/tau)*(fM1 - f1) 
                    sfx = sfx + at[i][q]*F1x
                    sfy = sfx + at[i][q]*F1y
                                
                if(q == 1):
                    sm = sm + (a[i][q]/tau)*(fM2 - f2) 
                    sfx = sfx + at[i][q]*F2x
                    sfy = sfy + at[i][q]*F2y
                        
            B = f - (dt/dx)*sfx - (dt/dy)*sfy + dt*sm
           
            if(i == 1):
                f2 = (tau*B + dt*a[i][i]*fM2)/(tau + dt*a[i][i])
                fs = f2
            if(i == 2):
                f3 = (tau*B + dt*a[i][i]*fM3)/(tau + dt*a[i][i])
                fs = f3
            
        for j in range(1,Nx):
            for k in range(0,Ny+1):
                for l in range(0,Nvx):
                    for m in range(0,Nvy):
                        #sigma[j,k] = (f1[j+1,k] - f1[j-1,k])/2
                        s1 = fs[j+1,k,l,m] - fs[j,k,l,m]
                        s2 = fs[j,k,l,m] - fs[j-1,k,l,m]
                        if (s1>0 and s2>0):
                            sigmax[j,k,l,m] = min(s1,s2)
                        elif (s1<0 and s2<0):
                            sigmax[j,k,l,m] = max(s1,s2)
                        else:
                            sigmax[j,k,l,m]= 0      
                    
        for j in range(0,Nx+1):             
            for k in range(0,Ny+1):
                for l in range(0,Nvx):             
                    for m in range(0,Nvy):
                            
                        if(j == 0):
                            Fp = max(vx[l],0)*fs[j,k,l,m] + min(vx[l],0)*fs[j+1,k,l,m]
                            Fn = vx[l]*fs[j,k,l,m]
        
                        elif(j == Nx):  
                            Fp = vx[l]*fs[j,k,l,m]
                            Fn = max(vx[l],0)*fs[j-1,k,l,m] + min(vx[l],0)*fs[j,k,l,m]
                    
                        else:
                            Fp = max(vx[l],0)*(fs[j,k,l,m]+0.5*sigmax[j,k,l,m]) + min(vx[l],0)*(fs[j+1,k,l,m]-0.5*sigmax[j+1,k,l,m])
                            Fn = max(vx[l],0)*(fs[j-1,k,l,m]+0.5*sigmax[j-1,k,l,m]) + min(vx[l],0)*(fs[j,k,l,m]-0.5*sigmax[j,k,l,m])
                            #Fn = max(v[k],0)*(f1[j-1,k]-0.5*sigma[j-1]) + min(v[k],0)*(f1[j,k]+0.5*sigma[j])
                    
                        if(i == 0):
                            F1x[j,k,l,m] = Fp - Fn
                        
                        elif(i == 1):
                            F2x[j,k,l,m] = Fp - Fn
                    
                        elif(i == 2):
                            F3x[j,k,l,m] = Fp - Fn   
        
        for j in range(0,Nx+1):
            for k in range(1,Ny):
                for l in range(0,Nvx):
                    for m in range(0,Nvy):
                        #sigma[j,k] = (f1[j+1,k] - f1[j-1,k])/2
                        s1 = fs[j,k+1,l,m] - fs[j,k,l,m]
                        s2 = fs[j,k,l,m] - fs[j,k-1,l,m]
                        if (s1>0 and s2>0):
                            sigmay[j,k,l,m] = min(s1,s2)
                        elif (s1<0 and s2<0):
                            sigmay[j,k,l,m] = max(s1,s2)
                        else:
                            sigmay[j,k,l,m]= 0  

        for j in range(0,Nx+1):             
            for k in range(0,Ny+1):
                for l in range(0,Nvx):             
                    for m in range(0,Nvy):
                            
                        if(k == 0):
                            Fp = max(vy[m],0)*fs[j,k,l,m] + min(vy[m],0)*fs[j,k+1,l,m]
                            Fn = vy[m]*fs[j,k,l,m]
        
                        elif(k == Nx):  
                            Fp = vy[m]*fs[j,k,l,m]
                            Fn = max(vy[m],0)*fs[j,k-1,l,m] + min(vy[m],0)*fs[j,k,l,m]
                    
                        else:
                            Fp = max(vy[m],0)*(fs[j,k,l,m]+0.5*sigmay[j,k,l,m]) + min(vy[m],0)*(fs[j,k+1,l,m]-0.5*sigmay[j,k+1,l,m])
                            Fn = max(vy[m],0)*(fs[j,k-1,l,m]+0.5*sigmay[j,k-1,l,m]) + min(vy[m],0)*(fs[j,k,l,m]-0.5*sigmay[j,k,l,m])
                            #Fn = max(v[k],0)*(f1[j-1,k]-0.5*sigma[j-1]) + min(v[k],0)*(f1[j,k]+0.5*sigma[j])
                
                        if(i == 0):
                            F1y[j,k,l,m] = Fp - Fn
                
                        elif(i == 1):
                            F2y[j,k,l,m] = Fp - Fn
                    
                        elif(i == 2):
                            F3y[j,k,l,m] = Fp - Fn                        
        
        if(i == 0):
            stx = stx + wt[i]*F1x
            sty = sty + wt[i]*F1y
            sc = sc + (w[i]/tau)*(fM1 - f1)
        if(i == 1):
            stx = stx + wt[i]*F2x
            sty = sty + wt[i]*F2y
            sc = sc + (w[i]/tau)*(fM2 - f2)
        if(i == 2):
            stx = stx + wt[i]*F3x
            sty = sty + wt[i]*F3y
            sc = sc + (w[i]/tau)*(fM3 - f3)
        
    f = f - (dt/dx)*stx - (dt/dy)*sty + dt*sc
  
    rho = trap2(vx,vy,f,phi1,phi1,dvx,dvy)
    mx = trap2(vx,vx,f,phi2,phi1,dvx,dvy)
    my = trap2(vy,vy,f,phi1,phi2,dvx,dvy)
    E11 = trap2(vx,vx,f,phi3,phi1,dvx,dvy)
    #E12 = 0.5*trap2(vx,vy,f,phi1,phi2,dvx,dvy)
    E22 = trap2(vy,vy,f,phi1,phi3,dvx,dvy)
    E = E11 + E22    
    ux = np.divide(mx,rho)
    uy = np.divide(my,rho)
    T = np.divide((2*E-np.multiply(rho,np.square(ux)+np.square(uy))),(rho*N*R))
    
    #time = time + dt
    #print(time)
'''
xm = np.zeros((Nx+1,Ny+1));
ym = np.zeros((Nx+1,Ny+1));

for i in range(0,Ny+1):
    xm[:,i] = x;
    
for j in range(0,Nx+1):
    ym[j,:] = y;'''


xm, ym = np.meshgrid(x,y)

'''    
plt.figure(1);
plt.clf();
plt.gca().set_aspect('auto');
plt.gca().set_xlim([x[0]-0.5*dx,x[Nx-1]+0.5*dx]);
plt.gca().set_ylim([y[0]-0.5*dy,y[Ny-1]+0.5*dy]);
plt.pcolor(x,y,rho);
plt.grid(True);
title = r"Density: $\rho(x,y)$";
plt.title(title);
plt.xlabel(r'$x$');
plt.ylabel(r'$y$');

plt.figure(2);
plt.clf();
plt.gca().set_aspect('auto');
plt.gca().set_xlim([x[0]-0.5*dx,x[Nx-1]+0.5*dx]);
plt.gca().set_ylim([y[0]-0.5*dy,y[Ny-1]+0.5*dy]);
plt.pcolor(x,y,ux);
plt.grid(True);
title = r"Velocity: $ux(x,y)$";
plt.title(title);
plt.xlabel(r'$x$');
plt.ylabel(r'$y$');

plt.figure(3);
plt.clf();
plt.gca().set_aspect('auto');
plt.gca().set_xlim([x[0]-0.5*dx,x[Nx-1]+0.5*dx]);
plt.gca().set_ylim([y[0]-0.5*dy,y[Ny-1]+0.5*dy]);
plt.pcolor(x,y,uy);
plt.grid(True);
title = r"Velocity: $uy(x,y)$";
plt.title(title);
plt.xlabel(r'$x$');
plt.ylabel(r'$y$');

plt.figure(4);
plt.clf();
plt.gca().set_aspect('auto');
plt.gca().set_xlim([x[0]-0.5*dx,x[Nx-1]+0.5*dx]);
plt.gca().set_ylim([y[0]-0.5*dy,y[Ny-1]+0.5*dy]);
plt.pcolor(x,y,T);
plt.grid(True);
title = r"Temperature: $T(x,y)$";
plt.title(title);
plt.xlabel(r'$x$');
plt.ylabel(r'$y$');
'''    

plt.figure(1)
plt.clf()
plt.grid()
#plt.plot(x,rho,label = 'Numerical');
#plt.plot(x,rhot,label = 'Exact');
ax = plt.axes(projection='3d')
ax.plot_surface(xm, ym, rho)
plt.title(r'Surface plot of rho vs x and y');
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('rho');
plt.show()

plt.figure(2)
plt.clf()
plt.grid()
ax = plt.axes(projection='3d')
ax.plot_surface(xm, ym, ux)
plt.title(r'Surface plot of ux vs x and y');
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('ux');
plt.show()

plt.figure(3)
plt.clf()
plt.grid()
ax = plt.axes(projection='3d')
ax.plot_surface(xm, ym, uy)
plt.title(r'Surface plot of uy vs x and y');
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('uy');
plt.show() 

plt.figure(4)
plt.clf()
plt.grid()
ax = plt.axes(projection='3d')
ax.plot_surface(xm, ym, T)
plt.title(r'Surface plot of T vs x and y');
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('T');
plt.show() 

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



