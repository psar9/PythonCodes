# import critical libraries
import numpy as np;
import matplotlib;
import matplotlib.pyplot as plt

# set some plotting options
matplotlib.rcParams.update({'font.size':14, 'font.family':'serif'});
matplotlib.rcParams.update({'text.usetex':'true'});
matplotlib.rcParams.update(matplotlib.rcParamsDefault)

'''
X = -7.5
Nx = 20
V = -3
Nv = 100
t = 15
dx = abs(2*X/Nx)
dv = abs((4-V)/(Nv-1))
n = 8
e = 3**-n'''

X = 0
Nx = 100
#V = -10
V = -4.5
Nv = 100
t = 0.16
dx = abs(1/Nx)
dv = abs((2*V)/(Nv-1))
dt = abs(0.9*dx/V)
nt = int(np.ceil(t/dt))
dt = t/nt

def trap(F,phi):
    sum = np.zeros(len(F))
    sum = sum + 0.5*F[:,0]*phi(v[0]) + 0.5*F[:,-1]*phi(v[-1])
    for k in range(1,len(v)-1):
        sum = sum + F[:,k]*phi(v[k])          
    return dv*sum  

def trap1(F,phi):
    sum = 0
    sum = sum + 0.5*F[0]*phi(v[0]) + 0.5*F[-1]*phi(v[-1])
    for k in range(1,len(v)-1):
        sum = sum + F[k]*phi(v[k])          
    return dv*sum

def m1(y):
    return 1

def m2(y):
    return y

def m3(y):
    return 0.5*y**2

def m4(y):
    return 0.5*y**3
    
def vanleer(r):
    psi = (r+abs(r))/(1+abs(r))
    #psi = 1
    return psi

n = 14
    
e = 2**-n
x = np.zeros(Nx+1)
v = np.zeros(Nv)
rho = np.ones(Nx+1)
T = np.ones(Nx+1)
u = np.ones(Nx+1)
rhom = np.ones(Nx+2)
Tm = np.ones(Nx+2)
um = np.ones(Nx+2)
Mav = np.zeros((Nx+2,Nv))
M = np.zeros((Nx+1,Nv))
g = np.zeros((Nx+2,Nv))
phi1 = np.zeros((Nx+2,Nv))
phi2 = np.zeros((Nx+2,Nv))
gd = np.zeros((Nx+1,Nv))
a1 = np.zeros((Nx+2,Nv))
b1 = np.zeros((Nx+2,Nv))
a2 = np.zeros((Nx+2,Nv))
b2 = np.zeros((Nx+2,Nv))
pi1 = np.zeros((Nx+2,Nv))
pi2 = np.zeros((Nx+2,Nv))
f = np.zeros((Nx+1,3))
s1 = np.zeros(Nx+2)
s2 = np.zeros(Nx+2)
s3 = np.zeros(Nx+2)
l1 = np.zeros(3)
l2 = np.zeros(3)
l3 = np.zeros(3)
r1 = np.zeros(3)
r2 = np.zeros(3)
r3 = np.zeros(3)
z1 = np.zeros((Nx+2,3))
z2 = np.zeros((Nx+2,3))
z3 = np.zeros((Nx+2,3))
theta1 = np.ones(Nx+2)
theta2 = np.ones(Nx+2)
theta3 = np.ones(Nx+2)
Ap = np.zeros((Nx+2,3))
An = np.zeros((Nx+2,3))
F = np.zeros((Nx+2,3))
Q = np.zeros((Nx+1,3))

'''for j in range(0,Nx+1):
    x[j] = X + j*dx
    if(x[j] < 0):
        rho[j] = 1
        u[j] = 1.2
        T[j] = 0.1
        
    else:
        rho[j] = 1.65
        u[j] = 0.72
        T[j] = 0.4
    U1[j] = rho[j]
    U2[j] = rho[j]*u[j]
    U3[j] = 0.5*rho[j]*u[j]**2 + 0.5*rho[j]*T[j]'''
    
for j in range(0,Nx+1):
    x[j] = X + j*dx
    if(x[j] <= 0.5):
        rho[j] = 1
        u[j] = 0
        T[j] = 1     
    else:
        rho[j] = 0.125
        u[j] = 0
        T[j] = 0.8
        
    Q[j,0] = rho[j]
    Q[j,1] = rho[j]*u[j]
    Q[j,2] = (0.5*rho[j]*u[j]**2) + (0.5*rho[j]*T[j])
 
for k in range(0,Nv):
    v[k] = V + k*dv   
      
for i in range(0,nt):
    
    for j in range(0,Nx+1):
        for k in range(0,Nv):
            M[j,k] = (rho[j]/((2*np.pi*T[j])**0.5))*np.exp((-(v[k]-u[j])**2)/(2*T[j]))
    
    for j in range(0,Nx+2):
        jm1 = (j > 0)*(j-1) + (j == 0)*0
        j0 = (j < Nx+1)*j + (j == Nx+1)*(Nx)
        rhom[j] = (rho[j0]+rho[jm1])/2
        um[j] = (u[j0]+u[jm1])/2
        Tm[j] = (T[j0]+T[jm1])/2
    
    for j in range(0,Nx+2):
        for k in range(0,Nv):
            Mav[j,k] = (rhom[j]/((2*np.pi*Tm[j])**0.5))*np.exp((-(v[k]-um[j])**2)/(2*Tm[j]))
          
    for j in range(0,Nx+2):
        jm1 = (j > 0)*(j-1) + (j == 0)*0
        j0 = (j < Nx+1)*j + (j == Nx+1)*(Nx)
        jp1 = (j < Nx)*(j+1) + (j == Nx)*(Nx)
        for k in range(0,Nv):        
            phi1[j,k] = (max(v[k],0)*((g[j,k] - g[jm1,k])/dx)) + (min(v[k],0)*((g[jp1,k] - g[j,k])/dx)) 
            phi2[j,k] = v[k]*((M[j0,k] - M[jm1,k])/dx)
   
    for j in range(0,Nx+2):
        for k in range(0,Nv): 
            a1[j,k] = (v[k] - um[j])*phi1[j,k]
            a2[j,k] = (v[k] - um[j])*phi2[j,k]
            b1[j,k] = ((((v[k] - um[j])**2)/(2*Tm[j])) - 0.5)*phi1[j,k]      
            b2[j,k] = ((((v[k] - um[j])**2)/(2*Tm[j])) - 0.5)*phi2[j,k]

    c1 = trap(a1,m1)
    c2 = trap(a2,m1)
    d1 = trap(b1,m1)
    d2 = trap(b2,m1)
    
    p1 = trap(phi1,m1)
    p2 = trap(phi2,m1)
    
    for j in range(0,Nx+2):
        for k in range(0,Nv): 
            pi1[j,k] = (1/rhom[j])*(p1[j] + (((v[k]-um[j])*c1[j])/Tm[j]) + (((((v[k]-um[j])**2)/(2*Tm[j]))-0.5)*2*d1[j]))*Mav[j,k]
            pi2[j,k] = (1/rhom[j])*(p2[j] + (((v[k]-um[j])*c2[j])/Tm[j]) + (((((v[k]-um[j])**2)/(2*Tm[j]))-0.5)*2*d2[j]))*Mav[j,k]
    
    g = ((e/(e+dt))*(g - (dt*(phi1 - pi1)))) - ((dt/(e+dt))*(phi2 - pi2))
    
    for j in range(0,Nx+1):
        for k in range(0,Nv):
            gd[j,k] = v[k]*((g[j+1,k] - g[j,k])/dx)

    for j in range(0,Nx+1):
        f[j,0] = rho[j]*u[j]
        f[j,1] = rho[j]*(u[j]**2 + T[j])
        f[j,2] = (0.5*u[j]**2 + 1.5*T[j])*u[j]*rho[j]

    for j in range(0,Nx+2):  
        jm1 = (j > 0)*(j-1) + (j == 0)*0
        j0 = (j < Nx+1)*j + (j == Nx+1)*(Nx)
        Tav = 0.5*(T[jm1]+T[j0])
        uav = 0.5*(u[jm1]+u[j0])
        c = np.sqrt(3*Tav)
        oc2 = 1/(3*Tav)
        s1[j] = uav - c
        s2[j] = uav
        s3[j] = uav + c
        l1[0] = 0.5*uav*(uav+c)*oc2
        l1[1] = -0.5*(c+2*uav)*oc2
        l1[2] = oc2
        l2[0] = (c*c-uav*uav)*oc2
        l2[1] = 2.0*uav*oc2
        l2[2] = -2.0*oc2
        l3[0] = 0.5*uav*(uav-c)*oc2
        l3[1] = 0.5*(c-2*uav)*oc2
        l3[2] = oc2
        r1[0] = 1.0
        r1[1] = uav-c
        r1[2] = 0.5*(uav-c)**2
        r2[0] = 1.0
        r2[1] = uav
        r2[2] = 0.5*uav**2
        r3[0] = 1.0
        r3[1] = uav+c
        r3[2] = 0.5*(uav+c)**2
        z1[j,:] = np.dot(l1,(f[j0,:] - f[jm1,:]))*r1
        z2[j,:] = np.dot(l2,(f[j0,:] - f[jm1,:]))*r2
        z3[j,:] = np.dot(l3,(f[j0,:] - f[jm1,:]))*r3

    for j in range(0,Nx+2):
        
        Ap[j,:] = 0
        An[j,:] = 0
        
        theta1[j] = 1
        theta2[j] = 1
        theta3[j] = 1
        
        jm1 = (j > 0)*(j-1) + (j == 0)*0
        jp1 = (j < Nx)*(j+1) + (j == Nx)*(Nx)
        
        z1sq = np.dot(z1[j,:],z1[j,:])
        if(s1[j]>0):
            Ap[j,:] = Ap[j,:] + z1[j,:]
            if(z1sq>0):
                theta1[j] = np.dot(z1[jm1,:],z1[j,:])/np.dot(z1[j,:],z1[j,:])        
                
        elif(s1[j]<0):
            An[j,:] = An[j,:] + z1[j,:]
            if(z1sq>0):
                theta1[j] = np.dot(z1[jp1,:],z1[j,:])/np.dot(z1[j,:],z1[j,:])
              
        else:
            Ap[j,:] = Ap[j,:] + 0.5*z1[j,:]
            An[j,:] = An[j,:] + 0.5*z1[j,:]
           
        z2sq = np.dot(z2[j,:],z2[j,:])
        if(s2[j]>0):
            Ap[j,:] = Ap[j,:] + z2[j,:]
            if(z2sq>0):
                theta2[j] = np.dot(z2[jm1,:],z2[j,:])/np.dot(z2[j,:],z2[j,:])
      
        elif(s2[j]<0):
            An[j,:] = An[j,:] + z2[j,:]
            if(z2sq>0):
                theta2[j] = np.dot(z2[jp1,:],z2[j,:])/np.dot(z2[j,:],z2[j,:])
                
        else:
            Ap[j,:] = Ap[j,:] + 0.5*z2[j,:]
            An[j,:] = An[j,:] + 0.5*z2[j,:]
           
        z3sq = np.dot(z3[j,:],z3[j,:])
        if(s3[j]>0):
            Ap[j,:] = Ap[j,:] + z3[j,:]
            if(z3sq>0):
                theta3[j] = np.dot(z3[jm1,:],z3[j,:])/np.dot(z3[j,:],z3[j,:])
          
        elif(s3[j]<0):
            An[j,:] = An[j,:] + z3[j,:]
            if(z3sq>0):
                theta3[j] = np.dot(z3[jp1,:],z3[j,:])/np.dot(z3[j,:],z3[j,:])
           
        else:
            Ap[j,:] = Ap[j,:] + 0.5*z3[j,:]
            An[j,:] = An[j,:] + 0.5*z3[j,:]
   
    for j in range(0,Nx+2):
        '''
        F[j,:] = np.sign(s1[j])*(1-((abs(s1[j])*dt)/dx))*z1[j,:]\
               + np.sign(s2[j])*(1-((abs(s2[j])*dt)/dx))*z2[j,:]\
               + np.sign(s3[j])*(1-((abs(s3[j])*dt)/dx))*z3[j,:]
        '''  
        F[j,:] = np.sign(s1[j])*(1-((abs(s1[j])*dt)/dx))*z1[j,:]*vanleer(theta1[j])\
               + np.sign(s2[j])*(1-((abs(s2[j])*dt)/dx))*z2[j,:]*vanleer(theta2[j])\
               + np.sign(s3[j])*(1-((abs(s3[j])*dt)/dx))*z3[j,:]*vanleer(theta3[j])
        F[j,:] = 0.5*F[j,:]
    
    for j in range(0,Nx+1):
        #Q[j,:] = Q[j,:] - ((dt/dx)*(An[j+1,:] + Ap[j,:])) - ((dt/dx)*(F[j+1,:]-F[j,:])) - (dt*e*trap1(gd[j,:],m1))
        Q[j,0] = Q[j,0] - ((dt/dx)*(An[j+1,0] + Ap[j,0])) - ((dt/dx)*(F[j+1,0]-F[j,0])) - (dt*e*trap1(gd[j,:],m1))
        Q[j,1] = Q[j,1] - ((dt/dx)*(An[j+1,1] + Ap[j,1])) - ((dt/dx)*(F[j+1,1]-F[j,1])) - (dt*e*trap1(gd[j,:],m2))
        Q[j,2] = Q[j,2] - ((dt/dx)*(An[j+1,2] + Ap[j,2])) - ((dt/dx)*(F[j+1,2]-F[j,2])) - (dt*e*trap1(gd[j,:],m3))
    
    rho = Q[:,0]
    u = np.divide(Q[:,1],Q[:,0])
    T = np.divide(((2*Q[:,2])-np.multiply(rho,np.square(u))),rho)
      
plt.figure(1)
plt.grid('True')
plt.plot(x,rho);
plt.title(r'Plot of rho vs x');
plt.xlabel(r'x');
plt.ylabel(r'rho');
  
plt.figure(2)
plt.grid('True')
plt.plot(x,u);
plt.title(r'Plot of u vs x');
plt.xlabel(r'x');
plt.ylabel(r'u');

plt.figure(3)
plt.grid('True')
plt.plot(x,T);
plt.title(r'Plot of T vs x');
plt.xlabel(r'x');
plt.ylabel(r'T');

x = np.append(x,x[-1]+dx)
plt.figure(4)
plt.clf()
plt.pcolor(x,v,np.transpose(g))
plt.colorbar()
plt.xlabel(r'x');
plt.ylabel(r'v');
plt.show(); 





