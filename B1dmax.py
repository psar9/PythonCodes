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
M = np.zeros((Nx+1,Nv))
G = np.zeros((Nx+1,Nv))
Gc = np.zeros((Nx+1,Nv))
Z = np.zeros((Nx+1,Nv))
Zc = np.zeros((Nx+1,Nv))
H = np.zeros(Nx+1)
A1 = np.zeros(Nx+1)
A2 = np.zeros(Nx+1)
A3 = np.zeros(Nx+1)
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
        jm1 = (j > 0)*(j-1) + (j == 0)*0
        jp1 = (j < Nx)*(j+1) + (j == Nx)*(Nx)
        for k in range(0,Nv):     
            M[j,k] = (rho[j]/((2*np.pi*T[j])**0.5))*np.exp((-(v[k]-u[j])**2)/(2*T[j]))
            Gc[j,k] = -(((v[k]-u[j])**3/(2*T[j]))-1.5*(v[k]-u[j]))*((T[jp1]-T[jm1])/(2*dx*T[j]))*M[j,k]
            Z[j,k] = (max(v[k],0)*((G[j,k] - G[jm1,k])/dx)) + (min(v[k],0)*((G[jp1,k] - G[j,k])/dx)) 
      
    for j in range(0,Nx+1):
        A1[j] = 0
        A2[j] = 0
        A3[j] = 0
        for k in range(0,Nv): 
            A1[j] = A1[j] + Z[j,k]
            A2[j] = A2[j] + (v[k] - u[j])*Z[j,k]
            A3[j] = A3[j] + ((((v[k] - u[j])**2)/(2*T[j])) - 0.5)*Z[j,k]   
  
    for j in range(0,Nx+1):
        for k in range(0,Nv): 
            Zc[j,k] = (dv/rho[j])*(A1[j] + (((v[k]-u[j])*A2[j])/T[j]) + (((((v[k]-u[j])**2)/(2*T[j]))-0.5)*2*A3[j]))*M[j,k]
    
    G = ((e/(e+dt))*(G - (dt*(Z - Zc)))) + ((dt/(e+dt))*Gc)
    
    for j in range(0,Nx+1):
        H[j] = 0
        for k in range(0,Nv):
            H[j] = H[j] + v[k]**3*G[j,k]
        
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
        jm1 = (j > 0)*(j-1) + (j == 0)*0
        jp1 = (j < Nx)*(j+1) + (j == Nx)*(Nx)
        #Q[j,:] = Q[j,:] - ((dt/dx)*(An[j+1,:] + Ap[j,:])) - ((dt/dx)*(F[j+1,:]-F[j,:])) - (dt*e*trap1(gd[j,:],m1))
        Q[j,0] = Q[j,0] - ((dt/dx)*(An[j+1,0] + Ap[j,0])) - ((dt/dx)*(F[j+1,0]-F[j,0])) 
        Q[j,1] = Q[j,1] - ((dt/dx)*(An[j+1,1] + Ap[j,1])) - ((dt/dx)*(F[j+1,1]-F[j,1]))
        Q[j,2] = Q[j,2] - ((dt/dx)*(An[j+1,2] + Ap[j,2])) - ((dt/dx)*(F[j+1,2]-F[j,2])) - (((dt*e*dv)/(4*dx))*(H[jp1]-H[jm1]))
    
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

plt.figure(4)
plt.clf()
plt.pcolor(x,v,np.transpose(G))
plt.colorbar()
plt.xlabel(r'x');
plt.ylabel(r'v');
plt.show(); 





