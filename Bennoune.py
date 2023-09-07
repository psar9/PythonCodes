# import critical libraries
import numpy as np;
import matplotlib;
import matplotlib.pyplot as plt

# set some plotting options
matplotlib.rcParams.update({'font.size':14, 'font.family':'serif'});
matplotlib.rcParams.update({'text.usetex':'true'});
matplotlib.rcParams.update(matplotlib.rcParamsDefault)

'''X = -7.5
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
V = -4.5#-4
Nv = 100#80
t = 0.16
dx = abs(1/Nx)
dv = abs((2*V)/(Nv-1))
#dv = abs((4-V)/(Nv-1))
n = 14
e = 2**-n

dt = abs(0.9*dx/V)
nt = int(np.ceil(t/dt))
dt = t/nt

x = np.zeros(Nx+1)
v = np.zeros(Nv)
rho = np.ones(Nx+1)
T = np.ones(Nx+1)
u = np.ones(Nx+1)
U1 = np.ones(Nx+1)
U2 = np.ones(Nx+1)
U3 = np.ones(Nx+1)
Fp1 = np.zeros(Nx+1)
Fp2 = np.zeros(Nx+1)
Fp3 = np.zeros(Nx+1)
Fn1 = np.zeros(Nx+1)
Fn2 = np.zeros(Nx+1)
Fn3 = np.zeros(Nx+1)
#Nvn = int(np.ceil(3/dv))
#Nvp = int(np.ceil(4/dv))
#vn = np.zeros(Nvn)
#vp = np.zeros(Nvp)
f = np.zeros((Nx+1,Nv))
dM = np.zeros((Nx+1,Nv))
M = np.zeros((Nx+1,Nv))
M0 = np.zeros((Nx+1,Nv))
M1 = np.zeros((Nx+1,Nv))
M2 = np.zeros((Nx+1,Nv))
g = np.zeros((Nx,Nv))
phi1 = np.zeros((Nx,Nv))
phi2 = np.zeros((Nx,Nv))
gd = np.zeros((Nx+1,Nv))
a1 = np.zeros((Nx,Nv))
b1 = np.zeros((Nx,Nv))
a2 = np.zeros((Nx,Nv))
b2 = np.zeros((Nx,Nv))
pi1 = np.zeros((Nx,Nv))
pi2 = np.zeros((Nx,Nv))
Fp = np.zeros((Nx+1,Nv))
Fn = np.zeros((Nx+1,Nv))

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
    U1[j] = rho[j]
    U2[j] = rho[j]*u[j]
    U3[j] = (0.5*rho[j]*u[j]**2) + (0.5*rho[j]*T[j])     
 
for k in range(0,Nv):
    v[k] = V + k*dv   

'''for k in range(0,Nvn):
    vn[k] = V + k*dv   
    
for k in range(0,Nvp):
    vp[k] = 0 + k*dv 
   
v = np.array(list(vn)+list(vp))'''
       
def trap(v,F,phi,dv):
    sum = np.zeros(len(F))
    sum = sum + 0.5*F[:,0]*phi(v[0]) + 0.5*F[:,-1]*phi(v[-1])
    for k in range(1,len(v)-1):
        sum = sum + F[:,k]*phi(v[k])          
    return dv*sum  

def m1(y):
    return 1

def m2(y):
    return y

def m3(y):
    return 0.5*y**2

'''def pi(v,u,phi,M,T,rho):
    return (1/rho) * (trap(phi) + (((v-u)*(trap((v-u)*phi)))/T) + \
                      (((((v - u)**2)/(2*T)) - 0.5)*2*trap(((((v - u)**2)/(2*T)) - (1/2))*phi)))*M'''

for i in range(0,nt):
    
    for j in range(0,Nx):
        
        rhom = (rho[j]+rho[j+1])/2
        um = (u[j]+u[j+1])/2
        Tm = (T[j]+T[j+1])/2
              
        for k in range(0,Nv):
            
            M[j,k] = (rhom/((2*np.pi*Tm)**0.5))*np.exp((-(v[k]-um)**2)/(2*Tm))
            dM[j,k] = (rho[j+1]/((2*np.pi*T[j+1])**0.5))*np.exp((-(v[k]-u[j+1])**2)/(2*T[j+1]))- (rho[j]/((2*np.pi*T[j])**0.5))*np.exp((-(v[k]-u[j])**2)/(2*T[j]))
            phi2[j,k] = v[k]*(dM[j,k]/dx)
            
            if(j == 0):
                phi1[j,k] = min(v[k],0)*((g[j+1,k] - g[j,k])/dx)
                
            elif(j == Nx-1): 
                phi1[j,k] = max(v[k],0)*((g[j,k] - g[j-1,k])/dx)
                
            else:
                phi1[j,k] = (max(v[k],0)*((g[j,k] - g[j-1,k])/dx)) + (min(v[k],0)*((g[j+1,k] - g[j,k])/dx)) 
            
            a1[j,k] = (v[k] - u[j])*phi1[j,k]
            a2[j,k] = (v[k] - u[j])*phi2[j,k]
            b1[j,k] = ((((v[k] - u[j])**2)/(2*T[j])) - 0.5)*phi1[j,k]      
            b2[j,k] = ((((v[k] - u[j])**2)/(2*T[j])) - 0.5)*phi2[j,k]
            
    c1 = trap(v,a1,m1,dv)
    c2 = trap(v,a2,m1,dv)
    d1 = trap(v,b1,m1,dv)
    d2 = trap(v,b2,m1,dv)
    
    p1 = trap(v,phi1,m1,dv)
    p2 = trap(v,phi2,m1,dv)
    
    for j in range(0,Nx):
          
        for k in range(0,Nv):
            
            pi1[j,k] = (1/rho[j])*(p1[j] + (((v[k]-u[j])*c1[j])/T[j]) + (((((v[k]-u[j])**2)/(2*T[j]))-0.5)*2*d1[j]))*M[j,k]
            pi2[j,k] = (1/rho[j])*(p2[j] + (((v[k]-u[j])*c2[j])/T[j]) + (((((v[k]-u[j])**2)/(2*T[j]))-0.5)*2*d2[j]))*M[j,k]
 
    g = (e/(e+dt))*(g - (dt*(phi1 - pi1))) - ((dt/(e+dt))*(phi2 - pi2)) 
    
    for j in range(0,Nx+1):
          
        for k in range(0,Nv):
            
            if(j == 0):
                
                M0[j,k] = (rho[j]/((2*np.pi*T[j])**0.5))*np.exp((-(v[k]-u[j])**2)/(2*T[j]))               
                M2[j,k] = (rho[j+1]/((2*np.pi*T[j+1])**0.5))*np.exp((-(v[k]-u[j+1])**2)/(2*T[j+1]))
                gd[j,k] = 0               
              
            elif(j == Nx): 
                
                M0[j,k] = (rho[j-1]/((2*np.pi*T[j-1])**0.5))*np.exp((-(v[k]-u[j-1])**2)/(2*T[j-1]))             
                M2[j,k] = (rho[j]/((2*np.pi*T[j])**0.5))*np.exp((-(v[k]-u[j])**2)/(2*T[j]))         
                gd[j,k] = 0
                    
            else:
                    
                M0[j,k] = (rho[j-1]/((2*np.pi*T[j-1])**0.5))*np.exp((-(v[k]-u[j-1])**2)/(2*T[j-1])) 
                M2[j,k] = (rho[j+1]/((2*np.pi*T[j+1])**0.5))*np.exp((-(v[k]-u[j+1])**2)/(2*T[j+1]))                              
                gd[j,k] = v[k]*((g[j,k] - g[j-1,k])/dx)
            
            M1[j,k] = (rho[j]/((2*np.pi*T[j])**0.5))*np.exp((-(v[k]-u[j])**2)/(2*T[j])) 
            Fp[j,k] = max(v[k],0)*M1[j,k] + min(v[k],0)*M2[j,k]
            Fn[j,k] = max(v[k],0)*M0[j,k] + min(v[k],0)*M1[j,k]
            
        
    U1 = U1 - ((dt/dx)*(trap(v,Fp,m1,dv) - trap(v,Fn,m1,dv))) - (dt*e*trap(v,gd,m1,dv))
    U2 = U2 - ((dt/dx)*(trap(v,Fp,m2,dv) - trap(v,Fn,m2,dv))) - (dt*e*trap(v,gd,m2,dv))
    U3 = U3 - ((dt/dx)*(trap(v,Fp,m3,dv) - trap(v,Fn,m3,dv))) - (dt*e*trap(v,gd,m3,dv))
    
    '''U1 = U1 - ((dt/dx)*(trap(v,Fp,m1,dv) - trap(v,Fn,m1,dv)))
    U2 = U2 - ((dt/dx)*(trap(v,Fp,m2,dv) - trap(v,Fn,m2,dv))) 
    U3 = U3 - ((dt/dx)*(trap(v,Fp,m3,dv) - trap(v,Fn,m3,dv)))'''
            
    '''f1 = U2
    f2 = 2*U3
    f3 = (np.divide(U2,U1))*((3*U3) - np.divide(np.square(U2),U1))
    
    for j in range(0,Nx+1):
        
        if(j == 0):
            
            F1[j] = (0.5*(f1[j]+f1[j+1])) - (0.5*abs(V)*(U1[j+1]-U1[j])) - f1[j]
            F2[j] = (0.5*(f2[j]+f2[j+1])) - (0.5*abs(V)*(U2[j+1]-U2[j])) - f2[j]
            F3[j] = (0.5*(f3[j]+f3[j+1])) - (0.5*abs(V)*(U3[j+1]-U3[j])) - f3[j]
          
        elif(j == Nx): 
            
            F1[j] = f1[j]- (0.5*(f1[j-1]+f1[j])) + (0.5*abs(V)*(U1[j]-U1[j-1]))
            F2[j] = f2[j]- (0.5*(f2[j-1]+f2[j])) + (0.5*abs(V)*(U2[j]-U2[j-1]))
            F3[j] = f3[j]- (0.5*(f3[j-1]+f3[j])) + (0.5*abs(V)*(U3[j]-U3[j-1]))         
                
        else:
                
            F1[j] = (0.5*(f1[j]+f1[j+1])) - (0.5*abs(V)*(U1[j+1]-U1[j])) - (0.5*(f1[j-1]+f1[j])) + (0.5*abs(V)*(U1[j]-U1[j-1]))
            F2[j] = (0.5*(f2[j]+f2[j+1])) - (0.5*abs(V)*(U2[j+1]-U2[j])) - (0.5*(f2[j-1]+f2[j])) + (0.5*abs(V)*(U2[j]-U2[j-1]))
            F3[j] = (0.5*(f3[j]+f3[j+1])) - (0.5*abs(V)*(U3[j+1]-U3[j])) - (0.5*(f3[j-1]+f3[j])) + (0.5*abs(V)*(U3[j]-U3[j-1]))                           
          
    U1 = U1 - (dt/dx)*F1
    U2 = U2 - (dt/dx)*F2
    U3 = U3 - (dt/dx)*F3'''
    
    rho = U1
    u = np.divide(U2,U1)
    T = np.divide(((2*U3)-np.multiply(rho,np.square(u))),rho)
    
plt.figure(1)
plt.clf()
plt.grid()
plt.plot(x,rho);
#plt.plot(x,rho,label = 'Numerical');
plt.title(r'Plot of rho vs x');
plt.xlabel(r'x');
plt.ylabel(r'rho');
#plt.legend();
plt.show();  

plt.figure(2)
plt.clf()
plt.grid()
plt.plot(x,u);
plt.title(r'Plot of u vs x');
plt.xlabel(r'x');
plt.ylabel(r'u');
plt.show(); 

plt.figure(3)
plt.clf()
plt.grid()
plt.plot(x,T,label = 'Numerical');
plt.title(r'Plot of T vs x');
plt.xlabel(r'x');
plt.ylabel(r'T');
plt.show(); 

plt.figure(4)
plt.clf()
plt.pcolor(np.transpose(g))
plt.colorbar()
plt.xlabel(r'x');
plt.ylabel(r'v');
plt.show(); 





