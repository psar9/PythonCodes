# import critical libraries
import numpy as np;
import matplotlib;
import matplotlib.pyplot as plt

# set some plotting options
matplotlib.rcParams.update({'font.size':14, 'font.family':'serif'});
matplotlib.rcParams.update({'text.usetex':'true'});
matplotlib.rcParams.update(matplotlib.rcParamsDefault)

X = 0
Y = 0
Nx = 20
Ny = 20
Vx = -4.5
Vy = -4.5
Nvx = 10
Nvy = 10
t = 0.16
dx = abs(1/Nx)
dy = abs(1/Ny)
dvx = abs((2*Vx)/(Nvx-1))
dvy = abs((2*Vy)/(Nvy-1))
dtx = abs(0.9*dx/Vx)
dty = abs(0.9*dy/Vy)
nt = int(np.ceil(t/dtx))
dt = t/nt

def trap2d(vx,vy,F,phix,phiy):
    
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

def trap2(vx,vy,F,phix,phiy):
    
    sum = 0
    
    for l in range(1,len(vx)-1):
        for m in range(1,len(vy)-1):
            sum = sum + 4*F[l,m]*phix(vx[l])*phiy(vy[m])
    
    for l in range(1,len(vx)-1):
        sum = sum + 2*(F[l,0]*phix(vx[l])*phiy(vy[0])+F[l,-1]*phix(vx[l])*phiy(vy[-1]))
        
    for m in range(1,len(vy)-1):
        sum = sum + 2*(F[0,m]*phix(vx[0])*phiy(vy[m])+F[-1,m]*phix(vx[-1])*phiy(vy[m]))
    
    sum = sum + F[0,0]*phix(vx[0])*phiy(vy[0]) + F[0,-1]*phix(vx[0])*phiy(vy[-1])\
    + F[-1,0]*phix(vx[-1])*phiy(vy[0]) + F[-1,-1]*phix(vx[-1])*phiy(vy[-1])
    
    return dvx*dvy*sum/4

def m1(y):
    return 1

def m2(y):
    return y

def m3(y):
    return 0.5*y**2
    
def vanleer(r):
    psi = (r+abs(r))/(1+abs(r))
    return psi

n = 14
    
e = 2**-n
x = np.zeros(Nx+1)
y = np.zeros(Ny+1)
vx = np.zeros(Nvx)
vy = np.zeros(Nvy)
rho = np.ones((Nx+1,Ny+1))
T = np.ones((Nx+1,Ny+1))
ux = np.ones((Nx+1,Ny+1))
uy = np.ones((Nx+1,Ny+1))
M = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
g = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
phi1 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
phi2 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
gdx = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
gdy = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
ax1 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
ay1 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
ax2 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
ay2 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
b1 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
b2 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
pi1 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
pi2 = np.zeros((Nx+1,Ny+1,Nvx,Nvy))
fx = np.zeros((Nx+1,Ny+1,4))
fy = np.zeros((Nx+1,Ny+1,4))
s1x = np.zeros((Nx+2,Ny+2))
s2x = np.zeros((Nx+2,Ny+2))
s3x = np.zeros((Nx+2,Ny+2))
s4x = np.zeros((Nx+2,Ny+2))
s1y = np.zeros((Nx+2,Ny+2))
s2y = np.zeros((Nx+2,Ny+2))
s3y = np.zeros((Nx+2,Ny+2))
s4y = np.zeros((Nx+2,Ny+2))
z1x = np.zeros((Nx+2,Ny+2,4))
z2x = np.zeros((Nx+2,Ny+2,4))
z3x = np.zeros((Nx+2,Ny+2,4))
z4x = np.zeros((Nx+2,Ny+2,4))
z1y = np.zeros((Nx+2,Ny+2,4))
z2y = np.zeros((Nx+2,Ny+2,4))
z3y = np.zeros((Nx+2,Ny+2,4))
z4y = np.zeros((Nx+2,Ny+2,4))
theta1x = np.ones((Nx+2,Ny+2))
theta2x = np.ones((Nx+2,Ny+2))
theta3x = np.ones((Nx+2,Ny+2))
theta4x = np.ones((Nx+2,Ny+2))
theta1y = np.ones((Nx+2,Ny+2))
theta2y = np.ones((Nx+2,Ny+2))
theta3y = np.ones((Nx+2,Ny+2))
theta4y = np.ones((Nx+2,Ny+2))
Apx = np.zeros((Nx+2,Ny+2,4))
Anx = np.zeros((Nx+2,Ny+2,4))
Apy = np.zeros((Nx+2,Ny+2,4))
Any = np.zeros((Nx+2,Ny+2,4))
Fx = np.zeros((Nx+2,Ny+2,4))
Fy = np.zeros((Nx+2,Ny+2,4))
Q1 = np.zeros((Nx+1,Ny+1))
Q2 = np.zeros((Nx+1,Ny+1))
Q3 = np.zeros((Nx+1,Ny+1))
Q4 = np.zeros((Nx+1,Ny+1))
    
for j in range(0,Nx+1):
    x[j] = X + j*dx
    if(x[j] <= 0.5):
        rho[j,:] = 1
        ux[j,:] = 0
        uy[j,:] = 0
        T[j,:] = 1     
    else:
        rho[j,:] = 0.125
        ux[j,:] = 0
        uy[j,:] = 0
        T[j,:] = 0.8
        
for k in range(0,Ny+1):
    y[k] = Y + k*dy
    '''
    if(y[k] <= 0.5):
        rho[:,k] = 1
        ux[:,k] = 0
        uy[:,k] = 0
        T[:,k] = 1     
    else:
        rho[:,k] = 0.125
        ux[:,k] = 0
        uy[:,k] = 0
        T[:,k] = 0.8'''

for j in range(0,Nx+1):
    for k in range(0,Ny+1):            
        Q1[j,k] = rho[j,k]
        Q2[j,k] = rho[j,k]*ux[j,k]
        Q3[j,k] = rho[j,k]*uy[j,k]
        Q4[j,k] = (0.5*rho[j,k]*(ux[j,k]**2+uy[j,k]**2)) + (rho[j,k]*T[j,k])
        
for l in range(0,Nvx):
    vx[l] = Vx + l*dvx

for m in range(0,Nvy):
    vy[m] = Vy + m*dvy
      
for i in range(0,nt):
    
    for j in range(0,Nx+1):
        for k in range(0,Ny+1):
            for l in range(0,Nvx):
                for m in range(0,Nvy):
                        M[j,k,l,m] = (rho[j,k]/(2*np.pi*T[j,k]))*np.exp((-((vx[l]-ux[j,k])**2+(vy[m]-uy[j,k])**2))/(2*T[j,k]))

    for j in range(0,Nx+1):
        for k in range(0,Ny+1):
            jm1 = (j > 0)*(j-1) + (j == 0)*0
            jp1 = (j < Nx)*(j+1) + (j == Nx)*(Nx)
            km1 = (k > 0)*(k-1) + (k == 0)*0
            kp1 = (k < Ny)*(k+1) + (k == Ny)*(Ny)
            for l in range(0,Nvx):
                for m in range(0,Nvy): 
                    phi1[j,k,l,m] = (max(vx[l],0)*((g[j,k,l,m] - g[jm1,k,l,m])/dx)) + (min(vx[l],0)*((g[jp1,k,l,m] - g[j,k,l,m])/dx))\
                                  + (max(vy[m],0)*((g[j,k,l,m] - g[j,km1,l,m])/dy)) + (min(vy[m],0)*((g[j,kp1,l,m] - g[j,k,l,m])/dy))
                    phi2[j,k,l,m] = vx[l]*((M[j,k,l,m] - M[jm1,k,l,m])/dx) + vy[m]*((M[j,k,l,m] - M[j,km1,l,m])/dy)
    
    for j in range(0,Nx+1):
        for k in range(0,Ny+1):
            for l in range(0,Nvx): 
                for m in range(0,Nvy): 
                    ax1[j,k,l,m] = (vx[l] - ux[j,k])*phi1[j,k,l,m]
                    ax2[j,k,l,m] = (vx[l] - ux[j,k])*phi2[j,k,l,m]
                    ay1[j,k,l,m] = (vy[m] - uy[j,k])*phi1[j,k,l,m]
                    ay2[j,k,l,m] = (vy[m] - uy[j,k])*phi2[j,k,l,m]
                    b1[j,k,l,m] = ((((vx[l] - ux[j,k])**2 + (vy[m] - uy[j,k])**2)/(2*T[j,k])) - 1)*phi1[j,k,l,m]      
                    b2[j,k,l,m] = ((((vx[l] - ux[j,k])**2 + (vy[m] - uy[j,k])**2)/(2*T[j,k])) - 1)*phi2[j,k,l,m]
    
    A11 = trap2d(vx,vy,phi1,m1,m1)
    A12 = trap2d(vx,vy,phi2,m1,m1)
    A21 = trap2d(vx,vy,ax1,m1,m1)
    A22 = trap2d(vx,vy,ax2,m1,m1)
    A31 = trap2d(vx,vy,ay1,m1,m1)
    A32 = trap2d(vx,vy,ay2,m1,m1)
    A41 = trap2d(vx,vy,b1,m1,m1)
    A42 = trap2d(vx,vy,b2,m1,m1)
    
    for j in range(0,Nx+1):
        for k in range(0,Ny+1):
            for l in range(0,Nvx): 
                for m in range(0,Nvy): 
                    pi1[j,k,l,m] = (1/rho[j,k])*(A11[j,k] + (((vx[l]-ux[j,k])*A21[j,k])/T[j,k]) + (((vy[m]-uy[j,k])*A31[j,k])/T[j,k])\
                                 + (((((vx[l]-ux[j,k])**2 + (vy[m]-uy[j,k])**2)/(2*T[j,k]))-1)*A41[j,k]))*M[j,k,l,m]
                    pi2[j,k,l,m] = (1/rho[j,k])*(A12[j,k] + (((vx[l]-ux[j,k])*A22[j,k])/T[j,k]) + (((vy[m]-uy[j,k])*A32[j,k])/T[j,k])\
                                 + (((((vx[l]-ux[j,k])**2 + (vy[m]-uy[j,k])**2)/(2*T[j,k]))-1)*A42[j,k]))*M[j,k,l,m]
 
    g = ((e/(e+dt))*(g - (dt*(phi1 - pi1)))) - ((dt/(e+dt))*(phi2 - pi2))
    
    for j in range(0,Nx+1):
        for k in range(0,Ny+1):
            for l in range(0,Nvx):
                for m in range(0,Nvy):
                    jm1 = (j > 0)*(j-1) + (j == 0)*0
                    jp1 = (j < Nx)*(j+1) + (j == Nx)*(Nx)
                    km1 = (k > 0)*(k-1) + (k == 0)*0
                    kp1 = (k < Ny)*(k+1) + (k == Ny)*(Ny)
                    gdx[j,k,l,m] = vx[l]*((g[jp1,k,l,m] - g[jm1,k,l,m])/(2*dx))
                    gdy[j,k,l,m] = vy[m]*((g[j,kp1,l,m] - g[j,km1,l,m])/(2*dy))

    for j in range(0,Nx+1):
        for k in range(0,Ny+1):
            fx[j,k,0] = rho[j,k]*ux[j,k]
            fx[j,k,1] = rho[j,k]*(ux[j,k]**2 + T[j,k])
            fx[j,k,2] = rho[j,k]*ux[j,k]*uy[j,k]
            fx[j,k,3] = (0.5*(ux[j,k]**2+uy[j,k]**2) + 2*T[j,k])*ux[j,k]*rho[j,k]

    for j in range(0,Nx+2):
        for k in range(0,Ny+2): 
            jm1 = (j > 0)*(j-1) + (j == 0)*0
            j0 = (j < Nx+1)*j + (j == Nx+1)*(Nx)
            k0 = (k < Ny+1)*k + (k == Ny+1)*(Ny)
            Tav = 0.5*(T[jm1,k0]+T[j0,k0])
            uxav = 0.5*(ux[jm1,k0]+ux[j0,k0])
            uyav = 0.5*(uy[jm1,k0]+uy[j0,k0])
            rhoav = 0.5*(rho[jm1,k0]+rho[j0,k0])
            Eav = 0.5*rhoav*(uxav**2+uyav**2) + (rhoav*Tav)
            cx = np.sqrt(2*Tav)
            Hx = (Eav+ rhoav*Tav)/rhoav
            s1x[j,k] = uxav - cx
            s2x[j,k] = uxav
            s3x[j,k] = uxav
            s4x[j,k] = uxav + cx
            r1x = np.array([1.0, uxav-cx, uyav, Hx - uxav*cx])
            r2x = np.array([1.0, uxav, uyav, 0.5*(uxav**2+uyav**2)])
            r3x = np.array([0, 0, 1, uyav])
            r4x = np.array([1.0, uxav + cx, uyav, Hx + uxav*cx])
            lex = np.linalg.inv(np.transpose(np.array([r1x,r2x,r3x,r4x])))
            l1x = lex[0]
            l2x = lex[1]
            l3x = lex[2]
            l4x = lex[3]
            z1x[j,k,:] = np.dot(l1x,(fx[j0,k0,:] - fx[jm1,k0,:]))*r1x
            z2x[j,k,:] = np.dot(l2x,(fx[j0,k0,:] - fx[jm1,k0,:]))*r2x
            z3x[j,k,:] = np.dot(l3x,(fx[j0,k0,:] - fx[jm1,k0,:]))*r3x
            z4x[j,k,:] = np.dot(l4x,(fx[j0,k0,:] - fx[jm1,k0,:]))*r4x  
            
    for j in range(0,Nx+2):
        for k in range(0,Ny+2):
            
            Apx[j,k,:] = 0
            Anx[j,k,:] = 0
            
            theta1x[j,k] = 1
            theta2x[j,k] = 1
            theta3x[j,k] = 1 
            
            jm1 = (j > 0)*(j-1) + (j == 0)*0
            jp1 = (j < Nx)*(j+1) + (j == Nx)*(Nx)
            
            z1xsq = np.dot(z1x[j,k,:],z1x[j,k,:])
            if(s1x[j,k]>0):
                Apx[j,k,:] = Apx[j,k,:] + z1x[j,k,:]
                if(z1xsq>0):
                    theta1x[j,k] = np.dot(z1x[jm1,k,:],z1x[j,k,:])/np.dot(z1x[j,k,:],z1x[j,k,:])        
                    
            elif(s1x[j,k]<0):
                Anx[j,k,:] = Anx[j,k,:] + z1x[j,k,:]
                if(z1xsq>0):
                    theta1x[j,k] = np.dot(z1x[jp1,k,:],z1x[j,k,:])/np.dot(z1x[j,k,:],z1x[j,k,:])
                  
            else:
                Apx[j,k,:] = Apx[j,k,:] + 0.5*z1x[j,k,:]
                Anx[j,k,:] = Anx[j,k,:] + 0.5*z1x[j,k,:]
               
            z2xsq = np.dot(z2x[j,k,:],z2x[j,k,:])
            if(s2x[j,k]>0):
                Apx[j,k,:] = Apx[j,k,:] + z2x[j,k,:]
                if(z2xsq>0):
                    theta2x[j,k] = np.dot(z2x[jm1,k,:],z2x[j,k,:])/np.dot(z2x[j,k,:],z2x[j,k,:])
          
            elif(s2x[j,k]<0):
                Anx[j,k,:] = Anx[j,k,:] + z2x[j,k,:]
                if(z2xsq>0):
                    theta2x[j,k] = np.dot(z2x[jp1,k,:],z2x[j,k,:])/np.dot(z2x[j,k,:],z2x[j,k,:])
                    
            else:
                Apx[j,k,:] = Apx[j,k,:] + 0.5*z2x[j,k,:]
                Anx[j,k,:] = Anx[j,k,:] + 0.5*z2x[j,k,:]
               
            z3xsq = np.dot(z3x[j,k,:],z3x[j,k,:])
            if(s3x[j,k]>0):
                Apx[j,k,:] = Apx[j,k,:] + z3x[j,k,:]
                if(z3xsq>0):
                    theta3x[j,k] = np.dot(z3x[jm1,k,:],z3x[j,k,:])/np.dot(z3x[j,k,:],z3x[j,k,:])
              
            elif(s3x[j,k]<0):
                Anx[j,k,:] = Anx[j,k,:] + z3x[j,k,:]
                if(z3xsq>0):
                    theta3x[j,k] = np.dot(z3x[jp1,k,:],z3x[j,k,:])/np.dot(z3x[j,k,:],z3x[j,k,:])   
            else:
                Apx[j,k,:] = Apx[j,k,:] + 0.5*z3x[j,k,:]
                Anx[j,k,:] = Anx[j,k,:] + 0.5*z3x[j,k,:]
                
            z4xsq = np.dot(z4x[j,k,:],z4x[j,k,:])
            if(s4x[j,k]>0):
                Apx[j,k,:] = Apx[j,k,:] + z4x[j,k,:]
                if(z4xsq>0):
                    theta4x[j,k] = np.dot(z4x[jm1,k,:],z4x[j,k,:])/np.dot(z4x[j,k,:],z4x[j,k,:])
              
            elif(s4x[j,k]<0):
                Anx[j,k,:] = Anx[j,k,:] + z4x[j,k,:]
                if(z4xsq>0):
                    theta4x[j,k] = np.dot(z4x[jp1,k,:],z4x[j,k,:])/np.dot(z4x[j,k,:],z4x[j,k,:])   
                    
            else:
                Apx[j,k,:] = Apx[j,k,:] + 0.5*z4x[j,k,:]
                Anx[j,k,:] = Anx[j,k,:] + 0.5*z4x[j,k,:]
   
    for j in range(0,Nx+2):
        for k in range(0,Ny+2):
     
            Fx[j,k,:] = np.sign(s1x[j,k])*(1-((abs(s1x[j,k])*dt)/dx))*z1x[j,k,:]*vanleer(theta1x[j,k])\
                      + np.sign(s2x[j,k])*(1-((abs(s2x[j,k])*dt)/dx))*z2x[j,k,:]*vanleer(theta2x[j,k])\
                      + np.sign(s3x[j,k])*(1-((abs(s3x[j,k])*dt)/dx))*z3x[j,k,:]*vanleer(theta3x[j,k])\
                      + np.sign(s4x[j,k])*(1-((abs(s4x[j,k])*dt)/dx))*z4x[j,k,:]*vanleer(theta4x[j,k])
            Fx[j,k,:] = 0.5*Fx[j,k,:]
           
    for j in range(0,Nx+1):
        for k in range(0,Ny+1):
            Q1[j,k] = Q1[j,k] - ((dt/dx)*(Anx[j+1,k,0] + Apx[j,k,0])) - ((dt/dx)*(Fx[j+1,k,0]-Fx[j,k,0])) 
            Q2[j,k] = Q2[j,k] - ((dt/dx)*(Anx[j+1,k,1] + Apx[j,k,1])) - ((dt/dx)*(Fx[j+1,k,1]-Fx[j,k,1])) 
            Q3[j,k] = Q3[j,k] - ((dt/dx)*(Anx[j+1,k,2] + Apx[j,k,2])) - ((dt/dx)*(Fx[j+1,k,2]-Fx[j,k,2])) 
            Q4[j,k] = Q4[j,k] - ((dt/dx)*(Anx[j+1,k,3] + Apx[j,k,3])) - ((dt/dx)*(Fx[j+1,k,3]-Fx[j,k,3]))
    
    rho = Q1
    ux = np.divide(Q2,Q1)
    uy = np.divide(Q3,Q1)
    T = np.divide((Q4 - 0.5*np.multiply(rho,np.square(ux)+np.square(uy))),rho)
    
    for j in range(0,Nx+1):
        for k in range(0,Ny+1):
            fy[j,k,0] = rho[j,k]*uy[j,k]
            fy[j,k,1] = rho[j,k]*ux[j,k]*uy[j,k]
            fy[j,k,2] = rho[j,k]*(uy[j,k]**2 + T[j,k])
            fy[j,k,3] = (0.5*(ux[j,k]**2+uy[j,k]**2) + 2*T[j,k])*uy[j,k]*rho[j,k]
    
    for j in range(0,Nx+2):  
        for k in range(0,Ny+2):
            j0 = (j < Nx+1)*j + (j == Nx+1)*(Nx)
            km1 = (k > 0)*(k-1) + (k == 0)*0
            k0 = (k < Ny+1)*k + (k == Ny+1)*(Ny)
            Tavy = 0.5*(T[j0,km1]+T[j0,k0])
            uavx = 0.5*(ux[j0,km1]+ux[j0,k0])
            uavy = 0.5*(uy[j0,km1]+uy[j0,k0])
            rhoy = 0.5*(rho[j0,km1]+rho[j0,k0])
            Eavy = 0.5*rhoy*(uavx**2+uavy**2) + (rhoy*Tavy)
            cy = np.sqrt(2*Tavy)
            Hy = (Eavy + rhoy*Tavy)/rhoy
            s1y[j,k] = uavy - cy
            s2y[j,k] = uavy
            s3y[j,k] = uavy
            s4y[j,k] = uavy + cy
            r1y = np.array([1.0, uavx, uavy-cy, Hy - uavy*cy])
            r2y = np.array([1.0, uavx, uavy, 0.5*(uavy**2+uavx**2)])
            r3y = np.array([0, 1, 0, uavx])
            r4y = np.array([1.0, uavx, uavy+cy, Hy + uavy*cy])
            ley = np.linalg.inv(np.transpose(np.array([r1y,r2y,r3y,r4y])))
            l1y = ley[0]
            l2y = ley[1]
            l3y = ley[2]
            l4y = ley[3]
            z1y[j,k,:] = np.dot(l1y,(fy[j0,k0,:] - fy[j0,km1,:]))*r1y
            z2y[j,k,:] = np.dot(l2y,(fy[j0,k0,:] - fy[j0,km1,:]))*r2y
            z3y[j,k,:] = np.dot(l3y,(fy[j0,k0,:] - fy[j0,km1,:]))*r3y
            z4y[j,k,:] = np.dot(l4y,(fy[j0,k0,:] - fy[j0,km1,:]))*r4y
    
    for j in range(0,Nx+2):
        for k in range(0,Ny+2):
                
            Apy[j,k,:] = 0
            Any[j,k,:] = 0
            
            theta1y[j,k] = 1
            theta2y[j,k] = 1
            theta3y[j,k] = 1
            
            km1 = (k > 0)*(k-1) + (k == 0)*0
            kp1 = (k < Ny)*(k+1) + (k == Ny)*(Ny)
            
            z1ysq = np.dot(z1y[j,k,:],z1y[j,k,:])
            if(s1y[j,k]>0):
                Apy[j,k,:] = Apy[j,k,:] + z1y[j,k,:]
                if(z1ysq>0):
                    theta1y[j,k] = np.dot(z1y[j,km1,:],z1y[j,k,:])/np.dot(z1y[j,k,:],z1y[j,k,:])        
                    
            elif(s1y[j,k]<0):
                Any[j,k,:] = Any[j,k,:] + z1y[j,k,:]
                if(z1ysq>0):
                    theta1y[j,k] = np.dot(z1y[j,kp1,:],z1y[j,k,:])/np.dot(z1y[j,k,:],z1y[j,k,:])
                  
            else:
                Apy[j,k,:] = Apy[j,k,:] + 0.5*z1y[j,k,:]
                Any[j,k,:] = Any[j,k,:] + 0.5*z1y[j,k,:]
               
            z2ysq = np.dot(z2y[j,k,:],z2y[j,k,:])
            if(s2y[j,k]>0):
                Apy[j,k,:] = Apy[j,k,:] + z2y[j,k,:]
                if(z2ysq>0):
                    theta2y[j,k] = np.dot(z2y[j,km1,:],z2y[j,k,:])/np.dot(z2y[j,k,:],z2y[j,k,:])
          
            elif(s2y[j,k]<0):
                Any[j,k,:] = Any[j,k,:] + z2y[j,k,:]
                if(z2ysq>0):
                    theta2y[j,k] = np.dot(z2y[j,kp1,:],z2y[j,k,:])/np.dot(z2y[j,k,:],z2y[j,k,:])
                    
            else:
                Apy[j,k,:] = Apy[j,k,:] + 0.5*z2y[j,k,:]
                Any[j,k,:] = Any[j,k,:] + 0.5*z2y[j,k,:]
               
            z3ysq = np.dot(z3y[j,k,:],z3y[j,k,:])
            if(s3y[j,k]>0):
                Apy[j,k,:] = Apy[j,k,:] + z3y[j,k,:]
                if(z3ysq>0):
                    theta3y[j,k] = np.dot(z3y[j,km1,:],z3y[j,k,:])/np.dot(z3y[j,k,:],z3y[j,k,:])
              
            elif(s3y[j,k]<0):
                Any[j,k,:] = Any[j,k,:] + z3y[j,k,:]
                if(z3ysq>0):
                    theta3y[j,k] = np.dot(z3y[j,kp1,:],z3y[j,k,:])/np.dot(z3y[j,k,:],z3y[j,k,:])   
            else:
                Apy[j,k,:] = Apy[j,k,:] + 0.5*z3y[j,k,:]
                Any[j,k,:] = Any[j,k,:] + 0.5*z3y[j,k,:]
                
            z4ysq = np.dot(z4y[j,k,:],z4y[j,k,:])
            if(s4y[j,k]>0):
                Apy[j,k,:] = Apy[j,k,:] + z4y[j,k,:]
                if(z4ysq>0):
                    theta4y[j,k] = np.dot(z4y[j,km1,:],z4y[j,k,:])/np.dot(z4y[j,k,:],z4y[j,k,:])
              
            elif(s4y[j,k]<0):
                Any[j,k,:] = Any[j,k,:] + z4y[j,k,:]
                if(z4ysq>0):
                    theta4y[j,k] = np.dot(z4y[j,kp1,:],z4y[j,k,:])/np.dot(z4y[j,k,:],z4y[j,k,:])   
                    
            else:
                Apy[j,k,:] = Apy[j,k,:] + 0.5*z4y[j,k,:]
                Any[j,k,:] = Any[j,k,:] + 0.5*z4y[j,k,:]
    
    for j in range(0,Nx+2):
        for k in range(0,Ny+2):
            
            Fy[j,k,:] = np.sign(s1y[j,k])*(1-((abs(s1y[j,k])*dt)/dy))*z1y[j,k,:]*vanleer(theta1y[j,k])\
                      + np.sign(s2y[j,k])*(1-((abs(s2y[j,k])*dt)/dy))*z2y[j,k,:]*vanleer(theta2y[j,k])\
                      + np.sign(s3y[j,k])*(1-((abs(s3y[j,k])*dt)/dy))*z3y[j,k,:]*vanleer(theta3y[j,k])\
                      + np.sign(s4y[j,k])*(1-((abs(s4y[j,k])*dt)/dy))*z4y[j,k,:]*vanleer(theta4y[j,k])
            Fy[j,k,:] = 0.5*Fy[j,k,:]
         
    for j in range(0,Nx+1):
        for k in range(0,Ny+1):
            Q1[j,k] = Q1[j,k] - ((dt/dy)*(Any[j,k+1,0] + Apy[j,k,0])) - (dt*e*trap2(vx,vy,gdx[j,k,:,:] + gdy[j,k,:,:],m1,m1)) - ((dt/dy)*(Fy[j,k+1,0]-Fy[j,k,0])) 
            Q2[j,k] = Q2[j,k] - ((dt/dy)*(Any[j,k+1,1] + Apy[j,k,1])) - (dt*e*trap2(vx,vx,gdx[j,k,:,:] + gdy[j,k,:,:],m2,m1)) - ((dt/dy)*(Fy[j,k+1,1]-Fy[j,k,1]))
            Q3[j,k] = Q3[j,k] - ((dt/dy)*(Any[j,k+1,2] + Apy[j,k,2])) - (dt*e*trap2(vy,vy,gdx[j,k,:,:] + gdy[j,k,:,:],m1,m2)) - ((dt/dy)*(Fy[j,k+1,2]-Fy[j,k,2])) 
            Q4[j,k] = Q4[j,k] - ((dt/dy)*(Any[j,k+1,3] + Apy[j,k,3])) - (dt*e*(trap2(vx,vx,gdx[j,k,:,:] + gdy[j,k,:,:],m3,m1)+(trap2(vy,vy,gdx[j,k,:,:] + gdy[j,k,:,:],m1,m3)))) - ((dt/dy)*(Fy[j,k+1,3]-Fy[j,k,3])) 
    '''
    for j in range(0,Nx+1):
        for k in range(0,Ny+1):
            Q1[j,k] = Q1[j,k] - ((dt/dy)*(Any[j,k+1,0] + Apy[j,k,0])) - ((dt/dy)*(Fy[j,k+1,0]-Fy[j,k,0])) 
            Q2[j,k] = Q2[j,k] - ((dt/dy)*(Any[j,k+1,1] + Apy[j,k,1])) - ((dt/dy)*(Fy[j,k+1,1]-Fy[j,k,1]))
            Q3[j,k] = Q3[j,k] - ((dt/dy)*(Any[j,k+1,2] + Apy[j,k,2])) - ((dt/dy)*(Fy[j,k+1,2]-Fy[j,k,2])) 
            Q4[j,k] = Q4[j,k] - ((dt/dy)*(Any[j,k+1,3] + Apy[j,k,3])) - ((dt/dy)*(Fy[j,k+1,3]-Fy[j,k,3]))              
    '''          
    rho = Q1
    ux = np.divide(Q2,Q1)
    uy = np.divide(Q3,Q1)
    T = np.divide((Q4 - 0.5*np.multiply(rho,np.square(ux)+np.square(uy))),rho)
'''
plt.figure(1)
plt.pcolor(x, y, np.transpose(rho))
plt.colorbar()
plt.title(r'Surface plot of rho vs x and y');
plt.xlabel('x')
plt.ylabel('y')
plt.show()

plt.figure(2)
plt.pcolor(x, y, np.transpose(ux))
plt.colorbar()
plt.title(r'Surface plot of ux vs x and y');
plt.xlabel('x')
plt.ylabel('y')
plt.show()

plt.figure(3)
plt.pcolor(x, y, np.transpose(uy))
plt.colorbar()
plt.title(r'Surface plot of uy vs x and y');
plt.xlabel('x')
plt.ylabel('y')
plt.show() 

plt.figure(4)
plt.pcolor(x, y, np.transpose(T))
plt.colorbar()
plt.title(r'Surface plot of T vs x and y');
plt.xlabel('x')
plt.ylabel('y')
plt.show() 
'''
plt.figure(1)
plt.grid('True')
plt.plot(x,rho[:,10]);
plt.title(r'Plot of rho vs x');
plt.xlabel(r'x');
plt.ylabel(r'rho');
  
plt.figure(2)
plt.grid('True')
plt.plot(x,ux[:,10]);
plt.title(r'Plot of ux vs x');
plt.xlabel(r'x');
plt.ylabel(r'u');

plt.figure(3)
plt.grid('True')
plt.plot(x,uy[:,10]);
plt.title(r'Plot of uy vs x');
plt.xlabel(r'x');
plt.ylabel(r'u');

plt.figure(4)
plt.grid('True')
plt.plot(x,T[:,10]);
plt.title(r'Plot of T vs x');
plt.xlabel(r'x');
plt.ylabel(r'T');

plt.show(); 





