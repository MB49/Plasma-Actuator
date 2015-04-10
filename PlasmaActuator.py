from __future__ import division
import numpy as np
import time
from pylab import *
from matplotlib import *
T = np.transpose
import matplotlib.animation as animation

nx = 61
ny = 61
nt = 100

Length=0.006#6 mm
dx = Length/(nx-1)
dy = Length/(ny-1)
x = np.linspace(-Length/2,Length/2,nx)
y = np.linspace(-Length/2,Length/2,ny)
#################################
phi = np.zeros((ny,nx)) #potential
eps = np.zeros((ny,nx)) #permitivity
roc=np.ones((ny,nx))#charge
ld=np.ones((nx,ny))#Debye length
# initial/boundary conditions

phi[30:32,11:14]=0
phi[27:29,13:24]=0
phi[30:32,23:26]=1000
phi[30:32,35:38]=0
phi[27:29,37:48]=0
phi[30:32,47:50]=1000
eps[0:30,:]=2.7
eps[30:,:]=1
eps[30,:]=(eps[31,:]*eps[29,:])/((eps[31,:]\
        *(dx/(2*dx)))+((eps[29,:]*(dx/(2*dx)))))
#########################################
# convergence parameters
tol = 1e-8

# Gauss-Seidel solver
tic = time.time()
max_phi_diff = 1;
while (max_phi_diff > tol):

    phi_old = phi.copy()

    

    phi[1:nx-1,1:ny-1]=((eps[1:nx-1,1:ny-1]*((phi[2:nx,1:ny-1]/(dx**2))+(phi[1:nx-1,2:ny]/dy**2)))+\
                              ((eps[0:nx-2,1:ny-1]*phi[0:nx-2,1:ny-1])/dx**2)+((eps[1:nx-1,0:ny-2]*phi[1:nx-1,0:ny-2])/dy**2))\
                          /(((eps[1:nx-1,1:ny-1]+eps[0:nx-2,1:ny-1])/dx**2)+((eps[1:nx-1,1:ny-1]+eps[1:nx-1,0:ny-2])/dy**2))

    phi[-1,:] =phi[-2,:]    ##dphi/dy = 0 at y = 3 mm
    phi[0,:] = phi[1,:]     ##dphi/dy = 0 at y = -3 mm
    phi[:,-1] =phi[:,-2]    ##dphi/dx = 0 at x = 3 mm
    phi[:,0] = phi[:,1]     ##dphi/dx = 0 at x = -3 mm



    
    phi[30:32,11:14]=0
    phi[27:29,13:24]=0
    phi[30:32,23:26]=1000
    phi[30:32,35:38]=0
    phi[27:29,37:48]=0
    phi[30:32,47:50]=1000
    
    # check for convergence
    phi_diff = phi - phi_old
    max_phi_diff = np.absolute(phi_diff).max()

toc = time.time() - tic
print(toc)


X,Y=np.meshgrid(x,y)
fig = plt.figure(figsize = (11,7), dpi=100)
plt.title('Electric Potential',fontsize=16)
plt.contourf(X*1000,Y*1000,phi,50,cmap=plt.cm.jet)
cbar = plt.colorbar()
ax = fig.add_subplot(111)
rect1 = matplotlib.patches.Rectangle((-1.9,0), 0.2, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect1)
rect2 = matplotlib.patches.Rectangle(((-1.7,-0.3)), 1, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect2)
rect3 = matplotlib.patches.Rectangle((-0.7,0), 0.2, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect3)
rect4 = matplotlib.patches.Rectangle((0.5,0), 0.2, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect4)
rect5 = matplotlib.patches.Rectangle((0.7,-0.3), 1, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect5)
rect6 = matplotlib.patches.Rectangle((1.7,0), 0.2, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect6)
plt.xlabel('X Axis (mm)', fontsize=16)
plt.ylabel('Y Axis (mm)', fontsize=16)
cbar.ax.set_ylabel('Potential (V)',fontsize=16)
show()


ld[30:,:]=0.00017 #charge
ld[0:30,:]=10000#Debye lenght

roc[27:29,13:24]=0.00750

roc[27:29,37:48]=0.00750

eps[0:30,:]=2.7
eps[30:,:]=1
eps[30,:]=(eps[31,:]*eps[29,:])/((eps[31,:]\
        *(dx/(2*dx)))+((eps[29,:]*(dx/(2*dx)))))
roc[:,-1] =0
roc[-1,:] = 0
roc[:,0] =0
roc[0,:] = 0

# convergence parameters
tol = 1e-8

# Gauss-Seidel solver
tic = time.time()
max_roc_diff = 1;
while (max_roc_diff > tol):

    roc_old = roc.copy()

    

    roc[1:nx-1,1:ny-1]=((eps[1:nx-1,1:ny-1]*((roc[2:nx,1:ny-1]/(dx**2))+(roc[1:nx-1,2:ny]/dy**2)))+\
                              ((eps[0:nx-2,1:ny-1]*roc[0:nx-2,1:ny-1])/dx**2)+((eps[1:nx-1,0:ny-2]*roc[1:nx-1,0:ny-2])/dy**2))\
                              /(((eps[1:nx-1,1:ny-1]+eps[0:nx-2,1:ny-1])/dx**2)+((eps[1:nx-1,1:ny-1]+eps[1:nx-1,0:ny-2])\
                                                                                 /dy**2)+(1/(ld[1:nx-1,1:ny-1]**2)))
    
    ld[30:,:]=0.00017 
    ld[0:30,:]=10000
    roc[27:29,13:24]=0.00750    
    roc[27:29,37:48]=0.00750    
    roc_diff = roc - roc_old
    max_roc_diff = np.absolute(roc_diff).max()

toc = time.time() - tic
print(toc)


X,Y=np.meshgrid(x,y)
fig = plt.figure(figsize = (11,7), dpi=100)
plt.title('Electric Charge',fontsize=16)
plt.contourf(X*1000,Y*1000,roc,50,cmap=plt.cm.jet)
cbar = plt.colorbar()
ax = fig.add_subplot(111)
rect1 = matplotlib.patches.Rectangle((-1.9,0), 0.2, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect1)
rect2 = matplotlib.patches.Rectangle(((-1.7,-0.3)), 1, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect2)
rect3 = matplotlib.patches.Rectangle((-0.7,0), 0.2, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect3)
rect4 = matplotlib.patches.Rectangle((0.5,0), 0.2, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect4)
rect5 = matplotlib.patches.Rectangle((0.7,-0.3), 1, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect5)
rect6 = matplotlib.patches.Rectangle((1.7,0), 0.2, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect6)
plt.xlabel('X Axis (mm)', fontsize=16)
plt.ylabel('Y Axis (mm)', fontsize=16)
cbar.ax.set_ylabel('Charge Density (C/m^3)',fontsize=16)
show()
F=np.ones((ny,nx))
Fx=np.ones((ny,nx))
Fy=np.ones((ny,nx))
F[1:nx-1,1:ny-1]=roc[1:nx-1,1:ny-1]*(-(((phi[1:nx-1,1:ny-1]-phi[0:nx-2,1:ny-1])/dx)+((phi[1:nx-1,1:ny-1]-phi[1:nx-1,0:ny-2])/dy)))
F[30:32,11:14]=0
F[27:29,13:24]=0
F[30:32,23:26]=0
F[30:32,35:38]=0
F[27:29,37:48]=0
F[30:32,47:50]=0
F[0:30,:]=0
fig = plt.figure(figsize = (11,7), dpi=100)
plt.title('Body Force',fontsize=16)
plt.contourf(X*1000,Y*1000,F,50,cmap=plt.cm.jet)
cbar = plt.colorbar()
ax = fig.add_subplot(111)
rect1 = matplotlib.patches.Rectangle((-1.9,0), 0.2, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect1)
rect2 = matplotlib.patches.Rectangle(((-1.7,-0.3)), 1, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect2)
rect3 = matplotlib.patches.Rectangle((-0.7,0), 0.2, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect3)
rect4 = matplotlib.patches.Rectangle((0.5,0), 0.2, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect4)
rect5 = matplotlib.patches.Rectangle((0.7,-0.3), 1, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect5)
rect6 = matplotlib.patches.Rectangle((1.7,0), 0.2, 0.1, facecolor='none',linewidth=2)
ax.add_patch(rect6)
plt.xlabel('X Axis (mm)', fontsize=16)
plt.ylabel('Y Axis (mm)', fontsize=16)
cbar.ax.set_ylabel('Body Force',fontsize=16)
show()
Fx[1:nx-1,1:ny-1]=roc[1:nx-1,1:ny-1]*(-(((phi[1:nx-1,1:ny-1]-phi[0:nx-2,1:ny-1])/dx)))
Fy[1:nx-1,1:ny-1]=roc[1:nx-1,1:ny-1]*(-((phi[1:nx-1,1:ny-1]-phi[1:nx-1,0:ny-2])/dy))
Fx=Fx
Fy=Fy





#Navier Stokes
# projection method in primitive variables on a colocated mesh

dt=0.0000025;#time step
nstep=2001;#number of steps
mu=0.000018;#0.1kinematic viscosity
maxit=100;
beta=1.2;
h=Length/(nx-1);
time=0;
u=np.zeros((nx,ny));
v=np.zeros((nx,ny));
p=np.zeros((nx,ny));
ut=np.zeros((nx,ny));
vt=np.zeros((nx,ny));
u[1:nx-1,ny-1]=0;
for iter in range (1,nstep):
  
  Fn=np.zeros((nx,ny))
  Fxn=np.zeros((nx,ny))
  Fyn=np.zeros((nx,ny))
  def Fn(F,iter):
      Fn=F*abs(sin(2*np.pi*5000*(iter*dt)))
      return Fn
  Fn=Fn(F,iter)
  def Fxn(Fx,iter):
      Fxn=Fx*abs(sin(2*np.pi*5000*(iter*dt)))
      return Fxn
  Fxn=Fxn(Fx,iter)
  def Fyn(Fy,iter):
      Fyn=Fy*abs(sin(2*np.pi*5000*(iter*dt)))
      return Fyn
  Fyn=Fyn(Fy,iter)
     
    
  
  
  
  for i in range (1,nx-1):
    for j in range (1,ny-1):
       
       
       ut[i,j]=u[i,j]+dt*(-(0.5/h)*(u[i,j]*(u[i+1,j]-u[i-1,j])+v[i,j]*(u[i,j+1]-u[i,j-1]))+\
                          (mu/h**2)*(u[i+1,j]+u[i-1,j]+u[i,j+1]+u[i,j-1]-4*u[i,j]))+dt*Fxn[i,j];
       vt[i,j]=v[i,j]+dt*(-(0.5/h)*(u[i,j]*(v[i+1,j]-v[i-1,j])+v[i,j]*(v[i,j+1]-v[i,j-1]))+\
                          (mu/h**2)*(v[i+1,j]+v[i-1,j]+v[i,j+1]+v[i,j-1]-4*v[i,j]))+dt*Fyn[i,j];

  vt[1:nx-1,0]=(mu*dt/h**2)*(v[1:nx-1,2]-2*v[1:nx-1,1]);
  vt[1:nx-1,ny-1]=(mu*dt/h**2)*(v[1:nx-1,ny-3]-2*v[1:nx-1,ny-2]);
  ut[0,1:ny-1]=(mu*dt/h**2)*(u[2,1:ny-1]-2*u[1,1:ny-1]);
  ut[nx-1,1:ny-1]=(mu*dt/h**2)*(u[nx-3,1:ny-1]-2*u[nx-2,1:ny-1]);
  for it in range (1,maxit): 
    p[1:nx-1,0]=p[1:nx-1,1]+(h/dt)*vt[1:nx-1,0];
    p[1:nx-1,ny-1]=p[1:nx-1,ny-2]-(h/dt)*vt[1:nx-1,ny-1];
    p[0,1:ny-1]=p[1,1:ny-1]+(h/dt)*ut[0,1:ny-1];
    p[nx-1,1:ny-1]=p[nx-2,1:ny-1]-(h/dt)*ut[nx-1,1:ny-1];
    for i in range (1,nx-1):
     for j in range (1,ny-1):
       p[i,j]=beta*0.25*(p[i+1,j]+p[i-1,j]+p[i,j+1]+p[i,j-1]-(0.5*h/dt)*(ut[i+1,j]-ut[i-1,j]+vt[i,j+1]-vt[i,j-1]))+(1-beta)*p[i,j];

    p[floor(nx/2),floor(ny/2)]=0.0; # set the pressure in the center. Needed since bc is not incorporated into eq

  u[1:nx-1,1:ny-1]=ut[1:nx-1,1:ny-1]-(0.5*dt/h)*(p[2:nx,1:ny-1]-p[0:nx-2,1:ny-1]);
  v[1:nx-1,1:ny-1]=vt[1:nx-1,1:ny-1]-(0.5*dt/h)*(p[1:nx-1,2:ny]-p[1:nx-1,0:ny-2]);
  u[0:30,:] = 0
  u[:,-1] = 0
  v[0:30,:] = 0
  v[:,-1]=0
  u[30:32,11:14]=0
  u[27:29,13:24]=0
  u[30:32,23:26]=0
  u[30:32,35:38]=0
  u[27:29,37:48]=0
  u[30:32,47:50]=0
      
  v[30:32,11:14]=0
  v[27:29,13:24]=0
  v[30:32,23:26]=0
  v[30:32,35:38]=0
  v[27:29,37:48]=0
  v[30:32,47:50]=0
  
  time=time+dt
  w=np.sqrt(v**2+u**2)
  if (np.mod(iter,20)==0):
   fig = plt.figure(figsize = (11,7), dpi=100)
   plt.title('U',fontsize=16)
   plt.contourf(x*1000,y*1000,u,55,cmap=plt.cm.jet)
   cbar = plt.colorbar()
   ax = fig.add_subplot(111)
   rect1 = matplotlib.patches.Rectangle((-1.9,0), 0.2, 0.1, facecolor='none',linewidth=2)
   ax.add_patch(rect1)
   rect2 = matplotlib.patches.Rectangle(((-1.7,-0.3)), 1, 0.1, facecolor='none',linewidth=2)
   ax.add_patch(rect2)
   rect3 = matplotlib.patches.Rectangle((-0.7,0), 0.2, 0.1, facecolor='none',linewidth=2)
   ax.add_patch(rect3)
   rect4 = matplotlib.patches.Rectangle((0.5,0), 0.2, 0.1, facecolor='none',linewidth=2)
   ax.add_patch(rect4)
   rect5 = matplotlib.patches.Rectangle((0.7,-0.3), 1, 0.1, facecolor='none',linewidth=2)
   ax.add_patch(rect5)
   rect6 = matplotlib.patches.Rectangle((1.7,0), 0.2, 0.1, facecolor='none',linewidth=2)
   ax.add_patch(rect6)
   plt.xlabel('X Axis (mm)', fontsize=16)
   plt.ylabel('Y Axis (mm)', fontsize=16)
   cbar.ax.set_ylabel('Velocity u (m/s)',fontsize=16)
   show()
   fig = plt.figure(figsize = (11,7), dpi=100)
   plt.title('Flow',fontsize=16)
   plt.quiver(x*1000,y*1000,v,u,w)
   cbar = plt.colorbar()
   ax = fig.add_subplot(111)
   rect1 = matplotlib.patches.Rectangle((-1.9,0), 0.2, 0.1, facecolor='none',linewidth=2)
   ax.add_patch(rect1)
   rect2 = matplotlib.patches.Rectangle(((-1.7,-0.3)), 1, 0.1, facecolor='none',linewidth=2)
   ax.add_patch(rect2)
   rect3 = matplotlib.patches.Rectangle((-0.7,0), 0.2, 0.1, facecolor='none',linewidth=2)
   ax.add_patch(rect3)
   rect4 = matplotlib.patches.Rectangle((0.5,0), 0.2, 0.1, facecolor='none',linewidth=2)
   ax.add_patch(rect4)
   rect5 = matplotlib.patches.Rectangle((0.7,-0.3), 1, 0.1, facecolor='none',linewidth=2)
   ax.add_patch(rect5)
   rect6 = matplotlib.patches.Rectangle((1.7,0), 0.2, 0.1, facecolor='none',linewidth=2)
   ax.add_patch(rect6)
   plt.xlabel('X Axis (mm)', fontsize=16)
   plt.ylabel('Y Axis (mm)', fontsize=16)
   cbar.ax.set_ylabel('Flow (m/s)',fontsize=16)
   show()
   #savefig('wflow='+str(iter)+'.png')
   #clf()

   
