# diffc.py
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#  
# Simulation of a wave going through media with different wave velocities
#
# This project is licensed under the GNU GPL v3. For more information about the 
# license go online to <http://www.gnu.org/licenses>. Any reproduction of code
# in this file MUST contain this header in its entirety.
#
from matplotlib.pylab import *
# import scipy.weave as weave (not used)
import time
Nx=513 # number of nodes
Ny=257
Ncxl = 300
Ncxh = 413
Ncyl = 230
Ncyh = 250
cred = 1.5
Nxwall=256 # x position of the dividing wall
w=256       # width of opening in dividing wall
Nbeg=(Ny-1-w)/2
Nend=(Ny-1+w)/2
Nt=2048 # number of time steps
Nstep=8 # number of time steps between saves.
ntau=200 # tau/dt - number of points per rise time.
Nsave=Nt/Nstep
nstart=-3*ntau # start three sigmas earlier than peak.
dt=1.0
alpha0=1/sqrt(2)
 # alpha=c*dt/dx is the courant number in 2-D
c=ones((Ny,Nx))     # speed of wave (no dispersion)
# Initialize the Electric and magnetic fields. For this, we need a dx
# that respects Courant condition. Leave it as a parameter.
c[Ncyl:Ncyh,Ncxl:Ncxh] = cred 
alpha=ones((Ny,Nx))/sqrt(2)
alpha[Ncyl:Ncyh,Ncxl:Ncxl] = cred/sqrt(2)
dx=c*dt/alpha
Lx=(Nx-1)*dx;Ly=(Ny-1)*dx
xwall=Nxwall*dx
xbeg=Nbeg*dx
xend=Nend*dx
x=zeros(Nx)
y=zeros(Ny)
x[:]=cumsum(dx[Ncyl+1,:])
y[:]=arange(0,Ly[0,0]+dx[0,0],dx[0,0])
xsrc=Lx/4.0
ysrc=Ly/2.0 # location of source is (xsrc,ysrc)
scalex=Nx/16.0 # size of source along x,y. 
scaley=Ny/8.0
Ntsave=4*ntau/Nstep # time when a snapshot is to be taken
# define coefficients
#alpha=120*pi*alpha  # just a scalar, since we assume a homogeneous region.
#alpha=alpha/(120*pi)
t = arange(nstart,Nt+nstart)*dt
Jsrc=exp(-(t/(ntau+0.0))**2/2.0)*sin(0.1*t)
X,Y=meshgrid(x,y)
Jshape=exp(-(X-xsrc)**2/scalex**2-(Y-ysrc)**2/scaley**2)
def getname(prompt,string):
  #str0="wall-w%d-Nx%d-Ny%d.dat" % (w,Nx,Ny)
  print "Enter name of save file for %s" % (prompt)
  str=raw_input(string+": ")
  if str=="":
    str=string
  return(str)
def ucalc(u,dudx,dudy,alpha):
  unew[1:-1,1:-1]=u[1:-1,1:-1]+alpha[1:,1:-1]*(dudy[1:,1:-1]-dudy[0:-1,1:-1] \
    -dudx[1:-1,1:]+dudx[1:-1,0:-1])
  unew[-1,:]=u[-2,:]+(1-alpha0)/(1+alpha0)*(u[-1,:]-unew[-2,:])
  unew[:,-1]=u[:,-2]+(1-alpha0)/(1+alpha0)*(u[:,-1]-unew[:,-2])
  unew[0,:]=u[1,:]+(1-alpha0)/(1+alpha0)*(u[0,:]-unew[1,:])
  unew[:,0]=u[:,1]+(1-alpha0)/(1+alpha0)*(u[:,0]-unew[:,1])
  unew[0:Nbeg,Nxwall]=0
  unew[Nend:,Nxwall]=0
  unew[0:Nbeg,Nxwall-1]=u[0:Nbeg,Nxwall-2]+(1-alpha0)/(1+alpha0)* \
                         (u[0:Nbeg,Nxwall-1]-unew[0:Nbeg,Nxwall-2])
  unew[Nend:,Nxwall-1]=u[Nend:,Nxwall-2]+(1-alpha0)/(1+alpha0)* \
                         (u[Nend:,Nxwall-1]-unew[Nend:,Nxwall-2])
  u[:]=unew[:]+(Jsrc[n]*Jshape)[:]
def dudycalc(u,dudy,alpha):
  dudy[:,1:-1]=alpha[1:,1:-1]*(u[1:,1:-1]-u[0:-1,1:-1])+dudy[:,1:-1]
def dudxcalc(u,dudx,alpha):
  dudx[1:-1,:]=-alpha[1:-1,1:]*(u[1:-1,1:]-u[1:-1,0:-1])+dudx[1:-1,:]
u=zeros((Ny,Nx)) # define the arrays holding the current time info
dudx=zeros((Ny,Nx-1))
dudy=zeros((Ny-1,Nx))
unew=zeros((Ny,Nx))
Uvals=zeros((Nsave,Ny,Nx))
Duxvals=zeros((Nsave,Ny,Nx-1))
Duyvals=zeros((Nsave,Ny-1,Nx))
t1=time.time()
nsave=0  # index that keeps track of save number
# this loop calculates the next time instant of the wave.
for n in range(Nt): # n indexes time
    ucalc(u,dudx,dudy,alpha)
    dudycalc(u,dudy,alpha)
    dudxcalc(u,dudx,alpha)
    # save data every Nstep iterations
    if n%Nstep==0:
      Uvals[nsave,:,:]=u
      Duyvals[nsave,:,:]=dudy
      Duxvals[nsave,:,:]=dudx
      nsave += 1
t2=time.time()
elapsed=(t2-t1)*1000.0
print "Time for iteration part in python-iter: %.3f msec" % elapsed
vals=linspace(-3,2,11)
fig2=figure(2)
set_cmap('jet')
contourf(x,y,log10(Uvals[0,:,:]),vals)
colorbar()
show()
hold(False)
Umag=zeros(Nsave)
Uwall=zeros(Ny)
# create the movie
for n in range(1,Nsave):
  if n*Nstep<Nt*0.5:
      contourf(x,y,log10(Uvals[n,:,:]),vals)
      hold(True)
      #plot([xwall,xwall],[0,xbeg])
      #plot([xwall,xwall],[xend,Ly])
      title(r"$E_z$ at time %.1f" % t[n*Nstep])
      draw()
      hold(False)
  if n==Ntsave:
      str=getname("Enter name of save file for contour", \
                  "contour-w%d-Nx%d-Ny%d.pdf" % (w,Nx,Ny))
      savefig(str,format="pdf")
  Umag[n]=norm(Uvals[n,:,:])**2
  Uwall+=Uvals[n,:,-1]**2
  # print n,Umag[n]
  time.sleep(0.05)

str=getname("Enter name of save file for Wall Intensity", \
            "wall-w%d-Nx%d-Ny%d.dat" % (w,Nx,Ny))
savefig(str,format="pdf")
savetxt(str,Uwall,header="w=%d,Nx=%d,Ny=%d,alpha=%.1f" % (w,Nx,Ny,alpha))
fname="energy-w%d-Nx%d-Ny%d.dat" % (Nend-Nbeg,Nx,Ny)
savetxt(fname,Umag,header="w=%d,Nx=%d,Ny=%d,alpha=%.1f" % (w,Nx,Ny,alpha))
figure(1)
plot(t,Jsrc)
grid()
title("Plot of the source as a function of time")
figure(4)
semilogy(t[::Nstep],Umag)
grid()
title("Plot of max Energy in wave as a function of time")
figure(3)
plot(y,Uwall/norm(Jsrc))
grid()
title("Time averaged wall intensity at $x=L_x$")
show()
