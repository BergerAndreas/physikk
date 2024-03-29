from math import sin,cos
from matplotlib import pyplot as plt
from matplotlib import pylab
from pylab import *
from math import log
import numpy as np

def iptrack(filename):
	data=np.loadtxt(filename,skiprows=2)
	return np.polyfit(data[:,1],data[:,2],15)


def trvalues(p,x):
	y=np.polyval(p,x)
	dp=np.polyder(p)
	dydx=np.polyval(dp,x)
	ddp=np.polyder(dp)
	d2ydx2=np.polyval(ddp,x)
	alpha=np.arctan(-dydx)
	R=(1.0+dydx**2)**1.5/d2ydx2
	return [y,dydx,d2ydx2,alpha,R]

polynomial = iptrack("mass_A.txt")

steps = 45000
xn = -5.794662415E-1 # initial position x = 0
vn = 0 # initial value y(0) = 1
dt = 0.001 # step size
m = 0.0051
g = 9.81
c = 0.0013

def dvdt(xn,vn): # right hand side of the ODE written on the form y' = f(x,y)

	f_l = c*abs(vn)
	I = (2/5)*m
	alpha_x =  trvalues(polynomial,xn)[-2] #Se på denne etterpå

	if(vn==0):
		return ((m*g*sin(alpha_x))/(I+m))

	return ((m*g*sin(alpha_x)-(f_l)*(vn/abs(vn)))/(I+m))


def dxdt(xn,vn):
	return vn*cos(trvalues(polynomial,xn)[-2])

vn_sol = []
xn_sol = []
yn_sol = []
kn_sol = []
pn_sol = []
en_sol = []

for n in range(0,steps):
    vn = vn + dt*dvdt(xn,vn)
    xn = xn + dt*dxdt(xn,vn)
    yn = trvalues(polynomial,xn)[0]

    xn_sol.append(xn)
    vn_sol.append(vn)
    yn_sol.append(yn)
    pn_sol.append(m*g*yn)
    kn_sol.append(0.5 * m * vn**2 + 0.5 * 0.4 * m * vn**2)

	#print(n+1,xn,vn)


#Fiix potensiell energi
pn_sol_min = min(pn_sol)
for i in range(len(pn_sol)):
    pn_sol[i] = pn_sol[i] - pn_sol_min

#bergen total energi
for i in range(len(pn_sol)):
    en_sol.append(pn_sol[i] + kn_sol[i])


x = np.linspace(0,45,steps)
def fitExponent(tList,yList,ySS=0):
   '''
   This function finds a 
       tList in sec
       yList - measurements
       ySS - the steady state value of y
   returns
       amplitude of exponent
       tau - the time constant
   '''
   bList = [log(max(y-ySS,1e-6)) for y in yList]
   b = matrix(bList).T
   rows = [ [1,t] for t in tList]
   A = matrix(rows)
   #w = (pinv(A)*b)
   (w,residuals,rank,sing_vals) = lstsq(A,b)
   tau = -1.0/w[1,0]
   amplitude = exp(w[0,0])
   return (amplitude,tau)
print(fitExponent(x,en_sol, min(en_sol)))


#plt.plot(x,bergenet_k)
#plt.plot(x,yn_sol)
#plt.plot(x,vnx_sol)
func = en_sol[0]*np.exp(- x/6.840358227028433)
plt.axhline(y=0, color = "black")
plt.axvline(x=0, color = "black")
plt.xlabel("Time [s]")
plt.ylabel("Energy [J]")
plt.title("Numerical total energy plot")
plt.plot(x, en_sol)
#func = 0.6*pn_sol[0]*np.exp(-0.071*x)
#plt.plot(x, func)
plt.show()
