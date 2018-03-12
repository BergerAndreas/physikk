from math import sin,cos
from matplotlib import pyplot as plt
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

xn = -5.794662415E-1 # initial position x = 0
vn = 0 # initial value y(0) = 1
dt = 0.001 # step size


def dvdt(xn,vn): # right hand side of the ODE written on the form y' = f(x,y)
	m = 0.0051
	g = 9.81
	c = 0.0015
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

for n in range(0,42000):
	vn = vn + dt*dvdt(xn,vn)
	xn = xn + dt*dxdt(xn,vn)
	xn_sol.append(trvalues(polynomial,xn)[0])
	vn_sol.append(vn)
	#print(n+1,xn,vn)

plt.plot(vn_sol)
plt.plot(xn_sol)
plt.axhline(y=0)
print(max(vn_sol))
plt.show()