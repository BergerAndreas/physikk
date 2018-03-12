from math import sin,cos,e
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
yn = -5.794662415E-1 # initial position x = 0
#yn = 0.25
vnx = 0 # initial value y(0) = 1
vn = 0
dt = 0.001 # step size


def dvxdt(yn,vnx): # right hand side of the ODE written on the form y' = f(x,y)
	m = 0.0051
	g = 9.81
	c = 0.0015
	f_l = c*abs(vnx)
	I = (2/5)*m
	alpha_x =  trvalues(polynomial,yn)[-2] #Se på denne etterpå
	
	if(vnx==0):
		return ((m*g*sin(alpha_x))/(I+m))

	return ((m*g*sin(alpha_x)-(f_l)*(vnx/abs(vnx)))/(I+m))*cos(trvalues(polynomial,yn)[-2])

def dvdt(yn,vn): # right hand side of the ODE written on the form y' = f(x,y)
	m = 0.0051
	g = 9.81
	c = 0.0015
	f_l = c*abs(vnx)
	I = (2/5)*m
	alpha_x =  trvalues(polynomial,yn)[-2] #Se på denne etterpå
	
	if(vnx==0):
		return ((m*g*sin(alpha_x))/(I+m))

	return ((m*g*sin(alpha_x)-(f_l)*(vnx/abs(vnx)))/(I+m))


def dxdt(yn,vn):
	return vn*cos(trvalues(polynomial,yn)[-2])

vnx_sol = []
vn_sol =[]
yn_sol = []

for n in range(0,42000):
	vnx = vnx + dt*dvxdt(yn,vnx)
	vn = vn + dt*dvdt(yn,vn)
	yn = yn + dt*dxdt(yn,vnx)
	vnx_sol.append(vnx)
	yn_sol.append(trvalues(polynomial,yn)[0])
	vn_sol.append(vn)
	#print(n+1,yn,vnx)

bergenet_k = []
for i in vn_sol:
	bergenet_k.append(((1/2)*0.0051*(i**2))+((0.4/2)*0.0051*(i**2)))

beregnet_p = []
for i in yn_sol:
	beregnet_p.append(0.0051*9.81*i)

min_p = min(beregnet_p)
for p in range(len(beregnet_p)):
	beregnet_p[p] =  beregnet_p[p]-min_p


total_e = []
for x in range(0,len(vnx_sol)):
	total_e.append(beregnet_p[x]+bergenet_k[x])
#plt.plot(vnx_sol)
"""
x_axis = [0.001*x for x in range(0,42000)]
plt.plot(x,yn_sol)
#t = [(0.31629073689431875*(e**(-0.001443169*x))*10**(-3)) for x in range(0,42000)]
#plt.plot(t)
print(max(yn_sol))
"""
x = np.linspace(0,42,42000)
#plt.plot(x,bergenet_k)
#plt.plot(x,yn_sol)
#plt.plot(x,vnx_sol)
#plt.axhline(y=0)
#plt.plot(x, beregnet_p)
plt.xlabel("Time s")
plt.ylabel("Energy [J]")
plt.title("Numerical total energy plot")
plt.plot(x,total_e)
#func = beregnet_p[0]*np.exp(-0.071*x)
#plt.plot(x, func)
plt.show()