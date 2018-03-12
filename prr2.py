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

steps = 45000
xn = -5.794662415E-1 # initial position x = 0
vn = 0 # initial value y(0) = 1
dt = 0.001 # step size
m = 0.0051
g = 9.81
c = 0.0015

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

#Load actual values
real_v = []
real_t = []
real_y = []
real_energy = []
with open("mass_A.txt") as file:
	next(file)
	next(file)
	for line in file:
		lines = line.split()
		real_t.append(float(lines[0]))
		real_y.append(m*g*float(lines[1]))
		real_v.append(0.5 * m * float(lines[2])**2 + 0.5 * 0.4 * m * float(lines[2])**2)
for i in real_v:
	real_energy.append(real_v+real_y)

x = np.linspace(0,45,steps)
#plt.plot(x,bergenet_k)
#plt.plot(x,yn_sol)
#plt.plot(x,vnx_sol)
plt.axhline(y=0, color = "black")
plt.axvline(x=0, color = "black")
plt.xlabel("Time s")
plt.ylabel("Energy [J]")
plt.title("Numerical total energy plot")
plt.plot(x, en_sol)
plt.plot(real_t,real_energy)
#func = 0.6*pn_sol[0]*np.exp(-0.071*x)
#plt.plot(x, func)
plt.show()
