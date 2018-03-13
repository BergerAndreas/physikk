from matplotlib import pyplot as plt
from matplotlib import pylab
from pylab import *
from math import log
import numpy as np
m = 0.0051
g = 9.81
c = 0.0015
#Load actual values
real_t = []
real_y = []
real_v = []
with open("mass_A_withspeed.txt") as file:
	next(file)
	next(file)
	for line in file:
		lines = line.split()
		real_t.append(float(lines[0]))
		real_y.append(float(lines[1]))
		real_v.append(float(lines[2]))


pn_sol = []
kn_sol = []
real_energy = []


for i in range(len(real_t)):

    pn_sol.append(m*g*real_y[i])
    kn_sol.append(0.5 * m * real_v[i]**2 + 0.5 * 0.4 * m * real_v[i]**2)

for i in range(len(pn_sol)):
	pn_sol[i] = pn_sol[i]-pn_sol[-1]
	if(pn_sol[i]<0):
		pn_sol[i]=pn_sol[i]/1
for i in range(len(real_t)):
	real_energy.append(pn_sol[i]+kn_sol[i])

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

print(fitExponent(real_t,real_energy, min(real_energy)))

x = np.linspace(0,40,45000)
func = 0.016679188680116492*np.exp(- x/6.840358227028433)
plt.plot(x, func)
#plt.plot(0,0103203543e-0,0006409612x)
plt.plot(real_t,real_energy)
plt.axhline(y=0, color = "black")
plt.axvline(x=0, color = "black")
plt.xlabel("Time [s]")
plt.ylabel("Energy [J]")
plt.title("Measured real energy plot")
plt.show()