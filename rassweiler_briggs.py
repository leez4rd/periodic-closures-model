#Rassweiler-Briggs model of reef dynamics 
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint 

phiC = 0.001 #open recruitment of coral
phiM = 0.0001 #open recruitment of macroalgae 

rM = 0.5 #production of vulnerable macroalgae from invulnerable stage

#growth rates
gTC = 0.1 #combined rate of growth and local recruitment of corals over free space 
gTV = 0.2 #growth rate of vulnerable macroalgae over free space 
gTI = 0.4 #growth rate of invulnerable macroalgae over free space 
rH = 0.49

gamma = 0.4 #growth of macroalgae over coral vs free space
omega = 2 #maturation rate of macroalgae from vulnerable to invulnerable class 

#death rates
dC = 0.5 #death rate of coral 
dI = 0.4 #death rate of invulnerable macroalgae
dV = 0.58 #death rate of vulnerable macroalgae per unit biomass of herbivores 

K = 20
Graze = 0.58 

#reference points  
C_max = .509  # coral cover with no fishing
P_max = 20  # parrotfish with no fishing
M_max = .466    # algal cover with really high fishing - note this is Mi only
f_crit = .275  # highest fishing at which coral wins (regardless of initial)

#IC's 
X0 = [0.1, 0.4, 0.4, 0.01]

N = 100 #total amount of time
steps = 1000 #number of timesteps
t = np.linspace(0, N, steps) #timestep array
period = 20
f = 0.3
p = 0.3

def system(X, t, period, f, p):
	P, C, Mv, M = X

	T = 1 - C - Mv - M 
	dPdt = rH*P*(1-P/K)-f*P
	dCdt = (phiC*T) + gTC*T*C - gamma*gTI*M*C - dC*C
	dMvdt = phiM*T + rM*T*M + gTV*T*Mv - dV*Mv - P*Mv*Graze - omega * Mv
	#why is that second term multiplied by M not Mv? 
	dMdt = omega*Mv + gTI*T*M + gamma*gTI*M*C - dI*M
	
	return [dPdt, dCdt, dMvdt, dMdt]

sol = odeint(system, X0, t, args = (period, f, p))
P_sol = sol[:, 0]
C_sol = sol[:, 1]
Mv_sol = sol[:, 2]
Mi_sol = sol[:, 3]
plt.figure()
plt.plot(t, C_sol, label='Coral')
plt.plot(t, Mi_sol, label='Invulnerable Macroalgae')
plt.plot(t, Mv_sol, label='vulnerable Macroalgae')
#plt.plot(t, P_sol, label = 'Parrotfish')
plt.plot(t, 1-Mi_sol-C_sol-Mv_sol, label = 'Algal Turfs')
plt.legend(loc=0)

plt.xlabel('Time')
plt.ylabel('Abundances')
plt.title('Rassweiler-Briggs Model')
plt.show()