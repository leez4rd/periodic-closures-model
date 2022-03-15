import numpy as np 
import scipy
from scipy.integrate import odeint 
import matplotlib.pyplot as plt
import seaborn as sb
from mpl_toolkits.mplot3d import Axes3D
import math


#v a r i a b l e s
#############################################################################33


#parameters
r = 0.3
i_C = 0.05
i_M = 0.05
ext_C = 0.0001
ext_P = 0.0001
gamma = 0.8
d = 0.1
g = 1
s = 1
sigma = .5 #strength of coral-herbivore feedback
eta = 2 #strength of algae-herbivore feedback
alpha = 0.5 #strength of algae-coral feedback 

#time
N = 2000 #total amount of time
steps = 1000 #number of timesteps
t = np.linspace(0, N, steps) #timestep array

#baseline reference points 
C_max = .509  # coral cover with no fishing
P_max = 20  # parrotfish with no fishing
M_max = .466    # algal cover with really high fishing - note this is Mi only

P0 = P_max*.5
C0H = C_max*.95
M0H = M_max*.8

C0L = C_max*.05
M0L = M_max*.01

P0 = 0.1
XHI = [0.1, C0H, M0L]
XLO = [0.1, C0L, M0H]
#high and low coral starting points are 0.48355, 0.02545

#######################################################



#f u n c t i o n s
###########################################################
def initialize_patch_model(n, frac_nomove):
	frac_dispersed = (1-frac_nomove)*(1/(n)) #fraction of fish that disperse to other patches symmetrically

	#transition matrix for dispersal: element [i,j] of kP describes influx of P from j to i
	#doesn't include fish staying in patch (just subtracts the ones that leave) but math should still be consistent 
	global kP, P, C, M, dPs, dCs, dMs 
	kP = np.empty((n,n))
	for i in range(n):
		for j in range(n):
			kP[i][j] = frac_dispersed
			if i == j:
				kP[i][j] = -frac_dispersed


	P_influx = np.empty(n)
	P = np.empty(n) 
	C = np.empty(n) 
	M = np.empty(n)
	dPs = np.empty(n)
	dCs = np.empty(n)
	dMs = np.empty(n)

	#concatenate initial condition arrays 
	global X1
	X1 = [P0]*n + [C0L]*n + [M0H]*n
	global X2
	X2 = [P0]*n + [C0H]*n + [M0L]*n


def K(sigma, C):
	return (1-sigma)+sigma*C

def square_signal(t, closure_length, region, m, n):
	start = int((t % (n*closure_length))/closure_length)
	if start+m-1 >= n:
		end = (start + m - 1)%n
	else:
		end = (start + m - 1)
	if region >= start and region <= end:
		return 0
	elif start + m - 1 >= n and (region >= start or region <= end):
		return 0
	else:
		return 1
	
#this signal function is not quite working yet 
def sigmoid_signal(t, period, p):
	if period == 0:
		return 0
	else:
		return 1.0 / (1 + math.exp(-(t % period - p * period)))


def patch_system(X, t, closure_length, f, m, n): 
	P_influx = [0]*n
	C_influx = [0]*n
	M_influx = [0]*n
	P,C,M = X.reshape(3, n) 

	for i in range(n):
		for j in range(n):
			P_influx[i] += (kP[i][j]) *P[j]  
		dPs[i] = P_influx[i]+ s*P[i]*(1 - (P[i] / K(sigma,C[i]))) - f*P[i] *(square_signal(t, closure_length, i, m, n)) 
		dCs[i] = (i_C + r*C[i])*(1-M[i]-C[i])*(1-alpha*M[i]) - d*C[i]
		dMs[i] = (i_M+gamma*M[i])*(1-M[i]-C[i])-g*M[i]*P[i]/(g*eta*M[i]+1) 
	#concatenate into 1D vector to pass to next step
	return np.concatenate((dPs, dCs, dMs), axis=0)


#general purpose graphing function 
def graph_sol(closure_length, f, m, n, frac_nomove, coral_high):

	initialize_patch_model(n, frac_nomove)

	if (coral_high):
		IC_set = X2
	else:
		IC_set = X1
	plt.figure()
	sol = odeint(patch_system, IC_set, t, args = (closure_length, f/(1-m/n), m, n))
	print(sol[:, n:2*n])
	plt.plot(sol[:, n:2*n],  label = 'coral')
	plt.xlabel('time')
	plt.ylabel('abundance')
	plt.title('spatial model')
	plt.legend(loc=0)
	txt="parameters" + "\nfishing when open: " + str(f/(1-m/n)) + "\npercent closure: " + str(m/n) +"\nclosure time: " + str(closure_length)
	plt.figtext(0.3, .31, txt, wrap=True, fontsize=8)
	plt.show()

#########################################################################################
'''
#some graphs for exploring the model

graph_sol(20, 0.25, 2, 5, 0.95, False)
#two patches 
graph_sol(50, 0.25, 0, 2, 0.95, False)
graph_sol(50, 0.25, 1, 2, 0.95, False)

#three patches
graph_sol(50, 0.25, 1, 3, 1, False)

#five patches

graph_sol(20, 0.25, 1, 5, 0.95, False)

#why are some patches recovering and not others for some of these?


graph_sol(40, 0.35, 2, 5, 0.95, False)


#ten patches!
graph_sol(10, 0.22, 1, 10, 0.95, False)

graph_sol(10, 0.22, 2, 10, 0.95, False)

#why does this not reduce to the 1 out of 2 case? overlap alters the scenario? 
graph_sol(10, 0.22, 5, 10, 0.95, False)

'''

initialize_patch_model(10, 0.9)

n = 10
M =  10

#plot average coral over all regions after convergence to equilibrium versus m, versus n, and with different dispersals
f = 0.25
closure_lengths = [5, 10, 20, 50, 100]
ms = np.empty(10)
coral_covers = np.empty(10)
for closure_length in closure_lengths: 
	avg = 0 
	for m in range(M):
		sol = odeint(patch_system, X1, t, args = (closure_length, f/(1-m/n), m, n))
		
		for year in range(999- (999 % closure_length*n) - 2*closure_length*n, 999 - (999 % closure_length*n)):
			avg += sol[year][1]
		avg = avg / (2*closure_length*n + 1)
		ms[m] = m
		coral_covers[m] = avg
	plt.plot(ms, coral_covers, label = 'coral starts low, CL = %d' % closure_length)
plt.xlabel('percentage time closed')
plt.ylabel('coral cover at end')
plt.legend(loc=0)
plt.show()

#now plotting 1/n closures as we vary n 
f = 0.25
m = 1
closure_lengths = [1, 5, 10, 20, 50, 100]
ns = np.empty(10)
coral_covers = np.empty(10)
for closure_length in closure_lengths: 
	avg = 0 
	for n in range(1,M+1):
		sol = odeint(patch_system, X1, t, args = (closure_length, f/(1-m/n), m, n))
		
		for year in range(999- (999 % closure_length*n) - 2*closure_length*n, 999 - (999 % closure_length*n)):
			avg += sol[year][1]
		avg = avg / (2*closure_length*n + 1)
		ns[n] = n
		coral_covers[n] = avg
	plt.plot(ns, coral_covers, label = 'coral starts low, CL = %d' % closure_length)
plt.xlabel('percentage time closed')
plt.ylabel('coral cover at end')
plt.legend(loc=0)
plt.show()










#heatmap of period vs m -- a bit messy right now 
coral_array_HI =  np.zeros(shape=(M,M))
coral_array_LO =  np.zeros(shape=(M,M))
coral_array_AVG =  np.zeros(shape=(M,M))
period_array = np.empty(M)
m_array = np.empty(M)
fishing = 0.25

for period in range(1,M+1):
		for m in range(M):
			displacer = 1/(1-float(m)/float(n))
			#final_coral_HI = odeint(patch_system, XHI, t, args = (period, displacer*np.float128(fishing) / np.float128(n) ,np.float128(p) / np.float128(n)))#full_output = 1)
			final_coral_LO = odeint(patch_system, X1, t, args = (period*20, displacer*float(fishing) ,m, 10))# full_output = 1)
			avg1 = 0
			avg2 = 0
			for year in range(999- (999 % (period)) - 2*(period), 999 - (999 % (period))):
				#avg1 += final_coral_HI[year][1]
				avg2 += final_coral_LO[year][1]
			avg1 = avg1 / (2*(period) + 1)
			avg2 = avg2 / (2*(period) + 1)
			#coral_array_HI[period-1][m] = avg1
			coral_array_LO[period-1][m] = avg2
			coral_array_AVG[period-1][m] = 0.5 * (avg1 + avg2)
			period_array[period-1] = period
			m_array[m] = m
			show_labels = False

plt.title('heatmap', fontsize = 20)
sb.heatmap(coral_array_LO, vmin = 0.0, vmax = 1.0, cmap = "viridis", annot = show_labels) #YlGnBu for original color scheme
plt.ylabel('period', fontsize = 10) # x-axis label with fontsize 15
plt.xlabel('number of closures', fontsize = 10) # y-axis label with fontsize 15
plt.yticks(rotation=0)
plt.show()





'''
#RHS's for each ODE
def parrot_eq(P, C, t, period, p):
	return s*P*(1 - (P / K(sigma,C))) - f*P *square_signal(t, period, p)
def coral_eq(C, M):
	return (i_C + r*C)*(1-M-C)*(1-alpha*M) - d*C #+ ext_C
def algae_eq(P, C, M):
	return (i_M+gamma*M)*(1-M-C)-g*M*P/(g*eta*M+1) #+ ext_C
'''