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
N = 1000 #total amount of time
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

#fishing effort density dependence
steepness = 25
shift = 0.0001

#poaching
poaching = 0


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

	if frac_nomove == 0:
		print(kP)

	#concatenate initial condition arrays 
	global X1
	X1 = [P0]*n + [C0L]*n + [M0H]*n
	global X2
	X2 = [P0]*n + [C0H]*n + [M0L]*n


def K(sigma, C):
	return (1-sigma)+sigma*C

def square_signal(t, closure_length, region, m, n, poaching):
	start = int((t % (n*closure_length))/closure_length)
	if start+m-1 >= n:
		end = (start + m - 1)%n

	else:
		end = (start + m - 1)
	if region >= start and region <= end:
		return poaching
	elif start + m - 1 >= n and (region >= start or region <= end):
		return poaching
	else:
		#determine whether to bake displacement into signal
		return (1-(m/n)*poaching)#/(1-(m/n))
	
#this signal function is not quite working yet 
def sigmoid_signal(t, period, p):
	if period == 0:
		return 0
	else:
		return 1.0 / (1 + math.exp(-(t % period - p * period)))

def fishing(parrotfish, f):
	return f/(1+math.exp(-steepness*(parrotfish-shift)))

def patch_system(X, t, closure_length, f, m, n, poaching): 
	P_influx = [0]*n
	C_influx = [0]*n
	M_influx = [0]*n
	P,C,M = X.reshape(3, n) 

	for i in range(n):
		for j in range(n):
			P_influx[i] += (kP[i][j]) *P[j]  
		#print(P_influx[i])
		dPs[i] = P_influx[i]+ s*P[i]*(1 - (P[i] / K(sigma,C[i]))) - fishing(P[i], f)*P[i] *(square_signal(t, closure_length, i, m, n, poaching))
		dCs[i] = (i_C + r*C[i])*(1-M[i]-C[i])*(1-alpha*M[i]) - d*C[i]
		dMs[i] = (i_M+gamma*M[i])*(1-M[i]-C[i])-g*M[i]*P[i]/(g*eta*M[i]+1) 
	
	#concatenate into 1D vector to pass to next step
	return np.concatenate((dPs, dCs, dMs), axis=0)

# some of these should go into a separate utility function file
#general purpose graphing function 
def graph_sol(closure_length, f, m, n, frac_nomove, coral_high, poaching):

	initialize_patch_model(n, frac_nomove)

	if (coral_high):
		IC_set = X2
	else:
		IC_set = X1
	plt.figure()
	sol = odeint(patch_system, IC_set, t, args = (closure_length, f/(1-m/n), m, n, poaching))
	print(sol[:, n:2*n])
	patch_num = [x for x in range(1,n+1)]
	for i in range(len(patch_num)):
		plt.plot(sol[:, n+i],  label = 'patch % d'% patch_num[i])
	plt.xlabel('time')
	plt.ylabel('abundance')
	plt.title('spatial model')
	plt.legend(loc=0)
	txt="parameters" + "\nfishing when open: " + str(f/(1-m/n)) + "\npercent closure: " + str(m/n) +"\nclosure time: " + str(closure_length)
	plt.figtext(0.3, .31, txt, wrap=True, fontsize=8)
	plt.show()

#plots C vs f 
def hysteresis_zone():

	closure_length = 10
	m = 1
	n = 10
	initialize_patch_model(n, 0.95)

	max_fishing = 100
	end_coral_low = np.empty(100)
	end_coral_high = np.empty(100)
	fshing = np.empty(100)
	for f in range(max_fishing):
		hi_sol = odeint(patch_system, X2, t, args = (closure_length,(f/(max_fishing))/(1-(m/n)),m, n, poaching))
		lo_sol = odeint(patch_system, X1, t, args = (closure_length, (f/(max_fishing))/(1-m/n), m, n, poaching))
		coral_at_end_low = lo_sol[999][n]
		coral_at_end_high = hi_sol[999][n]
		fshing[f] = f/max_fishing
		end_coral_low[f] = coral_at_end_low
		end_coral_high[f] = coral_at_end_high
	plt.figure()
	plt.plot(fshing, end_coral_low, label = 'algae starts high')
	plt.plot(fshing, end_coral_high, label = 'coral starts high')
	plt.xlabel('fishing')
	plt.ylabel('coral cover at end')
	plt.legend(loc=0)
	txt="parameters" + "\nsigma: " + str(sigma) + "\neta: " + str(eta) +"\nalpha: " + str(alpha)
	plt.figtext(0.2, 0.3, txt, wrap=True, fontsize=8)
	plt.show()



def heatmap():
	patches = 20

	#heatmap of period vs m -- a bit messy right now 
	coral_array_HI =  np.zeros(shape=(patches,patches))
	coral_array_LO =  np.zeros(shape=(patches,patches))
	coral_array_AVG =  np.zeros(shape=(patches,patches))
	period_array = np.empty(patches)
	m_array = np.empty(patches)
	fishin = 0.23
	n = 20

	frac_nomove = 1
	for period in range(1,patches+1):
		print(period)
		print(patches)
		for m in range(patches):

			initialize_patch_model(n, frac_nomove)
			print(n, m, sep = ' ')
			displacer = 1/(1-m/float(n))
			final_coral_LO = odeint(patch_system, X1, t, args = (period, displacer*float(fishin), m, n, 0))# full_output = 1)
			#graph_sol(period*4, displacer*float(fishin), m, n, 1, False, 0)
			#graph_sol(5, 0.45, 4, 10, 1, False, 0) -- INCONSISTENT
			#avg1 = 0
			avg2 = 0
			for year in range(999- (999 % (n*period)) - 2*(n*period), 999 - (999 % (n*period))):
				#avg1 += final_coral_HI[year][n]
				avg2 += final_coral_LO[year][n]
			#avg1 = avg1 / (2*(period) + 1)
			avg2 = avg2 / (2*(period*n) + 1)
			print(avg2)
			print("________________________")
			#coral_array_HI[period-1][m] = avg1
			coral_array_LO[period-1][m] = avg2
			#coral_array_AVG[period-1][m] = 0.5 * (avg1 + avg2)
			period_array[period-1] = period
			m_array[m] = m
			show_labels = False

	plt.title('heatmap', fontsize = 20)
	sb.heatmap(coral_array_LO, vmin = 0.0, vmax = 1.0, cmap = "viridis", annot = show_labels) #YlGnBu for original color scheme
	plt.ylabel('period', fontsize = 10) # x-axis label with fontsize 15
	plt.xlabel('number of closures', fontsize = 10) # y-axis label with fontsize 15
	plt.yticks(rotation=0)
	plt.show()


# just finds how long it takes to reach C_1 from C_2 with zero fishing
def time_to_eq(C_1, C_2):
	# solve ode over sufficient number of timesteps 
	# search through solution for point where coral is first greater than C_2
	# the index of this will be the time


#fast heatmap
#simplifies m over n to lowest term and does least number of patches necessary 
import numpy as np 

def fast_heatmap():
	patches = 10

	#heatmap of period vs m -- a bit messy right now 
	coral_array_HI =  np.zeros(shape=(patches,patches))
	coral_array_LO =  np.zeros(shape=(patches,patches))
	coral_array_AVG =  np.zeros(shape=(patches,patches))
	period_array = np.empty(patches)
	m_array = np.empty(patches)
	fishin = 0.35
	n = 10

	frac_nomove = 1
	for period in range(1,patches+1):
		for m in range(patches):
			GCD = np.gcd(m, n)
			m = m / GCD
			n = n / GCD
			initialize_patch_model(n, frac_nomove)
			displacer = 1/(1-m/float(n))
			final_coral_LO = odeint(patch_system, X1, t, args = (period*2, displacer*float(fishin), m, n, 0))
			avg2 = 0
			for year in range(999- (999 % (n*period*2)) - 2*(n*period*2), 999 - (999 % (n*period*2))):
				avg2 += final_coral_LO[year][n]
			avg2 = avg2 / (2*(period*n*2) + 1)
			print(avg2)
			print("________________________")
			coral_array_LO[period-1][m] = avg2
			period_array[period-1] = period
			m_array[m] = m
			show_labels = False

	plt.title('heatmap', fontsize = 20)
	sb.heatmap(coral_array_LO, vmin = 0.0, vmax = 1.0, cmap = "viridis", annot = show_labels) #YlGnBu for original color scheme
	plt.ylabel('period', fontsize = 10) # x-axis label with fontsize 15
	plt.xlabel('number of closures', fontsize = 10) # y-axis label with fontsize 15
	plt.yticks(rotation=0)
	plt.show()

