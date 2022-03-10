import numpy as np 
import scipy
from scipy.integrate import odeint 
import matplotlib.pyplot as plt
import seaborn as sb
import math

####### parameters and initial conditions #############
#baseline reference points (from R code)
C_max = .509  # coral cover with no fishing
P_max = 20  # parrotfish with no fishing
M_max = .466    # algal cover with really high fishing

P0 = 0.1

P0 = P_max*.5
C0H = C_max*.95
M0H = M_max*.8

C0L = C_max*.05
M0L = M_max*.01

#parameters
r = 0.3
i_C = 0.05
i_M = 0.05
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

#both corals start in the low state
X_opp = [0.1, C0H, M0L, 0.1, C0L, M0H]

#one region starts with high coral, one starts low 
X_bothlow = [0.1, C0L, M0H, 0.1, C0L, M0H]
#######################################################


#################### functions ########################
def K(sigma, C):
	return (1-sigma)+sigma*C


#signal functions -- same as original 
def square_signal(t, period, p):
	if (((np.float128(steps) / np.float128(N)) * t) % ((np.float128(steps)/np.float128(N))*period))/(np.float128(steps)*period/np.float128(N)) < p:
		print("region closed, fishing off")
		return 0
	else:
		print("region open, fishing on")
		return 1

def sigmoid_signal(t, period, p):
	#this one is still buggy
	if period == 0:
		return 0
	else:
		return 1.0 / (1 + math.exp(-(t % period - p * period)))

#system evolution
def deriv(X, t, period, f, p, frac_nomove): 

	#dispersal vector -- just adding / subtracting these values at each timestep
	#with assumption that if patch 1 loses k% of its fish to patch 2, patch 2 gains (patch 1 cover / patch 2 cover)*k in % fish
	#ie if patch 1 has C1 = 0.7 and loses 1% to patch 2 with C2 = 0.3, we end up with C1 = 0.693 and C2 = (0.7/0.3)*0.01 +0.3 = 0.323
	#above logic might not work 

	frac_dispersed = (1-frac_nomove)*(1/2)
	k = [frac_dispersed, frac_dispersed]

	P1 = X[0] 
	C1 = X[1]
	M1 = X[2]
	P2 = X[3] 
	C2 = X[4]
	M2 = X[5]

	P_deriv1 = s*P1*(1 - (P1 / K(sigma,C1))) - f*P1 *square_signal(t, period, p) - k[0]*P1 + k[1]*P2#/P1
	C_deriv1 = (i_C + r*C1)*(1-M1-C1)*(1-alpha*M1) - d*C1
	M_deriv1 = (i_M+gamma*M1)*(1-M1-C1)-g*M1*P1/(g*eta*M1+1)
	P_deriv2 = s*P2*(1 - (P2 / K(sigma,C2))) - f*P2 *(1-square_signal(t, period, p)) - k[1]*P2 + k[0]*P1#/P2
	C_deriv2 = (i_C + r*C2)*(1-M2-C2)*(1-alpha*M2) - d*C2 
	M_deriv2 = (i_M+gamma*M2)*(1-M2-C2)-g*M2*P2/(g*eta*M2+1) 

	return [P_deriv1, C_deriv1, M_deriv1, P_deriv2, C_deriv2, M_deriv2]

#general grapher function 
def graph_sol(period, f, p, diff_cond, frac_nomove):
	if (diff_cond):
		IC_set = X_opp
	else:
		IC_set = X_bothlow

	sol = odeint(deriv, IC_set, t, args = (period, f, p, frac_nomove))
	patch1 = sol[:, 1]
	patch2 = sol[:, 4]
	print(sol)
	plt.figure()
	#plt.plot(M_sol, C_sol, label = 'phase trajectory')
	#plt.show()
	plt.plot(t, patch1, label='Coral Cover, Patch 1', color = 'turquoise')
	plt.plot(t, patch2, label='Coral Cover, Patch 2', color = 'purple')
	plt.xlabel('Time')
	plt.ylabel('% cover')
	plt.title('Two Patch Model')
	plt.legend(loc=0)
	burst_fishing = f 
	txt="parameters" + "\nfishing when open: " + str(burst_fishing) + "\npercent closure: " + str(p) +"\nperiod: " + str(period) + "\nfrac no move: "+str(frac_nomove)
	plt.figtext(0.2, 0.2, txt, wrap=True, fontsize=8)
	plt.show()
 

 ###############################################################




############################ plots #############################
'''
#one fully closed, the other fully open, no movement 
graph_sol(100, 1, 0, False, 1)

#each closed for same amount of time, 200 yr periods, no movement
graph_sol(200, 0.3/(1-0.5), 0.5, False, 1)

#each closed for same amount of time, 50 yr periods but f is only 0.2; still no movement 
graph_sol(50, 0.2/(1-0.5), 0.5, False, 1)

#now we vary the dispersal
graph_sol(200, 0.3/(1-0.5), 0.5, False, 0.95)

#the dispersal appears to make the system more resilient, so where is its new collapse threshold?
graph_sol(200, 0.5/(1-0.5), 0.5, False, 0.95)
'''
#below are examples where periodic closures save both starting from low coral values 

graph_sol(100, 0.26/(1-0.5), 0.5, False, 0.95)
graph_sol(100, 0.25/(1-0), 0, False, 0.95)

graph_sol(100, 0.22/(1-0), 0, False, 1)
graph_sol(100, 0.22/(1-0.5), 0.5, False, 1)
#note that this only works because p = 0.5 is in the sweet spot of our heatmap for f = 0.22
#what does this suggest about the behavior of multi patch models for which p = m/n is or is not in the sweet spot? 

'''
#we can also explore temporal asymmetry like below but this will not be relevant in spatial model 
#now we close one for 30 / 70 while the other is closed 70 / 30
graph_sol(100, 0.3/(1-0.3), 0.3, False)

#when the closures are not symmetric in time, we only get simultaneous recovery if (p, period) is in the light
#blue part of the heatmap for BOTH regions -- which cannot occur for some values of f 
graph_sol(100, 0.4/(1-0.6),0.6, False)
graph_sol(100, 0.4/(1-0.15), 0.15, False)

#f = 0.3 with 25 percent closure 
graph_sol(100, 0.3/(1-0.25), 0.25, False)100
'''



#final coral versus fishing effort in opp vs same scenarios
p =  0
period = 20

max_fishing = 100
end_coral_low = np.empty(100)
end_coral_high = np.empty(100)
fshing = np.empty(100)
for f in range(max_fishing):
	hi_sol = odeint(deriv, X_opp, t, args = (period,(np.float128(f)/(np.float128(max_fishing))/(1-p)),p, 0.8))
	lo_sol = odeint(deriv, X_bothlow, t, args = (period,(np.float128(f)/(np.float128(max_fishing))/(1-p)),p, 0.8))
	coral_at_end_low = .5*(lo_sol[999][1]+lo_sol[999][4])
	coral_at_end_high = .5*(hi_sol[999][1]+hi_sol[999][4])
	fshing[f] = np.float128(f) / np.float128(max_fishing)
	end_coral_low[f] = coral_at_end_low
	end_coral_high[f] = coral_at_end_high
plt.figure()
plt.plot(fshing, end_coral_low, label = 'both start low')
plt.plot(fshing, end_coral_high, label = 'one low, one high')
plt.xlabel('fishing')
plt.ylabel('coral cover at end')
plt.legend(loc=0)
txt="parameters" + "\nsigma: " + str(sigma) + "\neta: " + str(eta) +"\nalpha: " + str(alpha)
plt.figtext(0.2, 0.3, txt, wrap=True, fontsize=8)
plt.show()


frac_nomove = 0.95
'''
#plot final coral cover versus percentage time closed 

sim_period = 50
percentages = np.empty(100)
coral_covers = np.empty(100)
is_coral_high = True

for frac_closed in range(1, 100):
	fishing = 0.6
	frac = frac_closed
	frac *= 0.01
	fishing = fishing / (1.0 - frac)
	hi_sol = odeint(deriv, X_opp, t, args = (sim_period,fishing,frac, frac_nomove))
	lo_sol = odeint(deriv, X_bothlow, t, args = (sim_period,fishing,frac, frac_nomove))
	percentages[frac_closed] = frac
	avg = 0.0
	if is_coral_high:
		soln = hi_sol
	else:
		soln = lo_sol
	for year in range(999- (999 % sim_period) - 2*sim_period, 999 - (999 % sim_period)):
		avg += soln[year][1]
	avg = avg / (2*sim_period + 1)
	coral_covers[frac_closed] = avg 

plt.figure()
plt.plot(percentages, coral_covers, label = 'avg coral over whole region')
plt.xlabel('percentage time closed')
plt.ylabel('coral cover at end')
plt.legend(loc=0)
plt.show()
'''


#heatmap -- code is a lil wonky sorry
M = 20 #mesh fineness 
coral_array_OPP = np.random.rand(M,M)
coral_array_LOW = np.random.rand(M,M)
per_array = np.empty(M)
fishing_array = np.empty(M)
p = 0.5 
for tau in range(0,M):
	for fishing in range(0,M):
		index = tau
		TAU = tau*20+1
		displacer = 1/(1-np.float128(p)/np.float128(M))
		#final_coral_OPP = odeint(deriv, X_opp, t, args = (TAU, displacer*float(fishing) ,float(p) / float(M)), frac_nomove)# full_output = 1)
		final_coral_LOW = odeint(deriv, X_bothlow, t, args = (TAU, 2*displacer*float(fishing)/float(M), p, frac_nomove))# full_output = 1)
		avg2 = 0

		#truncate incomplete period and average over two periods 
		for year in range(999- (999 % (TAU)) - 2*(TAU), 999 - (999 % (TAU))):
			#avg1 += final_coral_OPP[year][1]+final_coral_OPP[year][4])/2
			avg2 += (final_coral_LOW[year][1]+final_coral_LOW[year][4])/2
		#avg1 = avg1 / (2*(TAU) + 1)
		avg2 = avg2 / (2*(TAU) + 1)

		#coral_array_OPP[index][p] = avg1
		coral_array_LOW[index][fishing] = avg2
		per_array[index] = np.float128(TAU) / np.float128(M)
		fishing_array[fishing] = np.float128(fishing)/ np.float128(M)
		show_labels = False


#plt.title('Coral Start Opposite', fontsize = 20)
#sb.heatmap(coral_array_OPP, vmin = 0.0, vmax = 1.0, cmap = "mako", annot = show_labels)
#plt.show()
plt.title('both patches start low', fontsize = 20)
sb.heatmap(coral_array_LOW, vmin = 0.0, vmax = 1.0, cmap = "mako", annot = show_labels) #YlGnBu for original color scheme
plt.ylabel('period / 20', fontsize = 10)
plt.xlabel('fishing (scaled)', fontsize = 10)
plt.show()

