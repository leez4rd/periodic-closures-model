import numpy as np 
import scipy
from scipy.integrate import odeint 
import matplotlib.pyplot as plt
import seaborn as sb
from mpl_toolkits.mplot3d import Axes3D
import math
#from decimal import *
from PyPDF2 import PdfFileMerger


#baseline reference points
C_max = .509  # coral cover with no fishing
P_max = 20  # parrotfish with no fishing
M_max = .466    # algal cover with really high fishing - note this is Mi only

P0 = 0.1
'''
C0H = C_max*.8
M0H = M_max*.8

C0L = C_max*.05
M0L = M_max*.05
'''
P0 = P_max*.5
C0H = C_max*.95
M0H = M_max*.8

C0L = C_max*.05
M0L = M_max*.01

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

#initial conditions
CHI = 0.5
MHI = 0.5

MLO = 0.02
CLO = 0.02

P0 = 0.1
#note: setting parrotfish to 0.99 and 0.01 might alter bistable zone
XHI = [0.1, C0H, M0L]
XLO = [0.1, C0L, M0H]

#unstable equilibrium!!
Xeq = [0.1, 0.383235864, 1-0.383235864]

#compute these individually by embedding hysteresis
#sim code in function that returns equilibria
C_loeq = 0.01
C_hieq = 0.7


def K(sigma, C):
	return (1-sigma)+sigma*C

def square_signal(t, period, p):
	#for integer times, (t % period) / period = p is the switch point 
	#to account for non-integer periods and times, multiply numerator and denominator by steps / N 
	#modified to make everything np.float128
	if (((np.float128(steps) / np.float128(N)) * t) % ((np.float128(steps)/np.float128(N))*period))/(np.float128(steps)*period/np.float128(N)) < p:
		return 0
	else:
		return 1
def sigmoid_signal(t, period, p):
	if period == 0:
		return 0
	else:
		return 1.0 / (1 + math.exp(-(t % period - p * period)))
#def f(parrotfish):
	#return 1/(1+math.exp(-(parrotfish-0.1)/200))

def deriv(X, t, period, f, p): 
	P = X[0] 
	C = X[1]
	M = X[2]


	#variable for closing once then returning to original net fishing 
	CLOSE_ONCE = 1
	if t>period:
		CLOSE_ONCE = 0

	#band-aids
	if P<0:
		P *= -1
	if P <= float(0.0000000001):
		P = float(0.0000000001)
	if P > K(sigma, C):
		P = K(sigma, C)
	if C > 1:
		C = 1

	#very rough approx for collapse threshold
	slope_thresh = 10000
	if (1-p)*period != 0:
		slope_thresh = (C_loeq - C_hieq)/((1-p)*(period))

	P_deriv = s*P*(1 - (P / K(sigma,C))) - f*P *square_signal(t, period, p) #*CLOSE_ONCE - (1-p)*f*P*(1-CLOSE_ONCE)
	C_deriv = (i_C + r*C)*(1-M-C)*(1-alpha*M) - d*C + ext_C
	M_deriv = (i_M+gamma*M)*(1-M-C)-g*M*P/(g*eta*M+1) + ext_C

	if C_deriv < slope_thresh and t < period:
		print("time: ", t, "dC/dt:", C_deriv)
		print("threshold: ", slope_thresh)

	return [P_deriv, C_deriv, M_deriv]


def graph_sol(period, f, p, coral_high):
	if (coral_high):
		IC_set = XHI
	else:
		IC_set = XLO

	sol = odeint(deriv, IC_set, t, args = (period, f, p))
	P_sol = sol[:, 0]
	C_sol = sol[:, 1]
	M_sol = sol[:, 2]
	plt.figure()
	#plt.plot(M_sol, C_sol, label = 'phase trajectory')
	#plt.show()
	plt.plot(t, C_sol, label='Coral')
	plt.plot(t, M_sol, label='Macroalgae')
	plt.plot(t, P_sol, label = 'Parrotfish')
	plt.plot(t, 1-M_sol-C_sol, label = 'Algal Turfs')
	plt.xlabel('Time')
	plt.ylabel('Abundances')
	plt.title('van de Leemput Model')
	plt.legend(loc=0)
	burst_fishing = f / (1-p)
	txt="parameters" + "\nfishing when open: " + str(burst_fishing) + "\npercent closure: " + str(p) +"\nperiod: " + str(period)
	plt.figtext(0.2, 0.2, txt, wrap=True, fontsize=8)
	plt.show()
'''
parr = np.linspace(0,100,2)
y = [f(p) for p in parr]
plt.plot(parr, y)
plt.show()
#graph_sol(200, 0, 0, False)
#50.476
'''
graph_sol(40, 0.3/(1-0.4), 0.4, False)

x  = np.linspace(0,1,1000)
fx = lambda x: 50.476/(1-x+0.0000001)
y  = [fx(val) for val in x] 
ax= plt.plot(x, y, color = "red")
plt.title("tau as a function of p at threshold")
plt.xlabel("fraction closed")
plt.ylabel("period in years")
plt.xlim([0,1])
plt.ylim([100*5+1,0])
plt.show()


#testing for hysteresis over a range of values for fishing 
p = 0
period = 20

max_fishing = 100
end_coral_low = np.empty(100)
end_coral_high = np.empty(100)
fshing = np.empty(100)
for f in range(max_fishing):
	hi_sol = odeint(deriv, XHI, t, args = (period,(np.float128(f)/(np.float128(max_fishing))/(1-p)),p))
	lo_sol = odeint(deriv, XLO, t, args = (period,(np.float128(f)/(np.float128(max_fishing))/(1-p)),p))
	coral_at_end_low = lo_sol[999][1]
	coral_at_end_high = hi_sol[999][1]
	fshing[f] = np.float128(f) / np.float128(max_fishing)
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

'''
#finding the unstable equilibrium for f = 0.26
coral_covers = np.linspace(0, 1, 100)
final_coral_covers = np.empty(len(coral_covers))
count = 0		
for initial_coral in coral_covers:
	sol = odeint(deriv, [0.1, initial_coral, 1-initial_coral], t, args = (100, 0.26, 0))
	final_coral_covers[count] = sol[999][1]
	count += 1
plt.plot(coral_covers, final_coral_covers)
plt.show()
#about 0.3785?
'''


graph_sol(100, 0, 0, False)
graph_sol(50, 0.26/(1-0.7),0.7, False)
#HYSTERESIS BOUNDS - 0.232 to 0.482, approx

#analogous plots changing period instead -- note that
#f/(1-p) > fcrit
'''
graph_sol(40, 0.26/(1-0.52), 0.52, False)
graph_sol(80, 0.26/(1-0.52), 0.52, False)
graph_sol(120, 0.26/(1-0.52), 0.52, False)
graph_sol(200, 0.26/(1-0.52), 0.52, False)
'''

#these graphs best illustrate the behavior 
#as we increase p
graph_sol(150, 0.30/(1-0.20), 0.20, False)
graph_sol(150, 0.30/(1-0.25), 0.25, False)
graph_sol(150, 0.30/(1-0.30), 0.30, False)
graph_sol(150, 0.30/(1-0.35), 0.35, False)
graph_sol(150, 0.30/(1-0.40), 0.40, False)
graph_sol(150, 0.30/(1-0.82), 0.82, False)




'''
*OBSERVATIONS*

FIG 1
fishing when open puts us in the bistable range, 
but coral still gains more than it loses and eventually wins

FIG 2
stuck in a low cycle despite high start; coral doesn't have time to go high

FIG 3
winning on a fast cycle with a low start

FIG 4
which equilibrium coral goes to in long periods depends only on f/(1-p)
because all populations have time for dynamics to play out 
e.g. if period is 500 years, then all that matters for whether system collapses is
whether displacement fishing takes you past the hysteresis zone 
an empirical measure of this might be the ratio of coral to parrotfish decay rates (scaled appropriately)
as that ratio will be closer to one if the system is closer to collapse 

FIG 5 
comparing increase in p at 200 yrs; these only make sense if our hysteresis region is incomplete 


FOR FIG 3
closures are too short to allow coral recovery; f = 0.6 collapses system, f = 0.45 preserves it;
the coral is losing more than it is gaining up to the equilibrium, with that difference lessening with every passing period 

FIG 6
very fast cycling, coral starts high and remains high despite large fishing

FIG 7
even though the period is a generous 100 years, this isn't quite long enough to undo the plummet resulting
from fishing effort of 0.36 

FIG 8 
nice graph

'''
'''
graph_sol(20, 0, 0, False)
#FIG 1
#displaced fishing in hysteresis zone but coral still wins
graph_sol(25, 0.20/(1-0.3), 0.3, False)

#FIG 2
#stuck in a low cycle despite high start 
graph_sol(75, 0.5/(1-0.3), 0.3, True)

#FIG 3 - winning on a fast cycle with low start
graph_sol(5, 0.19/(1-0.5), 0.5, False)

#FIG 4i and 4ii - super long period
graph_sol(500, 0.26/(1-0.5), 0.5, False)
graph_sol(500, 0.22/(1-0.5), 0.5, False)

#FIG 5(i, ii, iii) - comparing three scenarios and only varying p
graph_sol(200, 0.27, 0, False)
graph_sol(200, 0.27/(1-0.6), 0.6, False)
graph_sol(200, 0.27/(1-0.3), 0.3, False)

#FIG 6 - very fast cycling, coral starts high and remains high despite large fishing
graph_sol(10, 0.22/(1-0.9), 0.9, True)

#FIG 7 - not open for long enough to go all the way down, but it wants to 
graph_sol(100, 0.26/(1-0.6), 0.6, False)

#FIG 8 - fifty year cycle -- good graph
graph_sol(50, 0.22/(1-0.5), 0.5, False)
''' 




#plot final coral cover versus percentage time closed 
sim_period = 200
percentages = np.empty(100)
coral_covers = np.empty(100)
is_coral_high = False

for frac_closed in range(1, 100):
	fishing = 0.35
	frac = frac_closed
	frac *= 0.01
	fishing = fishing / (1.0 - frac)
	hi_sol = odeint(deriv, XHI, t, args = (sim_period,fishing,frac))
	lo_sol = odeint(deriv, XLO, t, args = (sim_period,fishing,frac))
	percentages[frac_closed] = frac
	avg = 0.0
	if is_coral_high:
		soln = hi_sol
	else:
		soln = lo_sol
	for year in range(999- (999 % sim_period) - 2*sim_period, 999 - (999 % sim_period)):
		avg += soln[year][1]
	avg = avg / (2*sim_period + 1)
	coral_covers[frac_closed] = avg #lo_sol[999 - (999 % (period))][1]

	#graphing intermediate steps
	'''
	if frac_closed % 30 == 0:
		graph_sol(sim_period, fishing, frac, is_coral_high)
	'''
	
plt.figure()
plt.plot(percentages, coral_covers, label = 'coral starts low')
plt.xlabel('percentage time closed')
plt.ylabel('coral cover at end')
plt.legend(loc=0)
plt.show()

'''
#coral versus period 
periods = np.empty(1000)
coral_covers = np.empty(1000)
is_coral_high = False
frac = 0.7
fishing = 0.26
fishing = fishing / (1.0 - frac)
#net reduction in fishing 
#high_sol = odeint(deriv, XHI, t, args = (100, 0.5, 0))
#low_sol = odeint(deriv, XLO, t, args = (100, 0.5, 0))

for sim_period in range(1, 1000):
	hi_sol = odeint(deriv, XHI, t, args = (sim_period,fishing,frac))
	lo_sol = odeint(deriv, XLO, t, args = (sim_period,fishing,frac))
	periods[sim_period] = sim_period
	avg = 0.0
	if is_coral_high:
		soln = hi_sol
	else:
		soln = lo_sol
	for year in range(999- (999 % sim_period) - 2*sim_period, 999 - (999 % sim_period)):
		avg += soln[year][1]
	avg = avg / (2*sim_period + 1)
	coral_covers[sim_period] = avg #lo_sol[999 - (999 % (period))][1]

plt.figure()
plt.plot(periods, coral_covers, label = 'coral starts low')
plt.xlabel('period of closure')
plt.ylabel('coral cover at end')
plt.legend(loc=0)
plt.show()
'''


'''
#heatmap
n = 100 #mesh fineness 
m = 100
coral_array_HI = np.random.rand(n,n)
coral_array_LO = np.random.rand(n,n)
coral_array_AVG = np.random.rand(n,n)
per_array = np.empty(n)
p_array = np.empty(n)
pdata = []
taudata = []
fishing = 26
for tau 		in range(0,n):
		for p in range(0,n):
			index = tau
			TAU = tau*5+1
			displacer = 1/(1-np.float128(p)/np.float128(n))
			#final_coral_HI = odeint(deriv, XHI, t, args = (TAU, displacer*np.float128(fishing) / np.float128(n) ,np.float128(p) / np.float128(n)))#full_output = 1)
			final_coral_LO = odeint(deriv, XLO, t, args = (TAU, displacer*np.float128(fishing) / np.float128(n) ,np.float128(p) / np.float128(n)))# full_output = 1)
			#FCHI = final_coral_HI[999 - (999 % (5*tau+1))][1]
			#FCLO = final_coral_LO[999 - (999 % (5*tau+1))][1]
			#avg1 = 0
			avg2 = 0
			for year in range(999- (999 % (TAU)) - 2*(TAU), 999 - (999 % (TAU))):
				#avg1 += final_coral_HI[year][1]
				avg2 += final_coral_LO[year][1]
			#avg1 = avg1 / (2*(TAU) + 1)
			#avg2 = avg2 / (2*(TAU) + 1)

			if(avg2 > 0.6):
				plt.title("closure parameters for optimal coral cover")
				plt.xlabel("percent closed")
				plt.ylabel("period")
				pdata += [p]
				taudata += [TAU]
				plt.scatter(p, TAU)
			

			for i in range (tau, 100):
				FCHI += final_coral_HI[i][1]
				FCLO += final_coral_LO[i][1]
			FCHI = FCHI / np.float128(100)
			FCLO = FCLO / np.float128(100)
			#coral_array_HI[index][p] = avg1
			coral_array_LO[index][p] = avg2
			#coral_array_AVG[index][p] = 0.5 * (avg1 + avg2)
			per_array[index] = np.float128(TAU) / np.float128(n)
			p_array[p] = np.float128(p)/ np.float128(n)
			show_labels = False
#fit_function = scipy.optimize.curve_fit(guess_function, pdata, taudata)
#plt.plot(pdata, guess_function(pdata, *fit_function), 'g--', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(fit_function))
#plt.show()
#sb.heatmap(coral_array_AVG, vmin = 0.0, vmax = 1.0, cmap = "mako", annot = show_labels)

#sb.color_palette("viridis", as_cmap=True)
#plt.title('Heatmap of Coral Dominance', fontsize = 20) # title with fontsize 20
#plt.ylabel('period / 25', fontsize = 10) # x-axis label with fontsize 15
#plt.xlabel('percent time closed', fontsize = 10) # y-axis label with fontsize 15
#plt.yticks(rotation=0)
#plt.show()
#plt.title('Coral Starts High', fontsize = 20)
#sb.heatmap(coral_array_HI, vmin = 0.0, vmax = 1.0, cmap = "mako", annot = show_labels)
#plt.show()
plt.title('Coral Starts Low', fontsize = 20)
sb.heatmap(coral_array_LO, vmin = 0.0, vmax = 1.0, cmap = "mako", annot = show_labels) #YlGnBu for original color scheme
plt.ylabel('period / 5', fontsize = 10)
plt.xlabel('percent time closed', fontsize = 10)
plt.scatter(21, 38)
plt.scatter(8, 88)
plt.scatter(48, 66)
plt.scatter(72, 10)
plt.show()
'''

n = 25
#heatmap of f vs p
coral_array_HI = np.random.rand(n,n)
coral_array_LO = np.random.rand(n,n)
coral_array_AVG = np.random.rand(n,n)
f_array = np.empty(n)
p_array = np.empty(n)
period = 50
for fishing in range(1,n):
		for p in range(1,n):
			displacer = 1/(1-np.float128(p)/np.float128(n))
			#final_coral_HI = odeint(deriv, XHI, t, args = (period, displacer*np.float128(fishing) / np.float128(n) ,np.float128(p) / np.float128(n)))#full_output = 1)
			final_coral_LO = odeint(deriv, XLO, t, args = (period, displacer*np.float128(fishing) / np.float128(n) ,np.float128(p) / np.float128(n)))# full_output = 1)
			#FCHI = final_coral_HI[999 - (999 % (5*tau+1))][1]
			#FCLO = final_coral_LO[999 - (999 % (5*tau+1))][1]
			avg1 = 0
			avg2 = 0
			for year in range(999- (999 % (period)) - 2*(period), 999 - (999 % (period))):
				#avg1 += final_coral_HI[year][1]
				avg2 += final_coral_LO[year][1]
			avg1 = avg1 / (2*(period) + 1)
			avg2 = avg2 / (2*(period) + 1)
			'''
			for i in range (tau, 100):
				FCHI += final_coral_HI[i][1]
				FCLO += final_coral_LO[i][1]
			FCHI = FCHI / np.float128(100)
			FCLO = FCLO / np.float128(100)
			'''
			#coral_array_HI[fishing][p] = avg1
			coral_array_LO[fishing][p] = avg2
			coral_array_AVG[fishing][p] = 0.5 * (avg1 + avg2)
			f_array[fishing] = np.float128(fishing) / np.float128(n)
			p_array[p] = np.float128(p)/ np.float128(n)
			show_labels = False
sb.heatmap(coral_array_AVG, vmin = 0.0, vmax = 1.0, cmap = "viridis", annot = show_labels)

sb.color_palette("viridis", as_cmap=True)
plt.title('Heatmap of Coral Dominance', fontsize = 20) # title with fontsize 20
plt.ylabel('fishing', fontsize = 10) # x-axis label with fontsize 15
plt.xlabel('percent time closed', fontsize = 10) # y-axis label with fontsize 15
plt.yticks(rotation=0)
plt.show()
#plt.title('Coral Starts High', fontsize = 20)
#sb.heatmap(coral_array_HI, vmin = 0.0, vmax = 1.0, cmap = "viridis", annot = show_labels)
#plt.show()
plt.title('Coral Starts Low', fontsize = 20)
sb.heatmap(coral_array_LO, vmin = 0.0, vmax = 1.0, cmap = "viridis", annot = show_labels) #YlGnBu for original color scheme
plt.show()









#plt.savefig("g"+chr(tau)+".pdf")
#print("plot done...")


'''
#big simulation
#work in progress...
coral_array_HI = np.random.rand(n,n)
coral_array_LO = np.random.rand(n,n)
coral_array_AVG = np.random.rand(n,n)
fishing_array = np.empty(n)
p_array = np.empty(n)
for tau in range(n):
	for fishing in range(n):
		for p in range(n):
			final_coral_HI = odeint(deriv, XHI, t, args = (tau*10+1, np.float128(fishing) / np.float128(n) ,np.float128(p) / np.float128(n)))[999][1]
			final_coral_LO = odeint(deriv, XLO, t, args = (tau*10+1, np.float128(fishing) / np.float128(n) ,np.float128(p) / np.float128(n)))[999][1]
			coral_array_HI[fishing][p] = final_coral_HI
			coral_array_LO[fishing][p] = final_coral_LO
			coral_array_AVG[fishing][p] = 0.5 * (final_coral_LO + final_coral_HI)
			fishing_array[fishing] = fishing / n
			p_array[p] = p / n
			show_labels = False
	sb.heatmap(coral_array_AVG, vmin = 0.0, vmax = 1.0, cmap = "YlGnBu", annot = show_labels)
	sb.color_palette("viridis", as_cmap=True)
	plt.title('Heatmap of Coral Dominance', fontsize = 20) # title with fontsize 20
	plt.ylabel('Fishing Intensity (as %)', fontsize = 10) # x-axis label with fontsize 15
	plt.xlabel('Percentage Time Closed', fontsize = 10) # y-axis label with fontsize 15
	plt.yticks(rotation=0)
	g = ''
	g += chr(tau+1)
	plt.savefig(g+".pdf")
	print("plot done...")

#for contour plot
#fishing_array = np.reshape(fishing_array, 10)
#p_array = np.reshape(p_array, 10)

#plt.contour(fishing_array, p_array, coral_array)
#plt.show()
'''

'''
#heatmap (original)
show_labels = False 
sb.heatmap(coral_array, vmin = 0.0, vmax = 1.0, cmap="YlGnBu", annot = show_labels)
sb.color_palette("viridis", as_cmap=True)
plt.title('Heatmap of Coral Dominance', fontsize = 20) # title with fontsize 20
plt.ylabel('Fishing Intensity (as %)', fontsize = 10) # x-axis label with fontsize 15
plt.xlabel('Percentage Time Closed', fontsize = 10) # y-axis label with fontsize 15
plt.yticks(rotation=0)
plt.show()
'''

#heatmap (modifying to accent bistability region)

'''
pdfs = ['g1.pdf','g2.pdf','g3.pdf','g4.pdf','g5.pdf','g6.pdf','g7.pdf','g8.pdf','g9.pdf','g0.pdf']
merger = PdfFileMerger()
for pdf in pdfs:
	merger.append(pdf)

merger.write("BigSim.pdf")
merger.close()
'''


'''
more graphs
graph_sol(500, 0.35, 0.5, True)
graph_sol(500, 0.35, 0.5, False)
graph_sol(500, 0.45, 0.5, True)
graph_sol(500, 0.55, 0.5, False)
graph_sol(500, 0.55, 0.5, True)
graph_sol(500, 0.65, 0.5, True)
graph_sol(500, 0.75, 0.5, True)
graph_sol(500, 0.85, 0.5, True)
graph_sol(500, 0.95, 0.5, True)
graph_sol(500, 1.05, 0.5, True)
graph_sol(500, 1.15, 0.5, True)


graph_sol(100, 0, 0, True)
graph_sol(100, 0.6, 0, True)
graph_sol(100, 0.6/(1-0.5), 0.5, True)
graph_sol(100, 0.6/(1-0.75), 0.75, True)
graph_sol(4, 0.9/(1-0.5), 0.5, False)
'''
#187468

