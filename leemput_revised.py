
import numpy as np 
import scipy.integrate 
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import seaborn as sb
import math

#parameters
r = 0.3 #coral growth rate
i_C = 0.05 #coral recruitment
i_M = 0.05 #algae recruitment
gamma = 0.8 #algae growth rate
d = 0.1 #coral death rate
g = 1 #grazing rate
s = 1 #parrotfish growth rate
sigma = .5 #strength of coral-herbivore feedback
eta = 2 #strength of algae-herbivore feedback
alpha = 0.5 #strength of algae-coral feedback 

#time
N = 1000 #years
steps = 1000 #timesteps per year
t = np.linspace(0, N, steps) #time array

#initial conditions
#baseline reference points (from R code)
C_max = .509  # coral cover with no fishing
P_max = 20  # parrotfish with no fishing
M_max = .466    # algal cover with really high fishing - note this is Mi only

P0 = P_max*.5
C0H = C_max*.95
M0H = M_max*.8

C0L = C_max*.05
M0L = M_max*.01
'''
CHI = 0.5
MHI = 0.5

MLO = 0.02
CLO = 0.02
'''

P0 = 0.1
XHI = [P0, C0H, M0L]
XLO = [P0, C0L, M0H]


#parameters of interest
p = 0.5
f = 0.5
period = 100

#estimated equilibria values
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
	return 1.0 / (1 + math.exp(-(t % period - p * period)))

''''
#plot sigmoid
plt.figure()
for time in t:
	plt.scatter(time, sigmoid_signal(time, period, p))
plt.show()
'''

def system(t, X, period, f, p): 

	#set state variables to state vec
	P, C, M = X


	#variable for implementing one time closures
	#when burst term is included, we do one cycle
	#then return to the original net fishing value
	BURST = 1
	if t>period:
		BURST = 0

	#floating point error band-aids
	if P<0.0:
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

	#ODE system
	P_deriv = s*P*(1 - (P / K(sigma,C))) - f*P *square_signal(t, period, p) #*BURST - (1-p)*f*P*(1-BURST)
	C_deriv = (i_C + r*C)*(1-M-C)*(1-alpha*M) - d*C
	M_deriv = (i_M+gamma*M)*(1-M-C)-g*M*P/(g*eta*M+1)

	#looking under the hood
	'''
	if C_deriv < slope_thresh and t < period:
		print("time: ", t, "dC/dt:", C_deriv)
		print("threshold: ", slope_thresh)
	'''
	return [P_deriv, C_deriv, M_deriv]

def return_solution(t, period, f, p, coral_high):
	if (coral_high):
		IC_set = XHI
	else:
		IC_set = XLO
	sol = scipy.integrate.odeint(system, IC_set, t, args = (period, f, p))
	t = np.linspace(0, len(sol.y[0, :]), len(sol.y[0,:]))
	P_sol = sol.y[0, :]
	C_sol = sol.y[1, :]
	M_sol = sol.y[2, :]
	return sol.y


'''
fig = plt.figure(figsize = (9, 4))
Tslider_ax = plt.axes([0.1, 0.05, 0.8, 0.05])
fslider_ax = plt.axes([0.1, 0.15, 0.8, 0.05])
pslider_ax = plt.axes([0.1, 0.25, 0.8, 0.05])
ax = fig.add_subplot(111)
[graph] = ax.plot(t, return_solution(t, 100, 0, 0, False))
Tslider = Slider(Tslider_ax, 'period', 1, 500, valinit = 50)
fslider = Slider(fslider_ax, 'f', 0, 1.5, valinit = 0.26)
pslider = Slider(pslider_ax, 'p', 0, 1, valinit = 0.5)

def update_graph(val):
	line.set_ydata(return_solution(t, Tslider.val, fslider.val, pslider.val, False))
	fig.canvas.draw_idle()

fslider.on_changed(update_graph)
Tslider.on_changed(update_graph)
pslider.on_changed(update_graph)
plt.show()
'''

def graph_sol(period, f, p, coral_high):
	if (coral_high):
		IC_set = XHI
	else:
		IC_set = XLO
	sol = scipy.integrate.solve_ivp(system, [0,N], IC_set, method = "BDF", args = (period, f, p))
	t = np.linspace(0, len(sol.y[0, :]), len(sol.y[0,:]))
	P_sol = sol.y[0, :]
	C_sol = sol.y[1, :]
	M_sol = sol.y[2, :]
	plt.figure()
	plt.plot(t, C_sol, label='Coral')
	plt.plot(t, M_sol, label='Macroalgae')
	plt.plot(t, P_sol, label = 'Parrotfish')
	plt.plot(t, 1-M_sol-C_sol, label = 'Algal Turfs')
	plt.xlabel('Time')
	plt.ylabel('Abundances')
	plt.title('van de Leemput Model')
	plt.legend(loc=0)
	txt="parameters" + "\nfishing: " + str(f) + "\npercent closure: " + str(p) +"\nperiod: " + str(period)
	plt.figtext(0.2, 0.2, txt, wrap=True, fontsize=8)
	plt.show()


'''
#testing for hysteresis over a range of values for fishing 
max_fishing = 100
end_coral_low = np.empty(100)
end_coral_high = np.empty(100)
fshing = np.empty(100)
for f in range(max_fishing):
	hi_sol = scipy.integrate.solve_ivp(system, [0,10000], XHI, method = "BDF", args = (period,float(f)/(float(max_fishing)),p))
	lo_sol = scipy.integrate.solve_ivp(system, [0,10000], XLO, method = "BDF", args = (period,float(f)/(float(max_fishing)),p))
	coral_at_end_low = lo_sol.y[1][len(lo_sol)-1]
	coral_at_end_high = hi_sol.y[1][len(hi_sol)-1]
	fshing[f] = float(f) / float(max_fishing)
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

'''
#time series plots 
graph_sol(200, 0.4/(1-0.5), 0.5, False)
graph_sol(40, 0.4/(1-0.5), 0.5, True)
graph_sol(40, 0.3/(1-0.5), 0.5, True)

'''




























'''
def holling2_signal(t, period, p):
	return (t % period - p*period) / (1 + t % period - p*period)
'''