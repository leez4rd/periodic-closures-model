import numpy as np
from scipy.integrate import odeint  
import matplotlib.pyplot as plt 
import math

class Model:
	
	# variable arguments based on model type? 
	# need a function to tinker with individual parameters, and a function to list them 
	def __init__(self, model_type, n, frac_nomove): 
		self.model_type = model_type
		self.n = n 
		self.frac_nomove = frac_nomove

		self.f = 0
		self.closure_length = 0
		self.m = 0 
		self.poaching = 0


		'''
		#time -- maybe move this down to run model for dynamic stuff, same with ICs
		self.yrs = yrs #total amount of time
		self.t = np.linspace(0, yrs, yrs) #timestep array -- same number of timesteps as years 

		#initial conditions
		self.ICs = ICs

		initialize_patch_model()
		'''

	def initialize_patch_model(self, P0, C0L, C0H, M0L, M0H):
		# do variables defined in this function only exist in the scope of this function if it is called by constructor?

		frac_dispersed = (1-self.frac_nomove)*(1/(self.n)) # fraction of fish that disperse to other patches symmetrically

		# transition matrix for dispersal: element [i,j] of kP describes influx of P from j to i
		kP = np.empty((self.n,self.n)) 
		for i in range(self.n):
			for j in range(self.n):
				kP[i][j] = frac_dispersed
				if i == j:
					kP[i][j] = -frac_dispersed*(self.n - 1)

		setattr(self,'kP', kP)

		# P_influx = np.empty(n) -- I think this is unnecessary 
		setattr(self, 'P', np.empty(self.n))
		setattr(self, 'C' , np.empty(self.n))
		setattr(self,'M', np.empty(self.n))
		setattr(self,'dPs', np.empty(self.n))
		setattr(self,'dCs', np.empty(self.n))
		setattr(self,'dMs', np.empty(self.n))
		setattr(self, 'state_vec', [self.P,self.C,self.M])
		# concatenate initial condition arrays 
		# need to define baseline values somewhere
		# should these be attributes ?  
		setattr(self, 'X1', [P0]*self.n + [C0L]*self.n + [M0H]*self.n)
		setattr(self, 'X2', [P0]*self.n + [C0H]*self.n + [M0L]*self.n)

	# rass briggs model is only one with different state variables, so just have a separate init function 
	# better way might be to give initialize... an optional argument 
	def RB_initialize_patch_model(n, frac_nomove):

		frac_dispersed = (1-frac_nomove)*(1/(n)) # fraction of fish that disperse to other patches symmetrically
		# transition matrix for dispersal: element [i,j] of kP describes influx of P from j to i
		kP = np.empty((n,n))
		for i in range(n):
			for j in range(n):
				kP[i][j] = frac_dispersed
				if i == j:
					kP[i][j] = -frac_dispersed*(n - 1)
		setattr(self,'kP', kP)

		# P_influx = np.empty(n) -- I think this is unnecessary 
		setattr(self, 'P', np.empty(self.n))
		setattr(self, 'C' , np.empty(self.n))
		setattr(self,'Mi', np.empty(self.n))
		setattr(self,'Mv', np.empty(self.n))
		setattr(self,'dPs', np.empty(self.n))
		setattr(self,'dCs', np.empty(self.n))
		setattr(self,'dMis', np.empty(self.n))
		setattr(self,'dMvs', np.empty(self.n))
		setattr(self, 'state_vec', [P,C,Mv,Mi])
		# instead of global variables, can we pass these directly to patch_system? 
		# concatenate initial condition arrays -- amend to match state vec
		# X1 = [P0]*n + [C0L]*n + [M0H]*n 
		# X2 = [P0]*n + [C0H]*n + [M0L]*n


	# this defines the system conditionally based on user input -- this is the function we pass to odeint method
	# will class formatting mess up odeint?  
	# maybe put this outside class with Model as an arg? 

	'''
	def patch_system(self):

		P_influx = [0]*n

		for i in range(n):
			for j in range(n):
				P_influx[i] += (kP[i][j]) *P[j]  
	
			# better way to do this ?

			if self.model_type == 'RB':
				X = [self.P, self.C, self.Mv, Mi]
				return rass_briggs(X, i)
			elif self.model_type == 'BM':
				X = [self.P, self.C, self.M]
				return blackwood(X, i)
			elif self.model_type == 'vdL_PC':
				X = [self.P, self.C, self.M]
				return leemput(X, i)

			elif self.model_type == 'vdL_MP':
				X = [self.P, self.C, self.M]
				return leemput(X, i)
				
			elif self.model_type == 'vdL_MC':
				X = [self.P, self.C, self.M]
				return leemput(X, i)
				
			elif self.model_type == 'vdL': # all feedbacks active 
				X = [self.P, self.C, self.M]
				return leemput(X, i)
				
			else:
				print("Bad input, defaulting to Blackwood-Mumby!")
				X = [self.P, self.C, self.M]
				return blackwood(X, i)
	'''

	# returns the model run for a certain set of parameters 
	# maybe create a subclass for a model run to separate plotting jobs ? 
	def run_model(self, IC_set, t):
		# can we not pass self as an argument here ? 
		sol = odeint(patch_system, IC_set, t, args = (self, ))
		return sol 

	# make setting for show versus save (make these optional arguments) 
	# maybe add second way to use it by plugging in array from run_model 
	def time_series(self, save, show, IC_set, t, closure_length, f, m, n, poaching):
		if (coral_high):
			IC_set = X2
		else:
			IC_set = X1
		plt.figure()
		sol = odeint(self.patch_system(model_type), IC_set, t, args = (closure_length, f/(1-m/n), m, n, poaching))
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
		if save:
			plt.savefig('model_time_series1.jpg')
		if show:
			plt.show()

	# make a flag for fast or slow version 
	def coral_recovery_map(self):

		# slow version 

		coral_array =  np.zeros(shape=(self.n,self.n))
		
		period_array = np.empty(n)
		m_array = np.empty(n)
		
		for period in range(1,n+1):
			for m in range(n):
				displacer = 1/(1-m/float(n))
				final_coral = odeint(patch_system, X1, t, args = (period, displacer*float(self.f), m, self.n, 0))
				avg = 0
				for year in range(MAX_TIME - (MAX_TIME % (self.n*period)) - 2*(self.n*period), MAX_TIME - (MAX_TIME % (self.n*period))):
					avg += final_coral[year][n]
				avg = avg / (2*(period*self.n) + 1)
				coral_array[period-1][m] = avg2
				period_array[period-1] = period
				m_array[m] = m

		plt.title('coral end state', fontsize = 20)
		sb.heatmap(coral_array, vmin = 0.0, vmax = 1.0, cmap = "viridis", annot = show_labels) #YlGnBu for original color scheme
		plt.ylabel('period', fontsize = 10) # x-axis label with fontsize 15
		plt.xlabel('number of closures', fontsize = 10) # y-axis label with fontsize 15
		plt.yticks(rotation=0)
		plt.show()

		return None

	def bistable_zone():

		return None
	def find_unstable_equilibrium():
		return None
	def scenario_plot():

		return None 

	


	# modify this to take custom feedback parameters, or maybe custom anything? 
	# this doesn't work for some reason
	def load_parameters(self):
		if self.model_type == 'vdL':
			params = {
			"r": 0.3,
			"i_C" : 0.05,
			"i_M" : 0.05,
			"ext_C" : 0.0001,
			"ext_P" : 0.0001,
			"gamma" : 0.8,
			"d" : 0.1,
			"g" : 1,
			"s" : 1,
			"sigma" : .5, #strength of coral-herbivore feedback
			"eta" : 2, #strength of algae-herbivore feedback
			"alpha" : 0.5, #strength of algae-coral feedback 
			"P0" : 0.1,
			"C_HI" : .4,
			"M_LO" : .04,
			"C_LO" : .04,
			"M_HI" : .4
			}

		elif self.model_type == 'vdL_MC':
			params = {
			"r": 0.3,
			"i_C" : 0.05,
			"i_M" : 0.05,
			"ext_C" : 0.0001,
			"ext_P" : 0.0001,
			"gamma" : 0.8,
			"d" : 0.1,
			"g" : 1,
			"s" : 1,
			"sigma" : 0, #strength of coral-herbivore feedback
			"eta" : 0, #strength of algae-herbivore feedback
			"alpha" : 0.5, #strength of algae-coral feedback 
			"P0" : 0.1,
			"C_HI" : .4,
			"M_LO" : .04,
			"C_LO" : .04,
			"M_HI" : .4
			}

		elif self.model_type == 'vdL_MP':
			params = {
			"r": 0.3,
			"i_C" : 0.05,
			"i_M" : 0.05,
			"ext_C" : 0.0001,
			"ext_P" : 0.0001,
			"gamma" : 0.8,
			"d" : 0.1,
			"g" : 1,
			"s" : 1,
			"sigma" : 0, #strength of coral-herbivore feedback
			"eta" : 2, #strength of algae-herbivore feedback
			"alpha" : 0, #strength of algae-coral feedback 
			"P0" : 0.1,
			"C_HI" : .4,
			"M_LO" : .04,
			"C_LO" : .04,
			"M_HI" : .4
			}

		elif self.model_type == 'vdL_PC':
			params = {
			"r": 0.3,
			"i_C" : 0.05,
			"i_M" : 0.05,
			"ext_C" : 0.0001,
			"ext_P" : 0.0001,
			"gamma" : 0.8,
			"d" : 0.1,
			"g" : 1,
			"s" : 1,
			"sigma" : .5, #strength of coral-herbivore feedback
			"eta" : 0, #strength of algae-herbivore feedback
			"alpha" : 0, #strength of algae-coral feedback 
			"P0" : 0.1,
			"C_HI" : .4,
			"M_LO" : .04,
			"C_LO" : .04,
			"M_HI" : .4
			}

		elif self.model_type == 'BM':
			params = {
			"gamma" : 0.8,
			"beta" : 1,
			"alpha" : 1,
			"s" : 0.49,
			"r" : 1,
			"d" : 0.44,
			"a" : 0.1,
			"i_C" : 0.05,
			"i_M" : 0.05,
			"P0" : 0.1,
			"C_HI" : .4,
			"M_LO" : .04,
			"C_LO" : .04,
			"M_HI" : .4}

		elif self.model_type == 'RB':
			params = {
			"phiC" : 0.001, #open recruitment of coral
			"phiM" : 0.0001, #open recruitment of macroalgae

			"rM" : 0.5, #production of vulnerable macroalgae from invulnerable stage"

			#growth rates
			"gTC" : 0.1, #combined rate of growth and local recruitment of corals over free space 
			"gTV" : 0.2, #growth rate of vulnerable macroalgae over free space 
			"gTI" : 0.4, #growth rate of invulnerable macroalgae over free space 
			"rH" : 0.49,

			"gamma" : 0.4, #growth of macroalgae over coral vs free space
			"omega" : 2, #maturation rate of macroalgae from vulnerable to invulnerable class "

			#death rates
			"dC" : 0.5, #death rate of coral 
			"dI" : 0.4, #death rate of invulnerable macroalgae
			"dV" : 0.58, #death rate of vulnerable macroalgae per unit biomass of herbivores "

			"K" : 20,
			"Graze" : 0.58,

			#reference points  
			"C_max" : .509, # coral cover with no fishing
			"P_max" : 20, # parrotfish with no fishing
			"M_max" : .466, # algal cover with really high fishing - note this is Mi only
			"f_crit" : .275, # highest fishing at which coral wins (regardless of initial)
			"P0" : 0.1,
			"C_HI" : .4,
			"M_LO" : .04,
			"C_LO" : .04,
			"M_HI" : .4
			}

		for name, val in params.items():
			# initialize a new variable within class with that parameter's name and value 
			# there is probably a less hacky way to do this, but i aint got time to find it
			#exec(f"{name} = {val}", globals())
			setattr(self, name, val)

	# management parameter setter 
	def set_mgmt_params(self, closure_length, f, m, poaching):
		self.closure_length = closure_length
		self.f = f
		self.m = m
		self.poaching = poaching 

	# test accessor
	def get_M(self):
		return self.M



# puttin patch system outside the class so we don't have to pass self as first parameter 
# this way we can pass this function directly to odeint, and just call it from within the Model class 
# don't necessarily need outside class, just can't modify object...
# how could we bundle all of this into a library (ie import periodic_closures_modeltools) ? 

def patch_system(X, t, system_model):
		P_influx = [0]*system_model.n

		for i in range(system_model.n):
			for j in range(system_model.n):
				P_influx[i] += (system_model.kP[i][j]) * system_model.P[j]  
	
			# better way to do this ?

			if system_model.model_type == 'RB':
				
				return rass_briggs(X, i, system_model)
			elif system_model.model_type == 'BM':
				
				return blackwood(X, i, system_model)
			elif system_model.model_type == 'vdL_PC':
				
				return leemput(X, i, system_model)

			elif system_model.model_type == 'vdL_MP':
				
				return leemput(X, i, system_model)
				
			elif system_model.model_type == 'vdL_MC':
				
				return leemput(X, i, system_model)
				
			elif system_model.model_type == 'vdL': # all feedbacks active 
				
				return leemput(X, t, i, system_model, P_influx)
				
			else:
				print("Bad input, defaulting to Blackwood-Mumby!")
				return blackwood(X, i, system_model)

def rass_briggs(X, i, system_model):

		P, C, Mv, Mi = X.reshape(4, n)
		T = 1 - C - Mv - Mi 
		dPdt = rH*P*(1-P/K) - system_model.fishing(P[i], f)*P[i] *(system_model.square_signal(t, closure_length, i, m, n, poaching))
		dCdt = (phiC*T) + gTC*T*C - gamma*gTI*Mi*C - dC*C
		dMvdt = phiM*T + rM*T*Mi + gTV*T*Mv - dV*Mv - P*Mv*Graze - omega * Mv
		# conceptual question: why is that second term multiplied by M not Mv? 
		dMidt = omega*Mv + gTI*T*Mi + gamma*gTI*Mi*C - dI*Mi
		return [dPdt, dCdt, dMvdt, dMidt]

		# check input
		return None 

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

# idk
def fishing(parrotfish, f):
	steepness = 0.1
	shift = -10
	return f/(1+math.exp(-steepness*(parrotfish-shift)))

# is self keyword necessary here ? or can we call these functions regardless?
# is it necessary to send X as a parameter if it is embedded in the class? 
# should all of these be class methods which can change self.X? 
def blackwood(X, i, system_model):
	P, C, M = X.reshape(3, n)

	dPs[i] = s*P[i]*(1 - (P[i] / (beta*system_model.K(C[i])))) - system_model.fishing(P[i], f)*P[i]*system_model.square_signal(t, closure_length, i, m, n, poaching)
	dCs[i] = r*(1-M[i]-C[i])*C[i]-d*C[i] - a*M[i]*C[i] + i_C*(1-M[i]-C[i])
	# need to define g(P) before this model is used 
	dMs[i] = a*M[i]*C[i] - g(P[i])*M[i] *(1/(1-C[i])) + gamma*M[i]*(1-M[i]-C[i])+i_M*(1-M[i]-C[i])

	return np.concatenate((dPs, dCs, dMs), axis=0)

# will need to pass [self.P, self.C, self.M] to this for it to work 
def leemput(X, t, i, system_model, P_influx): # COPY THIS FORMAT FOR OTHER MODELS 
	# check input 
	P,C,M = X.reshape(3, system_model.n) # will reshaping work since we are passing arrays of length n? 
	system_model.dPs[i] = P_influx[i]+ system_model.s*P[i]*(1 - (P[i] / K(system_model.sigma,C[i]))) - fishing(P[i], system_model.f)*P[i] *(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching))
	system_model.dCs[i] = (system_model.i_C + system_model.r*C[i])*(1-M[i]-C[i])*(1-system_model.alpha*M[i]) - system_model.d*C[i]
	
	system_model.dMs[i] = (system_model.i_M+system_model.gamma*M[i])*(1-M[i]-C[i])-system_model.g*M[i]*P[i]/(system_model.g*system_model.eta*M[i]+1)

	#concatenate into 1D vector to pass to next step
	return np.concatenate((system_model.dPs, system_model.dCs, system_model.dMs), axis=0)


def main():

	
	yrs = 1000 #total amount of time
	t = np.linspace(0, yrs, yrs) #timestep array -- same number of timesteps as years 
	P0, C0L, C0H, M0L, M0H = 0.1, 0.04, 0.4, 0.04, 0.4

	x = Model('vdL', 2, 1)
	y = Model('vdL', 10, 0.97)

	print(x.model_type)
	print(y.model_type)

	x.load_parameters()
	y.load_parameters()

	x.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	y.initialize_patch_model(P0, C0L, C0H, M0L, M0H)

	print(y.gamma)
	print(y.get_M())
	x.set_mgmt_params(20, 0.25, 1, 0)
	print(x.run_model(x.X1, t)) # IT WORKED !!!!!!
	# print(x.run_model(parameter list))
	# x.time_series(parameter list, show = True)

	return 

if __name__ == '__main__':
	main()