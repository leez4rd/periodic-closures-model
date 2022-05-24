import numpy as np
from scipy.integrate import odeint  
import matplotlib.pyplot as plt 

class Model:
	
	# variable arguments based on model type? 
	# need a function to tinker with individual parameters, and a function to list them 
	def __init__(self, model_type, n, frac_nomove): 

		# load_parameters() -- hypothetical function to load all params of RB model into this object  
		if model_type == 'RB':
			self.load_parameters('RB')
			# self.RB_initialize_patch_model(n, frac_nomove)

		elif model_type == 'BM':
			# self.initialize_patch_model(n, frac_nomove)
			self.load_parameters('BM')
		elif model_type == 'vdL_PC':
			self.load_parameters('vdL_PC')
			# self.initialize_patch_model(n, frac_nomove)
		elif model_type == 'vdL_MP':
			self.load_parameters('vdL_MP')
			# self.initialize_patch_model(n, frac_nomove)
		elif model_type == 'vdL_MC':
			# self.initialize_patch_model(n, frac_nomove)
			self.load_parameters('vdL_MC')
		elif model_type == 'vdL': # all feedbacks active 
			self.initialize_patch_model(n, frac_nomove)
		else:
			print("Bad input, defaulting to Blackwood-Mumby!")
			self.initialize_patch_model(n, frac_nomove)
			# self.load_parameters('BM')

	
	def initialize_patch_model(self, n, frac_nomove):
		# do variables defined in this function only exist in the scope of this function if it is called by constructor?

		frac_dispersed = (1-frac_nomove)*(1/(n)) # fraction of fish that disperse to other patches symmetrically

		# transition matrix for dispersal: element [i,j] of kP describes influx of P from j to i
		kP = np.empty((n,n)) 
		for i in range(n):
			for j in range(n):
				kP[i][j] = frac_dispersed
				if i == j:
					kP[i][j] = -frac_dispersed*(n - 1)

		# instead of global variables here, can we initialize these within patch system? or pass them as parameters? 

		P_influx = np.empty(n)
		P = np.empty(n) 
		C = np.empty(n) 
		M = np.empty(n)
		dPs = np.empty(n)
		dCs = np.empty(n)
		dMs = np.empty(n)

		#concatenate initial condition arrays 
		# need to define baseline values somewhere 
		# X1 = [P0]*n + [C0L]*n + [M0H]*n
		# X2 = [P0]*n + [C0H]*n + [M0L]*n

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


		# P_influx = np.empty(n) -- I think this is unnecessary 
		P = np.empty(n) 
		C = np.empty(n) 
		Mi = np.empty(n)
		Mv = np.empty(n)
		dPs = np.empty(n)
		dCs = np.empty(n)
		dMis = np.empty(n)
		dMvs = np.empty(n)


		# instead of global variables, can we pass these directly to patch_system? 
		# concatenate initial condition arrays -- amend to match state vec
		# X1 = [P0]*n + [C0L]*n + [M0H]*n 
		# X2 = [P0]*n + [C0H]*n + [M0L]*n


	# this defines the system conditionally based on user input -- this is the function we pass to odeint method
	def patch_system(X, model_type):
		P_influx = [0]*n

		for i in range(n):
			for j in range(n):
				P_influx[i] += (kP[i][j]) *P[j]  
	
			# better way to do this ?

			if model_type == 'RB':
				return rass_briggs(X, i)
			elif model_type == 'BM':
				return blackwood(X, i)

			elif model_type == 'vdL_PC':
				return leemput(X, i)

			elif model_type == 'vdL_MP':
				return leemput(X, i)
				
			elif model_type == 'vdL_MC':
				return leemput(X, i)
				
			elif model_type == 'vdL': # all feedbacks active 
				return leemput(X, i)
				
			else:
				print("Bad input, defaulting to Blackwood-Mumby!")
				return blackwood(X, i)

	# returns the model run for a certain set of parameters 
	# maybe create a subclass for a model run to separate plotting jobs ? 
	def run_model(IC_set, t, closure_length, f, m, n, poaching):
		sol = odeint(patch_system, IC_set, t, args = (closure_length, f/(1-m/n), m, n, poaching))
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

		'''

		# slow version 

		coral_array_LO =  np.zeros(shape=(self.n,self.n))
		
		period_array = np.empty(n)
		m_array = np.empty(n)
		
		for period in range(1,n+1):
			for m in range(n):
				displacer = 1/(1-m/float(n))
				final_coral = odeint(patch_system, X1, t, args = (period, displacer*float(self.f), m, self.n, 0)))
				avg = 0
				for year in range(MAX_TIME - (MAX_TIME % (self.n*period)) - 2*(self.n*period), MAX_TIME - (MAX_TIME % (self.n*period))):
					avg += final_coral[year][n]
				avg = avg / (2*(period*self.n) + 1)
				coral_array_LO[period-1][m] = avg2
				period_array[period-1] = period
				m_array[m] = m

		plt.title('coral end state', fontsize = 20)
		sb.heatmap(coral_array, vmin = 0.0, vmax = 1.0, cmap = "viridis", annot = show_labels) #YlGnBu for original color scheme
		plt.ylabel('period', fontsize = 10) # x-axis label with fontsize 15
		plt.xlabel('number of closures', fontsize = 10) # y-axis label with fontsize 15
		plt.yticks(rotation=0)
		plt.show()
		'''


		return None

	def bistable_zone():
		return None
	def find_unstable_equilibrium():
		return None
	def scenario_plot():
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

	def fishing(parrotfish, f):
		return f/(1+math.exp(-steepness*(parrotfish-shift)))

	def blackwood(X, i):
		P, C, M = X.reshape(3, n)

		dPs[i] = s*P[i]*(1 - (P[i] / (beta*K(C[i])))) - f*P[i]*square_signal(t, closure_length, i, m, n, poaching)
		dCs[i] = r*(1-M[i]-C[i])*C[i]-d*C[i] - a*M[i]*C[i] + i_C*(1-M[i]-C[i])
		# need to define g(P)
		dMs[i] = a*M[i]*C[i] - g(P[i])*M[i] *(1/(1-C[i])) + gamma*M[i]*(1-M[i]-C[i])+i_M*(1-M[i]-C[i])

		return np.concatenate((dPs, dCs, dMs), axis=0)
	# how to pass dPs dCs and dMs to this function? 
	def leemput(X, i):
		# check input 
		P,C,M = X.reshape(3, n) # will reshaping work since we are passing arrays of length n? 
		dPs[i] = P_influx[i]+ s*P[i]*(1 - (P[i] / K(sigma,C[i]))) - fishing(P[i], f)*P[i] *(square_signal(t, closure_length, i, m, n, poaching))
		dCs[i] = (i_C + r*C[i])*(1-M[i]-C[i])*(1-alpha*M[i]) - d*C[i]
		
		dMs[i] = (i_M+gamma*M[i])*(1-M[i]-C[i])-g*M[i]*P[i]/(g*eta*M[i]+1)

		#concatenate into 1D vector to pass to next step
		return np.concatenate((dPs, dCs, dMs), axis=0)

	def rass_briggs(X, i):

		P, C, Mv, Mi = X.reshape(4, n)
		T = 1 - C - Mv - Mi 
		dPdt = rH*P*(1-P/K) - fishing(P[i], f)*P[i] *(square_signal(t, closure_length, i, m, n, poaching))
		dCdt = (phiC*T) + gTC*T*C - gamma*gTI*Mi*C - dC*C
		dMvdt = phiM*T + rM*T*Mi + gTV*T*Mv - dV*Mv - P*Mv*Graze - omega * Mv
		# conceptual question: why is that second term multiplied by M not Mv? 
		dMidt = omega*Mv + gTI*T*Mi + gamma*gTI*Mi*C - dI*Mi
		return [dPdt, dCdt, dMvdt, dMidt]

		# check input
		return None 

	def load_parameters(self, model_type):
		if model_type == 'vdL':
			# import vdL parameters from parameters folder 
			# maybe store as a dictionary? 
		str = "name"
		 
		exec("self.%s = %d" % (str, val))
 
		return None 



x = Model('vdL', 2, 1)
print(x.run_model(parameter list))
x.time_series(parameter list, show = True)