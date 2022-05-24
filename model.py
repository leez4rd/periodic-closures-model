import numpy as np
from scipy.integrate import odeint  
import matplotlib.pyplot as plt 

# how do we solve the global variable problem ? 
# how do we solve the parameter problem? 
class Model:
	# hmm
	# variable arguments based on model type? or just use string to specify parameters...
	def initialize_patch_model(n, frac_nomove):
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
		
		X1 = [P0]*n + [C0L]*n + [M0H]*n
		X2 = [P0]*n + [C0H]*n + [M0L]*n

	def __init__(self, model_type, n, frac_nomove): 

		# instead of global variables, we can initialize the arrays here
		# so that we could just do self.P[i] for the ith element of the parrotfish array 
		# how can we just load the parameters conditionally from another file 

		if model_type == 'RB':
			# load_parameters() -- hypothetical function to load all params of RB model into this object  
			load_parameters('RB')
			# RB_initialize_patch_model(n, frac_nomove)

		elif model_type == 'BM':
			# initialize_patch_model(n, frac_nomove)
			load_parameters('BM')
		elif model_type == 'vdL_PC':
			load_parameters('vdL_PC')
			# initialize_patch_model(n, frac_nomove)
		elif model_type == 'vdL_MP':
			load_parameters('vdL_MP')
			# initialize_patch_model(n, frac_nomove)
		elif model_type == 'vdL_MC':
			# initialize_patch_model(n, frac_nomove)
			load_parameters('vdL_MC')
		elif model_type == 'vdL': # all feedbacks active 
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

			#concatenate initial condition arrays -- we can put this somewhere else 
			
			X1 = [P0]*n + [C0L]*n + [M0H]*n
			X2 = [P0]*n + [C0H]*n + [M0L]*n

			# initialize_patch_model(n, frac_nomove)
			# load_parameters('vdL')
		else:
			print("Bad input, defaulting to Blackwood-Mumby!")
			# initialize_patch_model(n, frac_nomove)
			# load_parameters('BM')

	

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
	def time_series(save, show, IC_set, t, closure_length, f, m, n, poaching):
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
	def coral_recovery_map():
		'''
		patches = 20

		#heatmap of period vs m -- a bit messy right now 
		coral_array_HI =  np.zeros(shape=(patches,patches))
		coral_array_LO =  np.zeros(shape=(patches,patches))
		coral_array_AVG =  np.zeros(shape=(patches,patches))
		period_array = np.empty(patches)
		m_array = np.empty(patches)
		fishin = 0.35
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



x = Model('vDL', 2, 1)